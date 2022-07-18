clear 
clc
close all
%% User input
delta = 50;% in px Size of the ROI around particles detected(radius 50 = 100x100 pixel
nParticles = 2;%number of particles expected in the movie has to be exact
width = 5; %for fitting (3 for 200nm beads, 400 nm beads to be determined, 0 to let the code find the width)

pxSize = 95;%in nm
minDist = 6; %in pixels (min distant expected between particles
scaleBar = 2; %in um
tail = 20;%Length of the tail in frames, for plotting the traces on top of the movie
frameRate = 30; %for saving the movie
info.type = 'normal';%Transmission or normal movie
toAnalyze = '.ome.tif';%accepted: .mp4, .ome.tif, folder. (folder that contain folders of .ome.tif.
outputFolder = 'Results';%name of the folder to output the results
%% Loading
switch toAnalyze %switch depending on what user want to analyze
    case '.mp4'
        %load all .mp4 in the selected folder
        [folder2Mov,path,outDir] = Load.Folder(toAnalyze,outputFolder);
        mkdir(outDir);
    case '.ome.tif'
        %load all .ome.tif in the selected folder
        [folder2Mov,path,outDir] = Load.Folder(toAnalyze,outputFolder);
        mkdir(outDir);

    case 'folder'
        [path] = uigetdir();%allow user to select folder
        tmp = dir(path);%make it a directory in matlab (list of all thing contained inside)
        tmp = tmp(cell2mat({tmp.isdir}));%only keep the subfolders
        tmp = tmp(3:end);%remove the access to parent folders
        folder2Mov = [];
        for i = 1:size(tmp,1)
            path2File = [tmp(i).folder filesep tmp(i).name];%get path to subfolders
            file2Analyze = Core.Movie.getFileInPath(path2File,'.ome.tif');%get .ome.tif in subfolder
            folder2Mov = [folder2Mov file2Analyze];%store the info to the file if it found one (otherwise file2Analyze is empty)
        end
        %create the out directory
        outDir = [path filesep 'Results']; 
        mkdir(outDir)
    otherwise
        error('Unknown extension');%throw an error if none of the above possibilities is given

end
%check that folder2Mov is not empty
assert(~isempty(folder2Mov),'Error no file was found, check that you put the correct analysis type');
%% detection of the center of the beads
%preallocate memory for storing data of the different files
allData = struct('fileName',[],'locPos',[]);
allData(size(folder2Mov,2)).locPos = [];

for i =1: size(folder2Mov,2)
    
    file.path = folder2Mov(i).folder;
    file.ext  = '.ome.tif';
    switch toAnalyze
        case '.mp4'
            p2file = [folder2Mov(i).folder filesep folder2Mov(i).name];
            v = VideoReader(p2file);%Create VideoReader object
            nFrames = floor(v.Duration*v.FrameRate);%extract the number of frames
            fullStackIn = zeros(v.Height,v.Width,nFrames);%preallocate memory
            for j = 1:nFrames
                frame = readFrame(v);%read frames

                fullStackIn(:,:,j) = rgb2gray(frame);%extract the intensity data from frames
            end
        otherwise
            myMov = Core.Movie(file,info);%Create Movie Object
            fullStack = myMov.getFrame(1);%extract first frame
            frame = fullStack.Cam1;
            %check if cropping is necessary
            if size(frame,2) > 400
                myMov.cropIm;
            end
            %load full stack
            fullStack = myMov.getFrame;
            %revert the intensity scale
            if strcmpi(info.type,'transmission')
                fullStackIn = imcomplement(fullStack.Cam1);
            else
                fullStackIn = fullStack.Cam1;
            end
    end
    
    frame = 5;
    %get the n maxima where n is the number of particles expected and
    %minDist is the distance expected between them
    [pos] = goldProj.nMaxDetection (fullStackIn(:,:,frame),nParticles,minDist);
    
    x0 = pos(:,2);
    y0 = pos(:,1);

    pos = [y0 x0];
    %calculate center of mass of particles position for cropping later
    cropPos = round(mean(pos,1));
    %plot to check that the max detection worked
    figure
    imagesc(fullStackIn(:,:,frame))
    hold on
    plot(pos(:,2),pos(:,1),'r+')

    %% Cropping Movie
    %crop the desired region around the center of mass of particles located
    %using delta provided by user
    %fullStackIn = fullStackIn(cropPos(1)-delta:cropPos(1)+delta, cropPos(2)-delta:cropPos(2)+delta,:);
    %% Fitting
    
    nFrames = size(fullStackIn,3);
    %get the X Y dommain
    x = 1:size(fullStackIn,2);
    y = 1:size(fullStackIn,1);
    [domX,domY] = meshgrid(x,y);
    dom(:,:,1) = domX;
    dom(:,:,2) = domY;
    %preallocate memory
    data2Store = zeros(nFrames,2,nParticles);
    fitMov = zeros(100,100,size(fullStackIn,3));
    h = waitbar(0,'Fitting Data');%create waiting bar
    unSortedData =[];
    for j = 1:nFrames
        currentFrame = double(fullStackIn(:,:,j));
        %inital detection of particles on currentFrame
        [pos] = goldProj.nMaxDetection (currentFrame,nParticles,minDist);
%         
%         figure(1)
%         imagesc(currentFrame)
%         hold on
%         scatter(pos(:,2),pos(:,1));
%  
        x0 = pos(:,2);
        y0 = pos(:,1);
        %Multiple gaussian fitting occurs here
        [gPar,resnorm,res] = Localization.Gauss.MultipleFitting(currentFrame,x0,y0,dom,nParticles,width); 
        g = reshape(gPar(5:end),[2,nParticles]);
        g = g';
        unSortedData = [unSortedData; g zeros(size(g,1),1) ones(size(g,1),1)*j  ];
        
%         figure(1)
%         imagesc(currentFrame)
%         hold on
%         scatter(g(:,1),g(:,2));
%         
        
        %Generate an image of the Fit to be able to plot in case we want to
        %check
        xHRes = linspace(1,size(fullStackIn,2),100);
        yHRes = linspace(1,size(fullStackIn,1),100);
        [domHResX,domHResY] = meshgrid(xHRes,yHRes);
        domHRes(:,:,1) = domHResX;
        domHRes(:,:,2) = domHResY;
        F = Localization.Gauss.MultipleGauss(gPar, domHRes,nParticles);
        
        if j>1
            %Tracking based on MSD minimization 
            newOrder = goldProj.simpleTracking(gPar(5:end),prevPos,'2D');
            %reshaping to format the final data and sorting with new order
            gPos = reshape(gPar(5:end),[2,nParticles]);
            gPos = gPos(:,newOrder);
            prevPos = reshape(gPos,[1,nParticles*2]);
            gPos = reshape(gPos,[1,2,nParticles]);
            %store data
            data2Store(j,:,:) = gPos;
            clear gPos

        else
            %First frame we just reshape
            gPos = reshape(gPar(5:end),[2,nParticles]);
            gPos = reshape(gPos,[1,2,nParticles]);
            data2Store(j,:,:) = gPos;
            %store the gPar in the prev for checking particle order
            prevPos = gPar(5:end);
            clear gPos
        end
        %store fittings
        fitMov(:,:,j) = F;
        %update waitbar value
        waitbar(j/nFrames,h,'Fitting Data')
    end
    %save data to the current folder being analyze
    filename = [file.path filesep 'LocalizationData.mat'];
    save(filename,'data2Store');
    %store data in allData
    allData(i).locPos = data2Store;
    allData(i).fileName = file.path;
    %clear waitbar
    close(h);
    
    %% display figure
    %switch from pixel to nm
    data2plot = data2Store *pxSize;
    %calculate mean center of mass of all particles
    cm = [mean(mean(data2plot(:,1,:))) mean(mean(data2plot(:,2,:)))];
    %plot localization data
    Fig = figure;
    hold on
    for j = 1 : nParticles
        scatter(data2plot(:,1,j)-cm(1),data2plot(:,2,j)-cm(2),20,'filled')

    end
    axis image
    xlabel('X Position (nm)')
    ylabel('Y Position (nm)')
    title('All localized spot');
    %save the figure in the current folder path
    filename = [file.path filesep 'LocalizationDensity.fig'];
    saveas(Fig,filename);
    
%% MovieMaker
%save a movie where the traces is displayed on top of the image
filename = [file.path filesep 'TrackMovie.gif'];
goldProj.makeTraceMovie(data2Store,fullStackIn,filename,frameRate,scaleBar,tail);

end
%save all Data in the master folder
filename = [outDir filesep 'allLoc'];
save(filename,'allData');
