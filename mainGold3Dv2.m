%% GOLD TRACKiNG
clear 
clc
close all

%% User input

path2Cal  = 'E:\Users\Boris\Documents\TmpData - MuEn\data2Boris\2DCal';
file.path = 'E:\Users\Boris\Documents\TmpData - MuEn\data2Boris\ensembled data\data\mov1';
file.ext  = '.ome.tif';

focusPlane = 1;%=2 af
width.xy = 3; %for fitting (3 for 200nm beads, 400 nm beads to be determined, 0 to let the code find the width)
width.z  = 4;%see above
nParticles = 5;%number of particles expected in the movie has to be exact
minDist = 6; %in pixels (min distant expected between particles
pxSize = 95;%in nm
cropRadius = 40; %cut the frame to reduce the amount of data to be fitted
info.type = 'normal';%Transmission or normal movie
info.runMethod = 'load';
info.calibrate = true;
toAnalyze = '.ome.tif';%accepted: .mp4, .ome.tif, folder. (folder that contain folders of .ome.tif.
info.checkSync = true;
outputFolder = 'Results';%name of the folder to output the results
%% Loading
%we get the zCalibration directory
folder2Mov = dir(file.path);
folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
%loop through the content of the directory
mov = struct();
count = 0;
for i = 3:size(folder2Mov,1)
    %Check if the directory
    folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
    file2Analyze = Core.Movie.getFileInPath(folderPath,file.ext);

    if ~isempty(file2Analyze)
        count = count+1;
        filetmp.path = file2Analyze.folder;
        filetmp.ext  = file.ext;
        tmp = Core.MPMovie(filetmp , path2Cal,info);
        
        if count == 1
            tmp.giveInfo;
        else
            tmp.info = mov.g1.getInfo; 
        end
        tmp.calibrate;
        mov.(['g' num2str(i-2)]) = tmp;


    else

        warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);

    end

end
disp('=======> DONE ! <========')
%mov.g1.showFrame(1,5);

%% detection of the center of the beads

%preallocate memory for storing data of the different files
allData = struct('fileName',[],'traces',[]);
allData(length(mov)).locPos = [];
fields = fieldnames(mov);
frame = 35;

for i = 1:length(fields)
    currMov = mov.(fields{i});
    %check detection
    fullStack = currMov.getFrame(frame);
    %revert the intensity scale
    if strcmpi(info.type,'transmission')
        testFrame = imcomplement(fullStack);
    else
        testFrame = fullStack;
    end

    frame = 5;
    %get the n maxima where n is the number of particles expected and
    %minDist is the distance expected between them
    [pos] = goldProj.nMaxDetection (testFrame(:,:,focusPlane),nParticles,minDist);

    x0 = pos(:,2);
    y0 = pos(:,1);

    pos = [y0 x0];
    %calculate center of mass of particles position for cropping later
    cropPos = round(mean(pos,1));
    %plot to check that the max detection worked
    
    %crop and show the cropping so user can see if it is okay
    testFrame = testFrame(cropPos(:,1)-cropRadius:cropPos(:,1)+cropRadius,...
            cropPos(:,2)-cropRadius:cropPos(:,2)+cropRadius,:);
    figure
    imagesc(testFrame(:,:,focusPlane))
    hold on
    plot(pos(:,2),pos(:,1),'r+')
    
    
    %% Fitting
    nFrames = currMov.raw.movInfo.maxFrame(1);
    %get the X Y dommain
    x = (1:size(testFrame,2))*pxSize/1000;
    y = (1:size(testFrame,1))*pxSize/1000;
    z = currMov.calibrated.oRelZPos;
    [domX,domY,domZ] = meshgrid(x,y,z);
    dom(:,:,:,1) = domX;
    dom(:,:,:,2) = domY;
    dom(:,:,:,3) = domZ;
    %preallocate memory
       
    data2Store = zeros(nFrames,3,nParticles);
    fitMov = zeros(100,100,nFrames);
    h = waitbar(0,'Fitting Data');%create waiting bar
    %unSortedData =[];
    planePos = currMov.calibrated.oRelZPos;
    widths.xy = zeros(1,nFrames);
    widths.z  = zeros(1,nFrames);
    width.xy = width.xy*pxSize/1000;
    width.z  = width.z*mean(diff(planePos));
    
    for j = 1:nFrames
        currentFrame = double(currMov.getFrame(j));
        currentFrame = currentFrame(cropPos(:,1)-cropRadius:cropPos(:,1)+cropRadius,...
            cropPos(:,2)-cropRadius:cropPos(:,2)+cropRadius,:);
        %inital detection of particles on currentFrame
        [pos] = goldProj.nMaxDetection (currentFrame(:,:,focusPlane),nParticles,minDist);
 
        x0 = pos(:,2)*pxSize/1000;
        y0 = pos(:,1)*pxSize/1000;
       
        [~,idx] = max(max(max(currentFrame)));
        z0 = ones(size(x0))*currMov.calibrated.oRelZPos(focusPlane);
        
        %Multiple gaussian fitting occurs here
        [gPar,resnorm,res] = Localization.Gauss.MultipleGFit3D(currentFrame,x0,y0,z0,dom,nParticles,width); 
        
        widths.xy(j) = gPar(2);
        widths.z(j) = gPar(3);
        g = reshape(gPar(5:end),[3,nParticles]);
        g = g';
        
        nPart = size(x0,1);
        z = g(:,3);
        %Get Z position
        for k = 1:nPart
            
            %domain = 1:size(currentFrame,3);

%             if or(z(k)<0,z(k)>max(domain))
%                 z(k)   = NaN;                           
%             else
%                 tmpZ = floor(z(k));
%                 fracZ = z(k)-tmpZ;
%                 z(k) = planePos(tmpZ)+fracZ*(planePos(tmpZ+1) - planePos(tmpZ));
%                 z(k) = z(k)*1000;
%             end
            
        end
        
        if j>1
            %Tracking based on MSD minimization 
            newOrder = goldProj.simpleTracking(gPar(5:end),prevPos,'3D');
            %reshaping to format the final data and sorting with new order
            gPos = reshape(gPar(5:end),[3,nParticles]);
            gPos = gPos(:,newOrder);
            prevPos = reshape(gPos,[1,nParticles*3]);
            gPos = reshape(gPos,[1,3,nParticles]);
            gPos(:,3,:) = z;
            %store data
            data2Store(j,:,:) = gPos;
            clear gPos

        else
            %First frame we just reshape
            gPos = reshape(gPar(5:end),[3,nParticles]);
            gPos = reshape(gPos,[1,3,nParticles])*pxSize;
            gPos(:,3,:) = z;
            data2Store(j,:,:) = gPos;
            %store the gPar in the prev for checking particle order
            prevPos = gPar(5:end);
            clear gPos
        end
        %update waitbar value
        waitbar(j/nFrames,h,'Fitting Data')
    end
    
    %clear waitbar
    close(h);
    disp(['The Average fitted width in XY is: ' num2str(mean(widths.xy))])
    disp(['The Average fitted width in z is: ' num2str(mean(widths.z))])
    currentPath = currMov.raw.movInfo.Path;
    %save data to the current folder being analyze
    filename = [currentPath filesep 'LocalizationData.mat'];
    save(filename,'data2Store');
    %store data in allData
    allData(i).traces = data2Store*pxSize;
    allData(i).fileName = currentPath;  
  

end

%% convert Data to table
figure
hold on
for i = 1: size(data2Store,3)
    plot3(data2Store(:,1,i)/pxSize,data2Store(:,2,i)/pxSize,data2Store(:,3,i));
    
end
axis image
view(3)

%%
trackRes = save.convertData2TrackRes(allData,nParticles);

filename = [file.path filesep 'trackRes.mat'];
save(filename,'trackRes');
h = msgbox('Data succesfully saved');