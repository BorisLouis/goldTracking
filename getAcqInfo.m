%% getAcqInfo
% Extract exposure time, mean time between frames and fps from the OME-TIF
% metadata of every movie in a folder, and save them in a .csv.
% Only the metadata is read (no frame is loaded), so it runs in seconds.
%
% Expected structure: <selected folder>/<subfolder>/<movie>.ome.tif
% Run from (or add to the path) the goldTracking folder so +Core and +Load are found.
clear
clc

path = uigetdir();
assert(~isequal(path,0),'No folder selected');

sub = dir(path);
sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));   % subfolders only

rows = {};
for i = 1:numel(sub)
    f = Core.Movie.getFileInPath(fullfile(sub(i).folder,sub(i).name),'.ome.tif');
    if isempty(f)
        fprintf('No .ome.tif in %s, skipped\n',sub(i).name);
        continue
    end
    f = f(1);   % one movie per subfolder, as in mainTracking

    % metadata only, no camera sync check (same as info.checkSync = false)
    [fInfo,mInfo] = Load.Movie.ome.getInfo(fullfile(f.folder,f.name),false);
    if iscell(fInfo)                       % multi-image tiff
        fInfo = Load.Movie.ome.combineFrameInfo(fInfo);
    end

    cam1 = strcmp({fInfo.C},num2str(mInfo.Cam(1)));   % frames are interleaved between cameras
    t    = [fInfo(cam1).time];                        % ms, relative to first frame

    expT = str2double(string(mInfo.expT));            % exposure time (ms)
    dt   = (max(t)-min(t))/(numel(t)-1);              % mean time between frames (ms)
    fps  = 1000/dt;

    rows(end+1,:) = {sub(i).name,f.name,expT,dt,fps}; %#ok<SAGROW>
end

assert(~isempty(rows),'No .ome.tif found in any subfolder');
T = cell2table(rows,'VariableNames',{'folder','file','expT','dt','fps'});
writetable(T,fullfile(path,'AcquisitionInfo.csv'));
disp(T)
