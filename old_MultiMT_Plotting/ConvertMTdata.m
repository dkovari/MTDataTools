function [OutPath] = ConvertMTdata(FilePath)

%% GUI File Selection
persistent last_dir;

if nargin>0
    [Dir,File,ext] = fileparts(FilePath);
    File = [File,ext];
else
    %select file
    [File,Dir] = uigetfile(fullfile(last_dir,'*_ForceExtension_*.txt'),'Select Force Extension File');
    if File==0
        return
    end
end
if ~isempty(Dir)
    last_dir = Dir;
end

[~,name,~] = fileparts(File);

%% Load Processed Data Info
fid = fopen(fullfile(Dir,File),'r');

reading_comments = false;
Comments = [];
while ~feof(fid)
    tline = fgetl(fid);
    if strncmp(tline,'DATE',4)
        break;
    end
    if reading_comments
        Comments = [Comments,sprintf('\n'),tline];
        continue;
    end
    if strncmp(tline,'PxScale',7)
        PxScale = sscanf(tline,'PxScale:\t%f',1);
    end
    if strncmp(tline,'MagnetRotation:',15)
        MagRot = sscanf(tline,'MagnetRotation:\t%f',1);
    end
    
    if strncmp(tline,'FrameCount',10)
        FrameCount = sscanf(tline,'FrameCount:\t%f',1);
    end
    if strncmp(tline,'TotalTracks',11)
        num_tracks = sscanf(tline,'TotalTracks:\t%d',1);
    end
    if strncmp(tline,'RefenenceTracks',15)
        ref_tracks = str2num(tline(17:end));
    end
    if strncmp(tline,'Comments:',9)
        reading_comments = true;
        Comments = [Comments,tline(10:end)];
    end
    
        
end
if feof(fid)
    error('could not find any data in the file');
end

%loop over lines and collect data
DateVec = [];
MagH = [];
Data = [];
while ~feof(fid)
    tline = fgetl(fid);
    %                      yr- mm- dd   hh: mm: ss  mag
    ldata = sscanf(tline,'%4d-%2d-%2d\t%2d:%2d:%6f\t%5f',7);
    if numel(ldata)~=7
        disp('could not read line');
        break;
    end
    fdata = str2num(tline(31:end));
    fdata = reshape(fdata,1,5,[]); %rehape data to [[L1,Fx1,Fx2,dx1,dy1],[L2,Fx2,Fx2,dx2,dy2],...]
    DateVec = [DateVec;ldata(1:6)'];
    MagH = [MagH;ldata(7)];
    Data = [Data;fdata];
end
fclose(fid);

if size(Data,3)~=num_tracks
    error('size of data does not match number of tracks specified');
end

MeasuredTracks = 1:num_tracks; %make a list of measurement tracks
MeasuredTracks(ref_tracks) = [];
%% Load XY Data
file = fullfile(Dir,[name,'_XY.bin']);
fid=fopen(file,'r');
num_tracks = fread(fid,1,'uint8');
num_ref_tracks = fread(fid,1,'uint8');
if num_ref_tracks>0
    ref_tracks = fread(fid,num_ref_tracks,'uint8');
end
%data is save in row-order
%XYdata(frame1,:) = [time,X1,Y1,X2,Y2,...Xn,Yn]
idx = 1;
Time = [];
X = NaN(0,num_tracks);
Y = NaN(0,num_tracks);
while ~feof(fid)
     t = fread(fid,1,'double');
     if isempty(t)
         break;
     end
    Time(idx) = t;
    for trk=1:num_tracks
        X(idx,trk) = fread(fid,1,'double');
        Y(idx,trk) = fread(fid,1,'double');
    end
    idx=idx+1;
end
Time = Time';
fclose(fid);

%% Load Z data
file = fullfile(Dir,[name,'_Z.bin']);
fid=fopen(file,'r');
num_tracks = fread(fid,1,'uint8');
num_ref_tracks = fread(fid,1,'uint8');
if num_ref_tracks>0
    ref_tracks = fread(fid,num_ref_tracks,'uint8');
end
%get refID for each track (specifies which track was used for Z-referencing
refID = fread(fid,num_tracks,'uint8');
idx = 1;
TimeZ = [];
Z_rel = NaN(0,num_tracks);
Z_abs = NaN(0,num_tracks);
while ~feof(fid)
     t = fread(fid,1,'double');
     if isempty(t)
         break;
     end
    TimeZ(idx) = t;
    for trk=1:num_tracks
        Z_rel(idx,trk) = fread(fid,1,'double');
        Z_abs(idx,trk) = fread(fid,1,'double');
    end
    idx=idx+1;
end
TimeZ = TimeZ';
fclose(fid);

dZ = Z_rel(:,refID) - Z_rel;

%% Create Header
%config data
%prepare file
Config.FileType = 'Force Extension';
Config.CreationDate = datestr(Time(1),'yyyy-mm-dd HH:MM');
Config.Comments = Comments; %ADD COMMENTS CODE HERE
%Hardware Config
Config.InstrumentName = 'NA';
Config.CameraInterface = 'NA';
Config.PxScale = PxScale;
Config.PiezoController = 'NA';
Config.PiezoCOM = 'NA';
Config.PiezoBAUD = 'NA';
Config.MotorController = 'NA';
Config.MotorCOM = 'NA';
Config.MotorBAUD = 'NA';
Config.magztype = 'NA';
Config.magzaxis = 'NA';
Config.magrtype = 'NA';
Config.magraxis = 'NA';
Config.mag_rotscale = NaN;
Config.TemperatureController = 'NA';
Config.TemperatureUnits = 'C';
Config.Temperature = 25;
Config.LogFile = '';

Config.CalibrationMin = NaN;
Config.CalibrationMax = NaN;
Config.CalibrationStep = NaN;

Config.num_tracks = num_tracks;
%Track Info
Config.TrackingInfo = struct('Type',{},'Radius',{},'IsCalibrated',{},'ZRef',{},'Window',{});
for n=1:num_tracks
    if any(n==ref_tracks)
        Config.TrackingInfo(n).Type = 'Reference';
    else
        Config.TrackingInfo(n).Type = 'Measurement';
    end
    Config.TrackingInfo(n).Radius = NaN;
    Config.TrackingInfo(n).IsCalibrated = true;
    Config.TrackingInfo(n).ZRef = refID(n);
    Config.TrackingInfo(n).Window = NaN(1,4);
    Config.TrackingInfo(n).CalibrationPositions = NaN;
end

%Camera Settings
Config.CameraSettings.FrameRate = NaN;
Config.CameraSettings.UseFrameRate = false;
Config.CameraSettings.ResultingFrameRate = NaN;
Config.CameraSettings.Expsoure = NaN;
Config.CameraSettings.Gain = NaN;
Config.CameraSettings.ExposureAuto = NaN;
Config.CameraSettings.GainAuto = NaN;
Config.CameraSettings.TargetBrightness = NaN;
Config.CameraSettings.BlackLevel = NaN;
Config.CameraSettings.ROI = NaN;



Config.ForceExtensionScheme.MagnetRotation = MagRot;
Config.ForceExtensionScheme.MagnetHeightPositions = MagH;
Config.ForceExtensionScheme.FrameCount = FrameCount;

%% record structure
Record = struct('parameter',{},'format',{},'size',{});

Record(end+1).parameter = 'Date';
Record(end).format = 'double';
Record(end).size = [1,1];

Record(end+1).parameter = 'Step';
Record(end).format = 'uint32';
Record(end).size = [1,1];

Record(end+1).parameter = 'FrameCount';
Record(end).format = 'double';
Record(end).size = [1,1];

Record(end+1).parameter = 'ObjectivePosition';
Record(end).format = 'double';
Record(end).size = [1,1];

Record(end+1).parameter = 'MagnetHeight';
Record(end).format = 'double';
Record(end).size = [1,1];

Record(end+1).parameter = 'MagnetRotation';
Record(end).format = 'double';
Record(end).size = [1,1];

Record(end+1).parameter = 'X';
Record(end).format = 'double';
Record(end).size = [num_tracks,1];

Record(end+1).parameter = 'Y';
Record(end).format = 'double';
Record(end).size = [num_tracks,1];

Record(end+1).parameter = 'Z_REL';
Record(end).format = 'double';
Record(end).size = [num_tracks,1];

Record(end+1).parameter = 'Z_ABS';
Record(end).format = 'double';
Record(end).size = [num_tracks,1];

Record(end+1).parameter = 'dZ';
Record(end).format = 'double';
Record(end).size = [num_tracks,1];

Record(end+1).parameter = 'UsingTilt';
Record(end).format = 'int8';
Record(end).size = [1,1];


%% Create file
OutPath = fullfile(Dir,[name,'.mtdat']);
[mtfile,Record] = mtdat.writeheader(OutPath,Config,Record);

%% Assemble data into mtdat file
DateNum = datenum(DateVec);
nT = numel(DateNum);
for t=1:nT %loop over each step
    
    if t==nT
        nextT = Inf;
    else
        nextT = DateNum(t+1);
    end
    
    indT = find(Time>=DateNum(t)&Time<nextT);
    indTz = find(TimeZ>=DateNum(t)&TimeZ<nextT);
    nF = min(numel(indT),numel(indTz));
    indT = indT(1:nF);
    indTz = indTz(1:nF);
    
    for f=1:nF %loop over each frame
        %setup record
        thisRecord.Date = Time(indT(f));
        thisRecord.Step = t;
        thisRecord.FrameCount = f;
        thisRecord.ObjectivePosition = NaN;
        thisRecord.MagnetHeight = MagH(t);
        thisRecord.MagnetRotation = MagRot;
        
        thisRecord.X = X(indT(f),:)';
        thisRecord.Y = Y(indT(f),:)';
        thisRecord.Z_REL = Z_rel(indTz(f),:)';
        thisRecord.Z_ABS = Z_abs(indTz(f),:)';
        
        thisRecord.dZ = dZ(indTz(f),:)';
        
        thisRecord.UsingTilt = false;
        %write record to file
        mtdat.writerecord(mtfile,Record,thisRecord);
    end
    
end
%close the file
fclose(mtfile);


