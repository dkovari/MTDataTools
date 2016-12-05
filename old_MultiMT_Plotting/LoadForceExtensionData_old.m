function LoadForceExtensionData_old(FilePath)
%Load data saved from the Pre- September 2016 version of MultiMT

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

while ~feof(fid)
    tline = fgetl(fid);
    if strncmp(tline,'PxScale',7)
        PxScale = sscanf(tline,'PxScale:\t%f',1);
    end
    if strncmp(tline,'DATE',4)
        break;
    end
    if strncmp(tline,'TotalTracks',11)
        num_tracks = sscanf(tline,'TotalTracks:\t%d',1);
    end
    if strncmp(tline,'RefenenceTracks',15)
        ref_tracks = str2num(tline(17:end));
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

%% Recalculate Force with Errorbars
kBT=1.380648813e-23*300;
DateNum = datenum(DateVec);
nT = numel(DateNum);
Fx = NaN(nT,num_tracks);
FxLower95 = NaN(nT,num_tracks);
FxUpper95 = NaN(nT,num_tracks);
Lext = NaN(nT,num_tracks);
LextStdErr = NaN(nT,num_tracks);
Xvar = NaN(nT,num_tracks);
for t=1:nT
    
    if t==nT
        nextT = Inf;
    else
        nextT = DateNum(t+1);
    end
    
    thisX = PxScale*X( Time>=DateNum(t)&Time<nextT, MeasuredTracks);
    thisY = PxScale*Y( Time>=DateNum(t)&Time<nextT, MeasuredTracks);
    thisdZ = dZ( TimeZ>=DateNum(t)&TimeZ<nextT, MeasuredTracks);
    
    L = sqrt( ...
            bsxfun(@minus,thisX,nanmean(thisX,1)).^2 + ...
            bsxfun(@minus,thisY,nanmean(thisY,1)).^2 + ...
            thisdZ.^2 );
    Lext(t,MeasuredTracks) = nanmean(L,1);
    LextStdErr(t,MeasuredTracks) = nanstd(L,0,1)/sqrt(numel(thisX));
    
    varX = var(thisX,0,1,'omitnan');
    Fx(t,MeasuredTracks) = kBT*Lext(t,MeasuredTracks)./varX*10^6;
    
    FxLower95(t,MeasuredTracks) = ...
        kBT*10^6*...
        (Lext(t,MeasuredTracks)-2*LextStdErr(t,MeasuredTracks))./ ((numel(thisX)-1)*varX/chi2inv(0.05/2,(numel(thisX)-1)));
    FxUpper95(t,MeasuredTracks) = ...
        kBT*10^6*...
        (Lext(t,MeasuredTracks)+2*LextStdErr(t,MeasuredTracks))./ ((numel(thisX)-1)*varX/chi2inv(1-(.05/2),(numel(thisX)-1)));
    
    
    Xvar(t,MeasuredTracks) = varX/(PxScale.^2);
end




%% Export Data
uiextras.putvar(Data,MagH,DateVec,MeasuredTracks);
uiextras.putvar(Time,X,Y,TimeZ,Z_rel,Z_abs,dZ,num_tracks,refID,PxScale,FxLower95,FxUpper95,LextStdErr);

uiextras.putvar(Fx,Lext,Xvar);

%% Plot Raw Data

% x-y scatter
[~,hAx,~,hFig] = uiextras.plot_selectable(X,Y,uiextras.cell_sprintf('Track %d',1:num_tracks));
xlabel(hAx,'X [px]');
ylabel(hAx,'Y [px]');
title(hAx,'X-Y Scatter');
axis(hAx,'equal');
set(hFig,'NumberTitle','off',...
         'Name',[name,':XY Scatter']);

% x time series
[~,hAx,~,hFig] = uiextras.plot_selectable(repmat(Time,1,num_tracks), X, uiextras.cell_sprintf('Track %d',1:num_tracks));
xlabel(hAx,'Time [sec]');
ylabel(hAx,'X [px]');
title(hAx,'X Trace');
set(hFig,'NumberTitle','off',...
         'Name',[name,':X Trace']);

% y time series
[~,hAx,~,hFig] = uiextras.plot_selectable(repmat(Time,1,num_tracks), Y, uiextras.cell_sprintf('Track %d',1:num_tracks));
xlabel(hAx,'Time [sec]');
ylabel(hAx,'Y [px]');
title(hAx,'Y Trace');
set(hFig,'NumberTitle','off',...
         'Name',[name,':Y Trace']);

% z_rel time series
[~,hAx,~,hFig] = uiextras.plot_selectable(repmat(TimeZ,1,num_tracks), Z_rel, uiextras.cell_sprintf('Track %d',1:num_tracks));
xlabel(hAx,'Time [sec]');
ylabel(hAx,'Z, Relative to Reference [µm]');
title(hAx,'Z-Relative Trace');
set(hFig,'NumberTitle','off',...
         'Name',[name,':Zrel Trace']);

% z_abs time series
[~,hAx,~,hFig] = uiextras.plot_selectable(repmat(TimeZ,1,num_tracks), Z_abs, uiextras.cell_sprintf('Track %d',1:num_tracks));
xlabel(hAx,'Time [sec]');
ylabel(hAx,'Z, Relative to Self [µm]');
title(hAx,'Z-Absolute Trace');
set(hFig,'NumberTitle','off',...
         'Name',[name,':Zabs Trace']);

% dZ time series
[~,hAx,~,hFig] = uiextras.plot_selectable(repmat(TimeZ,1,num_tracks), dZ, uiextras.cell_sprintf('Track %d',1:num_tracks));
xlabel(hAx,'Time [sec]');
ylabel(hAx,'\DeltaZ, Relative to reference [µm]');
title(hAx,'\DeltaZ Trace');
set(hFig,'NumberTitle','off',...
         'Name',[name,':dZ Trace']);
     
%% Plot Force and Length Data

% L vs MagH
[~,hAx,~,hFig] = uiextras.plot_timeordered(repmat(MagH,1,numel(MeasuredTracks)),...
                            Lext(:,MeasuredTracks),...
                            [],[],...
                            2*LextStdErr(:,MeasuredTracks),...
                            2*LextStdErr(:,MeasuredTracks),...
                            uiextras.cell_sprintf('Track %d',MeasuredTracks));
xlabel(hAx,'Magnet Position [mm]');
ylabel(hAx,'L_{ext} [µm]');
title(hAx,'L_{ext} vs. MagH');
set(hFig,'NumberTitle','off',...
         'Name',[name,':L-MagH']);

% Fx vs MagH
[~,hAx,~,hFig] = uiextras.plot_timeordered(repmat(MagH,1,numel(MeasuredTracks)),...
                            Fx(:,MeasuredTracks)/1e-12,...
                            [],[],...
                            (Fx(:,MeasuredTracks)-FxLower95(:,MeasuredTracks))/1e-12,...
                            (FxUpper95(:,MeasuredTracks)-Fx(:,MeasuredTracks))/1e-12,...
                            uiextras.cell_sprintf('Track %d',MeasuredTracks));
xlabel(hAx,'Magnet Position [mm]');
ylabel(hAx,'Fx [pN]');
title(hAx,'Fx vs. MagH');
set(hFig,'NumberTitle','off',...
         'Name',[name,':Fx-MagH']);

% F vs L
[~,hAx,~,hFig] = uiextras.ForceExtension_timeordered(...
                    Lext(:,MeasuredTracks),...
                    Fx(:,MeasuredTracks)/1e-12,...
                    2*LextStdErr(:,MeasuredTracks),...
                    (Fx(:,MeasuredTracks)-FxLower95(:,MeasuredTracks))/1e-12,...
                    (FxUpper95(:,MeasuredTracks)-Fx(:,MeasuredTracks))/1e-12,...
                    uiextras.cell_sprintf('Track %d',MeasuredTracks));
xlabel(hAx,'Extension [µm]');
ylabel(hAx,'Fx [pN]');
title(hAx,'Fx vs. L');
set(hFig,'NumberTitle','off',...
         'Name',[name,':Fx-L']);
     
% sqrt(Xvar) vs L
[~,hAx,~,hFig] = uiextras.plot_timeordered(...
                            Lext(:,MeasuredTracks),...
                            sqrt(Xvar(:,MeasuredTracks)),...
                            [],[],...
                            [],[],...
                            uiextras.cell_sprintf('Track %d',MeasuredTracks));
xlabel(hAx,'Extension [µm]');
ylabel(hAx,'sqrt( var(X) ) [px]');
title(hAx,'sqrt(Xvar) vs L');
set(hFig,'NumberTitle','off',...
         'Name',[name,':Xvar']);
                    

