function PlotForceExtension(varargin)
%Plot Force Extension results from mtdat file
% Input:
%   PlotForceExtension(header,ExperimentData)
%       plot data in specified data structures header,ExpData are produced
%       by LoadExperimentData()
%   PlotForceExtension('path to your mtdat file')
%       load data from file and plot
%   PlotForceExtension()
%       prompt user to select file with uigetfile()
persistent last_dir;
if nargin>1
    header = varargin{1};
    ExperimentData = varargin{2};
    filename = '';
else
    if nargin<1
        %select file
        [File,Dir] = uigetfile({'*_ForceExtension_*.txt;*.mtdat','MT Experiment Data';...
                                '*_ForceExtension_*.txt','Legacy MT Data';...
                                '*.mtdat','MT-Dat File'},...
                                'Select Force Extension File',fullfile(last_dir,'*.*'));
        if File==0
            return
        end
        if ~isempty(Dir)
            last_dir = Dir;
        end
    end
    if nargin==1
        filepath = varargin{1};
    else
        filepath=fullfile(Dir,File);
    end
    [~,filename,ext] = fileparts(filepath);
    if strcmpi(ext,'.mtdat')
        [header,ExperimentData] = LoadExperimentData(filepath);
    else
        filepath = ConvertMTdata(filepath);
        [header,ExperimentData] = LoadExperimentData(filepath);
        [~,filename,~] = fileparts(filepath);
    end
    
end
if isempty(header)
    return;
end

if ~isfield(ExperimentData,'StepData') || isempty(ExperimentData(1).StepData)
    warning('Data could not be loaded from file. Experiment might have been canceled or data could be corrupted.');
    return;
end

%% Calculate various time-averaged parameters
if ~isfield(header,'Temperature')
    header.Temperature = 25;
end
kBT=1.380648813e-23*(273.15+header.Temperature)*10^6;
hWait = waitbar(0,'Calculating Step-averaged Metrics');

cMX = ~isfield(ExperimentData,'meanX');
cMY = ~isfield(ExperimentData,'meanY');
cVX = ~isfield(ExperimentData,'varX');
cVY = ~isfield(ExperimentData,'varY');
cVZ = ~isfield(ExperimentData,'varZ') || ~isfield(ExperimentData,'varZLower95') || ~isfield(ExperimentData,'varZUpper95');
cL = ~isfield(ExperimentData,'meanL');
cSEL = cL || ~isfield(ExperimentData,'seL');
cFX = cL || ~isfield(ExperimentData,'Fx');
cFxerr = cL || cFX || cSEL || ~isfield(ExperimentData,'FxLower95') || ~isfield(ExperimentData,'FxUpper95');
cAV = ~isfield(ExperimentData,'AllanVar') || ~isfield(ExperimentData.AllanVar,'X')|| ~isfield(ExperimentData.AllanVar,'M')|| ~isfield(ExperimentData.AllanVar,'Navg');


for n=1:numel(ExperimentData)
    ExperimentData(n).dT = 24*3600*nanmean(diff([ExperimentData(n).StepData.Date],1,2),2);
    %XY statistics
    if cMX
        ExperimentData(n).meanX = nanmean([ExperimentData(n).StepData.X],2);
    end
    if cMY
        ExperimentData(n).meanY = nanmean([ExperimentData(n).StepData.Y],2);
    end

    
    X = bsxfun(@minus,[ExperimentData(n).StepData.X],ExperimentData(n).meanX)*header.PxScale;
    Y = bsxfun(@minus,[ExperimentData(n).StepData.Y],ExperimentData(n).meanY)*header.PxScale;
    nx = sum(~isnan(X),2);
    
    if cVX || cVY
        ExperimentData(n).varX = nanvar([ExperimentData(n).StepData.X],0,2);
        ExperimentData(n).varY = nanvar([ExperimentData(n).StepData.Y],0,2);
    end
    
    %Z statiztics
    if cVZ
        thisZ = [ExperimentData(n).StepData.Z_ABS];
        indNaN = all(isnan(thisZ),2);
        Z_REL = [ExperimentData(n).StepData.Z_REL];
        thisZ(indNaN,:) = Z_REL(indNaN,:);
        nz = sum(~isnan(thisZ),2);
        ExperimentData(n).varZ = nanvar(thisZ,0,2);
        ExperimentData(n).varZLower95 = (nz-1).*ExperimentData(n).varZ./chi2inv(0.05/2,(nz-1));
        ExperimentData(n).varZUpper95 = (nz-1).*ExperimentData(n).varZ./chi2inv(1-0.05/2,(nz-1));
    end
    
    %L statistics
    if cL || cSEL
        for s=1:numel(ExperimentData(n).StepData)
            ExperimentData(n).StepData(s).L =sqrt(X(:,s).^2+Y(:,s).^2+ExperimentData(n).StepData(s).dZ.^2);
        end

        L = [ExperimentData(n).StepData.L];
        ExperimentData(n).meanL = nanmean(L,2);
        ExperimentData(n).seL = nanstd(L,0,2)./sqrt(size(L,2));
    end
    
    %Fx using variance
    if cFX || cFxerr
        ExperimentData(n).Fx = kBT*ExperimentData(n).meanL./(ExperimentData(n).varX * header.PxScale.^2);
        ExperimentData(n).FxLower95 = ...
                kBT*...
                    (ExperimentData(n).meanL-2*ExperimentData(n).seL)./ ((nx-1).*(ExperimentData(n).varX * header.PxScale.^2)./chi2inv(0.05/2,(nx-1)));
        ExperimentData(n).FxUpper95 = ...
                kBT*...
                    (ExperimentData(n).meanL+2*ExperimentData(n).seL)./ ((nx-1).*(ExperimentData(n).varX * header.PxScale.^2)./chi2inv(1-0.05/2,(nx-1)));
    end
    
    % Allan Variance
    if cAV
        [avX,M,Navg] = allanvariance(1000*header.PxScale*[ExperimentData(n).StepData.X]');
        ExperimentData(n).AllanVar.X = avX';
        ExperimentData(n).AllanVar.M = M';
        ExperimentData(n).AllanVar.Navg = Navg';
    end
    
    waitbar(n/numel(ExperimentData),hWait);
end
try
    delete(hWait);
catch
end
% %% Prompt to save data?
% if any([cMX,cMY,cVX,cVY,cVZ,cL,cSEL,cFX,cFxerr,cAV])
%     btn = questdlg('Save Calculations?','yes','no','yes');
%     if strcmpi(btn,'yes')
%         [D,F] = uiputfile(filepath,'Save Calculations?');
%         if F~=0
%             
%         end
%     end
% end

%% Putvars
uiextras.putvar(header,ExperimentData);

%% Calc other params to make plotting easier
num_tracks = numel(header.TrackingInfo);
MeasTrk = find(strcmpi({header.TrackingInfo.Type},'Measurement'));
RefTrk = find(strcmpi({header.TrackingInfo.Type},'Reference'));
AllNames = cell(num_tracks,1);
for n=1:num_tracks
    AllNames{n} = sprintf('Trk %d, %s',n,header.TrackingInfo(n).Type);
end

MeasNames = cell(numel(MeasTrk),1);
for n=1:numel(MeasTrk)
    MeasNames{n} = sprintf('Trk %d',MeasTrk(n));
end

%% Plots
mL = [ExperimentData.meanL]';
mL(:,RefTrk) = [];
seL = [ExperimentData.seL]';
seL(:,RefTrk) = [];

Fx = [ExperimentData.Fx]';
Fx(:,RefTrk) = [];
Fl95 = [ExperimentData.FxLower95]';
Fl95(:,RefTrk) = [];
Fu95 = [ExperimentData.FxUpper95]';
Fu95(:,RefTrk) = [];

% F vs L
[~,hAx,~,hFig] = uiextras.ForceExtension_timeordered(...
                    mL,...
                    Fx/1e-12,...
                    2*seL,...
                    (Fx-Fl95)/1e-12,...
                    (Fu95-Fx)/1e-12,...
                    uiextras.cell_sprintf('Track %d',MeasTrk));
xlabel(hAx,'Extension [µm]');
ylabel(hAx,'Fx [pN]');
title(hAx,'Fx vs. L');
set(hFig,'NumberTitle','off',...
         'Name',[filename,':Fx-L']);
 
% Allan Variance
if ~isfield(header,'ParticleRadius') || isnan(header.ParticleRadius)
    ploop = true;
    while ploop
        answer = inputdlg({'Particle Radius [nm]'},'Particle Radius',1,{'1400'});
        if isempty(answer)
            return
        end
        val = str2double(answer{1});
        if ~isnan(val)
            header.ParticleRadius = val;
            ploop=false;
        end
    end
end


a0 = 1.4e-5/500*header.ParticleRadius;
aLow = .5*a0;
aHigh = 1.5*a0;
k = NaN(numel(ExperimentData),num_tracks);
a = NaN(numel(ExperimentData),1);
for trk=reshape(MeasTrk,1,[])
    
    figure('NumberTitle','off',...
        'Name',[filename,sprintf('AV:Track %d',trk)]);
    gca;
    hold on;
    colors = lines(numel(ExperimentData));
    hLines = gobjects(numel(ExperimentData),1);
    for n=1:numel(ExperimentData)
       [k(n,trk),a(n,trk)] = FitAllanVariance_HO(ExperimentData(n).AllanVar.X(trk,2:end),...
           ExperimentData(n).dT*ExperimentData(n).AllanVar.M(trk,2:end),...
           ExperimentData(n).AllanVar.Navg(trk,2:end),...
           'KappaStart',2.5e-5,...
           'AlphaStart',a0,...
           'AlphaLower',aLow,...
           'AlphaUpper',aHigh);
    
        plot(ExperimentData(n).dT*ExperimentData(n).AllanVar.M(trk,:),ExperimentData(n).AllanVar.X(trk,:),'.','color',colors(n,:),'markersize',10);
        
        x = linspace(ExperimentData(n).dT*ExperimentData(n).AllanVar.M(trk,1),ExperimentData(n).dT*ExperimentData(n).AllanVar.M(trk,end),100);
        y = allanvar_HO(x,a(n,trk),k(n,trk));
        
        hLines(n) = plot(x,y,'--','color',colors(n,:));
    end
    set(gca,'yscale','log','xscale','log');
    axis tight;
    legend(hLines,uiextras.cell_sprintf('\\kappa=%g, \\alpha=%g',[k(:,trk),a(:,trk)]));
    xlabel('Duration \tau [s]');
    ylabel('Allan Variance [nm^2]');
end

%% Force Extension
[~,hAx,~,hFig] = uiextras.ForceExtension_timeordered(...
                    mL,...
                    k(:,MeasTrk).*mL*1000,...
                    2*seL,...
                    [],...
                    [],...
                    uiextras.cell_sprintf('Track %d',MeasTrk));
xlabel(hAx,'Extension [µm]');
ylabel(hAx,'F [pN]');
title(hAx,'Fav vs. L');
set(hFig,'NumberTitle','off',...
         'Name',[filename,':Fav-L']);

%%
x = linspace(0,2,100);
plot(hAx,x,Fwlc(2,50,x),'-r');
