function PlotForceExtension(varargin)

if nargin>1
    header = varargin{1};
    ExperimentData = varargin{2};
    filename = '';
else
    if nargin<1
        [header,ExperimentData,filepath] = LoadExperimentData();
    else
        [header,ExperimentData,filepath] = LoadExperimentData(varargin{1});
    end
    if isempty(header)
        return;
    end
    [~,filename,~] = fileparts(filepath); 
end

if ~isfield(ExperimentData,'StepData') || isempty(ExperimentData(1).StepData)
    warning('Data could not be loaded from file. Experiment might have been canceled or data could be corrupted.');
    return;
end

%% Calculate Tether Length & Force
kBT=1.380648813e-23*(273.15+header.Temperature)*10^6;
for n=1:numel(ExperimentData)
    X = bsxfun(@minus,[ExperimentData(n).StepData.X],ExperimentData(n).meanX)*header.PxScale;
    Y = bsxfun(@minus,[ExperimentData(n).StepData.Y],ExperimentData(n).meanY)*header.PxScale;
    for s=1:numel(ExperimentData(n).StepData)
        ExperimentData(n).StepData(s).L =sqrt(X(:,s).^2+Y(:,s).^2+ExperimentData(n).StepData(s).dZ.^2);
    end
    ExperimentData(n).meanL = nanmean([ExperimentData(n).StepData.L],2);
    ExperimentData(n).stdL = nanstd([ExperimentData(n).StepData.L],0,2);
    ExperimentData(n).Fx = kBT*ExperimentData(n).meanL./(ExperimentData(n).stdX.^2 * header.PxScale.^2);
    ExperimentData(n).FxErr = kBT*ExperimentData(n).stdL./(ExperimentData(n).stdX.^2 * header.PxScale.^2);
end
% meanL = ExperimentData(n).meanL
% stdL = ExperimentData(n).stdL
% stdX = ExperimentData(n).stdX
% Fx = ExperimentData(n).Fx


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
%% Plot Data
%=========================
%% L vs Mag
mL = [ExperimentData.meanL]';
mL(:,RefTrk) = [];
sL = [ExperimentData.stdL]';
sL(:,RefTrk) = [];

[~,hAx,~,hFig] = plot_timeordered(...
        repmat([ExperimentData.MagnetHeight]',1,numel(MeasTrk)),...
        mL,...
        [],[],...
        sL,sL,MeasNames);
hAx.Title.String = 'Length vs Magnet Height';
xlabel(hAx,'Magnet Height [mm]');
ylabel(hAx,'Avg. Tether Length [µm]');
hFig.Name = [filename,' L v. MagH'];
if ~isempty(filename)
    hFig.NumberTitle = 'off';
end

%% Fx v Mag
Fx = [ExperimentData.Fx]'*10^12;
Fx(:,RefTrk) = [];
FxErr = [ExperimentData.FxErr]'*10^12;
FxErr(:,RefTrk) = [];

[~,hAx,~,hFig] = plot_timeordered(...
    repmat([ExperimentData.MagnetHeight]',1,numel(MeasTrk)),...
    Fx,...
    [],[],...
    FxErr,FxErr,MeasNames);
hAx.Title.String = 'Force vs Magnet Height';
xlabel(hAx,'Magnet Height [mm]');
ylabel(hAx,'Force ( k_BTL/<dx^2>) [pN]');
%set(hAx,'yscale','log');
hFig.Name = [filename,' F v. MagH'];
if ~isempty(filename)
    hFig.NumberTitle = 'off';
end

%% Fx v L
[hEb_FvL,hAx_FvL,~,hFig_FvL] = plot_timeordered(...
    mL,...
    Fx,...
    sL,sL,...
    FxErr,FxErr,MeasNames);
hAx_FvL.Title.String = 'Force vs Length';
xlabel(hAx_FvL,'Avg. Tether Length [µm]');
ylabel(hAx_FvL,'Force ( k_BTL/<dx^2>) [pN]');
set(hAx_FvL,'yscale','log');

hFig_FvL.Name = [filename,' F v. L'];
if ~isempty(filename)
    hFig_FvL.NumberTitle = 'off';
end

%% Export Data to Workspace
putvar(header,ExperimentData)


%% Setup Fitting Menu
num_tracks = numel(MeasTrk);

hFitMenu = uimenu(hFig_FvL,'Label','Fit Data');
hShowFit = uimenu(hFitMenu,'Label','Show Fit');
hExclude = uimenu(hFitMenu,'Label','Exclude Data', 'Separator','on');
hInclude = uimenu(hFitMenu,'Label','Include Data');
hReset = uimenu(hFitMenu,'Label','Reset Data');

hFitLines = gobjects(num_tracks,1);
for n=1:num_tracks
    hFitLines(n) = line(hAx_FvL,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb_FvL(n).Color,'Visible','off');
end
for t=1:num_tracks
    %lbl = sprintf('Track %d',t);
    % show fit
    uimenu(hShowFit,'Label',MeasNames{t},...
        'Checked','off',...
        'Callback', @(h,e) ShowFit(h,e,hAx_FvL,hEb_FvL(t),hFitLines(t)));
%     if strcmpi(header.TrackingInfo(t).Type,'Reference')
%         hM.Enable = 'off';
%     end
    
    % Exclude Data
    uimenu(hExclude,'Label',MeasNames{t},...
        'Callback', @(h,e) ExcludeData(h,e,hAx_FvL,hEb_FvL(t),hFitLines(t)));
%     if strcmpi(header.TrackingInfo(t).Type,'Reference')
%         hM.Enable = 'off';
%     end
    
    % Include Data
    uimenu(hInclude,'Label',MeasNames{t},...
        'Callback', @(h,e) IncludeData(h,e,hAx_FvL,hEb_FvL(t),hFitLines(t)));
%     if strcmpi(header.TrackingInfo(t).Type,'Reference')
%         hM.Enable = 'off';
%     end
    
    % reset data
    uimenu(hReset,'Label',MeasNames{t},...
        'Callback', @(h,e) ResetData(h,e,hAx_FvL,hEb_FvL(t),hFitLines(t)));
%     if strcmpi(header.TrackingInfo(t).Type,'Reference')
%         hM.Enable = 'off';
%     end
    
end

function ShowFit(hMenu, ~, hAx, hEb, hFitLine)
if strcmpi(hMenu.Checked,'on')
    hMenu.Checked = 'off';
    if ishghandle(hFitLine)
        %excluded data
        if isappdata(hFitLine,'hExcLine') && ishghandle(getappdata(hFitLine,'hExcLine'))
            delete(getappdata(hFitLine,'hExcLine'))
            rmappdata(hFitLine,'hExcLine');
        end
        %legend update
        hLeg = findobj(hAx.Parent,'Type','legend');
        if ~isempty(hLeg)
            id = find(hLeg.PlotChildren==hFitLine);
            if ~isempty(id)
                hLeg.PlotChildren(id) = [];
            end
        end
        hFitLine.Visible = 'off';
    end
    return;
end
hMenu.Checked = 'on';

if ~ishghandle(hFitLine)
    'make new'
    hFitLine = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb.Color);
end
hFitLine.Visible = 'on';

UpdateFit(hAx,hEb,hFitLine);


function UpdateFit(hAx,hEb, hFitLine)
if ~ishghandle(hFitLine)
    'not handle'
    return;
end

if ~isappdata(hFitLine,'ExcludeIdx')
    ExcludeIdx = [];
    setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);
else
    ExcludeIdx = getappdata(hFitLine,'ExcludeIdx');
end

%calc fit
L = hEb.XData;
F = hEb.YData;
L(ExcludeIdx) = [];
F(ExcludeIdx) = [];
[Lo,P,LoCI, PCI, ~] = FitWLC(L,F);

%update plot
xl = hAx.XLim;
yl = hAx.YLim;
x = linspace(0,Lo,30);
y = Fwlc(Lo,P,x);
set(hFitLine,'XData',x,'YData',y)
hAx.YLim = yl;
hAx.XLim = xl;

%update legend
hLeg = findobj(hAx.Parent,'Type','legend');
if ~isempty(hLeg)
    hFitLine.DisplayName = sprintf('L_0=%0.2f[%0.2f,%0.2f]µm L_p=%0.1f[%0.1f,%0.1f]nm',Lo,LoCI(1),LoCI(2),P,PCI(1),PCI(2));
    if ~any(hLeg.PlotChildren==hFitLine)
        hLeg.PlotChildren = [hLeg.PlotChildren,hFitLine];
    end
end

%plot exclusion points
if isappdata(hFitLine,'hExcLine') && ishghandle(getappdata(hFitLine,'hExcLine'))
    hExcLine = getappdata(hFitLine,'hExcLine');
    set(hExcLine,'XData',hEb.XData(ExcludeIdx),'YData',hEb.YData(ExcludeIdx));
else
    if isempty(ExcludeIdx)
        hExcLine = line(hAx,NaN,NaN,'Marker','x','color','r','MarkerSize',12,'LineStyle','none');
    else
        hExcLine = line(hAx,hEb.XData(ExcludeIdx),hEb.YData(ExcludeIdx),'Marker','x','color','r','MarkerSize',12,'LineStyle','none');
    end
    setappdata(hFitLine,'hExcLine',hExcLine);
end



function ExcludeData(~, ~, hAx, hEb, hFitLine)
if ~ishghandle(hFitLine)
    'not hghandle'
    return;
    %hFitLine = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb.Color);
end
if ~isappdata(hFitLine,'ExcludeIdx')
    ExcludeIdx = [];
    setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);
else
    ExcludeIdx = getappdata(hFitLine,'ExcludeIdx');
end

rect = getrect(hAx);
thisExclude = find( hEb.XData>rect(1) & hEb.XData<(rect(1)+rect(3)) & hEb.YData>rect(2) & hEb.YData<(rect(2)+rect(4)));
for n=numel(thisExclude):-1:1
    if any(thisExclude(n) == ExcludeIdx)
        thisExclude(n) = [];
    end
end
ExcludeIdx = [ExcludeIdx;reshape(thisExclude,[],1)];

setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);

UpdateFit(hAx,hEb,hFitLine);


function IncludeData(~, ~, hAx, hEb, hFitLine)
if ~ishghandle(hFitLine)
    'not hghandle'
    return;
    %hFitLine = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb.Color);
end
if ~isappdata(hFitLine,'ExcludeIdx')
    ExcludeIdx = [];
    setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);
else
    ExcludeIdx = getappdata(hFitLine,'ExcludeIdx');
end

rect = getrect(hAx);
thisInclude = find( hEb.XData>rect(1) & hEb.XData<(rect(1)+rect(3)) & hEb.YData>rect(2) & hEb.YData<(rect(2)+rect(4)));
for n=numel(ExcludeIdx):-1:1
    if any(ExcludeIdx(n) == thisInclude)
        ExcludeIdx(n) = [];
    end
end
setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);

UpdateFit(hAx,hEb,hFitLine);

function ResetData(~, ~, hAx, hEb, hFitLine)
if ~ishghandle(hFitLine)
    'not hghandle'
    return;
    %hFitLine = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb.Color);
end
ExcludeIdx = [];
setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);

UpdateFit(hAx,hEb,hFitLine);

function [Lo,P,LoCI, PCI, fo] = FitWLC(L,Fx)
%Fit Fx vs L to WLC model
% Fx: pN
% L: same units as Lo (um a good choice)
% Lo = contour length (units of L)
% P = persistence length in nm
% fo = fitobject

ft = fittype('log10(4.11/P*(1/4*(1-x/Lo)^(-2)-1/4+x/Lo))');
%ft = fittype( @(Lo,P,x) log10( 4.11./P.*(1/4*(1- lessthan1(x./Lo) ).^(-2)-1/4+lessthan1(x./Lo)) ));

fo = fit(L,log10(Fx),ft,...
        'StartPoint',[100,1],...
        'Lower',[0,0],...
        'Upper',[Inf,Inf]);
    
coef = coeffvalues(fo);
coefint = confint(fo);
Lo = coef(1);
P = coef(2);
LoCI = coefint(:,1);
PCI = coefint(:,2);


function F = Fwlc(Lo,P,x)
x(x>=Lo) = NaN;
x(x<=0) = NaN;
F = 4.11./P.*(1/4*(1-x./Lo).^(-2)-1/4+x./Lo);

function x = lessthan1(x)
x(x>1) = Inf;