function [hErr,hAx,hPnl,hFig] = fitplot_selectable(X,Y,varargin)


import uiextras.*
%% parse inputs
p = inputParser;
p.CaseSensitive = false;
%p.KeepUnmatched = true;

addParameter(p,'Xlower',[]);
addParameter(p,'Xupper',[]);
addParameter(p,'Ylower',[]);
addParameter(p,'Yupper',[]);

addParameter(p,'FigureHandle',[]);
addParameter(p,'DataNames',[]);

addParameter(p,'FitType',[]);
addParameter(p,'FitParameters',{}); %cell array of name,value pairs to pass to fit
addParameter(p,'FitLHS',[],@(x) isempty(x)||isa(x,'function_handle'));
addParameter(p,'FitLHSinv',[],@(x) isempty(x)||isa(x,'function_handle'));

addParameter(p,'FitStringGenerator',@defaultStrGen,@(x) isa(x,'function_handle'));

parse(p,varargin{:});


%Save unmatched in a simple cell array which we will pass to set the axis
%properties
%UnMatchedParams = {reshape(fieldnames(p.Unmatched),1,[]);reshape(struct2cell(p.Unmatched),1,[])};
p = p.Results;

%% Make plot
[hErr,hAx,hPnl,hFig] = errorbar_selectable(...
    X,...
    Y,...
    p.Xlower,p.Xupper,...
    p.Ylower,p.Yupper,p.DataNames,p.FigureHandle);

%set appdata for fitting system
setappdata(hAx,'FitType',p.FitType);
setappdata(hAx,'FitParameters',p.FitParameters);
setappdata(hAx,'FitLHS',p.FitLHS);
setappdata(hAx,'FitLHSinv',p.FitLHSinv);
setappdata(hAx,'FitStringGenerator',p.FitStringGenerator);

%% Setup Fitting Menu
if iscell(X)
    num_tracks = numel(X);
else
    num_tracks = size(X,2);
end

if numel(p.DataNames)<num_tracks
    for n=numel(p.DataNames)+1:num_tracks
        p.DataNames{n} = sprintf('Data %d',n);
    end
end

hFitMenu = uimenu(hFig,'Label','Fit Data');
hShowFit = uimenu(hFitMenu,'Label','Show Fit');
hExclude = uimenu(hFitMenu,'Label','Exclude Data', 'Separator','on');
hInclude = uimenu(hFitMenu,'Label','Include Data');
hReset = uimenu(hFitMenu,'Label','Reset Data');

hFitLines = gobjects(num_tracks,1);
for n=1:num_tracks
    hFitLines(n) = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hErr(n).Color,'Visible','off');
end
for t=1:num_tracks
    % show fit
    uimenu(hShowFit,'Label',p.DataNames{t},...
        'Checked','off',...
        'Interruptible','off',...
        'Callback', @(h,e) ShowFit(h,e,hAx,hErr(t),hFitLines(t)));
    
    % Exclude Data
    uimenu(hExclude,'Label',p.DataNames{t},...
        'Interruptible','off',...
        'Callback', @(h,e) ExcludeData(h,e,hAx,hErr(t),hFitLines(t)));
    
    % Include Data
    uimenu(hInclude,'Label',p.DataNames{t},...
        'Interruptible','off',...
        'Callback', @(h,e) IncludeData(h,e,hAx,hErr(t),hFitLines(t)));
    
    % reset data
    uimenu(hReset,'Label',p.DataNames{t},...
        'Interruptible','off',...
        'Callback', @(h,e) ResetData(h,e,hAx,hErr(t),hFitLines(t)));
    
end

function ShowFit(hMenu, ~, hAx, hEb, hFitLine)
import uiextras.*
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
    hFitLine = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb.Color);
end
hFitLine.Visible = 'on';

UpdateFit(hAx,hEb,hFitLine);

function UpdateFit(hAx,hEb, hFitLine)
import uiextras.*
if ~ishghandle(hFitLine)
    return;
end

if ~isappdata(hFitLine,'ExcludeIdx')
    ExcludeIdx = [];
    setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);
else
    ExcludeIdx = getappdata(hFitLine,'ExcludeIdx');
end

if ~isappdata(hAx,'FitType')
    return;
end

%calc fit
X = hEb.XData;
Y = hEb.YData;
X(ExcludeIdx) = [];
Y(ExcludeIdx) = [];

FT = getappdata(hAx,'FitType');
FO = getappdata(hAx,'FitParameters');
FLHS = getappdata(hAx,'FitLHS');
FLHSinv = getappdata(hAx,'FitLHSinv');
StrGen = getappdata(hAx,'FitStringGenerator');

if isempty(FT)
    return;
end

if isa(FLHS,'function_handle')
    Y = FLHS(Y);
end

fitobj = fit(X,Y,FT,FO{:});

%update plot
xl = hAx.XLim;
yl = hAx.YLim;
x = linspace(min(min(X),xl(1)),max(max(X),xl(2)),100);
y = fitobj(x);

if isa(FLHSinv,'function_handle')
    y = FLHSinv(y);
end

set(hFitLine,'XData',x,'YData',y)
hAx.YLim = yl;
hAx.XLim = xl;

%update legend
hLeg = findobj(hAx.Parent,'Type','legend');
if ~isempty(hLeg)
    hFitLine.DisplayName = StrGen(fitobj);
    if ~any(hLeg.PlotChildren==hFitLine)
        hLeg.PlotChildren = [reshape(hLeg.PlotChildren,1,[]),reshape(hFitLine,1,[])];
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

function S = defaultStrGen(fitobj)
%default fit string generator
% output:
%   v1:###[###,###]\n
%   v2:###[###,###]\n
%    ...
names = coeffnames(fitobj);
vals = coeffvalues(fitobj);
CI = confint(fitobj);
S = [];
for n=1:numel(names)
    S = [S,sprintf('%s:%g[%g,%g]\n',names{n},vals(n),CI(:,n))];
end

function ExcludeData(~, ~, hAx, hEb, hFitLine)
import uiextras.*
if ~ishghandle(hFitLine)
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
import uiextras.*
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
import uiextras.*
if ~ishghandle(hFitLine)
    'not hghandle'
    return;
    %hFitLine = line(hAx,NaN,NaN,'Marker','none','LineStyle','--','Color',hEb.Color);
end
ExcludeIdx = [];
setappdata(hFitLine,'ExcludeIdx',ExcludeIdx);

UpdateFit(hAx,hEb,hFitLine);
