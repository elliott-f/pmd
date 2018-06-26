function plot_aquatech(FILENAME)

%-----------------------------------------------------------------------------------------
% [] = plot_aquatech(FILENAME)
%
% Part of the NIWA PMD toolbox
%
% Plots Aquatech turbidity aqualogger turbidity/temp/depth on 3 subplots
%
% HISTORY
% 01-Nov-2014 FE created
% 16-Dec-2016 FE fixed bugs
%
% NIWA moorings
% Nov 2014
%-----------------------------------------------------------------------------------------

%% --- Choose file ---
if nargin == 0
    [filename,pathname] = uigetfile('tur*.mat','Choose the aquatech .mat file to plot','multiselect','on');
    if pathname==0, disp('No file chosen, reading cancelled'), return, end
else
    [pathstr, name, ext] = fileparts(FILENAME);
    % Mod BG May 2013 - gives full working directory if run from pwd
    if isempty(pathstr)
        pathstr = pwd;
    end
    pathname = [pathstr '\'];
    filename = {[name ext]};
end
if ischar(filename)
    filename = {filename}; % convert filename to cell for consistency with multifile case
end
% Sort alphabetically
filenamelist = fliplr(sort(filename));

% create colour defaults for plots depending on number of units
col = linspecer(numel(filenamelist));

%% --- PLOT DATA ---
% Loop through plotting all files
for fileno = 1:length(filenamelist)
    filename = char(filenamelist(fileno));
    FILENAME = [pathname filename];
    
    % Load the data
    load(FILENAME)
    
    % Set up figure and plot
    mcatfig = findobj('tag','turbFig');
    
    if fileno == 1 % Setup the figure (no figure open)
        n = 4;
        if exist('SSC','var') 
            n = n+1;
        end
        
        [fh,ah] = subplot_setup(n);
        set(fh,'tag','turbFig')
    end
    
    ind = 0;
    % Depth
    if exist('dep','var')
        ind = ind + 1;
        axes(ah(ind))
        hold on
        plot(time,dep,'col',col(fileno,:),'Marker','.');
        ylabel('Depth (m)')
        axis auto
        set(gca,'ydir','reverse')
    elseif exist('pre','var')
        ind = ind + 1;
        axes(ah(ind))
        hold on
        plot(time,pre,'col',col(fileno,:),'Marker','.');
        ylabel('Pre (dbar)')
        axis auto
        set(gca,'ydir','reverse')
    end
    
    % Turbidity
    if exist('tur','var')
        ind = ind + 1;
        axes(ah(ind))
        hold on
        plot(time,tur,'col',col(fileno,:),'Marker','.');
        ylabel('Tur (FTU)')
        axis auto
    end
    
    % Turbidity raw
    if exist('tur_raw','var')
        ind = ind + 1;
        axes(ah(ind))
        hold on
        plot(time,tur_raw,'col',col(fileno,:),'Marker','.');
        ylabel('Tur raw (counts)')
        axis auto
    end
    
    % Temp
    if exist('tem','var')
        ind = ind + 1;
        axes(ah(ind))
        hold on
        plot(time,tem,'col',col(fileno,:),'Marker','.');
        ylabel('Temp ( ^oC)')
        axis auto
    end
    
        % SCC
    if exist('SSC','var')
        ind = ind + 1;
        axes(ah(ind))
        hold on
        plot(time,SSC,'col',col(fileno,:),'Marker','.');
        ylabel('SSC ( mgl^-^l)')
        axis auto
    end
end

% Label axes and add legend
axes(ah(1))
% Add title
filenames = [];
for ind = 1:length(filenamelist);
    thisfile = filenamelist{ind};
    thisfile = thisfile(1:end-4);
    filenames = [filenames thisfile ', '];
end
filenames = filenames(1:end-2);
if length(filename)>1
    titletxt = (['Aquatech logger, file: ' filenames]);
else
    titletxt = (['Aquatech logger, files: ' filenames ]);
end
title(titletxt,'interpreter','none','fontweight','normal','fontsize',10)


% Add timestamps
pmd_timestamps(ah,time);

% Tidy!
pmd_style(fh)
set(ah,'Xlim',[time(1) time(end)])

% Add legend to second axis if we're plotting the last
axes(ah(2))
if fileno > 1
    th = findobj('parent',gca); % child list for temp axis
    legendh = [];
    legendtxt = {};
    for ind = 1:length(th);
        thistag = get(th(ind),'tag');
        if ~isempty(thistag) && ~strcmp(thistag,'monthgl') && ~strcmp(thistag,'monthname') && ~strcmp(thistag,'daynum') % ignore month gridlines
            legendh(end+1) = th(ind);
            legendtxt(end+1) = {thistag};
        end
    end
    legend(legendh,legendtxt,'interpreter','none')
end

% set figure tag for print_plots
% Mod May 2013 BG with quotes in case file name has a space in it!
if length(filenamelist)==1
    set(fh,'tag',[filename(1:end-4) '-uncropped']); % include inst #
else
    set(fh,'tag',['''' filename(1:end-5) '''']); % exclude inst # for multiple inst
end

print_fig

% Crop dates 
x = ginput;
set(ah,'xlim',[x(1,1) x(2,1)])
set(fh,'tag',filename(1:end-4))
print_fig


