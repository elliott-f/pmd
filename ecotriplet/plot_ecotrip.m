function plot_ecotrip(FILENAME)
%-----------------------------------------------------------------------------------------
% [] = plot_ecotrip(FILENAME)
%
% Part of the NIWA PMD toolbox
%
% Plots Wetlabs ECO-Triplet data
%
% HISTORY
% 01-Jul-2008   Craig Stewart
% 22-Jul-2015   Fiona Elliott   formatting
% 06-Apr-2017   Fiona Elliott   Adding multiselect option
%
% NIWA moorings
%-----------------------------------------------------------------------------------------


%% --- CHOOSE FILE ---
if nargin == 0
    [filename,pathname] = uigetfile('*.mat','Choose the ECO .mat file to plot','multiselect','on');
    if pathname==0, disp('No file chosen, reading cancelled'), return, end
else
    [pathname, name, ext] = fileparts(FILENAME);
    if isempty(pathname)
        pathname = pwd;
    end
    filename = {[name ext]};
end
if ischar(filename)
    filename = {filename}; % convert filename to cell for consistency with multifile case
end
% Sort alphabetically
filenamelist = sort(filename);

% Create colour map for plots based on number of instruments
col = linspecer(numel(filenamelist)* 2);


%% --- PLOT DATA ---
% Loop through plotting all files
for fileno = 1:length(filenamelist)
    filename = char(filenamelist(fileno));
    FILENAME = fullfile(pathname,filename);
    ncol = fileno*2-1;
    
    % Load the data
    eval(['load ''' FILENAME ''''])
    
    % Set up figure and plot
    ecofig = findobj('tag','ECO-Triplet'); 
    
    if fileno == 1 % Setup the figure (no figure open)
        [fh,ah] = subplot_setup(4);
        set(fh,'tag','ECO-Triplet','name','ECO-Triplet','NumberTitle','off')
    end
    
    % Chlorophyll
    axes(ah(1))
    lh(fileno) = plot(time,chl,'col',col(ncol,:),'tag',filename(1:end-4));
    ylabel('Chl (mg/m^3)')
    hold on
    
    % CDOM
    axes(ah(2))
    plot(time,cdom,'col',col(ncol,:))
    ylabel('CDOM (ppb)')
    hold on
    
    % Total Volume scattering
    axes(ah(3))
    plot(time,beta650,'col',col(ncol,:))
    ylabel('TV scatter (/m/sr)')
    set(gca,'tag','no3ax')
    hold on
    
    % Backscatter
    axes(ah(4))
%     Hb_bp(fileno) = plot(time,bp650,'col',col(ncol,:),'tag',[filename(1:end-4) ' particulate']);
%     hold on
    Hb_b(fileno) = plot(time,bt650,'col',col(ncol,:),'tag',[filename(1:end-4) ' total']);
%     legend([Hb_bp Hb_b],{'Particulate','Total'})
    ylabel('Total Backscatter (/m)')
    hold on
end


%% --- TIDY PLOTS ---
% Title
if numel(filenamelist)>1
    titletxt = (['Outer FoT Mooring (' upper(filenamelist{1}(2:8)) ') ECO-Triplet Data ']);
else
    titletxt = (['ECO-Triplet Data, file: ' upper(filenamelist{1})]);
end
title(ah(1),titletxt,'interpreter','none')

% Legend, add if more than 1 inst
if numel(filenamelist) > 1
    %th = findobj('parent',gca); % child list for temp axis
    legendh = [];
    legendtxt = {};
    for ind = 1:length(lh);
        legendh(end+1) = lh(ind);
        legendtxt{end+1} = get(lh(ind),'tag');
    end
    legend(ah(1),legendh,legendtxt,'interpreter','none');
end

% % Backscatter legend
% legendh = [];
% legendtxt = {};
% for ind = 1:numel(Hb_bp);
%     legendh(end+1)      = Hb_bp(ind);
%     legendh(end+1)      = Hb_b(ind);
%     if numel(Hb_bp) == 1
%         legendtxt    = {'particulate' 'total'};
%     else
%         legendtxt{end+1}    = get(Hb_bp(ind),'tag');
%         legendtxt{end+1}    = get(Hb_b(ind),'tag');
%     end
% end
% legend(ah(4),legendh,legendtxt,'interpreter','none');


% Add timestamps
pmd_timestamps(ah,time);

% Tidy!
pmd_style(fh)

% Set figure tag for print_plots
if numel(filenamelist)==1
    set(fh,'tag',filenamelist{1}(1:end-4),'name',filenamelist{1}(1:end-4));
else
    set(fh,'tag',[filenamelist{1}(1:8) '_comparison'],'name',[filenamelist{1}(1:8) '_comparison']);
end

% Crop dates 
x = ginput;
set(ah,'xlim',[x(1,1) x(2,1)])
set(fh,'tag',filename(1:end-4))
print_fig