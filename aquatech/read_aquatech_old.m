function [imd,time,param] = read_aquatech(FILENAME,latitude)

%-----------------------------------------------------------------------------------------
% [imd,time,param] = read_aquatech(FILENAME)
%
% Reads Aquatech aqualogger data from .csv file
% 
% Returns:  imd (instrument metadata - structures)
%           param (parameter data - structure)
%           time (time vector - double)
%
% NOTE: Only tested with files using single gain of 100x, may need to
% change script for autogaining
%
% History:
% 01-Oct-2014   FE Created for WIAL project
% 11-Nov-2016   FE Fixed bugs, less hardcoding
%
% NIWA moorings
% Oct 2014
%-----------------------------------------------------------------------------------------
% Defaults
if nargin < 2
    latitude = -41.301568; % Greta Point
    disp(['No latitude supplied, using Greta Point position, setting latitude to ' num2str(latitude)])
end

if nargin == 0
    [filename,pathname] = uigetfile('*.csv');
    FILENAME = [pathname filename];
end

% Define ATM pressure
% Aqualogger pressure is absolute, substract the atmospheric pressure to get the gauge pressure
% http://cliflo.niwa.co.nz/ is a good place to get the mean monthly msl pressure in hPa
% 1 hectopascal (hPa) is equal to 0.01 decibar (dbar)
atm_pre = 1013.3;
disp(['Using an atmospheric pressure of ' num2str(atm_pre) ])
imd.notes = {['Using an atmospheric pressure of ' num2str(atm_pre)]};

% Load file
fid = fopen(FILENAME);

% Settings
timeformat = 'dd/mm/yyyy HH:MM:SS PM';
tline = '';
lineno = 0;
imd.model = 'aquatech';

% Verify this is an aquatech data file
tline1 = fgetl(fid);
tline2 = fgetl(fid);
if ~(strcmp(tline1(1:5),'FILE,') && strcmp(tline2(1:8),'VERSION,'))
    msgbox('This doesn''t look like an aquatech datafile, cancelling processing')
    uiwait
    return
end

frewind(fid);
% Read Header, header end with line 'BURSTSTART'
while isempty(strfind(tline,'BURSTSTART'))
    
    % Get next header line
    tline = fgetl(fid);
    lineno = lineno + 1;
    
    % Check software version
    if ~isempty(strfind(tline,'VERSION'))
        loc = strfind(tline,',');
        Version = tline(loc(1)+1:end);
        disp(['VERSION: ' Version])
%         % Give warning about file format
%         if ~strcmp(Version,'4.0')
%             warning(['read_aquatech.m was written for software version 4.0, this file is ' Version ' data file format may be different'])
%         end
    end
    
        % Get logger type
    if ~isempty(strfind(tline,'LOGGER TYPE'))
        loc = strfind(tline,',');
        disp(['LOGGER TYPE: ' tline(loc + 1:end)])
    end
    
    % Get serial number
    if ~isempty(strfind(tline,'LOGGER,'))
        loc = strfind(tline,',');
        imd.snstr = tline(loc(1)+1:loc(2)-1);
        imd.snstr = imd.snstr(strfind(imd.snstr,'-')+1:end);
        disp(['SERIAL NUMBER: ' tline(loc(1)+1:loc(2)-1)]);
    end
    
    % Save regime
    if ~isempty(strfind(tline,'REGIME'))
        loc = strfind(tline,',');
        if strcmp(tline(1:loc(1)-1),'REGIME')
            imd.notes{end+1,1} = ['Mode:' tline(loc(1)+1:loc(2)-1)];
            imd.notes{end+1,1} = ['Interval:' tline(loc(2)+1:loc(3)-1)];
            imd.notes{end+1,1} = ['Samples Per Burst:' tline(loc(3)+1:loc(4)-1)];
            imd.notes{end+1,1} = ['Sample Rate:' tline(loc(4)+1:loc(5)-1)];
            imd.notes{end+1,1} = ['Samples Per Average:'  tline(loc(5)+1:loc(6)-1)];
        end
    end
end

fclose(fid)

vars = 'burststart,time,tem_raw,tem,pres_raw,pre,tur_raw,tur,gain,status,bat_raw,bat'; % all variables in .csv file
my_vars = {'tem_raw','tem','pres_raw','pre','tur_raw','tur','gain','bat'}; % the variables we are interested in
formatstring = '%s%s%f%f%f%f%f%f%f%s%f%f';

% Get data
disp(['Reading file: ' FILENAME '...'])
eval(['[' vars '] = textread(FILENAME,''' formatstring ''',''delimiter'','','',''headerlines'',lineno);'])

% Check status, status = 'SAT', this usually only effects turbidity
% dummy = sum(~ismember(status,'OK'));
% if dummy
%     disp([num2str(dummy) ' data points flagged by instrument, setting turbidity at these points to NaN'])
%     badinds = find(~ismember(status,'OK'));
%     tur_raw(badinds) = NaN;
%     tur(badinds) = NaN;
% end

% How many gains were used?
num_gains = unique(gain);
gains = '';
if numel(num_gains) == 1
    disp(['Only one gain used: ' num2str(num_gains) 'x'])
else
    m = size(num_gains,1);
    for i = 1:m;
        gains = [gains num2str(num_gains(i)) 'x,'];
    end
    disp([ num2str(numel(num_gains)) ' gains used: ' gains(1:end-1)])
end

% Create datenum array
time = datenum(strrep(time,'.',''),timeformat);

% Averaging bursts
% Get the index of all burst starts by using the var 'burststart'
burststart = [1 ; find(ismember(burststart,'BURSTSTART'))];
disp([num2str(numel(burststart)) ' bursts found'])

% Get times at burststart
time = time(burststart);

% Now loop through and calculate averages .... not the most effiecient way
% of doing this.....
disp('Averaging bursts...')
for i = 1:numel(my_vars)
    
    disp(['Calculating the average over burst of ' my_vars{i}])
    
    % Create arrays
    this_var = eval(my_vars{i});
    this_var_avg = zeros(size(time));
    
    % Loop to do averages
    for ind = 1:numel(this_var_avg)
        if ind ~= numel(this_var_avg)
            this_var_avg(ind) = nanmean(this_var(burststart(ind):burststart(ind+1)-1));
        else
            this_var_avg(ind) = nanmean(this_var(burststart(ind):end));
        end
    end
    
    % Save to actual variables
    eval([my_vars{i} '_avg = this_var_avg;']);
    
end

% Pressure
pre_avg = (pre_avg*10) - (atm_pre* 0.01);

% Convert pressure to depth
dep_avg = sw_dpth(pre_avg,latitude);

%% Put data into format required by write_dat and define units etc
param(1).name = 'tem';
param(1).longname = 'Temperature';
param(1).unit = 'DegC';
param(1).data = tem_avg;
param(1).fmat = '%5.2f';

param(2).name = 'dep';
param(2).longname = 'Depth';
param(2).unit = 'm';
param(2).data = dep_avg;
param(2).fmat = '%5.2f';

param(3).name = 'pre';
param(3).longname = 'Pressure';
param(3).unit = 'dbar';
param(3).data = pre_avg;
param(3).fmat = '%5.2f';

param(4).name = 'tur_raw';
param(4).longname = 'Turbidity raw';
param(4).unit = 'counts';
param(4).data = tur_raw_avg;
param(4).fmat = '%5.2f';

param(5).name = 'tur';
param(5).longname = 'Turbidity';
param(5).unit = 'FTU';
param(5).data = tur_avg;
param(5).fmat = '%5.2f';

param(6).name = 'gain';
param(6).longname = 'Gain Setting';
param(6).unit = 'Gain no.';
param(6).data = gain_avg;
param(6).fmat = '%4.2f';

param(7).name = 'bat';
param(7).longname = 'Battery Voltage';
param(7).unit = 'V';
param(7).data = bat_avg;
param(7).fmat = '%4.2f';

