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
% 2018-May-15 FE included ver 5.1 format
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

%% --- Set up -------------------------------------------------------------
% Prepopulate variables
timeformat = 'dd/mm/yyyy HH:MM:SS PM';
tline       = '';
lineno      = 0;
meta        = {};
time        = [];
param       = [];
timeAvg     = [];

% Define atm pressure
% Aqualogger pressure is absolute, substract the atmospheric pressure to get the gauge pressure
% http://cliflo.niwa.co.nz/ is a good place to get the mean monthly msl pressure in hPa
% 1 hectopascal (hPa) is equal to 0.01 decibar (dbar)
atm_pre = 1013.3;

disp('INSTRUMENT SETUP:')

% Set instrument metadata
disp('Instrument: Aquatech aqualogger');
imd.model = 'aquatech';

% Load file
fid = fopen(FILENAME);

% Verify this is an aquatech data file and what version we have
tline1 = fgetl(fid);
tline2 = fgetl(fid);
if ~(contains(tline1,'VERSION,') || contains(tline2,'VERSION,'))
    close(findall(0,'Tag','TMWWaitbar')); % Close waitbar
    msgbox('This doesn''t look like an aquatech datafile, cancelling processing')
    uiwait
    return
else
    % It probably is an aquatech file, find out what version
    if contains(tline1,'VERSION,')
        fields = textscan(tline1,'%s%f%s%s','delimiter',',');
        [version,versionNo,softwareNo,osNo] = deal(fields{:});
    elseif ~contains(tline1,'VERSION,') && contains(tline2,'VERSION,')
        fields = textscan(tline2,'%s%f','delimiter',',');
        [version,versionNo] = deal(fields{:});
    end
    if numel(versionNo) > 1
        versionNo = versionNo(1);
    end
    % Check firmware version
    disp(['Firmware version: ' num2str(versionNo)]);
end

% Version 4.0 and 5.1 have different header, parse accordingly
frewind(fid);
% Read Header, header end with line 'BURSTSTART'
while ~(contains(tline,'BurstStart') || contains(tline,'BURSTSTART')  || contains(tline,'Burst Start') || contains(tline,'Burst_Start'))
    
    % Get next header line
    tline = fgetl(fid);
    lineno = lineno + 1;
    
    % Form firmware version less than 5.0
    % Get logger type
    if contains(tline,'LOGGER TYPE')
        loc = strfind(tline,',');
        disp(['Logger type: ' tline(loc + 1:end)]);
    end
    
    % Get serial number
    if contains(tline,'LOGGER,') && ~contains(tline,'type')
        loc = strfind(tline,',');
        lac = strfind(tline,'-');
        if versionNo < 5
            disp(['Serial number: ' tline(loc(1)+1:loc(2)-1)]);
            imd.snstr = tline(lac(1)+1:end);
        else
            disp(['Serial number: ' tline(loc(1)+1:end)]);
            imd.snstr = tline(lac(1)+1:end); % proabably wrong
        end
    end
    
    % Save regime
    if contains(tline,'REGIME')
        loc = strfind(tline,',');
        if versionNo < 5
            if strcmp(tline(1:loc(1)-1),'REGIME')
                disp([' ']);
                disp(['SAMPLING SETUP:']);
                disp(['Mode: ' tline(loc(1)+1:loc(2)-1)]);
                disp(['Interval: ' tline(loc(2)+1:loc(3)-1)]);
                disp(['Samples per burst: ' tline(loc(3)+1:loc(4)-1)]);
                disp(['Sample rate: ' tline(loc(4)+1:loc(5)-1)]);
                disp(['Samples per average: '  tline(loc(5)+1:loc(6)-1)]);
                disp([' ']);
                disp(['SENSORS:']);
            end
        else
            if isempty(loc)
                disp([' ']);
                disp(['SAMPLING SETUP:']);
                tline = fgetl(fid);
                lineno = lineno + 1;
                fields = textscan(tline,'%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
                tline = fgetl(fid);
                lineno = lineno + 1;
                fields2 = textscan(tline,'%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
                disp([char(fields{2}) ': ' char(fields2{2})]);
                disp([char(fields{3}) ': ' char(fields2{3})]);
                disp([char(fields{4}) ': ' char(fields2{4})]);
                disp([char(fields{5}) ': ' char(fields2{5})]);
                disp([char(fields{6}) ': ' char(fields2{6})]);
                disp([char(fields{7}) ': ' char(fields2{7})]);
                disp([char(fields{8}) ': ' char(fields2{8})]);
                disp([char(fields{9}) ': ' char(fields2{9})]);
                disp([' ']);
                disp(['SENSORS:']);
            end
        end
    end
    
    % Channel names
    if contains(tline,'[LOGGER CHANNEL')
        loc = strfind(tline,',');
        if numel(loc) < 2
            x = strtrim(tline(loc(1)+1:end));
        else
            x = strtrim(tline(loc(1)+1:loc(2)-1));
        end
        disp([tline(2:loc(1)-2) ': ' x]);
    end
    
    % Variable names
    if contains(tline,'HEADING')
        % Get heading names
        numFields = numel(strfind(tline,','));
        pres = repmat('%s',1,numFields);
        headings = textscan(tline,pres,'delimiter',',');
        
        % Get units names
        tline = fgetl(fid);
        lineno = lineno + 1;
        numFields = numel(strfind(tline,','));
        pres = repmat('%s',1,numFields);
        units = textscan(tline,pres,'delimiter',',');
        
        % Join them to get variable names
        varnames = cell(numFields,1);
        for i = 1:numel(headings)
            if isempty(char(headings{i}))
                headings{i} = headings{i-1};
            end
            myHeading = strtrim(char(headings{i}));
            myUnit = strtrim(char(units{i}));
            varnames{i} = [myHeading '_' myUnit];
            if contains(varnames{i},'HEADING_UNITS')
                varnames{i} = 'Burst_Start'; 
            elseif contains(varnames{i},'°')
                varnames{i} = strrep(varnames{i},'°','deg');
                varnames{i} = strrep(varnames{i},' ','_');
            elseif contains(varnames{i},'[dd/MM/yyyy HH:mm:ss]')
                varnames{i} = 'time';
            elseif contains(varnames{i},'[ms]')
                varnames{i} = strrep(varnames{i},'[ms]','ms');
            elseif contains(varnames{i},' ')
                varnames{i} = strrep(varnames{i},' ','_');
            elseif contains(varnames{i},'/')
                varnames{i} = strrep(varnames{i},'/','');
            end
        end
    end
end

%% --- Read in Data--------------------------------------------------------
startRow = lineno;

while ~feof(fid)
    fgetl(fid);
    lineno = lineno + 1;
end
endRow = lineno;
frewind(fid);

% Get format specifications
formatSpec = '';
stringVars = {'Turbidity_Status','Burst_Start','time'};
for i = 1:numel(varnames)
    if ismember(varnames{i},stringVars)
        formatSpec = [formatSpec '%s'];
    else
        formatSpec = [formatSpec '%f'];
    end
end
formatSpec = [formatSpec '%[^\n\r]'];  

% Get data
disp(' ');
disp('PROCESSING INFORMATION: ');
disp(['Using an atmospheric pressure of ' num2str(atm_pre)]);
disp(['Reading file: ' FILENAME '...']);

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fid, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fid, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', ',', 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fid);
    textscan(fid, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fid, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fid);

%% --- Clean Data ---------------------------------------------------------
dataArray(end) = [];
for i = 1:numel(varnames)
    eval([varnames{i} ' = dataArray{:,i};']);
end

% Check status, status = 'SAT', this usually only effects turbidity
% dummy = sum(~ismember(status,'OK'));
% if dummy
%     disp([num2str(dummy) ' data points flagged by instrument, setting turbidity at these points to NaN'];
%     badinds = find(~ismember(status,'OK'));
%     tur_raw(badinds) = NaN;
%     tur(badinds) = NaN;
% end

% How many gains were used?
if exist('gains','var')
    num_gains = unique(gain);
    gains = '';
    if numel(num_gains) == 1
        disp(['Only one gain used: ' num2str(num_gains) 'x']);
    else
        m = size(num_gains,1);
        for i = 1:m;
            gains = [gains num2str(num_gains(i)) 'x,'];
        end
        disp([ num2str(numel(num_gains)) ' gains used: ' gains(1:end-1)]);
    end
end

% Create datenum array
time = strrep(time,'.','');
if length(time{1})==15
    % No seconds
    timeformat = 'dd/mm/yyyy HH:MM';
end
time = datenum(time,timeformat);

% Averaging bursts
% Get the index of all burst starts by using the var 'burststart'
% if versionNo == 4
burststart = [find(ismember(Burst_Start,{'BURSTSTART','Burst Start'}))];
if isempty(burststart)
    % Theres weird differences in file structure in versionNo == 5.1
    Burst_Start = cellfun(@(x) x(2:end),Burst_Start,'UniformOutput', false);
    burststart = find(ismember(Burst_Start,'Burst Start'));
end

disp([num2str(numel(burststart)) ' bursts found']);

% Get times at burststart
timeAvg = time(burststart);

% Now loop through and calculate averages
disp(('Averaging bursts...'));
for i = 2:numel(varnames)
    if ~strcmp(varnames{i},'Turbidity_Status')
        disp(['Calculating the average over burst of ' varnames{i}]);
        disp(varnames{i})
        % Create arrays
        this_var = eval(varnames{i});
        this_var_avg = zeros(size(timeAvg));
        
        % Loop to do averages
        for ind = 1:numel(this_var_avg)
            if ind ~= numel(this_var_avg)
                this_var_avg(ind) = nanmean(this_var(burststart(ind):burststart(ind+1)-1));
            else
                this_var_avg(ind) = nanmean(this_var(burststart(ind):end));
            end
        end
        
        % Save to actual variables
        eval([varnames{i} '_avg = this_var_avg;']);
    end
end

% Pressure
Pressure_avg = (Pressure_bar_avg*10) - (atm_pre* 0.01);

% save to standard var name
if exist('Temperature_degC_avg','var')
    tem_avg = Temperature_degC_avg;
elseif exist('Ext_temperature_degC_avg','var')
    tem_avg = Ext_temperature_degC_avg;
end
if exist('Turbidity_Raw1_avg','var')
    tur_raw_avg = Turbidity_Raw1_avg;
elseif exist('Turbidity_Raw_avg','var')
    tur_raw_avg = Turbidity_Raw_avg;
end
if exist('Battery_V_avg','var')
    bat_avg = Battery_V_avg;
elseif exist('Battery_V_avg','var')
    bat_avg = Battery_V_avg;
end

timeAll = time;
time = timeAvg;

% Convert pressure to depth
%dep_avg = sw_dpth(pre_avg,latitude);

%% -- Put data into format required by write_dat and define units etc -----
ind = 1;
if exist('tem_avg','var')
    ind = ind + 1;
    param(1).name = 'tem';
    param(1).longname = 'Temperature';
    param(1).unit = 'DegC';
    param(1).data = tem_avg;
    param(1).fmat = '%5.2f';
end

if exist('Pressure_avg','var')
    ind = ind + 1;
    param(2).name = 'pre';
    param(2).longname = 'Pressure';
    param(2).unit = 'dbar';
    param(2).data = Pressure_avg;
    param(2).fmat = '%5.2f';
end

if exist('tur_raw_avg','var')
    ind = ind + 1;
    param(3).name = 'tur_raw';
    param(3).longname = 'Turbidity raw';
    param(3).unit = 'counts';
    param(3).data = tur_raw_avg;
    param(3).fmat = '%5.2f';
end

if exist('Turbidity_FTU_avg','var')
    ind = ind + 1;
    param(4).name = 'tur';
    param(4).longname = 'Turbidity';
    param(4).unit = 'FTU';
    param(4).data = Turbidity_FTU_avg;
    param(4).fmat = '%5.2f';
end

if exist('bat_avg','var')
    ind = ind + 1;
    param(5).name = 'bat';
    param(5).longname = 'Battery Voltage';
    param(5).unit = 'V';
    param(5).data = bat_avg;
    param(5).fmat = '%4.2f';
end

if exist('SSC_mgl_avg','var')
    ind = ind + 1;
    param(5).name = 'SSC';
    param(5).longname = 'SSC';
    param(5).unit = 'mg/l';
    param(5).data = SSC_mgl_avg;
    param(5).fmat = '%4.2f';
end
% param(6).name = 'gain';
% param(6).longname = 'Gain Setting';
% param(6).unit = 'Gain no.';
% param(6).data = gain_avg;
% param(6).fmat = '%4.2f';

