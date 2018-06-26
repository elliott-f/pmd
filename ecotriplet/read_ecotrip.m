function [imd,time,param] = read_ecotrip(FILENAME)
%-----------------------------------------------------------------------------------------
% [param,time,meta] = read_ecotriplet(FILENAME)
%
% Part of the PMD toolbox
%
% Reads raw ECO-triplet datafiles, files are usually corrupted so script reads through every 
% line. Notes bad lines. Coverts to eng units
%
% NOTE: - NEEDS A LOT MORE WORKS WITH FETCHING CAL COEFF'S
%         Include wavelength info can change 650/660 between units - is in dev file
%         Read cal coeffs directly from dev file?
%         create table/db with all cal data?
%
% HISTORY:
% 2011-Apr-15 CS
% 2016-Dec-08 FE Rewrote orginal script
%
% NIWA moorings
%-----------------------------------------------------------------------------------------

%% --- Set up ---

% Settings
applyFlag = 0;

% Choose file
if nargin == 0
    [filename,pathname] = uigetfile({'*.*'},'Choose a raw ECO-Triplet ascii file to load');
    if pathname==0, disp('No file chosen, reading cancelled'), return, end
    FILENAME = [pathname filename];
end

[pathstr, filename, ext] = fileparts(FILENAME);

% Set instrument metadata
imd.model   = 'ecotrip';
% loc         = strfind(filename,'_');
ecoCalSn    = {'476','477','1316','1227','1670'};
sn          = [];

disp('Instrument:       ECO-Triplet')
 
% Try get serial number if using niwa file naming convention
% if ~isempty(loc) && numel(loc)==1
%     sn = filename(loc+1:end);
% end

% If no seriral number is supplied, or its not numeric user selects
if isempty(sn) || ~isnumeric(str2num(sn))
 [s,v] = listdlg('PromptString','Select the serial number of this ECO-Triplet',...
        'SelectionMode','single',...
        'ListString',ecoCalSn);
    if v
        sn = ecoCalSn{s};
    end
end

imd.snstr = sn;
disp(['Serial number:   ' sn])

% Prepopulate variables
startRow    = 2;
endRow      = inf;
delimiter   = '\t';

switch sn
    case '476'
        % Define calibration coefficients (28-May-2014 SN 476)
        % 650	650
        beta650_scale  = 3.599e-06;
        beta650_dark   = 61;
        chl_scale      = 0.0112;
        chl_dark       = 55;
        cdom_scale     = 0.0861;
        cdom_dark      = 56;
    case '477'
        % Define calibration coefficients (29-Jan-08 SN 477)
        % 650	650
        beta650_scale  = 3.511e-06;
        beta650_dark   = 63;
        chl_scale      = 0.0122;
        chl_dark       = 48;
        cdom_scale     = 0.0964;
        cdom_dark      = 42;
    case '480'
        % Define calibration coefficients (26-Oct-15 SN 480)
        % 650	650
        beta650_scale  = 4.095e-06;
        beta650_dark   = 59;
        chl_scale      = 0.0113;
        chl_dark       = 54;
        cdom_scale     = 0.1122;
        cdom_dark      = 52;
    case '1316'
        % Define calibration coefficients (22-Jun-2015 SN1316)
        % 650	650
        beta650_scale  = 3.567e-6;
        beta650_dark   = 47;
        chl_scale      = 0.0122;
        chl_dark       = 58;
        cdom_scale     = 0.0899;
        cdom_dark      = 44;
    case '1227'
        % Define calibration coefficients (ECO BBFL2WB-1227 20 Oct 2015 SN1227) 
        % 650	650
        beta650_scale  = 3.698e-06;
        beta650_dark   = 47;
        chl_scale      = 0.0122;
        chl_dark       = 51;
        cdom_scale     = 0.0938;
        cdom_dark      = 47;
     case '1670'
        % Define calibration coefficients (ECO BBFL2WB-1227 20 Oct 2015 SN1227) 
        % 650	650
        beta650_scale  = 1.922e-06; %actually refers to wavelength 700nm
        beta650_dark   = 48;
        chl_scale      = 0.0072;
        chl_dark       = 44;
        cdom_scale     = 0.0908;
        cdom_dark      = 48;
end  


%% --- Load file ---
% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(FILENAME,'r');

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[3,4,5,6,7,8,9]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% Remove last row of data - usually 'etx' or corrupted
raw(end,:) = [];

% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [3,4,5,6,7,8,9]);
rawCellColumns = raw(:, [1,2,10,11,12,13,14,15,16,17]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

disp(['Data lines: ' num2str(length(R))])

% Flags
flag = zeros(length(R),9);

% Allocate imported array to column variable names
% Time
dateString = rawCellColumns(:, 1);
thisFlag = cellfun(@(x) length(x) ,dateString);
flag(:,1) = thisFlag ~= 8;

timeString = rawCellColumns(:, 2);
thisFlag = cellfun(@(x) length(x) ,timeString);
flag(:,2) = thisFlag ~= 8;

% Scattering
rawBeta650Ref = cell2mat(rawNumericColumns(:, 1));
flag(:,3) = rawBeta650Ref ~= 650;

rawBeta650 = cell2mat(rawNumericColumns(:, 2));
flag(:,4) = rawBeta650 > 4119 | isnan(rawBeta650);
disp(['rawBeta650 flagged data: ' num2str(sum(flag(:,4))) ' of ' num2str(numel(rawBeta650)) ' (' num2str(sum(flag(:,4))/numel(rawBeta650)*100,'%3.1f') '%)'])

% Chl-A
rawChlRef = cell2mat(rawNumericColumns(:, 3)); 
flag(:,5) = rawChlRef ~= 695;

rawChl = cell2mat(rawNumericColumns(:, 4)); 
flag(:,6) = rawChl > 4119 | isnan(rawChl);
disp(['rawChl flagged data: ' num2str(sum(flag(:,6))) ' of ' num2str(numel(rawChl)) ' (' num2str(sum(flag(:,6))/numel(rawChl)*100,'%3.1f') '%)'])

% CDOM
rawCdomRef = cell2mat(rawNumericColumns(:, 5)); 
flag(:,7) = rawCdomRef ~= 460;

rawCdom = cell2mat(rawNumericColumns(:, 6)); 
flag(:,8) = rawCdom > 4119 | isnan(rawCdom);
disp(['rawCdom flagged data: ' num2str(sum(flag(:,8))) ' of ' num2str(numel(rawCdom)) ' (' num2str(sum(flag(:,8))/numel(rawCdom)*100,'%3.1f') '%)'])

% Ref
ref = cell2mat(rawNumericColumns(:, 7)); 
flag(:,9) = ref > 600 | isnan(ref);

% Organise data
data = [rawBeta650Ref,rawBeta650,rawChlRef,rawChl,rawCdomRef,rawCdom,ref];
masterFlag = sum(flag,2) ~= 0;
if applyFlag
    disp(['Removing flagged data: ' num2str(sum(masterFlag)) ' of ' num2str(numel(masterFlag)) ' (' num2str(sum(masterFlag)/numel(masterFlag)*100,'%3.1f') '%)'])
    data(masterFlag,:) = [];
    imd.notes = {['Flagged ' num2str(sum(masterFlag)) ' of ' num2str(length(R)) ' lines (' num2str(sprintf('%2.1f',(sum(masterFlag)/length(R))*100)) '%) as corrupted']};
else
    disp(['NOT removing flagged data: ' num2str(sum(masterFlag)) ' of ' num2str(numel(masterFlag)) ' lines (' num2str(sum(masterFlag)/numel(masterFlag)*100,'%3.1f') '%)'])
end

% disp(['Flagged data lines: ' num2str(sum(masterFlag))]);
% disp(['Percentage of flagged data: ' num2str(sprintf('%2.1f',(sum(masterFlag)/length(R))*100)) '%']);

% Time
% dateString(masterFlag)  = {'01/01/16'};
% timeString(masterFlag)  = {'12:00:00'};
TIME                    = [char(dateString) repmat(' ',length(dateString),1) char(timeString)];
time                    = datenum(TIME,'mm/dd/yy HH:MM:SS');
if applyFlag
    time(masterFlag)    = [] ;
end

raw_beta650             = data(:,2);
raw_chl                 = data(:,4);
raw_cdom                = data(:,6);

% Convert to Eng Units by applying cal coefficients beta650_scale
beta650 = beta650_scale .* (raw_beta650 - beta650_dark);    % Beta 650 (/m/sr) (scattering)
chl     = chl_scale .* (raw_chl - chl_dark);                % Chlorophyll (?g/L = mg/m3)
cdom    = cdom_scale .* (raw_cdom - cdom_dark);             % CDOM concentration (ppb)

% Corrections
% Ideally should correct Scattering (labmda) for attenuation within the path. 
% See EcoTrip manual section 5.1 - scattering data corrections
% However for apg650 but this is only a 4%  error at apg660 = 1/m, so small 
% for most measures. Marine typically << 0.1 - 0.4 %.

% Calculate other derived parameters

%Ecotrip manual 5.2.1 Volume Scattering of Particles
%The corrected volume scattering of particles, ?(117°,lambda) values represent total
%volume scattering, i.e., scattering from particles and molecular scattering from water.
%To obtain the volume scattering of particles only, subtract the volume
%scattering of water,?w(117°,lambda)

% Volume scattering of Particles - MARINE (Salinity estimated to be 36 PSU) %
%beta650water = Derive[,1,]-(1.38*(660/500)^(-4.32)*(1+0.3*36/37)*(10^-4)*(1+(0.5+0.5*cos(2*117))*(1-0.09)/(1+0.09))) % Betap660 = Beta660 - Betaw660(/m/sr)
theta_deg   = 117;
theta       = (pi/180)*theta_deg;
lambda      = 650;
delta       = 0.09;
S           = 35; % estimated ocean salinity (PSU)
beta650w    = 1.38*(lambda/500)^(-4.32) * (1 + 0.3*S/37)*10^(-4) * (1+(cos(theta))^2*(1-delta)/(1+delta)); % Water volume scattering
beta650p    = beta650 - beta650w; % Particle volume scattering

bt_sw       = 0.0029308*(lambda/500)^(-4.24); % Total scattering of sea water (35-39 PSU)
bb_sw       = bt_sw/2; % Backscattering of seawater (half of total scattering)

% Calculate Particulate backscatter coefficient from Particulate volume scattering
X           = 1.1; % scale factor from Boss and Pegau 2001 (ecotrip manual 5.2.2)
b_bp        = 2*pi*X*beta650p; % Particulate backscatter coefficient (unit = m^-1)
b_b         = b_bp + bb_sw; % Total backscatter coefficient (units = m^-1)


%% Save variables to structure
%(1) Total volume scattering
param(1).name = 'beta650';
param(1).longname = 'Total Volume Scattering (117deg,650nm)';
param(1).unit = '/m/sr';
param(1).data = beta650;
param(1).fmat = '%7.6f';

%(2) Chlorophyll concentration
param(2).name = 'chl';
param(2).longname = 'Chlorophyll concentration';
param(2).unit = 'mg/m^3';
param(2).data = chl;
param(2).fmat = '%5.2f';

%(3) CDOM
param(3).name = 'cdom';
param(3).longname = 'Coloured Dissolved Organic Material concentration';
param(3).unit = 'ppb';
param(3).data = cdom;
param(3).fmat = '%5.1f';

%(4) Particulate volume scattering
param(4).name = 'bp650';
param(4).longname = 'Particulate backscatter coefficient (650nm)';
param(4).unit = '/m';
param(4).data = b_bp;
param(4).fmat = '%5.3f';

%(5) Total volume scattering
param(5).name = 'bt650';
param(5).longname = 'Total backscatter coefficient (650nm)';
param(5).unit = '/m';
param(5).data = b_b;
param(5).fmat = '%5.3f';


