function import_csv_file
%-----------------------------------------------------------------------------------------
% IMPORT_CSV_FILE Import current profile data from a csv AWAC file as column vectors.
%
% Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% NOTE: THIS SCRIPT ASSUMES 25 CELLS IN USE
%
% Example of .csv file format:
%  1        DateTime                         
%  2        Battery voltage                  (V)
%  3        Heading                          (degrees)
%  4        Pitch                            (degrees)
%  5        Roll                             (degrees)
%  6        Pressure                         (dbar)
%  7        Temperature                      (degrees C)
%  8        Analog input 1
%  9        Analog input 2
%  10       Speed cell 1                     (m/s)
%  11       Direction cell 1                 (degrees)
%  10-11    Repeated for all cells in profile.
%
% See also TEXTSCAN.
%
% History:
% Auto-generated by MATLAB on 2017/05/24 09:17:16
% 24-May-2017 FE Edited for use with process_AWAC.m
%
% NIWA Moorings
%-----------------------------------------------------------------------------------------

%% Initialize variables.
delimiter = ';';
startRow = 2;
endRow = inf;
[filename, pathname] = uigetfile('*.csv','Select a .csv file to import');
FILENAME = fullfile(pathname, filename);


%% Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: double (%f)
%   column45: double (%f)
%	column46: double (%f)
%   column47: double (%f)
%	column48: double (%f)
%   column49: double (%f)
%	column50: double (%f)
%   column51: double (%f)
%	column52: double (%f)
%   column53: double (%f)
%	column54: double (%f)
%   column55: double (%f)
%	column56: double (%f)
%   column57: double (%f)
%	column58: double (%f)
%   column59: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
DateTime = dataArray{:, 1};
Battery = dataArray{:, 2};
Heading = dataArray{:, 3};
Pitch = dataArray{:, 4};
Roll = dataArray{:, 5};
Pressure = dataArray{:, 6};
Temperature = dataArray{:, 7};
AnalogIn1 = dataArray{:, 8};
AnalogIn2 = dataArray{:, 9};
Speed1 = dataArray{:, 10};
Dir1 = dataArray{:, 11};
Speed2 = dataArray{:, 12};
Dir2 = dataArray{:, 13};
Speed3 = dataArray{:, 14};
Dir3 = dataArray{:, 15};
Speed4 = dataArray{:, 16};
Dir4 = dataArray{:, 17};
Speed5 = dataArray{:, 18};
Dir5 = dataArray{:, 19};
Speed6 = dataArray{:, 20};
Dir6 = dataArray{:, 21};
Speed7 = dataArray{:, 22};
Dir7 = dataArray{:, 23};
Speed8 = dataArray{:, 24};
Dir8 = dataArray{:, 25};
Speed9 = dataArray{:, 26};
Dir9 = dataArray{:, 27};
Speed10 = dataArray{:, 28};
Dir10 = dataArray{:, 29};
Speed11 = dataArray{:, 30};
Dir11 = dataArray{:, 31};
Speed12 = dataArray{:, 32};
Dir12 = dataArray{:, 33};
Speed13 = dataArray{:, 34};
Dir13 = dataArray{:, 35};
Speed14 = dataArray{:, 36};
Dir14 = dataArray{:, 37};
Speed15 = dataArray{:, 38};
Dir15 = dataArray{:, 39};
Speed16 = dataArray{:, 40};
Dir16 = dataArray{:, 41};
Speed17 = dataArray{:, 42};
Dir17 = dataArray{:, 43};
Speed18 = dataArray{:, 44};
Dir18 = dataArray{:, 45};
Speed19 = dataArray{:, 46};
Dir19 = dataArray{:, 47};
Speed20 = dataArray{:, 48};
Dir20 = dataArray{:, 49};
Speed21 = dataArray{:, 50};
Dir21 = dataArray{:, 51};
Speed22 = dataArray{:, 52};
Dir22 = dataArray{:, 53};
Speed23 = dataArray{:, 54};
Dir23 = dataArray{:, 55};
Speed24 = dataArray{:, 56};
Dir24 = dataArray{:, 57};
Speed25 = dataArray{:, 58};
Dir25 = dataArray{:, 59};

save(fullfile(pathname,[filename(1:end-4) '_sen.mat']),'month','day','year','hour',...
    'minute','second','errorCode','statusCode','batteryVoltage','soundSpeed','heading',...
    'pitch','roll','pressure','temperature','analog1','analog2')

save(fullfile(pathname,[filename(1:end-4) '_sen.mat']),'month','day','year','hour',...
    'minute','second','errorCode','statusCode','batteryVoltage','soundSpeed','heading',...
    'pitch','roll','pressure','temperature','analog1','analog2')

disp(['Saved ' filename ' data to matlab data file ' filename(1:end-4) '_sen.mat'])