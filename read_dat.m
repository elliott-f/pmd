function [time,param,moormeta,instmeta] = read_dat(FILENAME)

% [time,param,moormeta,instmeta] = read_dat(FILENAME) 
%
% Reads data from Malcolm Greig's oceanographic data format ascii file.
% i.e. 20 headerlines of metdata then 1 line per sample of
% homi damo dat dat dat ...
% where:
% homi = integer hour/minute
% damo = integer day/month
% and dat = space seperate ascii data
%
% Args:
% FILENAME = full path\filename of file
%
% Output:
% time = matlab datenum time vector (one time for each sample)
% data = structure containing variables 
% param = structure containing variable metadata and data
%    .name (variable name string)
%    .unit (variable unit string)
%    .fmat (variable output format string, e.g. '%5.2f')
%    .data (numeric data of the variable)
% various variables - output under the name used in the headerline
% 
% Optional outputs:
% moormeta = standard format mooring metadata (see ??)
% instmeta = standard format instrument metadata (see ??)
%
% Craig Stewart
% 11/11/2007
% Modified 22 Feb 2009
% 2010 Nov 2 : added param to output - for re-writing dat files in qcmd

% Get file name
if nargin==0
    [filename,pathname] = uigetfile('*.dat*','Select a DAT file to read');
    FILENAME = [pathname filename];
    if pathname==0, disp('No file chosen, reading cancelled'), return, end
else
    [pathstr, name, ext] = fileparts(FILENAME); 
    pathname = [pathstr '\']; 
    filename = [name ext];
end

%============================================================
% Read the header
 [moormeta,instmeta] = read_dat_header(FILENAME);

%============================================================
% Open the file
fid = fopen(FILENAME,'rt');
if fid == -1
    error(['Error> cannot open: ' FILENAME])
end
%disp(['Reading: ' FILENAME ])

try
% Skip the header
n = 18;
for ind = 1:n
    tline = fgetl(fid);
end

% Read the variable names
%% Changing to lower causes grief later on - removed BG Jan 2013
% varline = lower(fgetl(fid));
varline = fgetl(fid);
varline = strrep(varline,' dir ',' dirt '); % replaces dir with dirt as a 
% variable name as dir causes confusion with function dir
varnames = strread(varline,'%s');
varnames = varnames(2:end)'; % 1 is time name
numvar = length(varnames);

% Read the units
unitline = lower(fgetl(fid));
if ~isfield(instmeta,'datefmt')
    disp('Warning ''datefmt'' field not found in instmeta')
    keyboard
end
if instmeta.datefmt == 1 % Old style (homi damo) date format
    firstdatacol = 11;
elseif instmeta.datefmt == 2 % New Style (yyyy-mm-dd HH:MM:SS) date format
    firstdatacol = 21;
else
    disp(['trouble reading file: ' FILENAME]);
    disp(['unrecognised time units in unit line: ' unitline]);
    keyboard
    error('could not read file')
end
unitline = unitline(firstdatacol:end); 
unitnames = strread(unitline,'%s')';

% Read one line of data and get format
dataline = fgetl(fid);
fmat = getfmt(dataline(firstdatacol:end));

% Close file for reading
fclose(fid);

catch ME % close file then throw error
    fclose(fid);
    rethrow(ME)
end


% Read the data
if instmeta.datefmt == 1 % Old style (homi damo) date format
    varstr = ['homi damo ' varline(9:end)];
    fmtstr = ['%u%u' repmat('%f',1,numvar)];
    %eval(['[' eval('varstr') '] = textread(FILENAME,''' fmtstr ''',''headerlines'',20);'])
    eval(['[' varstr '] = textread(FILENAME,''' fmtstr ''',''headerlines'',20);'])

    % Extract time info
    ho = floor(homi/100);
    mi = round(100*mod(homi/100,1));
    da = floor(damo/100);
    mo = round(100*mod(damo/100,1));
    % Create year vector from header info
    t1 = datevec(instmeta.start_time);
    ye = repmat(t1(1),size(ho));
    nyind = find(diff(mo)<0);
    if ~isempty(nyind)
        for yearno = 1:length(nyind);
            ye(nyind(yearno)+1:end) = ye(nyind(yearno)+1:end) +1;
        end
    end
    time = datenum(ye,mo,da,ho,mi,0);

elseif instmeta.datefmt == 2 % New Style (yyyy-mm-dd HH:MM:SS) date format
    time_fmt = 'yyyy-mm-dd HH:MM:SS';
    varstr = ['timestr ' varline(11:end)];
    fmtstr = ['%19c' repmat('%f',1,numvar)];
    try
        eval(['[' varstr '] = textread(FILENAME,''' fmtstr ''',''headerlines'',20);'])
    catch
        disp(['Trouble reading file: ' filename]) 
        disp(['using: [' varstr '] = textread(FILENAME,''' fmtstr ''',''headerlines'',20);'])
        keyboard
    end
    time = datenum(timestr,time_fmt);
end

% Put all data in structure to export
for ind = 1:numvar;
    thisvarname = char(varnames(ind));
    param(ind).name = thisvarname;
    param(ind).unit = char(unitnames(ind));
    param(ind).fmat = fmat{ind}; % could try to detect format???
    param(ind).data = eval(char(varnames(ind)));
end


%============================================================
% Quality control - check for consistency
% Check start and end times
% if time(1) ~= moormeta.deploy_start_time;
%     disp('Warning: Start time is not that indicated in header')
%     disp(['Header: ' datestr(moormeta.deploy_start_time)])
%     disp(['Data:   ' datestr(time(1))])
%     uiwait(msgbox('Warning: Start time is not that indicated in header','modal'))
%     %eval(['!notepad ' FILENAME])
% end
% if time(end) ~= moormeta.deploy_stop_time;
%     disp(['Header: ' datestr(moormeta.deploy_stop_time)])
%     disp(['Data:   ' datestr(time(end))])
%     disp('Warning: End time is not that indicated in header')
%     uiwait(msgbox('Warning: End time is not that indicated in header','modal'))
%     %eval(['!notepad ' FILENAME])
%     %error(' ')
% end

% check time is always increasing
if any(diff(time)<0)
    disp(['Warning: Time not monotonically increasing in: ' FILENAME])
    uiwait(msgbox(['Warning: Time not monotonically increasing in: ' FILENAME],'modal'))
end

if 0 % save data to mat file
    newfilename = [FILENAME(1:end-3) 'mat'];
    disp(['Saving data to ' newfilename])
    eval(['save ' newfilename ' instmeta moormeta time data ' varline(9:end)])
end