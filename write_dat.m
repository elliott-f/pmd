function write_dat(filename,time,param,moormeta,instmeta,timefmt)

% write_dat(FILENAME,time,param,moormeta,instmeta)
%
% Writes data to Malcolm Greig's oceanographic data format ascii file.
% i.e. 20 headerlines of metdata then 1 line per sample of
% timestamp dat dat dat ...
% where:
% timestamp = data timestamp in format specified below (yyy-mm-dd HH:MM:SS)
% and dat = space seperate ascii data
%
% Args:
% FILENAME = full path\filename of output file
% time = matlab datenum time vector (one time for each sample)
% param = structure containing variable data and metadata
%    .name (variable name string)
%    .unit (variable unit string)
%    .fmat (variable output format string, e.g. '%5.2f')
%    .data (numeric data of the variable)
% 
% moormeta = standard format mooring metadata
% instmeta = standard format instrument metadata
%
% Optional args:
% timefmt: timeformat 
%          1 = homi damo
%          2 = yyyy-mm-dd HH:MM:SS
%          3 = yyyy-mm-dd HH:MM:SS.FFF 
%
% Craig Stewart
% 11/11/2007
% 18/May/08 - modified to new time format
% 5 August 2011 - added option time format in arg list

if nargin<5
    error('Not enough input arguments')
elseif nargin == 5
    timefmt = 2;
end

% Set time format
if timefmt==1
    timeformat = 'HHMM ddmm';
    hdrtimeformat = 'HHMM dd mm yy';
elseif timefmt == 2
    timeformat = 'yyyy-mm-dd HH:MM:SS';
    hdrtimeformat = timeformat;
elseif timefmt == 3
    timeformat = 'yyyy-mm-dd HH:MM:SS.FFF';
    hdrtimeformat = timeformat;    
else
    error('Unrecognised time format')
end

% Put all variable data into matrix
for ind = 1:length(param);
    dat = param(ind).data;
    DATA(:,ind) = dat(:);
end

fast = 0; % Condition for fastcat file FE 06/13

%=================================================
% extract path and filename
%[pathstr, name, ext] = fileparts(FILENAME); pathname = [pathstr '\']; filename = [name ext];

%=================================================
% Prepare metadata for headerlines
linetxt(1) = {filename};
if (exist('moormeta','var') || exist('instmeta','var'))
    try 
        linetxt(2) = {instmeta.instid}; 
        if strcmp(instmeta.instid,'fastcat 166') % Condition for fastcat file
            fast = 1;
        end
    end
    try 
        if ~any(isnan([moormeta.lat moormeta.lon]))
            [latstr,lonstr] = llnum2str(moormeta.lat,moormeta.lon);
            linetxt(3) = {[latstr ' ' lonstr]}; 
        end
    end
    try 
        linetxt(4) = {moormeta.location}; 
    end
    try 
        wd = num2str(moormeta.water_depth); 
    catch
        wd = '??'; 
    end
    try 
        mh = num2str(instmeta.height); 
    catch
        mh = '??'; 
    end
    if fast
        linetxt(5) = {[wd ' ~']};
    else
    linetxt(5) = {[wd ' ' mh]};
    end
    if ~isfield(instmeta,'utc_offset')
        instmeta.utc_offset = nan;
    end
    linetxt(6) = {[datestr(time(1),'yyyy-mm-dd HH:MM:SS') ', ' datestr(time(end),'yyyy-mm-dd HH:MM:SS') ', ' sprintf('%+2d',instmeta.utc_offset)]};
    if fast
        linetxt(7) = {'16Hz'};
    else
        intervalvec = datevec(round((time(2)-time(1))*(3600*24))./(3600*24));% interval (nearest second); 
        linetxt(7) = {int2str(intervalvec(3:end))};
    end  
    lineno = 8;
    if isfield(instmeta,'notes') 
        if iscell(instmeta.notes)
            inotes = instmeta.notes;
        else
            inotes = {instmeta.notes};
        end
        for ind = 1:length(inotes)
            linetxt(lineno) = inotes(ind);
            lineno = lineno + 1;
        end
    end
    if isfield(moormeta,'notes')
        if iscell(moormeta.notes)
            mnotes = moormeta.notes;
        else
            mnotes = {moormeta.notes};
        end
        for ind = 1:length(mnotes)
            linetxt(lineno) = mnotes(ind);
            lineno = lineno + 1;
        end
    end
else
    disp('Warning: mooring or instrument metadata not available - header data not writen to file.')
end


%=================================================
% Create variable name and units lines (header lines 19 and 20
varnametxt = ['time' repmat(' ',1,length(timeformat)-4)];
varunittxt = timeformat;
dataformat = [];
for varno = 1:length(param)
    % Check fieldwidth of each variable (by creating datastring and measuring)
    eval(['datstr = sprintf(''' char(param(varno).fmat) ''',DATA(1,varno));'])
    varwidth = length(datstr);
    varname = param(varno).name;
    varunit = param(varno).unit;
    blankvarname = repmat(' ',1,varwidth);

    % Make varname same length as fieldwidth
    if length(varname)>varwidth
        thisvarname = varname(1:varwidth); % crop variable name to field width
    else
        thisvarname = [blankvarname(1:varwidth-length(varname)) varname]; % pad front of short name
    end
    varnametxt = [varnametxt ' ' thisvarname]; % add param name to header line 19

    % Make param unit same length as fieldwidth
    if length(varunit)>varwidth
        thisvarunit = varunit(1:varwidth); % crop variable name to field width
    else
        thisvarunit = [blankvarname(1:varwidth-length(varunit)) varunit]; % pad front of short name
    end
    varunittxt = [varunittxt ' ' thisvarunit]; % add param unit to header line 20
    dataformat = [dataformat ' ' param(varno).fmat];
end
%dataformat = [dataformat '\n'];
linetxt(19) = {varnametxt};
linetxt(20) = {varunittxt};


%=================================================
% Open header template to read and new data file to write
headerfile = 'blank_header.txt';
fidh = fopen(headerfile,'rt');
if fidh == -1, error(['Could not open header file: ' headerfile]), end
fid = fopen(filename,'wt');
if fid == -1, error(['Could not open new data file: ' filename]), end
disp(['    Writing data to file ' filename])

%=================================================
% Write the headerlines
for lineno = 1:20;
    filler = fgetl(fidh);
    try
        hdrtxt = char(linetxt(lineno));
        hline = [hdrtxt filler(length(hdrtxt)+1:end)];
    catch
        hline = filler;
    end
    fprintf(fid,'%s\n',hline);
end
fclose(fidh);

%=================================================
% Write data
for ind = 1:size(DATA,1)
    eval(['datastr = sprintf(''' dataformat ''',DATA(ind,:));']);
    fprintf(fid,'%s\n',[datestr(time(ind),timeformat) datastr]);
end

% Close file
fclose(fid);
disp(['    ' int2str(ind) ' lines written to ' filename])

% Open file to check and edit with notepad
%mb = msgbox('Check file and edit if necessary.','WindowStyle','modal');
%uiwait(mb)
%disp('Close notepad when finished editing')
%eval(['!notepad ' FILENAME])
