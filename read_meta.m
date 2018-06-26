function meta = read_meta(FILENAME)
%-----------------------------------------------------------------------------------------
% function moormeta = read_meta(FILENAME)
%
% Part of the NIWA PMD toolbox
%
% Reads meta data for a mooring or instrument from a formatted text file
% and returns the data as a structure of character cells
%
% Each line contains the parameter name then a colon and the parameter
% value (as text), e.g.
% location     :Northern Biophysical
%
% Spaces before the colon are ignored - the string following the colon is taken as text
% Blank lines and lines begining with % are ignored
%
% HISTORY
% 01-Sep-2007 Craig Stewart
% 01-Jan-2009 Craig Stewart interpret lat lon and date strings
% 10-Jun-2009 Craig Stewart read instrument metadata
% 01-Apr-2016 Brett Grant	handle the case where there are multiple lines of notes in the meta file
%
% NIWA moorings
%-----------------------------------------------------------------------------------------


fid = fopen(FILENAME,'rt');
if fid == -1
    error(['Cannot open file ' FILENAME])
end

line_no = 0;
inststr = [];
instno = 0;
hardno = 0;
while ~feof(fid)
    line_no = line_no+1;
    tline = fgetl(fid);
    % Clean up line
    tline = strrep(tline,char(9),' '); % replace tabs with spaces
    if tline == -1
        tline = ' '; % Copes with empty last line (fgetl returns -1 for end of file)
    end
    if ~isempty(strrep(tline,' ','')) && ~strcmp(tline(1),'%')
        colonind = findstr(tline,':');
        if isempty(colonind)
            keyboard
            error(['colon not found in line ' int2str(line_no) ': ' tline])
        elseif length(colonind)>1
            %disp(['Warning: more than one colon found in line '
            %int2str(line_no) ': ' tline]) - don't warn - date time stamps
            %have this
            colonind = colonind(1);
        end
        
        param = tline(1:colonind-1);
        param = strrep(param,' ',''); % despace parameter name
        val = tline(colonind+1:end);
        if ~isempty(val)
            %param
            if strcmp(param,'model,sn,height')% Special case (1) - read instrument data for upright mooring
                instno = instno + 1;
                [meta.inst(instno).model,meta.inst(instno).snstr,meta.inst(instno).height] = strread(val,'%s%s%f','delimiter',',');
            elseif strcmp(param,'model,sn,height,notes')% Special case (2) - read instrument data with notes
                instno = instno + 1;
                [meta.inst(instno).model,meta.inst(instno).snstr,meta.inst(instno).height,meta.inst(instno).notes] = strread(val,'%s%s%f%s','delimiter',',');
            elseif strcmp(param,'model,sn,depth')% Special case (3) - read instrument data for hanging mooring
                instno = instno + 1;
                [meta.inst(instno).model,meta.inst(instno).snstr,meta.inst(instno).depth] = strread(val,'%s%s%f','delimiter',',');
%                 meta.inst(instno).height=str2num(meta.water_depth) - meta.inst(instno).depth;
            elseif strcmp(param,'model,sn,depth,notes')% Special case (4) - read instrument data with notes
                instno = instno + 1;
                [meta.inst(instno).model,meta.inst(instno).snstr,meta.inst(instno).depth,meta.inst(instno).notes] = strread(val,'%s%s%f%s','delimiter',',');
                meta.inst(instno).height=str2num(meta.water_depth)-meta.inst(instno).depth;
            elseif strcmp(param,'channel_name') % Special case (5) - read channel_name into a cell array of char - this is for read_rcm (rcm metadata)
                val = strrep(val,',',''',''');
                eval(['meta.' inststr param ' = {''' val '''};']);
            elseif strcmp(param,'type,unit,qty')% Special case (6) - read hardware data
                hardno = hardno + 1;
                [meta.hardware(hardno).type,meta.hardware(hardno).unit,meta.hardware(hardno).qty] = strread(val,'%s%s%f','delimiter',',');
            else
                % Generic input
                try
                    % first deal with notes (can be multiples)
                    if strcmp(param,'notes')
                        if isfield(meta,'notes')
                            meta.notes(end+1) = {val};
                        else
                            meta.notes = {val};
                        end
                    
                    else
                        %eval(['meta.' inststr param ' = ''' val ''';']);
                        eval(['meta.' param ' = ''' val ''';']);
                    end
                catch
                    disp(['Error reading meta file: ' FILENAME])
                    disp(['line ' int2str(line_no) ': ' tline])
                    disp(['param: ' param])
                    disp(['val  : ' val])
                    fclose(fid);
                    keyboard
                    error('Error reading meta file')
                end
            end
        end
    end
end
fclose(fid);

%set the water depth
%if strcmp(param,'water_depth')
    

% Convert known numeric parameters from string to double
if isfield(meta,'lat') && isfield(meta,'lon')
    [meta.lat,meta.lon] = llstr2num(meta.lat,meta.lon);
end
% Set the time format
if isfield(meta,'time_format')
    time_format = meta.time_format;
else
    time_format = 'yyyy/mm/dd HH:MM:SS'; % (default if not specified)
end
if isfield(meta,'start_time')
    try
        meta.start_dn = datenum(meta.start_time,time_format);
    catch
        %rmfield(meta,'start_time')
        disp(['Could not read start time ' meta.start_time ' in ' FILENAME])
    end
else
    disp(['start time not available in ' FILENAME])
end
if isfield(meta,'stop_time')
    try
        meta.stop_dn = datenum(meta.stop_time,time_format);
    catch
        %rmfield(meta,'stop_time')
        disp(['Could not read stop time ' meta.stop_time ' in ' FILENAME])
    end
else
    disp(['stop time not available in ' FILENAME])
end
% BG April 2016 - commented next three lines
% if isfield(meta,'notes')
%     meta.notes = {meta.notes};
% end

% Special cases - convert str2num
if isfield(meta,'utc_offset')
    meta.utc_offset = str2num(meta.utc_offset);
end
if isfield(meta,'water_depth')
    meta.water_depth = str2num(meta.water_depth);
end
if isfield(meta,'declination')
    meta.declination = str2num(meta.declination);
end
if isfield(meta,'interval_min')
    meta.interval_min = str2num(meta.interval_min);
end
if isfield(meta,'sn')
    meta.snstr = num2str(meta.sn);
end
if isfield(meta,'dsu_lag_s')
    meta.dsu_lag_s = str2num(meta.dsu_lag_s);
end