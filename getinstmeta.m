function instmeta = getinstmeta(meta,model,snstr)

% instmeta = getinstmeta(meta,model,snstr)
% Extracts instrument metadata from mooring metadata (files created after 29 March 2009)
% Craig Stewart
% 4 June 2009
%
% args: meta = structure of mooring and instrument metadata
%       model = string of instrument model name
%       snstr = integer or string (if containing letters)

% Get instrument metadata
if isfield(meta,'inst')
    for ind = 1:numel({meta.inst.snstr})
        if strcmpi(snstr,meta.inst(ind).snstr)
            disp(['ind: ' int2str(ind) ' sn match'])
            %keyboard
            if strcmpi(meta.inst(ind).model,model)
                instind = ind; % model and serial number match
            end
        end
    end
    if ~exist('instind','var')
        disp(['Error while searching for instrument metadata: Match not found for Model: ' model ', SN: ' snstr])
        keyboard
        error(['Error while searching for instrument metadata: Match not found for Model: ' model ', SN: ' snstr])
    end
    instmeta = meta.inst(instind); 
    instmeta.instid = [char(instmeta.model) ' ' char(instmeta.snstr)];
else
    error('Error while searching for instrument metadata: Field "inst" not found in metadata structure')
end


if ~isfield(meta,'time_format')
    meta.time_format = 'dd/mm/yyyy HH:MM:SS'; % (default if not specified)
end

if isfield(instmeta,'start_time')
    if ~isempty(instmeta.start_time)
        instmeta.start_dn = datenum(instmeta.start_time,meta.time_format);
    end
end
if isfield(instmeta,'stop_time')
    if ~isempty(instmeta.stop_time)
        instmeta.stop_dn = datenum(instmeta.stop_time,meta.time_format);
    end
end
if isfield(instmeta,'notes')
    if ~isa(instmeta.notes,'cell')
        instmeta.notes = {instmeta.notes};
    end
end

% Get instrument height
% if ~isfield(instmeta,'height')
%     error('Error: instrument height not defined in metadata')
% end