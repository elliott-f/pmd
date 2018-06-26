function join_dat(filelist,varlist,outputfilename)

% joins data from dat files
%
% Craig Stewart 20 Oct 2009
% 2010 Nov 7-10 Total re-write 

% Check args
if ~exist('outputfilename','var')
    outputfilename = input('Enter an output file name (no ext): ','s');
end

if isstr(filelist)
    if strcmp(filelist,'all')
        fl = dir('*.dat*');
        filelist = [{fl.name}]';
        filelist(1:2) = []; % directories
    end
end

% Preallocate
DATAMAT = []; % Data array
nv = length(varlist);
nf = length(filelist);

% Loop through files loading
for fn = 1:nf;
    
    % Load data
    %FILENAME = [filelist(fileno).path filelist(fileno).name];
    filename = char(filelist(fn));
    data = read_datx(filename);
    time = data.time;
    thisnomdep = data.moormeta.water_depth-data.instmeta.height;
    nomdepvec = repmat(thisnomdep,size(data.time));
    filenovec = repmat(fn,size(data.time));
    if ismember('sal',varlist) && isfield(data,'pre'); % then its a microcat: pre = dep
        data.dep = data.pre; % put pre into dep
        data = rmfield(data,'pre');
    end
    
    % Join data into single matrix    
    % create blank matrix and fill with data from this file
    datamat = repmat(nan*data.time,1,nv);
    for varno = 1:nv
        if isfield(data,varlist(varno))
            datamat(:,varno) = eval(['data.' char(varlist(varno))]);
        end
    end
    
    % Append to master data matrix - with new variables nomdepth and file#
    DATAMAT = [DATAMAT; datamat nomdepvec filenovec];
end
% Update variable names
varlist = [varlist {'nomdep','fileno'}];
nv = length(varlist);

% Sort chronologically
DATAMAT = sortrows(DATAMAT);
% Check time is always increasing
dt = diff(DATAMAT(:,1));
if min(dt<=0)
    badinds = find(dt<=0);
    error(['Warning: DT<=0 at ' datestr(DATAMAT(badinds(1),1))])
end

%% Save data as mat file
% Unpack data from matrix
for varno = 1:nv
    eval([char(varlist(varno)) ' = DATAMAT(:,varno);']);
end
outvarlist = ['time'];
for ind = 1:numel(varlist)
    outvarlist = [outvarlist ' ' char(varlist(ind))];
end
eval(['save ' outputfilename ' ' outvarlist])