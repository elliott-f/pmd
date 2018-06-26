%-----------------------------------------------------------------------------------------
% process_AWAC
%
% A new set of ancilary data generated at the start of every burst, this is
% contained in whd. The transformation matrix must be recalculated every time
% the orientation, heading, pitch or roll changes. - So calculate a new one
% at the start of every burst
% 
% ancillary = whd
% burst = wad
%
% Check .hdr file for configuration data and sensor descriptions for all
% files
% Read in hdr ??? T matrix, wave settings
%
% Refer to Transform.m from nortek
%
% WAVES: may want to use .stg file instead/as well???
%  .whd - Waves header file
%  1   Month                            (1-12)
%  2   Day                              (1-31)
%  3   Year
%  4   Hour                             (0-23)
%  5   Minute                           (0-59)
%  6   Second                           (0-59)
%  7   Burst counter
%  8   No of wave data records
%  9   Cell position                    (m)
% 10   Battery voltage                  (V)
% 11   Soundspeed                       (m/s)
% 12   Heading                          (degrees)
% 13   Pitch                            (degrees)
% 14   Roll                             (degrees)
% 15   Minimum pressure                 (dbar)
% 16   Maximum pressure                 (dbar)
% 17   Temperature                      (degrees C)
% 18   CellSize                         (m)
% 19   Noise amplitude beam 1           (counts)
% 20   Noise amplitude beam 2           (counts)
% 21   Noise amplitude beam 3           (counts)
% 22   Noise amplitude beam 4           (counts)
% 23   AST window start                 (m)
% 24   AST window size                  (m)
% 25   AST window offset                (m)
% 
% .wad - Waves data file (in beam coordinates)
%  1   Burst counter
%  2   Ensemble counter
%  3   Pressure                         (dbar)
%  4   AST Distance1 (Beam4)            (m)
%  5   AST Distance2 (Beam4)            (m)
%  6   AST Quality (Beam4)              (counts)
%  7   Analog input
%  8   Velocity (Beam1)                 (m/s)
%  9   Velocity (Beam2)                 (m/s)
% 10   Velocity (Beam3)                 (m/s)
% 11   Amplitude (Beam1)                (counts)
% 12   Amplitude (Beam2)                (counts)
% 13   Amplitude (Beam3)                (counts)
%
% History:
% 01-Sep-2016 FE Adapted from Transform.m
% 22-Jan-2017 FE Using .whd instead of .sen file
%
% NIWA Moorings
%-----------------------------------------------------------------------------------------
clear
dbstop if error

%% --- Constants ----
% Transformation matrix for beam to xyz coordinates,
% the values can be found from the header file that is generated in the
% conversion of data to ASCII format
% Transformation matrix for WPR 1685
T = [1.5774 -0.7891 -0.7891
     0.0000 -1.3662 1.3662
     0.3677 0.3677  0.3677];

% Instrument orientation
statusbit0 = 0; % Up

% ???? T = T/4096;   % Scale the transformation matrix correctly to floating point numbers

% If instrument is pointing down (bit 0 in status equal to 1)
% rows 2 and 3 must change sign
% NOTE: For the Vector the instrument is defined to be in
%       the UP orientation when the communication cable
%       end of the canister is on the top of the canister, ie when
%       the probe is pointing down.
%       For the other instruments that are used vertically, the UP
%       orientation is defined to be when the head is on the upper
%       side of the canister.
if statusbit0 == 1,
    T(2,:) = -T(2,:);
    T(3,:) = -T(3,:);
end

latitude = -41.2548;


%% --- Data Files ---
% WHD
[filename, pathname] = uigetfile('*_whd.mat','Select a waves header .whd mat file');
whd = load(fullfile(pathname,filename));

% WAD
[filename, pathname] = uigetfile('*_wad.mat','Select a waves data .wad mat file');
wad = load(fullfile(pathname,filename));


%% --- Process Data - Convert Beam to ENU ---
% First do some checks

% Number of bursts - this is usually different
disp(['Number of bursts in wave header (WHD): ' num2str(max(whd.burstCounter))])
disp(['Number of bursts in wave data (WAD)  : ' num2str(max(wad.burstCounter))])
% if max(wad.burstCounter) ~= max(whd.burstCounter)
%     error('Total number of bursts in WAD and WHD files is not the same, cancelling processing')
% end

% noWaveDataRecords is constant
if min(whd.noWaveDataRecords) ~= max(whd.noWaveDataRecords)
    error('Total number of samples in each burst is not the constant, cancelling processing') 
end

nBursts = whd.burstCounter(end);

% Create timestamps
whdTime = datenum(whd.year,whd.month,whd.day,whd.hour,whd.minute,whd.second);
whd.time = whdTime;
wadTime = ones(size(wad.burstCounter)); % Prepopulate variable
wadTime = (whdTime(wad.burstCounter) + (wad.ensembleCounter/86400)); % Use start of burst time + 1sec sampling rate 
wad.time = wadTime;
% Note - 1 second out?

disp(['Total number of bursts: ' num2str(nBursts)])
disp(['Samples per burst: ' num2str(numel(find(wad.burstCounter==1)))])
disp(['Start time: ' datestr(whdTime(1))])
disp(['End time: ' datestr(whdTime(end))])
disp('This script assumes a 1Hz sampling rate')

% Heading, pitch and roll are the angles output in the data in degrees
headingRad  = pi*(whd.heading-90)/180;
pitchRad    = pi*whd.pitch/180;
rollRad     = pi*whd.roll/180;

% Make heading matrix
% H = [cos(headingRad)    sin(headingRad) 0; 
%     -sin(headingRad)    cos(headingRad) 0; 
%     0                   0               1];
H = zeros(3,3,nBursts);
H(1,1,:) = cos(headingRad);
H(1,2,:) = sin(headingRad);
H(2,1,:) = -sin(headingRad);
H(2,2,:) = cos(headingRad);
H(3,3,:) = 1;

% Make tilt matrix
% P = [cos(pitchRad) -sin(pitchRad)*sin(rollRad) -cos(rollRad)*sin(pitchRad);
%      0             cos(rollRad)                -sin(rollRad);  
%      sin(pitchRad) sin(rollRad)*cos(pitchRad)  cos(pitchRad)*cos(rollRad)];
P = zeros(3,3,nBursts);
P(1,1,:) = cos(pitchRad);
P(1,2,:) = -sin(pitchRad).*sin(rollRad);
P(1,3,:) = -cos(rollRad).*sin(pitchRad);
P(2,2,:) = cos(rollRad);
P(2,3,:) = -sin(rollRad);
P(3,1,:) = sin(pitchRad);
P(3,2,:) = sin(rollRad).*cos(pitchRad);
P(3,3,:) = cos(pitchRad).*cos(rollRad);
  
% Beam is beam coordinates, for example beam = [0.23 ; -0.52 ; 0.12]
% enu is ENU coordinates
% enu = R*beam
beam = [wad.velocityBeam1 wad.velocityBeam2 wad.velocityBeam3]';
enu = zeros(size(beam));

% Make resulting transformation matrix 
% R = H*P*T
% Figure out how to do this via array multiplication
R = zeros(size(H));
for i = 1:nBursts
    R(:,:,i) = H(:,:,i)*P(:,:,i)*T;
end

% Given beam velocities, ENU coordinates are calculated as
% enu = R * beam;
disp('Calculating ENU from beam coordinates...')
h = waitbar(0,'Calculating ENU from beam coordinates...');
for i = 1:nBursts
    enu(:,wad.burstCounter==i) = R(:,:,i) * beam(:,wad.burstCounter==i);
    waitbar(i/nBursts,h)
end
close(h)

% Transform
enu = enu';

wad.velocityEast    = enu(:,1);
wad.velocityNorth   = enu(:,2);
wad.velocityUp      = enu(:,3);

figDim = [2 2 17.5 8];
axesDim = [ 1.5 1 15 6.2];


%% Plot to get an idea of data
% Plot 1 - AST and Pressure
disp('Create plot to get an idea of data')
disp('Plot AST and Depth (m) from all ensembles')
col = linspecer(7);
h1 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz - Waves AST and Depth','Tag','AWAC1MHz-AST');
ah1 = axes('Units','cen','Position',axesDim);
plot(wadTime,wad.ASTDistance1Beam4,'col',col(2,:),'marker','.','linestyle','none')
hold on
plot(wadTime,wad.ASTDistance2Beam4,'col',col(3,:),'marker','.','linestyle','none')
plot(wadTime,sw_dpth(wad.pressure,latitude),'col',col(4,:),'marker','.','linestyle','none') % convert to dep with sw toolbox
legend({'ASTDistance1Beam4 (m)','ASTDistance2Beam4 (m)','Depth (m)'})
title('AWAC AST and Depth - Wellington Harbour Monitoring')
ylabel('metres')
% Add timestamps
pmd_timestamps(gca,wadTime);
% Tidy!
pmd_style(h1)

% Plot 2 - AST and Pressure with AST quality - TODO - plot on different axes
disp('Plot AST quality from all ensembles')
h2 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz - Waves AST quality','Tag','AWAC1MHz-AST_quality');
ah2 = axes('Units','cen','Position',axesDim);
plot(wadTime,wad.ASTQualityBeam4,'col',col(1,:),'marker','.','linestyle','none')
hold on
% plot(wadTime,wad.ASTDistance1Beam4,'col',col(2,:),'marker','.','linestyle','none')
% plot(wadTime,wad.ASTDistance2Beam4,'col',col(3,:),'marker','.','linestyle','none')
% plot(wadTime,sw_dpth(wad.pressure,latitude),'col',col(4,:),'marker','.','linestyle','none') % convert to dep with sw toolbox
% legend({'ASTQualityBeam4 (counts)','ASTDistance1Beam4 (m)','ASTDistance2Beam4 (m)','Depth (m)'})
title('AWAC AST quality - Wellington Harbour Monitoring')
ylabel('counts')
% Add timestamps
pmd_timestamps(gca,wadTime);
% Tidy!
pmd_style(h2)

% Plot 3 - Velocity
disp('Plot east, north and vertical velocities (m/s) from all ensembles')
h3 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz - Waves velocity','Tag','AWAC1MHz-Velocity_all');
ah3 = axes('Units','cen','Position',axesDim);
plot(wadTime,enu(:,1),'col',col(1,:),'marker','.','linestyle','none')
hold on
plot(wadTime,enu(:,2),'col',col(2,:),'marker','.','linestyle','none')
plot(wadTime,enu(:,3),'col',col(3,:),'marker','.','linestyle','none')
set(gca,'YLim',[-1 1])
legend({'East','North','Up'})
title('AWAC waves velocity data - Wellington Harbour Monitoring')
ylabel('m/s')
% Add timestamps
pmd_timestamps(gca,wadTime);
% Tidy!
pmd_style(h3)

% Plot 4 - Amplitude
disp('Plot beam amplitudes (counts) from all ensembles')
h4 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz - Waves amplitude','Tag','AWAC1MHz-Amplitude');
ah4 = axes('Units','cen','Position',axesDim);
plot(wadTime,wad.amplitudeBeam1,'col',col(1,:),'marker','.','linestyle','none')
hold on
plot(wadTime,wad.amplitudeBeam2,'col',col(2,:),'marker','.','linestyle','none')
plot(wadTime,wad.amplitudeBeam3,'col',col(3,:),'marker','.','linestyle','none')
legend({'amplitudeBeam1','amplitudeBeam2','amplitudeBeam3'})
title('AWAC beam amplitudes - Wellington Harbour Monitoring')
ylabel('Counts')
% Add timestamps
pmd_timestamps(gca,wadTime);
% Tidy!
pmd_style(h4)

% Plot 5 - cellPosition
disp('Plot various depths (m) from burst data')
h5 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz','Tag','AWAC1MHz-depths');
ah5 = axes('Units','cen','Position',axesDim);
plot(whdTime,whd.cellPosition,'col',col(1,:),'marker','.');
hold on
plot(whdTime,sw_dpth(whd.minimumPressure,latitude),'col',col(2,:),'marker','.'); % convert to dep with sw toolbox
plot(whdTime,sw_dpth(whd.maximumPressure,latitude),'col',col(3,:),'marker','.'); % convert to dep with sw toolbox
plot(whdTime,whd.cellSize,'col',col(4,:),'marker','.'); 
plot(whdTime,whd.astWindowStart,'col',col(5,:),'marker','.'); 
plot(whdTime,whd.astWindowSize,'col',col(6,:),'marker','.'); 
plot(whdTime,whd.astWindowOffset,'col',col(7,:),'marker','.'); 
legend({'cellPosition','minimumDepth','maximumDepth','cellSize','astWindowStart','astWindowSize','astWindowOffset'})
title('AWAC cell heights - Wellington Harbour Monitoring')
ylabel('m')
% Add timestamps
pmd_timestamps(gca,whdTime);
% Tidy!
pmd_style(h5)

% Plot 6 - noiseAmplitudeBeam
disp('Plot noise amplitude beam (counts) from burst data')
h6 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz - noiseAmplitudeBeam','Tag','AWAC1MHz-noiseAmplitudeBeam');
ah6 = axes('Units','cen','Position',axesDim);
plot(whdTime,whd.noiseAmplitudeBeam1,'col',col(1,:),'marker','.');
hold on
plot(whdTime,whd.noiseAmplitudeBeam2,'col',col(2,:),'marker','.'); 
plot(whdTime,whd.noiseAmplitudeBeam3,'col',col(3,:),'marker','.'); 
plot(whdTime,whd.noiseAmplitudeBeam4,'col',col(4,:),'marker','.');
legend({'noiseAmplitudeBeam1','noiseAmplitudeBeam2','noiseAmplitudeBeam3','noiseAmplitudeBeam4'})
title('AWAC noise amplitude beam - Wellington Harbour Monitoring')
ylabel('counts')
% Add timestamps
pmd_timestamps(gca,whdTime);
% Tidy!
pmd_style(h6)

% Plot 7 - Attitide sensors
disp('Plot heading, pitch and roll from burst data')
h7 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz Attitude','Tag','AWAC1MHz-Attitude');
ah7 = axes('Units','cen','Position',axesDim);
plot(whdTime,whd.heading,'col',col(1,:),'marker','.');
hold on
plot(whdTime,whd.pitch,'col',col(2,:),'marker','.');
plot(whdTime,whd.roll,'col',col(3,:),'marker','.');
legend({'heading','pitch','roll'})
title('AWAC ancillary data - Wellington Harbour Monitoring')
ylabel('degrees')
% Add timestamps
pmd_timestamps(gca,whdTime);
% Tidy!
pmd_style(h7)

% Plot 8 - Temperature and sound speed - TODO - plot on different axes
disp('Plot temperature and sound speed from burst data')
h8 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz Temperature and sound speed ','Tag','AWAC1MHz-TSS');
ah8 = axes('Units','cen','Position',axesDim);
plot(whdTime,whd.temperature,'col',col(1,:),'marker','.');
hold on
plot(whdTime,whd.soundSpeed,'col',col(2,:),'marker','.');
legend({'temperature','soundSpeed'})
title('AWAC temperature and sound speed - Wellington Harbour Monitoring')
ylabel('degrees')
% Add timestamps
pmd_timestamps(gca,whdTime);
% Tidy!
pmd_style(h8)


%% --- Process Data - Average Burst ---
disp('Averaging ensembles')
% Get burst averages
wadAvg.time                         = whdTime; % or should we take time from halfway through the burst??
wadAvg.burstCounter                 = zeros(nBursts,1);
wadAvg.pressure                     = zeros(nBursts,1);
wadAvg.ASTDistance1Beam4            = zeros(nBursts,1);
wadAvg.ASTDistance2Beam4            = zeros(nBursts,1);
wadAvg.ASTQualityBeam4              = zeros(nBursts,1);
wadAvg.velocityBeam1                = zeros(nBursts,1);
wadAvg.velocityBeam2                = zeros(nBursts,1);
wadAvg.velocityBeam3                = zeros(nBursts,1);
wadAvg.velocityEast                 = zeros(nBursts,1);
wadAvg.velocityNorth                = zeros(nBursts,1);
wadAvg.velocityUp                   = zeros(nBursts,1);
wadAvg.amplitudeBeam1               = zeros(nBursts,1);
wadAvg.amplitudeBeam2               = zeros(nBursts,1);
wadAvg.amplitudeBeam3               = zeros(nBursts,1);

h = waitbar(0,'Calculating burst averages...');
for i = 1:nBursts
    index                           = wad.burstCounter == i;    
    wadAvg.burstCounter(i,:)        = i;
    wadAvg.pressure(i,:)            = mean(wad.pressure(index));
    wadAvg.ASTDistance1Beam4(i,:)   = mean(wad.ASTDistance1Beam4(index));
    wadAvg.ASTDistance2Beam4(i,:)   = mean(wad.ASTDistance2Beam4(index));
    wadAvg.ASTQualityBeam4(i,:)     = mean(wad.ASTQualityBeam4(index));
    wadAvg.velocityBeam1(i,:)       = mean(wad.velocityBeam1(index));
    wadAvg.velocityBeam2(i,:)       = mean(wad.velocityBeam2(index));
    wadAvg.velocityBeam3(i,:)       = mean(wad.velocityBeam3(index));
    wadAvg.velocityEast(i,:)        = mean(enu(index,1));
    wadAvg.velocityNorth(i,:)       = mean(enu(index,2));
    wadAvg.velocityUp(i,:)          = mean(enu(index,3));
    wadAvg.amplitudeBeam1(i,:)      = mean(wad.amplitudeBeam1(index));
    wadAvg.amplitudeBeam2(i,:)      = mean(wad.amplitudeBeam2(index));
    wadAvg.amplitudeBeam3(i,:)      = mean(wad.amplitudeBeam3(index));
    waitbar(i/nBursts,h)
end
close(h)


%% --- More plots ---
% Plot 9 - Velocity averaged
disp('Plot east, north and vertical velocities (m/s) from averaged ensembles')
h9 = figure('color','white','numberTitle','off',...
    'Unit', 'cen','Position', figDim,...
    'Name','AWAC 1MHz - Waves velocity','Tag','AWAC1MHz-Velocity_avg');
ah9 = axes('Units','cen','Position',axesDim);
plot(whdTime,wadAvg.velocityEast,'col',col(1,:),'marker','.')
hold on
plot(whdTime,wadAvg.velocityNorth,'col',col(2,:),'marker','.')
plot(whdTime,wadAvg.velocityUp,'col',col(3,:),'marker','.')
legend({'East','North','Up'})
title('AWAC waves velocity data - Wellington Harbour Monitoring')
ylabel('m/s')
set(gca,'YLim',[-0.5 0.5])
% Add timestamps
pmd_timestamps(gca,whdTime);
% Tidy!
pmd_style(h9)


% disp('Plot pressure, east, north and vertical velocities (m/s) from averaged ensembles')
% h10 = figure('color','white','numberTitle','off',...
%     'Unit', 'cen','Position', figDim,...
%     'Name','AWAC 1MHz - Waves velocity','Tag','AWAC1MHz-Velocity_subplot');
% ax(1) = subplot(311);
% plot(whdTime,wadAvg.pressure,'col',col(1,:),'marker','.')
% ylabel('Pre (dBar)')
% title('AWAC waves velocity data - Wellington Harbour Monitoring')
% 
% ax(2) = subplot(312);
% plot(whdTime,wadAvg.velocityEast,'col',col(2,:),'marker','.');
% hold on
% plot(whdTime,wadAvg.velocityNorth,'col',col(3,:),'marker','.');
% set(gca,'YLim',[-0.3 0.3])
% legend('Ve','Vn')
% ylabel('V (m/s)')
% 
% ax(3) = subplot(313);
% plot(whdTime,wadAvg.velocityUp,'col',col(4,:),'marker','.');
% set(gca,'YLim',[-0.1 0.1])
% ylabel('V_n (m/s)')
% linkaxes([ax(1),ax(2),ax(3)],'x')
% 
% subplotDim = get(,'Position',axesDim);
% 
% 
% % Add timestamps
% pmd_timestamps(ax,whdTime);
% % Tidy!
% pmd_style(h10)

