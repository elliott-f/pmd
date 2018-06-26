function [imd,time,param] = read_ecotrip(FILENAME)

% [imd,time,param] = read_ecotrip(FILENAME)
% Read data from ecotriplet, convert to eng units 
% output format for PMD
%
% Craig Stewart
% 2011 April 15 - 
% - loosely based on from Mark Gall ECOtriplet_Master.r 



%===================================================================
% Choose file
if nargin == 0
    [filename,pathname] = uigetfile({'eco*.raw'},'Choose an ecotriplet ascii file to load');
    if pathname==0, disp('No file chosen, reading cancelled'), return, end
    FILENAME = [pathname filename];
    ext = FILENAME(end-3:end);
else
    [pathstr, name, ext, versn] = fileparts(FILENAME); pathname = [pathstr '\']; filename = [name ext];
end

% Set instrument metadata (imd) model
imd.model = 'ecotrip';
loc = strfind(filename,'_');
imd.snstr = filename(loc+1:end-4); % get the serial number from the filename (its not in the file)

% example data from start of file July 2010
% 22567 records to read
% 07/18/10	12:00:17	650	495	695	399	460	96	568
% 07/18/10	12:00:20	650	520	695	297	460	99	567
% 07/18/10	12:00:23	650	526	695	239	460	98	566

% Read in datafile
n = textread(FILENAME,'%u%*[^\n]',1); % check hom many lines to read (so we don't read etx line at end)
[timestr,raw_beta650,raw_chl,raw_cdom] = textread(FILENAME,'%18c%*d%d%*d%d%*d%d%*d',n,'headerlines',1);
time = datenum(timestr,'mm/dd/yy HH:MM:SS');

% Because eco files are soooo dirty filter them first and note all the bad
% line



% Process for ecotrip 480 channel setup and cal coefficients
if strcmp(imd.snstr,'480')
    
    
    % Define calibration coefficients (28 Jan 2008 SN 480)
    beta650_scale = 3.520e-06;
    beta650_dark = 67;
    chl_scale = 0.0125
    chl_dark = 41;
    cdom_scale = 0.0985
    cdom_dark = 31;
    
elseif
    
else
    error(['Cal data not available for ecotriplet serial number ' imd.snstr])
end
    
% Convert to Eng Units by applying cal coefficients
beta650 = beta650_scale .* (raw_beta650 - beta650_dark); % % Beta 650 (/m/sr) (scattering)
chl = chl_scale .* (raw_chl - chl_dark); % % Chlorophyll (?g/L = mg/m3)
cdom = cdom_scale .* (raw_cdom - cdom_dark); % % CDOM concentration (ppb)

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
theta_deg = 117;
theta = (pi/180)*theta_deg;
lambda = 650;
delta = 0.09;
S = 35; % estimated ocean salinity (PSU)
beta650w = 1.38*(lambda/500)^(-4.32) * (1 + 0.3*S/37)*10^(-4) * (1+(cos(theta))^2*(1-delta)/(1+delta)); % Water volume scattering
beta650p = beta650 - beta650w; % Particle volume scattering

bt_sw = 0.0029308*(lambda/500)^(-4.24); % Total scattering of sea water (35-39 PSU)
bb_sw = bt_sw/2; % Backscattering of seawater (half of total scattering)
% Calculate Particulate backscatter coefficient from Particulate volume scattering
X = 1.1; % scale factor from Boss and Pegau 2001 (ecotrip manual 5.2.2)
b_bp = 2*pi*X*beta650p; % Particulate backscatter coefficient (unit = m^-1)
b_b = b_bp + bb_sw; % Total backscatter coefficient (units = m^-1)


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

% unit notes:
for pn = 1:length(param)
    imd.notes(pn) = {[param(pn).name ' = ' param(pn).longname '(' param(pn).unit ')']};
end

%**************************************************************************
% Program ends
%**************************************************************************


% Mark's R code below
% ## For reading in data from ECOTRIPLET ##
% 
% ## Required packages ##
% require(stats)
% library(lattice)
% library(stats)
% library(utils)
% library(nls2)
% library(gregmisc)
% library(rggobi)
% library(playwith)
% 
% ## Reading in datafiles ##
% setwd("Q:/Laboratories/Methods/Ecotriplet/2Checked")
% filenames <- dir(pattern=glob2rx("*.raw"))
% numsample <- length(filenames)
% a <- sapply(filenames, read.delim, header = FALSE, sep = "\t", skip=0, nrows=60, simplify = FALSE, col.names=c("Date", "Time", "TimeNU", "beta650", "beta650NU", "CHL", "CHLNU", "CDOM", "CDOMNU")) # same as lapply except gives control over adding names
% a[[1]]
% summary(a)
% ## Deriving parameters from calibration file (#268) and known equations ##
% rnames <- 1:60
% cnames <- c("Beta660", "CHL","CDOM","Betap660","bbp660","bp660")
% Derive <- array (, c(60, 6, numsample), dimnames=list(rnames, cnames, filenames)) # mapping an empty array for filling
% 
% for (i in 1:numsample)
%   {
%   n <- length(a[[i]][,1])
%   Derive[1:n,1,i] <- 3.905E-06*(a[[i]][,"beta650"]-53) # Beta660 (/m/sr)
%   Derive[1:n,2,i] <- 0.0128*(a[[i]][,"CHL"]-43) # Chl (ug/L - mg/m3)
%   Derive[1:n,3,i] <- 0.0976*(a[[i]][,"CDOM"]-49) # CDOM (ppb)
%   }
% summary(Derive)
% 
% # Ideally should correct for apg660 but this is only a 4 % error at apg660 = 1/m, so small for most measures. Marine typically << 0.1 - 0.4 %.
% 
% # Volume scattering of Particles - FRESHWATER (Salinity estimated to be 5 PSU) #
% # Also (0.5 + 0.5*cos(2*117) = cos^2(117) #
% Derive[,4,] <- Derive[,1,]-(1.38*(660/500)^(-4.32)*(1+0.3*5/37)*(10^-4)*(1+(0.5+0.5*cos(2*117))*(1-0.09)/(1+0.09))) # Betap660 (/m/sr)
% Derive[,5,] <- 2*pi*1.1*Derive[,4,] # bbp660 (/m). Total particulate backscattering coefficient 
% Derive[,6,] <- Derive[,5,]+0.000348 # bp660 (/m). Total backscattering coefficient. Adding backscatter of pure freshwater is salinity dependant: bw = 0.0022533*(660/500)^(-4.23)/2 = 0.000348.   The 2 is backscatter, half of total scattering. 
% 
% # Volume scattering of Particles - MARINE (Salinity estimated to be 36 PSU) #
% # Also (0.5 + 0.5cos(2*117) = cos^2(117) #
% Derive[,4,] <- Derive[,1,]-(1.38*(660/500)^(-4.32)*(1+0.3*36/37)*(10^-4)*(1+(0.5+0.5*cos(2*117))*(1-0.09)/(1+0.09))) # Betap660 = Beta660 - Betaw660(/m/sr)
% Derive[,5,] <- 2*pi*1.1*Derive[,4,] # bbp660 (/m). Total particulate backscattering coefficient 
% Derive[,6,] <- Derive[,5,]+0.0004515 # bp660 (/m). Total backscattering coefficient. Adding backscatter of saltwater is salinity dependant: bsw = 0.0029308*(660/500)^(-4.24)/2 = 0.0004515. The 2 is backscatter, half of total scattering.
% 
% ## Summary ##
% Median <-  data.frame(array (0, c(numsample, 6), dimnames=list
% (filenames,cnames))) # mapping an empty array for filling
% SD <- Median
% N <- Median
% for (i in 1:numsample)
%   {
%   Median[i,] <- sapply(data.frame(Derive[,,i]),median,na.rm=T)
%   SD[i,] <- sapply(data.frame(Derive[,,i]),sd,na.rm=T)
%   }
% CV <- SD/sqrt(60)/Median*100
% 
% # Plotting first three derived parameters with time to check data #
% for (i in 1:numsample)
%   {
%   png(filename=paste(filenames[i],".png",sep = ""))
%   par(mfrow=c(1,3))
%   #par(ask=T)
%   plot(Derive[,1,i],type='o',ylim=c(0,max(Derive[,1,i],na.rm=T)))
%   plot(Derive[,2,i],type='o',ylim=c(0,max(Derive[,2,i],na.rm=T)), main=filenames[i])
%   plot(Derive[,3,i],type='o',ylim=c(0,max(Derive[,3,i],na.rm=T)))
%   dev.off()
%   }




