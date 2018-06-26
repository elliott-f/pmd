## For reading in data from ECOTRIPLET ##

## Required packages ##
require(stats)
library(lattice)
library(stats)
library(utils)
library(nls2)
library(gregmisc)
library(rggobi)
library(playwith)

## Reading in datafiles ##
setwd("Q:/Laboratories/Methods/Ecotriplet/2Checked")
filenames <- dir(pattern=glob2rx("*.raw"))
numsample <- length(filenames)
a <- sapply(filenames, read.delim, header = FALSE, sep = "\t", skip=0, nrows=60, simplify = FALSE, col.names=c("Date", "Time", "TimeNU", "LAMBDA", "LAMBDANU", "CHL", "CHLNU", "CDOM", "CDOMNU")) # same as lapply except gives control over adding names
a[[1]]
summary(a)
## Deriving parameters from calibration file (#268) and known equations ##
rnames <- 1:60
cnames <- c("Beta660", "CHL","CDOM","Betap660","bbp660","bp660")
Derive <- array (, c(60, 6, numsample), dimnames=list(rnames, cnames, filenames)) # mapping an empty array for filling

for (i in 1:numsample)
  {
  n <- length(a[[i]][,1])
  Derive[1:n,1,i] <- 3.905E-06*(a[[i]][,"LAMBDA"]-53) # Beta660 (/m/sr)
  Derive[1:n,2,i] <- 0.0128*(a[[i]][,"CHL"]-43) # Chl (ug/L - mg/m3)
  Derive[1:n,3,i] <- 0.0976*(a[[i]][,"CDOM"]-49) # CDOM (ppb)
  }
summary(Derive)

# Ideally should correct for apg660 but this is only a 4 % error at apg660 = 1/m, so small for most measures. Marine typically << 0.1 - 0.4 %.

# Volume scattering of Particles - FRESHWATER (Salinity estimated to be 5 PSU) #
# Also (0.5 + 0.5*cos(2*117) = cos^2(117) #
Derive[,4,] <- Derive[,1,]-(1.38*(660/500)^(-4.32)*(1+0.3*5/37)*(10^-4)*(1+(0.5+0.5*cos(2*117))*(1-0.09)/(1+0.09))) # Betap660 (/m/sr)
Derive[,5,] <- 2*pi*1.1*Derive[,4,] # bbp660 (/m). Total particulate backscattering coefficient 
Derive[,6,] <- Derive[,5,]+0.000348 # bp660 (/m). Total backscattering coefficient. Adding backscatter of pure freshwater is salinity dependant: bw = 0.0022533*(660/500)^(-4.23)/2 = 0.000348.   The 2 is backscatter, half of total scattering. 

# Volume scattering of Particles - MARINE (Salinity estimated to be 36 PSU) #
# Also (0.5 + 0.5cos(2*117) = cos^2(117) #
Derive[,4,] <- Derive[,1,]-(1.38*(660/500)^(-4.32)*(1+0.3*36/37)*(10^-4)*(1+(0.5+0.5*cos(2*117))*(1-0.09)/(1+0.09))) # Betap660 = Beta660 - Betaw660(/m/sr)
Derive[,5,] <- 2*pi*1.1*Derive[,4,] # bbp660 (/m). Total particulate backscattering coefficient 
Derive[,6,] <- Derive[,5,]+0.0004515 # bp660 (/m). Total backscattering coefficient. Adding backscatter of saltwater is salinity dependant: bsw = 0.0029308*(660/500)^(-4.24)/2 = 0.0004515. The 2 is backscatter, half of total scattering.

## Summary ##
Median <-  data.frame(array (0, c(numsample, 6), dimnames=list
(filenames,cnames))) # mapping an empty array for filling
SD <- Median
N <- Median
for (i in 1:numsample)
  {
  Median[i,] <- sapply(data.frame(Derive[,,i]),median,na.rm=T)
  SD[i,] <- sapply(data.frame(Derive[,,i]),sd,na.rm=T)
  }
CV <- SD/sqrt(60)/Median*100

# Plotting first three derived parameters with time to check data #
for (i in 1:numsample)
  {
  png(filename=paste(filenames[i],".png",sep = ""))
  par(mfrow=c(1,3))
  #par(ask=T)
  plot(Derive[,1,i],type='o',ylim=c(0,max(Derive[,1,i],na.rm=T)))
  plot(Derive[,2,i],type='o',ylim=c(0,max(Derive[,2,i],na.rm=T)), main=filenames[i])
  plot(Derive[,3,i],type='o',ylim=c(0,max(Derive[,3,i],na.rm=T)))
  dev.off()
  }



