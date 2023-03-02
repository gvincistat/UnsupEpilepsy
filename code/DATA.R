# Description: This script produces the data sets EEG.Rdata and SURVEY.Rdata, which are used in 
# Giuseppe Vinci (2023). Unsupervised Learning in Epilepsy. In Statistical methods in epilepsy. Chapman & Hall/CRC. 
# Author: Giuseppe Vinci, Ph.D., Assistant Professor, Department of Applied and Computational Mathematics and Statistics, University of Notre Dame.


rm(list=ls())

######### SPECIFYING DIRECTORY WHERE ALL RESULTS WILL BE SAVED
FOLDER = '~/EPILEPSY'




######### INSTALLING R PACKAGES #########
# NOTE: install XQuartz from https://www.xquartz.org
if (!require(eegkit)) install.packages('eegkit')
if (!require(devtools)) install.packages('devtools')
library(devtools)
install_github("pisca46/edfReader", build_vignettes = TRUE)




######### DOWNLOAD EEG DATA ######### 
dir.create(paste0(FOLDER,'data'))
options(timeout=10000)
EDFURL = "https://physionet.org/files/siena-scalp-eeg/1.0.0/PN10/PN10-4.5.6.edf?download"
EDFPATH = paste0(FOLDER,'data/PN10-4.5.6.edf')
print(paste0('downloading EEG data from ', EDFURL, '...'))
download.file(url=EDFURL,destfile=EDFPATH)




######### LOAD EEG DATA ######### 
print(paste0('loading EDF file ', EDFPATH))
x=edfReader::readEdfHeader(EDFPATH)
dataset = edfReader::readEdfSignals(x)
dataset = dataset[c(1:8,14:23,34)] # selecting only active EEG channels: see https://physionet.org/content/siena-scalp-eeg/1.0.0/PN10/Seizures-list-PN10.txt
names(dataset)
print(paste0('Total length (sec): ',dataset[[1]]$total))




######### FEATURE EXTRACTION (POWER DENSITIES) #########
Signal = function(dataset,channel=1,start.sec=0,end.sec=10){
	rate = dataset[[channel]]$sRate
	dataset[[channel]]$signal[ceiling(start.sec*rate):ceiling(end.sec*rate)]
}

POWER = function(dataset,channel=1,start.sec=0,end.sec=1,lower=0, upper=30,Plot=FALSE){
	rate = dataset[[channel]]$sRate
	x = dataset[[channel]]$signal[ceiling(start.sec*rate):ceiling(end.sec*rate)]
	Power = eegkit::eegfft(x, Fs=rate, lower=lower, upper=upper)
	if(Plot==TRUE){
		plot(Power$frequency,10*log10(Power$strength^2),type='l',xlab='Frequency (Hz)',ylab='Power (dB)')
		abline(v=c(3,8,15,30),lty=2)
	}
	return(Power)
}

# This function converts a EEG waveform into a matrix with average power per freq bins for each time bin of length leng in the interval (start.sec, end.sec)
POWER.BINS = function(dataset,channel=1,start.sec=0,end.sec=10,nbins=5,freqs = FREQBINS){
	nfreq = length(freqs)-1
	M = matrix(NA,nc=nfreq,nr=nbins)
	bins = seq(start.sec,end.sec,len=nbins+1)
	for(i in 1:nbins){
		Power = POWER(dataset=dataset,channel=channel,start.sec=bins[i],end.sec=bins[i+1],lower=min(freqs), upper=max(freqs))
		FREQ = Power$frequency
		dB = 10*log10(Power$strength^2)
		for(j in 1:nfreq){
			M[i,j] = mean(dB[FREQ>=freqs[j] & FREQ<freqs[j+1]])
		}
	}
return(M)
}

Nw=length(dataset); bin.secs=10; FREQBINS=c(0,4,8,16,32,64,128)
START=1; END=16000
NBINS=END/bin.secs
print(paste0('EXTRACTING POWER DENSITIES from second ',START,' to ',END))
print(paste0('Number of EEG channels = ',Nw))
print(paste0('Frequency bins = ',paste0(FREQBINS,collapse='-')))
print(paste0('Time bin length = ',bin.secs, ' seconds'))
print(paste0('EEG recording portion = second ',START,' to ',END))
X = POWER.BINS(dataset,channel=1,start.sec=START,end.sec=END,nbins=NBINS,freqs = FREQBINS)
RAW = Signal(dataset,channel=1,start.sec=START,end.sec=END)
for(i in 2:Nw){
X = cbind(X,POWER.BINS(dataset,channel=i,start.sec=START,end.sec=END,nbins=NBINS))
RAW = cbind(RAW,Signal(dataset,channel=i,start.sec=START,end.sec=END))
}
print(paste0('Final data matrix X has n=',nrow(X),' rows (time bins) and p=',ncol(X),' columns (spectral features)'))




######### SEIZURES TIMES #########
# Seizures in PN10-4.5.6.edf
# N1. start 38min 29sec = 2309 sec || ends 38min 34sec = 2314 sec
# N2. start 6544 sec || end 6563 sec
# N3. start 11225  sec || end 11282 sec
seizures = rbind(c(2309,2314),c(6544,6563),c(11225,11282))
print(paste0('seizure times:'))
print(seizures)
seizure.bins.list = list(); seizure.rows = which(seizures[,2]<END & seizures[,2]>START)
for(j in 1:length(seizure.rows)){
	i = seizure.rows[j]
	seizure.bins.list[[j]] = floor((seizures[i,1]-START)/bin.secs):ceiling((seizures[i,2]-START)/bin.secs)
}
seizure.bins = unlist(seizure.bins.list)




######### SAVING EEG.Rdata #########
sRate = dataset[[1]]$sRate
EEG = list(X=X,RAW=RAW[1:(300*sRate),],seizures=seizure.bins.list,sRate=sRate)
print('saving the data')
save(EEG,file=paste0(FOLDER,'data/EEG.Rdata'))




######### GENERATING SURVEY DATA #########
gendata = function(n,m,catvec=c(2,3),mu.base){
	CAT = matrix(NA,nc=length(catvec),nr=m*100)
	for(j in 1:ncol(CAT)){
		CAT[,j] = sample(1:catvec[j],size=nrow(CAT),replace=TRUE)
	}
	CAT=unique(CAT)
	if(nrow(CAT)<m){
		print(paste0('m = ',nrow(CAT)))
		m = nrow(CAT)
	}
	CAT = CAT[1:m,]
	DATA = numeric()
	for(i in 1:nrow(CAT)){
		DATAj = matrix(NA,nr=n,nc=ncol(CAT)+length(mu.base))
		for(j in 1:ncol(CAT)){
			DATAj[1:n,j] = rep(CAT[i,j],n)
		}
		DATAj[,(ncol(CAT)+1):ncol(DATAj)] = MASS::mvrnorm(n,mu=mu.base+i,Sigma=diag(length(mu.base)))
		DATA = rbind(DATA,DATAj)
	}
return(list(DATA=DATA,m=m))
}

set.seed(1234)
n.each=20; m = 10; CATVEC = c(2,2,3,4); pnum=16
GENDATA = gendata(n=n.each,m=m,catvec=CATVEC,mu.base=rep(0,pnum))




######### SAVING SURVEY.Rdata #########
SURVEY = GENDATA$DATA
save(SURVEY,file=paste0(FOLDER,'/data/SURVEY.Rdata'))
