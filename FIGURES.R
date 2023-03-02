# Description: This script produces the figures 2-8 presented in
# Giuseppe Vinci (2023). Unsupervised Learning in Epilepsy. In Statistical methods in epilepsy. Chapman & Hall/CRC. 
# The following code uses the R package unSuper available at https://github.com/gvincistat/unSuper 
# Author: Giuseppe Vinci, Ph.D., Assistant Professor, Department of Applied and Computational Mathematics and Statistics, University of Notre Dame.


rm(list=ls())

######### SPECIFYING DIRECTORY CONTAINING THE FOLDER "data" CONTAINING EEG.Rdata and SURVEY.Rdata
FOLDER = '~/EPILEPSY'




######### LOADING R PACKAGE AND DATA #########
library(unSuper)
load(paste0(FOLDER,'data/EEG.Rdata'))
load(paste0(FOLDER,'data/SURVEY.Rdata'))
attach(EEG); Nw = ncol(RAW)
RESCALE = FALSE




######### CREATE FOLDER WHERE ALL FIGURES WILL BE STORED #########
dir.create(paste0(FOLDER,'figures'))




######### FUNCTIONS FOR VISUALISATION #########
PLOT.waves = function(RAW,channels=1:2,start.sec=0,end.sec=10, rate=512){
 par(mar=c(0,0,0,0))
 for(channel in channels){
	plot(RAW[ceiling(start.sec*rate):ceiling(end.sec*rate),channel],type='l',axes=FALSE,xlab='',ylab='',lwd=0.5)
 }
}

plot.seizure2D = function(Y,seizure.bins.list,legendpos='topleft',Legend=TRUE,...){
	for(j in 1:length(seizure.bins.list)){
	lines(Y[seizure.bins.list[[j]],1],Y[seizure.bins.list[[j]],2],lwd=1,col='blue')
	points(Y[seizure.bins.list[[j]][1],1],Y[seizure.bins.list[[j]][1],2],pch=15,col='red',cex=1)
	points(Y[seizure.bins.list[[j]][length(seizure.bins.list[[j]])],1],Y[seizure.bins.list[[j]][length(seizure.bins.list[[j]])],2],pch=17,col='seagreen',cex=1)
}
if(Legend==TRUE){
legend(legendpos,c('start seizure','end seizure'),pch=c(15,17),col=c('red','seagreen'),bty='n',...)
}
}




######### FIGURE 2: PCA #########
set.seed(123)
CPVEmin = 0.9
PCA = reduction(X,q=2,method='pca',rescale=RESCALE,cpve=CPVEmin,Plot=FALSE)
PC = PCA$Z

pdf(paste0(FOLDER,"figures/FIG2-PCA.pdf"),height=6,width=10)
par(oma=c(2,0,2,0))
COL = rgb(0,0,0,0.2)
layout(matrix(c(rep(1:Nw,10),rep(Nw+1,10*Nw),rep(c(rep(Nw+2,floor(Nw/2)),rep(Nw+3,Nw-floor(Nw/2))),10)),nc=30))
PLOT.waves(RAW,channels=1:Nw,100,160)
par(mar=c(2,2,1,2))
image(X,axes=FALSE,xlab='',ylab='')
mtext("time bins",side=1,line=1)
mtext("features",side=2,line=1)
par(mar=c(4,4.5,2,1))
plot(PCA$pca.extra$CPVE,ylim=c(0,1),pch=20,cex.lab=1.5,xlab='number of PCs (q)',ylab='CPVE')
abline(h=CPVEmin,col='blue')
abline(v=PCA$pca.extra$q.star)
plot(PC[,1],PC[,2],xlab="PC1",ylab="PC2",cex.lab=1.5,col=COL,pch=16,cex=1)
plot.seizure2D(PC,seizures,legendpos='topleft')
mtext('A',side=3,outer=TRUE,at=0.01,cex=1.5)
mtext('B',side=3,outer=TRUE,at=0.35,cex=1.5)
mtext('C',side=3,outer=TRUE,at=0.67,cex=1.5)
mtext('D',side=3,outer=TRUE,at=0.67,cex=1.5,line=-22)
dev.off()




######### FIGURE 3: MDS, ISOMAP, t-SNE, UMAP #########
MDS = reduction(X, method='mds', q=2)

K.neigh = 10

if (!require("BiocManager", quietly = TRUE)){
   install.packages("BiocManager")
   BiocManager::install("RDRToolbox")}
ISOMAP = reduction(X, method='isomap',q=2, neighbors=K.neigh, rescale=RESCALE)

set.seed(123)
TSNE = reduction(X,method='tsne',q=2,neighbors=K.neigh, rescale=RESCALE, verbose=TRUE, max_iter=5000)

set.seed(123)
UMAP = reduction(X, method='umap', q=2, neighbors=K.neigh, min_dist = 0.1, rescale=RESCALE, metric='euclidean')

pdf(paste0(FOLDER,"figures/FIG3-MDStSNEISOMAPUMAP.pdf"),height=10,width=12)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4.5,1,1))
COL = rgb(0,0,0,0.3)
plot(MDS$Z, xlab=expression(MDS[1]), ylab=expression(MDS[2]),main='',col=COL,pch=16,cex.lab=1.6)
plot.seizure2D(MDS$Z,seizures,legendpos='topleft',cex=1.5)
plot(ISOMAP$Z, xlab=expression(ISOMAP[1]), ylab=expression(ISOMAP[2]),col=COL,pch=16,cex.lab=1.6)
plot.seizure2D(ISOMAP$Z,seizures,Legend=FALSE)
plot(TSNE$Z, xlab=expression(t-SNE[1]), ylab=expression(t-SNE[2]),col=COL,pch=16,cex.lab=1.6)
plot.seizure2D(TSNE$Z,seizures,Legend=FALSE)
plot(UMAP$Z, xlab=expression(UMAP[1]), ylab=expression(UMAP[2]),col=COL,pch=16,cex.lab=1.6)
plot.seizure2D(UMAP$Z,seizures,Legend=FALSE)
dev.off()




######### FIGURE 4: ICA #########
set.seed(123)
n.ica = 6
nsec = 60
ICA = reduction(RAW[1:(nsec*sRate),], method='ica',rescale=RESCALE, q=n.ica, tol=0.0000001, verbose=TRUE, maxit=200)

pdf(paste0(FOLDER,"figures/FIG4-ICA.pdf"),height=6,width=10)
par(oma=c(2,0,2,0),mar=c(0,0,0,0))
layout(matrix(c(rep(1:Nw,each=n.ica),rep(1:n.ica,each=Nw)+Nw),nc=2))
for(i in 1:Nw){
	plot(RAW[ceiling(0*sRate):ceiling(nsec*sRate),i],type='l',axes=FALSE,xlab='',ylab='',lwd=0.5)
}
for(i in 1:n.ica){
	plot(ICA$Z[ceiling(0*sRate):ceiling(nsec*sRate),i],type='l',axes=FALSE,xlab='',ylab='',lwd=0.5)
}
mtext('A',side=3,outer=TRUE,at=0.01,cex=1.5)
mtext('B',side=3,outer=TRUE,at=0.5,cex=1.5)
dev.off()




######### FIGURE 5: K-MEANS CLUSTERING #########
set.seed(123)
K.clusters = c(6,4,2)
KMEANS=lapply(K.clusters,function(k) clustering(X, K=k, rescale=RESCALE, iter.max=10000, nstart=500))
PC = reduction(X,q=2,method='pca',rescale=RESCALE,Plot=FALSE)$Z

pdf(paste0(FOLDER,"figures/FIG5-KMEANS.pdf"),height=4,width=10)
layout(matrix(c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,
	3,3,3,3,4,3,3,3,3,4,3,3,3,3,4,3,3,3,3,4,
	5,5,5,5,6,5,5,5,5,6,5,5,5,5,6,5,5,5,5,6),nc=12))

for(i in 1:length(K.clusters)){
	par(mar=c(4,4.5,4,1))
	plot(PC[,1],PC[,2],col=KMEANS[[i]]$cluster+3,pch=20,xlab='PC1',ylab='PC2',cex.lab=1.5,cex=0.5,main=paste0("K = ",K.clusters[i]))
	plot.seizure2D(PC,seizures,legendpos='topleft',Legend=c(TRUE,FALSE,FALSE)[i])
	par(mar=c(2,4.5,2,1))
	image(1:nrow(X),1:ncol(X),X*0,axes=FALSE,xlab='',ylab='')
	for(j in 1:K.clusters[i]){
		abline(v=which(KMEANS[[i]]$cluster==j),col=j+3,lwd=.5)
	}
	mtext('time',side=1,outer=FALSE,cex=1,line=.5)
}
dev.off()

set.seed(123)
GAP = clustering(X, K=10, rescale=RESCALE, iter.max=100, nstart=50, gap.boot=500)




######### FIGURE 6: HIERARCHICAL CLUSTERING #########
DISS = dist(X)
HClust.complete = clustering(X,K=3,method='hierarchical',linkage='complete',Plot=FALSE)
HClust = HClust.complete$full
cuts = numeric()
n.clusters = numeric()
`%!in%` = Negate(`%in%`)
for(h in seq(0.001,max(DISS),len=10000)){
	n.clust = length(unique(cutree(HClust,h=h)))
	if(n.clust %!in% n.clusters){
		cuts = c(cuts,h)
		n.clusters = c(n.clusters,n.clust)
	}
}

pdf(paste0(FOLDER,"figures/FIG6-HCLUSTcomplete.pdf"),height=8,width=10)
layout(matrix(c(1,1,1,1,2,2,2,5,1,1,1,1,3,3,3,6,1,1,1,1,4,4,4,7),nc=3))
par(mar=c(4,4,2,1),oma=c(0,1,2,0))
hs = cuts[which(n.clusters %in% c(20,10,3))]+2
plot(HClust,main="Complete linkage", xlab="", sub="",cex=1,labels=FALSE,lwd=.2,cex.main=1.5,cex.lab=1.5)
abline(h=hs,col='red')
ll=1
for(h in hs){
 CS = cutree(HClust,h=h)
 par(mar=c(4,4.5,4,1))
	plot(PC[,1],PC[,2],col=CS+3,xlab='PC1',ylab='PC2',main=paste0('Height=',round(h,2),' (',max(CS),' clusters)'),cex.lab=1.5,cex=0.5,pch=20)
 plot.seizure2D(PC,seizures,legendpos='topleft',Legend=c(TRUE,FALSE,FALSE)[ll])
 ll=ll+1
}
for(h in hs){
	CS = cutree(HClust,h=h)
	par(mar=c(2,4.5,2,1))
	image(1:nrow(X),1:ncol(X),X*0,axes=FALSE,xlab='',ylab='')
	#TIME = seq(START,END,600)
	#TIME.BINS = TIME/bin.secs
	#TIME.LABELS = TIME/60
	#axis(1,at=TIME.BINS,labels=TIME.LABELS)
	for(j in 1:K.clusters[i]){
		abline(v=which(CS==j),col=j+3,lwd=.5)
	}
	mtext('time',side=1,outer=FALSE,cex=1,line=.5)
}
mtext('A',side=3,outer=TRUE,at=0.01,cex=1.5)
mtext('B',side=3,outer=TRUE,at=0.01,cex=1.5,line=-33)
dev.off()




######### FIGURE 7: KNN MATRIX COMPLETION
X=SURVEY
set.seed(123)
propmiss = 0.1
n = nrow(X); p=ncol(X)
XO = X; XO[sample(1:(n*p),size=round(propmiss*n*p))]=NA
KS = seq(1,80,2); REPS=50; NOTQUANT=1:4
COMPL.KNN = completion(XO,method='KNN',notquant=NOTQUANT,neighbors=KS,reps=REPS,Plot=TRUE)
Z = COMPL.KNN$Z
RISK = COMPL.KNN$risk
K.OPT = COMPL.KNN$neighbors.opt

COL=gray(level=seq(0,1,len=1000))
pdf(paste0(FOLDER,"figures/FIG7-KNNcompletion.pdf"),height=7,width=10)
par(mfrow=c(2,3),mar=c(4,4,4,1))
image(1:n,1:p,X,col=COL,xlab='samples',ylab='features',main='complete data')
image(1:n,1:p,XO,col=COL,xlab='samples',ylab='features',main=paste0(propmiss*100,'% missingness'))
XOmiss = matrix(NA,nc=p,nr=n)
XOmiss[is.na(XO)] = 1
image(1:n,1:p,XOmiss,col='red',add=TRUE)
SDS = sapply(1:length(KS),function(i) sd(COMPL.KNN$LOSS[,i]))/sqrt(REPS)
plot(KS,RISK,type='l',xlab='K',ylab='completion risk',main=paste0('Optimal K = ', K.OPT),lwd=3,ylim=range(c(RISK-2*SDS,RISK+2*SDS)))
polygon(c(KS,KS[length(KS):1]), c(RISK-2*SDS,(RISK+2*SDS)[length(KS):1]),border=NA,col=rgb(0,0,0,.2))
abline(v=K.OPT,lty=2)
image(1:n,1:p,Z,col=COL,xlab='samples',ylab='features',main=paste0('KNN completion (K = ',K.OPT,')'))
CATCORRECT = sapply(NOTQUANT,function(j){
	x=X[is.na(XO[,j]),j]
	z=Z[is.na(XO[,j]),j]
	sum(diag(table(x,z)))/sum(table(x,z))
	})
BARPLOT = barplot(CATCORRECT, xaxt="n",ylab='proportion correct recovery',ylim=c(0,1),border=NA)
labs = paste0('cat. ',NOTQUANT)
text(cex=1, x=BARPLOT, y=-0.15, labs, xpd=TRUE, srt=45)
X.NUM=X[,-NOTQUANT][is.na(XO[,-NOTQUANT])]; Z.NUM=Z[,-NOTQUANT][is.na(XO[,-NOTQUANT])]; LIM=range(c(X.NUM,Z.NUM))
plot(X.NUM,Z.NUM,xlab='missing numerical values',
	ylab='recovered numerical values',pch=16,col=rgb(0,0,0,0.5),xlim=LIM,ylim=LIM)
abline(a=0,b=1)
mtext('A',side=3,outer=TRUE,at=0.01,line=-2,cex=1.1)
mtext('B',side=3,outer=TRUE,at=0.34,line=-2,cex=1.1)
mtext('C',side=3,outer=TRUE,at=0.67,line=-2,cex=1.1)
mtext('D',side=3,outer=TRUE,at=0.01,line=-28,cex=1.1)
mtext('E',side=3,outer=TRUE,at=0.34,line=-28,cex=1.1)
mtext('F',side=3,outer=TRUE,at=0.67,line=-28,cex=1.1)
dev.off()




######### FIGURE 8: LOW-RANK MATRIX COMPLETION
X=EEG$X[1:200,]
set.seed(123)
propmiss = 0.1
n = nrow(X); p=ncol(X)
XO = X; XO[sample(1:(n*p),size=round(propmiss*n*p))]=NA

# LOW-RANK COMPLETION - ITERATIVE PCA
REPS=100; q.seq=10:25; TOL=10^-7; maxit=1000
COMPL.LR = completion(XO,method='lowrank',q=q.seq,max.iters=maxit,tol=TOL,reps=REPS,Plot=TRUE)
Z = COMPL.LR$Z
RISK = COMPL.LR$risk
q.opt = COMPL.LR$q.opt

# LOW-RANK - NUCLEAR NORM MINIMIZATION (FOR COMPARISON)
if (!require(filling)) install.packages('filling')
Z.nuclear <- filling::fill.nuclear(XO)$X

COL=gray(level=seq(0,1,len=1000))
M = is.na(XO)
pdf(paste0(FOLDER,"figures/FIG8-LOWRANKcompletion.pdf"),height=7,width=10)
par(mfrow=c(2,3),mar=c(4,4,4,1))
image(1:n,1:p,X,col=COL,xlab='samples',ylab='time bins',main='complete data')
image(1:n,1:p,XO,col=COL,xlab='samples',ylab='time bins',main=paste0(propmiss*100,'% missingness'))
XOmiss = matrix(NA,nc=p,nr=n)
XOmiss[is.na(XO)] = 1
image(1:n,1:p,XOmiss,col='red',add=TRUE)
SDS = sapply(1:length(q.seq),function(i) sd(COMPL.LR$LOSS[,i]))/sqrt(REPS)
plot(q.seq,RISK,type='l',xlab='q',ylab='completion risk',main=paste0('Optimal q = ', q.opt),lwd=3,ylim=range(c(RISK-2*SDS,RISK+2*SDS)))
polygon(c(q.seq,q.seq[length(q.seq):1]), c(RISK-2*SDS,(RISK+2*SDS)[length(q.seq):1]),border=NA,col=rgb(0,0,0,.2))
abline(v=q.opt,lty=2)
image(1:n,1:p,Z,col=COL,xlab='time bins',ylab='features',main=paste0('PCA Low-Rank completion (q=',q.opt,')'))
LIM = range(c(X[M==1],Z[M==1]))
plot(X[M==1],Z[M==1],xlab='missing values',
	ylab='PCA Low-Rank recovery',pch=19,col=rgb(0,0,0,0.5),xlim=LIM,ylim=LIM)
abline(a=0,b=1)
plot(Z[M==1],Z.nuclear[M==1],xlab='PCA Low-Rank recovery',
	ylab='Nuclearn Norm Low-Rank recovery',pch=19,col=rgb(0,0,0,0.5))
abline(a=0,b=1)
mtext('A',side=3,outer=TRUE,at=0.01,line=-2,cex=1.1)
mtext('B',side=3,outer=TRUE,at=0.34,line=-2,cex=1.1)
mtext('C',side=3,outer=TRUE,at=0.67,line=-2,cex=1.1)
mtext('D',side=3,outer=TRUE,at=0.01,line=-28,cex=1.1)
mtext('E',side=3,outer=TRUE,at=0.34,line=-28,cex=1.1)
mtext('F',side=3,outer=TRUE,at=0.67,line=-28,cex=1.1)
dev.off()

