library(ggplot2)
library(gplots)    # heatmap.2
library(vegan)
library(car)
library(MVN)
library(MASS)
rm(list=ls())

#######   Read in data  - - - - - - - - - - - - - - - -
data2 <- read.csv("Hoadley2019-data.csv",header=TRUE,sep=",")

########MDS plots for Figure 3 and 4 (panels A-D) - - - - - - - - - - -
########
perc<-data2
trenchii<-subset(perc, X4=="trenchii" & X2=="Nikko")
C40<-subset(perc, X4=="C40-sp" & X2=="Offshore")
C40goni<-subset(C40, X3=="Goniastrea" & X2=="Offshore")
C40pachy<-subset(C40, X3=="Pachyseris" & X2=="Offshore")
C21<-subset(perc, X4=="C21-sp" & X2=="Offshore")
C3u<-subset(perc, X4=="C3u-sp" & X2=="Offshore")
offshore<-rbind(C40goni, C40pachy, C21, C3u)
inshore<-trenchii

dataA<-list(inshore, offshore)
colors<-list(c("#fb9a99", "#e31a1c"),c("#a6cee3", "#1f78b4"))
namerss<-c("Inshore", "Offshore")

for(L in 1:2)
{
	data<-dataA[[L]]
	period<-cbind(data[,1:9], τPSII=data$P1tauAve, τPQ=data$P1MTtauAve, ETR=data$P1STETR, NPQ=data$P1mtNPQ, FvFm=data$P1MTfvfm, P=data$P1STp, σPSIIsig=data$P1STsigma, Lipids=data$P1symlip, Proteins=data$P1symprot, Carbs=data$P1symcarb, Volume=data$P1volume, Photo=data$P1grossPhotocell, Density=data$P1density, Chla=data$P1chlacell)
	
	perc<-period
	percs<-na.omit(perc)
	dim(percs)
	summary(percs)
	
	Acro1<-percs
	har<-t(Acro1[,10:23])
	colnames(har)<-Acro1$local
	stat<-har
	rs<-rownames(stat)
	cls<-colnames(stat)
	Boo<-dim(stat)[[1]]
	datum <- list()
	atum <- list()
	for (i in 1:Boo)
	{		H<-rownames(stat)[i]
			datum[[H]]<-((as.numeric(stat[i,])/range(as.numeric(stat[i,]), na.rm=TRUE)[2]))
			atum[[H]]<-t.test(stat[i,]~Acro1$X1)[[3]]
	}
	allerie2<-data.matrix(do.call(rbind, datum))
	ler2<-data.matrix(do.call(rbind, atum))
	print(ler2)
	colnames(allerie2)<-cls
	
	acfac<-Acro1[,1:9]
	acdata.mds<-metaMDS(log(t(allerie2)+1), distance='euclidean', k=2, trymax=40, autotransform=FALSE)
	acdata.mds2<-metaMDS(log(t(allerie2)+1), previous.best = acdata.mds, distance='euclidean', k=2, trymax=40, autotransform=FALSE)
	acord <- acdata.mds2$points
	Acroneword <- cbind(acfac, acord)
	acsym <- c(1,19)
	acroda<-data.frame(t(allerie2))
	acrovec<-envfit(acdata.mds2$points, acroda, perm = 9999, na.rm=TRUE)
	scores(acrovec, "vectors")
	acrovecT<-acrovec
	Acro<-Acroneword
	
	suber1<-subset(Acro, X3=="Goniastrea")
	suber2<-subset(Acro, X3=="Acropora")
	suber3<-subset(Acro, X3=="Cyphastrea")
	suber4<-subset(Acro, X3=="Pachyseris")
	suber5<-list(suber2, suber3, suber1, suber4)
	
	bal<-paste("Figures/4coralMDSplot", namerss[L], ".pdf", sep="")
	pdf(file = bal, width = 7, height = 5, bg="transparent")
	
	par(mfrow=c(2,2), mar = c(0.1, 0.1, 0.1, 0.1), oma = c(2, 4.25, 0.1, 0.1))
	for(I in 1:4)
	{	
		suber<-suber5[[I]]
		plot(Acro[,10], Acro[,11], col = c("grey", "grey")[as.numeric(Acro$X1)], bg = c("grey", "grey")[as.numeric(Acro$X1)], pch = c(24,21,22,23)[as.numeric(Acro$X3)], cex = 1.5, xlab='', ylab='', ylim=c(-0.6,0.6), xlim=c(-0.8, 0.8),axes=FALSE)
		par(new=TRUE)
		coll<-colors[[L]]
		plot(suber[,10], suber[,11], col = c("black", "black")[as.numeric(suber$X1)], bg = coll[as.numeric(suber$X1)], pch = c(24,21,22,23)[as.numeric(suber$X3)], cex = 2.25, xlab='', ylab='', ylim=c(-0.6,0.6), xlim=c(-0.8, 0.8),axes=FALSE)
		if(I==1)
		{
			par(new=TRUE)
			plot(acrovecT, col="black", p.max=0.001, cex=0.000001)
		}
		if(I==1)
		{
			axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
			#axis(1, at=, col.axis="black", las=1, cex.axis=2)
		}
		if(I==3)
		{
			axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
			axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
		}
		if(I==4)
		{
			#axis(2, at=, col.axis="black", las=2, cex.axis=2)
			axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
		}
		grid(col = "grey")
		box(col = "black")
	}
	dev.off()
}




##########heatmaps for panels e in figure 3 and 4 - - - - - - -
##########

perc<-data2
trenchii<-subset(perc, X4=="trenchii" & X2=="Nikko")
C40<-subset(perc, X4=="C40-sp" & X2=="Offshore")
C40goni<-subset(C40, X3=="Goniastrea" & X2=="Offshore")
C40pachy<-subset(C40, X3=="Pachyseris" & X2=="Offshore")
C21<-subset(perc, X4=="C21-sp" & X2=="Offshore")
C3u<-subset(perc, X4=="C3u-sp" & X2=="Offshore")
TAcro<-subset(trenchii, X3=="Acropora")
TGonia<-subset(trenchii, X3=="Goniastrea")
TCyph<-subset(trenchii, X3=="Cyphastrea")
TPachy<-subset(trenchii, X3=="Pachyseris")

dataA<-list(TAcro, TCyph, TGonia, TPachy, C21, C3u, C40goni, C40pachy)
names<-c("TAcro", "TCyph", "TGonia", "TPachy", "C21", "C3u", "C40goni", "C40pachy")
ler2<-list()
for(TG in 1:8)
{
	data<-dataA[[TG]]
	period<-cbind(data[,1:9], FvFm=data$P1MTfvfm, ETR=data$P1STETR, NPQ=data$P1mtNPQ, σPSIIsig=data$P1STsigma, P=data$P1STp, τPSII=data$P1tauAve, τPQ=data$P1MTtauAve, Chla=data$P1chlacell, Photo=data$P1grossPhotocell, Proteins=data$P1symprot, Lipids=data$P1symlip, Carbs=data$P1symcarb, Volume=data$P1volume, Density=data$P1density)
	
	perc<-period
	percs<-na.omit(perc)
	dim(percs)
	summary(percs)
	Acro1<-percs
	har<-t(Acro1[,10:23])
	colnames(har)<-Acro1$local
	stat<-har
	rs<-rownames(stat)
	cls<-colnames(stat)
	Boo<-dim(stat)[[1]]
	datum <- list()
	atum <- list()
	for (i in 1:Boo)
	{		
			if (as.matrix((shapiro.test(stat[i,]))[[2]])>0.05)
			{
				H<-rownames(stat)[i]
				atum[[H]]<-signif(t.test(stat[i,]~Acro1$X1)[[3]],2)
			}
			else if (as.matrix((shapiro.test(log(stat[i,]+1)))[[2]])>0.05)
			{
				H<-rownames(stat)[i]
				atum[[H]]<-signif(t.test(log(stat[i,]+1)~Acro1$X1)[[3]],2)
			}
			else if (as.matrix((shapiro.test(sqrt(stat[i,]+1)))[[2]])>0.05)
			{
				H<-rownames(stat)[i]
				atum[[H]]<-signif(t.test(sqrt(stat[i,]+1)~Acro1$X1)[[3]],2)
			}
			else
			{	
				H<-rownames(stat)[i]
				atum[[H]]<-signif(wilcox.test(log(stat[i,]+1)~factor(Acro1$X1), p.adjust="none", paired=FALSE)[[3]],2)
			}
	}
	ler2[[names[TG]]]<-(data.matrix(do.call(rbind, atum)))
}




offshore<-rbind(C40goni, C40pachy, C21, C3u)
inshore<-trenchii

dataB<-list(inshore, offshore)
data<-dataB[[1]]
period<-cbind(data[,1:9], FvFm=data$P1MTfvfm, ETR=data$P1STETR, NPQ=data$P1mtNPQ, σPSIIsig=data$P1STsigma, P=data$P1STp, τPSII=data$P1tauAve, τPQ=data$P1MTtauAve, Chla=data$P1chlacell, Photo=data$P1grossPhotocell, Proteins=data$P1symprot, Lipids=data$P1symlip, Carbs=data$P1symcarb, Volume=data$P1volume, Density=data$P1density)
perc<-period
percs<-na.omit(perc)
Nikki<-aggregate(percs[,10:23], by=list(Coral=percs$X3, Temp=percs$X1), FUN=mean)
Nikk<-as.matrix(Nikki[,3:16])
allerie2<-log2(Nikk[5:8,1:14]/Nikk[1:4,1:14])

for(L in 1:4)
{
	ler<-data.frame(ler2[L])
	for(I in 1:nrow(ler))
	{
		if(ler[I,1]>0.05)
		{
			allerie2[L,I]<-NA
		}
	}
}

pdf(file = "Figures/heatmap-trenchii.pdf", width = 9.5, height = 3.5, bg="transparent")
col_breaks = c(seq(-2,-0.3,length.out=115), seq(-0.2,0.2,length.out=275), seq(0.3,2,length.out=115))
my_palette <- colorRampPalette(c("#d73027","black","#2166ac"))(length(col_breaks)-1)
par(oma=c(0.1,0.1,0.1,0.1))
heatmap.2(allerie2, col=my_palette, breaks=col_breaks, margins = c(1,1), density.info="none", trace="none", dendrogram = "none", symm=FALSE, symkey=FALSE, symbreaks=FALSE, sepwidth=F,na.color = "grey", Colv=F, Rowv=F, scale="none", key=F, key.xlab="Z-score", key.title=NA, cexRow=0.000002, offsetRow=0, srtRow=0, cexCol=0.000000075, srtCol=90, lmat=rbind(c(4,3), c(2,1)), lhei=c(0.1, 3.5), lwid=c(0.5, 8.5))
dev.off()



dataB<-list(inshore, offshore)
data<-dataB[[2]]
period<-cbind(data[,1:9], FvFm=data$P1MTfvfm, ETR=data$P1STETR, NPQ=data$P1mtNPQ, σPSIIsig=data$P1STsigma, P=data$P1STp, τPSII=data$P1tauAve, τPQ=data$P1MTtauAve, Chla=data$P1chlacell, Photo=data$P1grossPhotocell, Proteins=data$P1symprot, Lipids=data$P1symlip, Carbs=data$P1symcarb, Volume=data$P1volume, Density=data$P1density)
perc<-period
percs<-na.omit(perc)
Nikki<-aggregate(percs[,10:23], by=list(Coral=percs$X3, Temp=percs$X1), FUN=mean)
Nikk<-as.matrix(Nikki[,3:16])
allerie2<-log2(Nikk[5:8,1:14]/Nikk[1:4,1:14])

for(L in 1:4)
{
	ler<-data.frame(ler2[L+4])
	for(I in 1:nrow(ler))
	{
		if(ler[I,1]>0.05)
		{
			allerie2[L,I]<-NA
		}
	}
}

pdf(file = "Figures/heatmap-offshore.pdf", width = 9.5, height = 3.5, bg="transparent")
col_breaks = c(seq(-2,-0.3,length.out=115), seq(-0.2,0.2,length.out=275), seq(0.3,2,length.out=115))
my_palette <- colorRampPalette(c("#d73027","black","#2166ac"))(length(col_breaks)-1)
par(oma=c(0.1,0.1,0.1,0.1))
heatmap.2(allerie2, col=my_palette, breaks=col_breaks, margins = c(1,1), density.info="none", trace="none", dendrogram = "none", symm=FALSE, symkey=FALSE, symbreaks=FALSE, sepwidth=F,na.color = "grey", Colv=F, Rowv=F, scale="none", key=F, key.xlab="Z-score", key.title=NA, cexRow=0.000002, offsetRow=0, srtRow=0, cexCol=0.000000075, srtCol=90, lmat=rbind(c(4,3), c(2,1)), lhei=c(0.1, 3.5), lwid=c(0.5, 8.5))
dev.off()






#########Regression plot in figure 5 - - - - - 



D10<-subset(data2, X5=="10")
offshore<-subset(D10, X2=="Offshore")
nikko<-subset(D10, X2=="Nikko")

samp<-paste("Figures/volume-sigma", ".pdf", sep="")
pdf(file = samp, width = 6, height = 6, bg="transparent")

par(mfrow=c(1,2), mar = c(2.0, 2.0, 0.1, 0.1), oma =c(0.1, 0.1, 0.1, 0.1))

plot(offshore$P1volume, offshore$P1STsigma, col = c("black", "black")[as.numeric(offshore$X1)], bg = c("#a6cee3", "#1f78b4")[as.numeric(offshore$X1)], pch = c(24,21,22,23)[as.numeric(offshore$X3)], cex = 1.25, ylim=c(150, 325), xlim=c(200,1700), xlab=c(""), ylab=c(""))
grid(col = "grey")

alli.mod1 = lm(offshore$P1STsigma ~ offshore$P1volume)
summary(alli.mod1)
alli.mod1[[1]]
abline(138.87491103, 0.08768918, col = "grey")

plot(nikko$P1volume, nikko$P1STsigma, col = c("black", "black")[as.numeric(nikko$X1)], bg = c("#fb9a99", "#e31a1c")[as.numeric(nikko$X1)], pch = c(24,21,22,23)[as.numeric(nikko$X3)], cex = 1.25, ylim=c(150, 325), xlim=c(200,1700), xlab=c(""), ylab=c(""))
grid(col = "grey")

dev.off()
