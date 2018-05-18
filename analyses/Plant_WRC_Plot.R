#                                                                              #
#	Long-term fert and mowing at WRC:  Community Charactorization                #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Ariane Peralta (2016/05/13)                                      #
# Modified by:  Mario Muscarella(2016/07/27)                                   #
#                                                                              #
#                                                                              #
################################################################################

#Code for multivariate analyses - AP work in progress - used file WRC_Importance_final.csv b/c had to hand enter 4 samples

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/WRC_FertMowing")
#setwd("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology")
opar <- par(no.readonly = TRUE)  # Saves plot defaults

# Add Summary Functions
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}
ci <- function(x, ...){1.96 * sd(x,na.rm = TRUE)}

# Code Dependencies
library(MASS)
library(nlme)
library(reshape2)
library(vegan)
library(reshape)
library(lme4)
library(ggplot2)
require("png")
require("grid")



# Import Data
PCC <- read.csv ("./data/WRC_Importance_final.csv", header=TRUE)
labels(PCC)
#soil <- read.csv("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology/data_original/WRC_Soil_Data.csv", header=TRUE)
#labels(soil)


treatments <- PCC$treatment
levels(treatments) <- c("UM/UF", "UM/F", "M/UF", "M/F")

###########################
# Simple Hypothesis Testing
###########################

#incorporating strata to restrict permutation within like treatments only (?? - double check with McCoy)
adonis = adonis(PCC[,-c(1:9)] ~ Fertilizer*Mowing*Year+(1|BLOCK/QUADRAT..), strata = PCC$treatment, method = "bray", data = PCC, perm=1000)
adonis

adonis2 = adonis(PCC[,-c(1:9)] ~ Fertilizer*Mowing*Year+(1|BLOCK/QUADRAT..), method = "bray", data = PCC, perm=1000)
adonis2


str(PCC)
WL.simper <- simper(PCC[,-c(1:9)], group = PCC$Fertilizer
simper <- summary(WL.simper)
simper

#SIMPER is somewhat frowned upon - Warton et al. 2012 - so trying to sort out which anlaysis to run to ask "Which plant species are responsible for group differences?")

#adonis(formula = PCC[, -c(1:9)] ~ Fertilizer * Mowing * Year,      data = PCC, permutations = 1000, method = "bray") 

#Permutation: free
#Number of permutations: 1000

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
#Fertilizer               1    10.988 10.9883  56.034 0.04454 0.000999 ***
#Mowing                   1    22.378 22.3785 114.118 0.09070 0.000999 ***
#Year                     1    16.649 16.6492  84.902 0.06748 0.000999 ***
#Fertilizer:Mowing        1     2.778  2.7775  14.164 0.01126 0.000999 ***
#Fertilizer:Year          1     1.345  1.3446   6.857 0.00545 0.000999 ***
#Mowing:Year              1     5.316  5.3162  27.110 0.02155 0.000999 ***
#Fertilizer:Mowing:Year   1     0.589  0.5894   3.006 0.00239 0.003996 ** 
#Residuals              952   186.687  0.1961         0.75664             
#Total                  959   246.731                 1.00000             
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#visualize plant community composition based on PCoA by year x treatment ([centroid of 8 replicate block x 3 quadrats] for each year x treatment)


#new.data <- read.csv("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology/WRC_Importance_final.csv", header=TRUE)


sampleREL.dist <- vegdist(decostand((PCC[,-c(1:9)]),method="log"),method="bray")
WRC_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE)
explainvar1 <- round(WRC_pcoa$eig[1]/sum(WRC_pcoa$eig)*100,1)
explainvar2 <- round(WRC_pcoa$eig[2]/sum(WRC_pcoa$eig)*100,1)
explainvar3 <- round(WRC_pcoa$eig[3]/sum(WRC_pcoa$eig)*100,1)
explainvar1
explainvar2
explainvar3
pcoap <- merge(as.data.frame(WRC_pcoa$points),PCC$treatment, by=0,all.x=T)
rownames(pcoap) <- rownames(WRC_pcoa$points)
pcoap <- merge(pcoap[,-1],PCC$Year, by=0,all.x=T)
rownames(pcoap) <- rownames(WRC_pcoa$points)
treatments <- PCC$treatment
year <- PCC$Year
levels(treatments) <- c("UM/UF", "UM/F", "M/UF", "M/F")
points <- cbind(as.data.frame(WRC_pcoa$points), treatments, year)
L.centroids <- melt(points, id=c("treatments", "year"), measure.vars = c("V1", "V2", "V3"))
centroids <- cast(L.centroids, ... ~ variable, mean)
centroids <- cast(L.centroids, ... ~ variable, fun.aggregate=c(mean,se))

cex.yr <- 1.2 + (as.numeric(unique(centroids$year)) - 2000) * 0.1



#write.csv(centroids, file="newdataClassified_genus_OTUsP.centroids.csv")


####  
  
#png(filename="./figures/Plant_PCoA.png",
#      width = 1800, height = 600, res = 96*2, bg = "white")

pdf(file="./figures/Plant_PCoA.pdf",
    width = 9, height = 3, bg = "white")


par(opar)
layout(matrix(1:4, 1))
  
par(mar=c(0.25,0.5,1,0), oma=c(5,5,1,1))
x.dim <- c((min(centroids$V1_mean)-(max(centroids$V1_mean)*0.15)) ,
           (max(centroids$V1_mean)+(max(centroids$V1_mean)*0.15)))

y.dim <- c((min(centroids$V2_mean)-(max(centroids$V2_mean*0.15))),   
           (max(centroids$V2_mean)+(max(centroids$V2_mean)*0.15)+0.025))
           
           
plot(pcoap$V1, pcoap$V2, xlab="", 
     ylab="", 
     xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",xaxt="n",yaxt="n", 
     cex.lab=1.5, cex.axis=1.2) 
axis(side=1, labels = T, las=1, cex = 0.75)
axis(side=2, las=1, cex = 0.75)
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
box(lwd=2)
arrows(centroids[which(centroids$treatments == "UM/UF"), ]$V1_mean, 
       y1 = centroids[which(centroids$treatments == "UM/UF"), ]$V2_mean - 
         centroids[which(centroids$treatments == "UM/UF"), ]$V2_se, 
       y0 = centroids[which(centroids$treatments == "UM/UF"), ]$V2_mean + 
         centroids[which(centroids$treatments == "UM/UF"), ]$V2_se,
       angle = 90,length=0.05, lwd = 2, code = 3)
arrows(centroids[which(centroids$treatments == "UM/UF"), ]$V2_mean, 
       x1 = centroids[which(centroids$treatments == "UM/UF"), ]$V1_mean - 
         centroids[which(centroids$treatments == "UM/UF"), ]$V1_se, 
       x0 = centroids[which(centroids$treatments == "UM/UF"), ]$V1_mean + 
         centroids[which(centroids$treatments == "UM/UF"), ]$V1_se,
       angle = 90, length=0.05, lwd = 2, code = 3)
points(centroids[which(centroids$treatments == "UM/UF"), ]$V1_mean, 
       centroids[which(centroids$treatments == "UM/UF"), ]$V2_mean, 
       pch=21, cex=cex.yr, col="gray10", bg="gray90")
text(centroids[which(centroids$treatments == "UM/UF"), ]$V1_mean, 
     centroids[which(centroids$treatments == "UM/UF"), ]$V2_mean + 
       centroids[which(centroids$treatments == "UM/UF"), ]$V2_se, 
     labels=centroids[which(centroids$treatments == "UM/UF"), ]$year, 
     pos=3, cex = 0.9, srt = 45, offset = 0.75)
rect(-0.15, 0.145, 0.22, 0.18, col = "white", border = NA)
text(-0.15, 0.165, "Unmowed/Unfertilized", adj = 0)


plot(pcoap$V1, pcoap$V2, xlab="", 
    ylab="", 
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",xaxt="n",yaxt="n", 
    cex.lab=1.5, cex.axis=1.2)
axis(side=1, labels = T, las=1, cex = 0.75)
axis(side=2, labels = F, las=1, cex = 0.75)
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
box(lwd=2)
arrows(centroids[which(centroids$treatments == "UM/F"), ]$V1_mean, 
      y1 = centroids[which(centroids$treatments == "UM/F"), ]$V2_mean - 
        centroids[which(centroids$treatments == "UM/F"), ]$V2_se, 
      y0 = centroids[which(centroids$treatments == "UM/F"), ]$V2_mean + 
        centroids[which(centroids$treatments == "UM/F"), ]$V2_se,
      angle = 90,length=0.05, lwd = 2, code = 3)
arrows(centroids[which(centroids$treatments == "UM/F"), ]$V2_mean, 
      x1 = centroids[which(centroids$treatments == "UM/F"), ]$V1_mean - 
        centroids[which(centroids$treatments == "UM/F"), ]$V1_se, 
      x0 = centroids[which(centroids$treatments == "UM/F"), ]$V1_mean + 
        centroids[which(centroids$treatments == "UM/F"), ]$V1_se,
      angle = 90, length=0.05, lwd = 2, code = 3)
points(centroids[which(centroids$treatments == "UM/F"), ]$V1_mean, 
      centroids[which(centroids$treatments == "UM/F"), ]$V2_mean, 
      pch=21, cex=cex.yr, col="gray10", bg="forestgreen")
text(centroids[which(centroids$treatments == "UM/F"), ]$V1_mean, 
    centroids[which(centroids$treatments == "UM/F"), ]$V2_mean + 
      centroids[which(centroids$treatments == "UM/F"), ]$V2_se, 
    labels=centroids[which(centroids$treatments == "UM/F"), ]$year, 
    pos=3, cex = 0.9, srt = 45, offset = 0.75)
rect(-0.13, 0.14, 0.22, 0.18, col = "white", border = NA)
text(-0.13, 0.165, "Unmowed/Fertilized", adj = 0)




plot(pcoap$V1, pcoap$V2, xlab="", 
    ylab="", 
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",xaxt="n",yaxt="n", 
    cex.lab=1.5, cex.axis=1.2)
axis(side=1, labels = T, las=1, cex = 0.75)
axis(side=2, labels = F, las=1, cex = 0.75)
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
box(lwd=2)
arrows(centroids[which(centroids$treatments == "M/UF"), ]$V1_mean, 
      y1 = centroids[which(centroids$treatments == "M/UF"), ]$V2_mean - 
        centroids[which(centroids$treatments == "M/UF"), ]$V2_se, 
      y0 = centroids[which(centroids$treatments == "M/UF"), ]$V2_mean + 
        centroids[which(centroids$treatments == "M/UF"), ]$V2_se,
      angle = 90,length=0.05, lwd = 2, code = 3)
arrows(centroids[which(centroids$treatments == "M/UF"), ]$V2_mean, 
      x1 = centroids[which(centroids$treatments == "M/UF"), ]$V1_mean - 
        centroids[which(centroids$treatments == "M/UF"), ]$V1_se, 
      x0 = centroids[which(centroids$treatments == "M/UF"), ]$V1_mean + 
        centroids[which(centroids$treatments == "M/UF"), ]$V1_se,
      angle = 90, length=0.05, lwd = 2, code = 3)
points(centroids[which(centroids$treatments == "M/UF"), ]$V1_mean, 
      centroids[which(centroids$treatments == "M/UF"), ]$V2_mean, 
      pch=21, cex=cex.yr, col="gray10", bg="gray90")
text(centroids[which(centroids$treatments == "M/UF"), ]$V1_mean, 
    centroids[which(centroids$treatments == "M/UF"), ]$V2_mean + 
      centroids[which(centroids$treatments == "M/UF"), ]$V2_se, 
    labels=centroids[which(centroids$treatments == "M/UF"), ]$year, 
    pos=3, cex = 0.9, srt = 45, offset = 0.75)
rect(-0.12, 0.14, 0.22, 0.18, col = "white", border = NA)
text(-0.12, 0.165, "Mowed/Unfertilized", adj = 0)
           


plot(pcoap$V1, pcoap$V2, xlab="", 
    ylab="", 
    xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",xaxt="n",yaxt="n", 
    cex.lab=1.5, cex.axis=1.2)
axis(side=1, las=1, cex = 0.75)
axis(side=2, labels = F, las=1, cex = 0.75)
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
box(lwd=2)
arrows(centroids[which(centroids$treatments == "M/F"), ]$V1_mean, 
      y1 = centroids[which(centroids$treatments == "M/F"), ]$V2_mean - 
        centroids[which(centroids$treatments == "M/F"), ]$V2_se, 
      y0 = centroids[which(centroids$treatments == "M/F"), ]$V2_mean + 
        centroids[which(centroids$treatments == "M/F"), ]$V2_se,
      angle = 90,length=0.05, lwd = 2, code = 3)
arrows(centroids[which(centroids$treatments == "M/F"), ]$V2_mean, 
      x1 = centroids[which(centroids$treatments == "M/F"), ]$V1_mean - 
        centroids[which(centroids$treatments == "M/F"), ]$V1_se, 
      x0 = centroids[which(centroids$treatments == "M/F"), ]$V1_mean + 
        centroids[which(centroids$treatments == "M/F"), ]$V1_se,
      angle = 90, length=0.05, lwd = 2, code = 3)
points(centroids[which(centroids$treatments == "M/F"), ]$V1_mean, 
      centroids[which(centroids$treatments == "M/F"), ]$V2_mean, 
      pch=21, cex=cex.yr, col="gray10", bg="forestgreen")
text(centroids[which(centroids$treatments == "M/F"), ]$V1_mean, 
    centroids[which(centroids$treatments == "M/F"), ]$V2_mean + 
      centroids[which(centroids$treatments == "M/F"), ]$V2_se, 
    labels=centroids[which(centroids$treatments == "M/F"), ]$year, 
    pos=3, cex = 0.9, srt = 45, offset = 0.75)
rect(-0.10, 0.14, 0.22, 0.18, col = "white", border = NA)
text(-0.10, 0.165, "Mowed/Fertilized", adj = 0)
           



mtext(paste("PCoA Axis 1 (",explainvar1, "%)", sep=""), side = 1, 
      line = 2.75, outer = T, cex = 1.25)
      
mtext(paste("PCoA Axis 2 (",explainvar2, "%)", sep=""), side = 2, 
      line = 2.75, outer = T, cex = 1.25)

dev.off() # this writes plot to folder
graphics.off() # shuts down open devices

# img <- readPNG("./figures/Plant_PCoA.png")
# grid.raster(img)


#####


# ordiellipse(cbind(pcoap$V1, pcoap$V2), pcoap$y, kind="se", conf=0.95, lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE, cex=2)
# levels(treatments) <- c("UM/UF", "UM/F", "m/UF", "MF")
# myColors <- c("#FFF000", "#CCFF00", "#33CC33", "#339933")
# names(myColors) <- levels(treatments)
