setwd("C:/Users/VDelahaut/Dropbox/PhD/01 Experiments/TSS-2017-EXP-02-GWAS sampling campaign/03 Results/05 Muscle samples/R/PCAs/")

##load packages
library(ggplot2)
library(gridExtra)




DAT_Metal <- read.csv("GWAS sampling campaign - Muscle samples_ResultsMetalAnalysis_NOOUTLIERS.csv", header = TRUE)
names(DAT_Metal)


##Rows are individual samples
##Columns are differen variables =  13 metals measured in muscle tissue

#subsets for regions
DAT_Maas <- subset.data.frame(DAT_Metal,Region=="Maas")
DAT_SE <- subset.data.frame(DAT_Metal,Region=="Scheldt-EastBank")
DAT_SW <- subset.data.frame(DAT_Metal,Region=="Scheldt-WestBank")

##MAAS
pca_Maas <- prcomp(~ng.Hg.g.DW+
                     ng.Cd.g.DW+
                     ng.Pb.g.DW+
                     ng.Al.g.DW+
                     ng.Cr.g.DW+
                     ng.Mn.g.DW+
                     ng.Fe.g.DW+
                     ng.Co.g.DW+
                     ng.Ni.g.DW+
                     ng.Cu.g.DW+
                     ng.Zn.g.DW+
                     ng.As.g.DW+
                     ng.Se.g.DW,
                   DAT_Maas,
                   na.action=na.exclude,
                   scale=TRUE)
summary(pca_Maas)

dtMaas <- data.frame('River' = DAT_Maas$Riv, pca_Maas$x[,1:2]) 
PCAMetalMaas<-ggplot(data = dtMaas)+
  geom_point(aes(x = PC1, y = PC2, col = River), size=6)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  ggtitle("Maas")+
  xlab("PC1 (44.02%)")+
  ylab("PC2 (11.51%)")+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 17),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size=25))




##SE
pca_SE <- prcomp(~ng.Hg.g.DW+
                     ng.Cd.g.DW+
                     ng.Pb.g.DW+
                     ng.Al.g.DW+
                     ng.Cr.g.DW+
                     ng.Mn.g.DW+
                     ng.Fe.g.DW+
                     ng.Co.g.DW+
                     ng.Ni.g.DW+
                     ng.Cu.g.DW+
                     ng.Zn.g.DW+
                     ng.As.g.DW+
                     ng.Se.g.DW,
                   DAT_SE,
                   na.action=na.exclude,
                   scale=TRUE)
summary(pca_SE)

dtSE <- data.frame('River' = DAT_SE$Riv, pca_SE$x[,1:2]) 
PCAMetalSE <- ggplot(data = dtSE)+
  geom_point(aes(x = PC1, y = PC2, col = River), size=6)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  ggtitle("Schelde East")+
  xlab("PC1 (34.69%)")+
  ylab("PC2 (15.94%)")+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 17),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size=25))



##SW
pca_SW <- prcomp(~ng.Hg.g.DW+
                   ng.Cd.g.DW+
                   ng.Pb.g.DW+
                   ng.Al.g.DW+
                   ng.Cr.g.DW+
                   ng.Mn.g.DW+
                   ng.Fe.g.DW+
                   ng.Co.g.DW+
                   ng.Ni.g.DW+
                   ng.Cu.g.DW+
                   ng.Zn.g.DW+
                   ng.As.g.DW+
                   ng.Se.g.DW,
                 DAT_SW,
                 na.action=na.exclude,
                 scale=TRUE)
summary(pca_SW)

dtSW <- data.frame('River' = DAT_SW$Riv, pca_SE$x[,1:2]) 
PCAMetalSW<-ggplot(data = dtSW)+
  geom_point(aes(x = PC1, y = PC2, col = River), size=8)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  ggtitle("Schelde West")+
  xlab("PC1 (45.50%)")+
  ylab("PC2 (11.53%)")+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 17),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size=25))



tiff(file="C:/Users/VDelahaut/Dropbox/PhD/01 Experiments/TSS-2017-EXP-02-GWAS sampling campaign/03 Results/05 Muscle samples/R/PCAs/PCAsmetal.tiff", 
     res=300, width = 8000, height = 3000, units='px')
grid.arrange(PCAMetalMaas, PCAMetalSE,PCAMetalSW , ncol =3, nrow = 1)
dev.off()


