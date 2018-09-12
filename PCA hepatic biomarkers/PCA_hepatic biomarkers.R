setwd("C:/Users/VDelahaut/Dropbox/PhD/01 Experiments/TSS-2017-EXP-02-GWAS sampling campaign/03 Results/04 Liver samples/R/PCA hepatic biomarkers/")

##load packages
library(ggplot2)
library(gridExtra)


DAT_hepbiomark <- read.csv("BiomarkersLiverRaw_APXDHARMDAexcluded_formultivariatestatistics.csv", header = TRUE)
names(DAT_hepbiomark)


##Rows are individual samples
##Columns are differen variables =  10 biomarkers determined in liver tissue

#subsets for regions
DAT_Hep_Maas <- subset.data.frame(DAT_hepbiomark,Region=="Maas")
names(DAT_Hep_Maas)
DAT_Hep_SE <- subset.data.frame(DAT_hepbiomark,Region=="Scheldt-East")
names(DAT_Hep_SE)
DAT_Hep_SW <- subset.data.frame(DAT_hepbiomark,Region=="Scheldt-West")


##MAAS
pca_HepMaas <- prcomp(~XO+
                     POX+
                     CAT+
                     SOD+
                     GR+
                     GPX+
                     GST+
                     FRAP+
                     Polyphenol+
                     ProteinContent,
                     DAT_Hep_Maas,
                   na.action=na.exclude,
                   scale=TRUE)
summary(pca_HepMaas)

dtHepMaas <- data.frame('River' = DAT_Hep_Maas$Riv, pca_HepMaas$x[,1:2]) 
PCAHepMaas<-ggplot(data = dtHepMaas)+
  geom_point(aes(x = PC1, y = PC2, col = River), size=6)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  ggtitle("Maas")+
  xlab("PC1 (37.1%)")+
  ylab("PC2 (18.26%)")+
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

pca_HepSE <- prcomp(~XO+
                        POX+
                        CAT+
                        SOD+
                        GR+
                        GPX+
                        GST+
                        FRAP+
                        Polyphenol+
                        ProteinContent,
                      DAT_Hep_SE,
                      na.action=na.exclude,
                      scale=TRUE)



summary(pca_HepSE)

dtHepSE <- data.frame('River' = DAT_Hep_SE$Riv, pca_HepSE$x[,1:2]) 
PCAHepSE<-ggplot(data = dtHepSE)+
  geom_point(aes(x = PC1, y = PC2, col = River), size=6)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  ggtitle("Schelde-oost")+
  xlab("PC1 (42.18%)")+
  ylab("PC2 (14.01%)")+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 17),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size=25))


##SW

pca_HepSW <- prcomp(~XO+
                      POX+
                      CAT+
                      SOD+
                      GR+
                      GPX+
                      GST+
                      FRAP+
                      Polyphenol+
                      ProteinContent,
                    DAT_Hep_SW,
                    na.action=na.exclude,
                    scale=TRUE)



summary(pca_HepSW)

dtHepSW <- data.frame('River' = DAT_Hep_SW$Riv, pca_HepSW$x[,1:2]) 
PCAHepSW<-ggplot(data = dtHepSW)+
  geom_point(aes(x = PC1, y = PC2, col = River), size=6)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  ggtitle("Schelde-west")+
  xlab("PC1 (47.06%)")+
  ylab("PC2 (22.94%)")+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 17),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size=25))


