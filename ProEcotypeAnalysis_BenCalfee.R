#######################################################################################
#                            Loading in required packages                             #
#######################################################################################
library(tidyverse)
library(data.table)
library(ape)
library(vegan)
library(ggvegan)
library(psych)
library(missForest)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(car)

#######################################################################################
#                            Read in csv files                                        #
#######################################################################################
#Data used for CCA
ALL<- fread("data/olddata/POWOW_2_and _3_Data_for _TJ.csv")
ALL <-ALL[-41,]
ALL <- ALL %>% column_to_rownames(., var = "V4")

#######################################################################################
#                  Seperate out required data and modify                              #
#######################################################################################
#Seperate out all species data
spec <- ALL[,c(2,10:13)] 

#Season specific species
summerspec <- spec[spec$Season == "Summer",] %>% .[-18,]
summerspec <-summerspec[!(row.names(summerspec) %in% c("14I","31G", "38I", "38G", "38E", "38B")),]
winterspec <- spec[spec$Season == "Winter",]
winterspec <-winterspec[!(row.names(winterspec) %in% c("19A", "19B", "19C", "19D", "19E", "19F", "19G", "19H", "19I", "41H", "41B", "41G", "41E", "10I", "27C", "37A")),]

#Hellinger transfermations Season specific species
summerspec.hell <- decostand(summerspec[,-1], "hellinger")
winterspec.hell <- decostand(winterspec[,-1], "hellinger")

#Environmental data
environment <- ALL %>% subset(., select = c("Season", "Temp. (C)", "Salinity (", "Density", "Lat.", "NH4 (nM)", "PO4 uM", "Depth (Assumed)", "chlA total", 
                                            "light (diffuse attenuation coefficient, kd)", "Syn FCM", "Peuk FCM"))
# Renaming colums for smoother down stream analysis
names(environment)[names(environment) == "Temp. (C)"] <- "Temp"
names(environment)[names(environment) == "Salinity ("] <- "Salinity"
names(environment)[names(environment) == "NH4 (nM)"] <- "NH4"
names(environment)[names(environment) == "Depth (Assumed)"] <- "Depth"
names(environment)[names(environment) == "chlA total"] <- "chlA"
names(environment)[names(environment) == "light (diffuse attenuation coefficient, kd)"] <- "Light"
names(environment)[names(environment) == "PO4 uM"] <- "PO4"
names(environment)[names(environment) == "Peuk FCM"] <- "Peuk_FCM"
names(environment)[names(environment) == "Syn FCM"] <- "Syn_FCM"

#Season specific environmental data
summerenvironment <- environment[environment$Season == "Summer",]
summerenvironment2 <- summerenvironment[!(row.names(summerenvironment) %in% c("14I","31G", "38I", "38G", "38E", "38B")),]
winterenvironment <- environment[environment$Season == "Winter",]
winterenvironment2 <- winterenvironment[!(row.names(winterenvironment) %in% c("19A", "19B", "19C", "19D", "19E", "19F", "19G", "19H", "19I", "41H", "41B", "41G", "41E", "10I", "27C", "37A")),]
which(winterenvironment2$Depth == 155)
winterenvironment2["22B", "Depth"] <- 150

#Z-score transfermation of enviromental data
SumEvocols <- colnames(summerenvironment[,-1])
SumEvo.Z <- summerenvironment
SumEvo.Z[SumEvocols] <- scale(SumEvo.Z[SumEvocols])

SumEvo.Z <- SumEvo.Z[-18,]
SumEvo.Z <-SumEvo.Z[!(row.names(SumEvo.Z) %in% c("14I","31G", "38I", "38G", "38E", "38B")),]

WinEvocols <- colnames(winterenvironment[,-1])
WinEvo.Z <- winterenvironment
WinEvo.Z[WinEvocols] <- scale(WinEvo.Z[WinEvocols])
unique(row.names(which(is.na(WinEvo.Z), arr.ind=TRUE)))
WinEvo.Z <-WinEvo.Z[!(row.names(WinEvo.Z) %in% c("19A", "19B", "19C", "19D", "19E", "19F", "19G", "19H", "19I", "41H", "41B", "41G", "41E", "10I", "27C", "37A")),]



#Plot vars against vars to look for co-correlation
pairs.panels(SumEvo.Z[,-c(1)], scale=T)
pairs.panels(WinEvo.Z[,-c(1)], scale=T)
#Temp and salinity co-cor so Temp is proxy for salinity

#######################################################################################
#                            Now to run the CCA models                                #
#######################################################################################
#Summer Mod cca
SumCCA <- cca(summerspec.hell ~ Lat. + Temp + Density + Light +   NH4 + Depth + Peuk_FCM + Syn_FCM, data =SumEvo.Z)
vif.cca(SumCCA)
#Temp neg co-cor with density. so we can use temp as a neg proxy for density. Temp is also co-cor with Lat. So we will use it as a proxy for lat
anova(SumCCA, permutations = how(nperm = 9999))
anova(SumCCA, permutations = how(nperm = 9999),by="margin")
ordistep(SumCCA)

FinalCCABen <- cca(formula = summerspec.hell ~ Lat. + Temp + Light + NH4 + Depth + Peuk_FCM, data = SumEvo.Z)
anova(FinalCCABen, permutations = how(nperm = 9999),by="margin")

#Winter
WinCCA <- cca(winterspec.hell ~ Lat. + Temp + Density + Light +  NH4 + Depth + Peuk_FCM + Syn_FCM, data =WinEvo.Z)
vif.cca(WinCCA)
#Temp is co-corr with Density and neg co-corr with lat. So it will be a proxy for both
WinCCA2 <- cca(winterspec.hell ~  Temp +  Light +  NH4 + Depth + Peuk_FCM + Syn_FCM, data =WinEvo.Z)
vif.cca(WinCCA2)

anova(WinCCA2, permutations = how(nperm = 9999))
anova(WinCCA2, permutations = how(nperm = 9999),by="margin")
ordistep(WinCCA2)
WinterFinalcca <- cca(formula = winterspec.hell ~ Temp + Light + Depth + Peuk_FCM + Syn_FCM, data = WinEvo.Z)
summary(WinterFinalcca)

depthcols <- brewer.pal(n = 9, name = "Spectral")
#######################################################################################
#                   Create dataframes needed for plots                                #
#######################################################################################
#Winter df for plot: Species
WinCCASpeciesDF_SpecScal <- data.frame(scores(WinterFinalcca, scaling = "species")$species) %>% rownames_to_column(., var = "ecotype")
WinCCASiteDF_SpecScal <- data.frame(scores(WinterFinalcca, scaling = "species")$sites) %>% rownames_to_column(., var = "sites")
WinCCASiteDF_SpecScal$Depth <- as.factor(winterenvironment2$Depth)
names(depthcols) <- unique(as.factor(winterenvironment2$Depth))
depthcolsdf <- data.frame(depthcols) %>% rownames_to_column(., var = "Depth")
WinCCASiteDF_SpecScal <- WinCCASiteDF_SpecScal %>% inner_join(., depthcolsdf, by = "Depth")

#Winter df for plot: Sites
WinCCASpeciesDF_SiteScal <- data.frame(scores(WinterFinalcca, scaling = "sites")$species) %>% rownames_to_column(., var = "ecotype")
WinCCASiteDF_SiteScal <- data.frame(scores(WinterFinalcca, scaling = "sites")$sites) %>% rownames_to_column(., var = "sites")
WinCCASiteDF_SiteScal$Depth <- as.factor(winterenvironment2$Depth)
names(depthcols) <- unique(as.factor(winterenvironment2$Depth))
depthcolsdf <- data.frame(depthcols) %>% rownames_to_column(., var = "Depth")
WinCCASiteDF_SiteScal <- WinCCASiteDF_SiteScal %>% inner_join(., depthcolsdf, by = "Depth")

summerenvironment3 <-summerenvironment2[-18,]
#Summer df for plot: Species
SumCCASpeciesDF_SpecScal <- data.frame(scores(FinalCCABen, scaling = "species")$species) %>% rownames_to_column(., var = "ecotype")
SumCCASiteDF_SpecScal <- data.frame(scores(FinalCCABen, scaling = "species")$sites) %>% rownames_to_column(., var = "sites")
SumCCASiteDF_SpecScal$Depth <- as.factor(summerenvironment3$Depth)
names(depthcols) <- unique(as.factor(summerenvironment3$Depth))
depthcolsdf <- data.frame(depthcols) %>% rownames_to_column(., var = "Depth")
SumCCASiteDF_SpecScal <- SumCCASiteDF_SpecScal %>% inner_join(., depthcolsdf, by = "Depth")

#Summer df for plot: Sites
SumCCASpeciesDF_SiteScal <- data.frame(scores(FinalCCABen, scaling = "sites")$species) %>% rownames_to_column(., var = "ecotype")
SumCCASiteDF_SiteScal <- data.frame(scores(FinalCCABen, scaling = "sites")$sites) %>% rownames_to_column(., var = "sites")
SumCCASiteDF_SiteScal$Depth <- as.factor(summerenvironment3$Depth)
names(depthcols) <- unique(as.factor(summerenvironment3$Depth))
depthcolsdf <- data.frame(depthcols) %>% rownames_to_column(., var = "Depth")
SumCCASiteDF_SiteScal <- SumCCASiteDF_SiteScal %>% inner_join(., depthcolsdf, by = "Depth")


#######################################################################################
#                                All CCA Plots                                        #
#######################################################################################
WinSpePlot <- autoplot(WinterFinalcca, scaling = "species") + theme_bw() + 
  ggtitle("Winter: Ecotypes Against Variables") +
  geom_vline(xintercept=c(0), linetype="dashed", color = "gray") +
  geom_hline(yintercept=c(0), linetype="dashed", color = "gray") +
  geom_point(data = WinCCASpeciesDF_SpecScal, mapping = aes(x = CCA1, y = CCA2),
             size = 6, shape=18, color = c("blue", "green", "darkgreen", "purple")) +
  geom_text(data = WinCCASpeciesDF_SpecScal, mapping = aes(x = CCA1, y = CCA2, label = ecotype), hjust=.5, vjust=-1) +
  geom_point(data = WinCCASiteDF_SpecScal, mapping = aes(x = CCA1, y = CCA2),
             size = 3, shape=19, color = WinCCASiteDF_SpecScal$depthcols)+
  theme(legend.position = "none") 
                      
WinSitePlot <- autoplot(WinterFinalcca, scaling = "sites") + theme_bw()+ 
  ggtitle("Winter: Sites Against Variables") +
  geom_vline(xintercept=c(0), linetype="dashed", color = "gray") +
  geom_hline(yintercept=c(0), linetype="dashed", color = "gray") +
  geom_point(data = WinCCASpeciesDF_SiteScal, mapping = aes(x = CCA1, y = CCA2),
             size = 6, shape=18, color = c("blue", "green", "darkgreen", "purple")) +
  geom_text(data = WinCCASpeciesDF_SiteScal, mapping = aes(x = CCA1, y = CCA2, label = ecotype), hjust=.5, vjust=-1) +
  geom_point(data = WinCCASiteDF_SiteScal, mapping = aes(x = CCA1, y = CCA2),
             size = 3, shape=19, color = WinCCASiteDF_SiteScal$depthcols) +
  theme(legend.position = "none")                        
                         
SumSpePlot <- autoplot(FinalCCABen, scaling = "species") + theme_bw()+ 
  ggtitle("Summer: Ecotypes Against Variables")+
  geom_vline(xintercept=c(0), linetype="dashed", color = "gray") +
  geom_hline(yintercept=c(0), linetype="dashed", color = "gray") +
  geom_point(data = SumCCASpeciesDF_SpecScal, mapping = aes(x = CCA1, y = CCA2),
             size = 6, shape=18, color = c("blue", "green", "darkgreen", "purple")) +
  geom_text(data = SumCCASpeciesDF_SpecScal, mapping = aes(x = CCA1, y = CCA2, label = ecotype), hjust=.5, vjust=-1) +
  geom_point(data = SumCCASiteDF_SpecScal, mapping = aes(x = CCA1, y = CCA2),
             size = 3, shape=19, color = SumCCASiteDF_SpecScal$depthcols)+
  theme(legend.position = "none")                       
                        
                         
SumSitePlot <- autoplot(FinalCCABen, scaling = "sites") + theme_bw() + 
  ggtitle("Summer: Sites Against Variables") +
  geom_vline(xintercept=c(0), linetype="dashed", color = "gray") +
  geom_hline(yintercept=c(0), linetype="dashed", color = "gray") +
  geom_point(data = SumCCASpeciesDF_SiteScal, mapping = aes(x = CCA1, y = CCA2),
             size = 6, shape=18, color = c("blue", "green", "darkgreen", "purple")) +
  geom_text(data = SumCCASpeciesDF_SiteScal, mapping = aes(x = CCA1, y = CCA2, label = ecotype), hjust=.5, vjust=-1) +
  geom_point(data = SumCCASiteDF_SiteScal, mapping = aes(x = CCA1, y = CCA2),
             size = 3, shape=19, color = SumCCASiteDF_SiteScal$depthcols) +
  theme(legend.position = "none")                         
                                                  
SumCCASiteDF_SiteScal$Depth <- factor(SumCCASiteDF_SiteScal$Depth, levels = c(5, 15, 25, 50, 75, 100, 125, 150, 200))
BensSpecies <-ggplot(SumCCASiteDF_SiteScal, aes(CCA1, CCA2, col = Depth)) +      # ggplot with legend
  geom_point() + scale_color_manual(values =depthcols) + theme_bw()
BensSpecies_alone <- get_legend(BensSpecies)
BensLegend<- ggdraw(BensSpecies_alone)

# ggsave("legend.pdf", plot =BensLegend)
# ggsave("WinSpePlot.pdf", plot =WinSpePlot, width=10, height=10, units = "in")
# ggsave("WinSitePlot.pdf", plot =WinSitePlot, width=10, height=10, units = "in")
# ggsave("SumSpePlot.pdf", plot =SumSpePlot, width=10, height=10, units = "in")
# ggsave("SumSitePlot.pdf", plot =SumSitePlot, width=10, height=10, units = "in")


















#######################################################################################
#  Construct linear models to look for correlations between Ecotypes and 16S phylum   #
#######################################################################################
top100phyvEco <- read_csv("data/Top100Phyla_vs_Eco.csv") %>% column_to_rownames(., var = "X1")
top100phyvEco.hell <- decostand(top100phyvEco, "hellinger")

top100phyvSum <- read.csv("data/NewSumWinData/POWOW_RawCounts_EcotypePhyla_Summer.csv") %>% 
  column_to_rownames(., var = "Sample") %>% 
  subset(., select = eMED4:Verrucomicrobia)
top100phyvSum.hell <- decostand(top100phyvSum, "hellinger")

top100phyvWin <- read.csv("data/NewSumWinData/POWOW_RawCounts_EcotypePhyla_Winter.csv") %>% 
  column_to_rownames(., var = "Sample") %>% 
  subset(., select = eMED4:Verrucomicrobia)
top100phyvWin.hell <- decostand(top100phyvWin, "hellinger")


#######################################################################################
#                 Linear models for all data regardless of Season                     #
#######################################################################################
eMIT9312lm <- lm(eMIT9312 ~Actinobacteria+Bacteria_unclassified+Bacteroidetes+ Chlamydiae+   
              Chloroflexi+  Cyanobacteria+Euryarchaeota+Firmicutes+Fusobacteria+  Gemmatimonadetes+Gracilibacteria +          
              Lentisphaerae+Planctomycetes+ Proteobacteria+Thaumarchaeota+Verrucomicrobia  , data = top100phyvEco.hell)
step(eMIT9312lm)
summary(eMIT9312lm)
pdf("data/eMIT9312lm.pdf")
avPlots(eMIT9312lm)
dev.off()
eNATL2Alm <- lm(eNATL2A ~Actinobacteria+Bacteria_unclassified+Bacteroidetes+ Chlamydiae+   
              Chloroflexi+  Cyanobacteria+Euryarchaeota+Firmicutes+Fusobacteria+  Gemmatimonadetes+Gracilibacteria +          
              Lentisphaerae+Planctomycetes+ Proteobacteria+Thaumarchaeota+Verrucomicrobia  , data = top100phyvEco.hell)
summary(eNATL2Alm)
avPlots(eNATL2Alm)
pdf("data/eNATL2Alm.pdf")
avPlots(eNATL2Alm)
dev.off()

eMIT9313lm <- lm(eMIT9313 ~Actinobacteria+Bacteria_unclassified+Bacteroidetes+ Chlamydiae+   
              Chloroflexi+  Cyanobacteria+Euryarchaeota+Firmicutes+Fusobacteria+  Gemmatimonadetes+Gracilibacteria +          
              Lentisphaerae+Planctomycetes+ Proteobacteria+Thaumarchaeota+Verrucomicrobia  , data = top100phyvEco.hell)
summary(eMIT9313lm)
pdf("data/eMIT9313lm.pdf")
avPlots(eMIT9313lm)
dev.off()

#######################################################################################
#            Linear models for Season specific ecotype/Phylum correlations            #
#######################################################################################

