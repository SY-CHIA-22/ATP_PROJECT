###########################Script for the figures for P. xylostella abundance data
#Load the required libraries

library(tidyverse)
library(RColorBrewer)
library(colortools)
library("GGally")
library(lattice)
library(ggpubr)

##Working directory
getwd()
setwd("D:/OneDrive - Wageningen University & Research/Entomology_PhD/Data/R-WD/Diversity")

##Load functions
sem <- function(x, na.rm = TRUE) {sd(x, na.rm = TRUE)/sqrt(length(x))} #Standard error around the mean


###Load an prepare the data
##Data loading
Abund0 <- read.csv("Abundance_Diversity_TC.csv", header=T, sep=",", 
                   na.strings=c("","NA"), stringsAsFactors = FALSE)

str(Abund0) #Check the data strucuture
summary(Abund0) #Check data stat summary

##Set the variable characteristics correctly for Abund0
Abund0$Block<- as.factor(Abund0$Block)
Abund0$Time.point <- as.factor(Abund0$Time.point)
Abund0$Feed_guild <- as.factor(Abund0$Feed_guild)
Abund0$Diversity <- as.factor(Abund0$Diversity)
Abund0$Treatments <- as.factor(Abund0$Treatments)
Abund0$Damage <- as.integer(Abund0$Damage)
Abund0$Plutella_n <- as.integer(Abund0$Plutella_n)

#Create a new variable
Abund0$FD <- as.factor(paste(Abund0$Feed_guild, Abund0$Diversity, sep = ""))

str(Abund0)
summary(Abund0) 

###Prepare the data for plotting 


Df <- Abund0%>%
  group_by(Feed_guild, Diversity, Time.point)%>%
  dplyr::summarise(Px = mean(Plutella_n), Px.se = sem(Plutella_n, na.rm = TRUE))%>%
  as.data.frame()

str(Df)
summary(Df)

#Order factor feeding guild correctly
Df$Feed_guild <- factor(Df$Feed_guild, levels = c("Ctrl", "Chewer", "Sap_feeder", "Mixture"))
levels(Df$Feed_guild)


#Add the compact letters from the posthoc test (manually done!!!!!!!!)

Df$Sig <- 0

for (i in 1:14) {
  
  if(Df$Time.point[[i]] == 1) {
    Df$Sig[[i]] <- as.character("ns")
  } else if (Df$Time.point[[i]] == 2 & ((Df$Feed_guild[[i]] == "Ctrl" & Df$Diversity[[i]] == "Ctrl") |
                                        (Df$Feed_guild[[i]] == "Chewer" & Df$Diversity[[i]] == "Div4"))) {
    Df$Sig[[i]] <- as.character("a")
  } else if (Df$Time.point[[i]] == 2 & ((Df$Feed_guild[[i]] == "Chewer" & Df$Diversity[[i]] == "Div1") |
                                        (Df$Feed_guild[[i]] == "Sap_feeder" & Df$Diversity[[i]] == "Div1") |
                                        (Df$Feed_guild[[i]] == "Mixture" & Df$Diversity[[i]] == "Div2") |
                                        (Df$Feed_guild[[i]] == "Sap_feeder" & Df$Diversity[[i]] == "Div4"))){
    Df$Sig[[i]] <- as.character("b")
  } else {
    Df$Sig[[i]] <- as.character("ab")
  }
  
}

Df

#Prepare colours for the graphics
MyColor <- c("#66A61E",  "#666666", "#A6761D",  "#E6AB02")
names(MyColor) <- levels(factor(levels(Df$Feed_guild))) #Asign them to specific levels of feeding guilds
MyColor

#Create label vector to change the Diversity levels into specie richness levels
Rich <- c(Ctrl = "Control", Div1 = "Richness: 1", Div2 = "Richness: 2",
          Div4 = "Richness: 4")

Rich <- c(Ctrl = "Control", Div1 = "R1", Div2 = "R2",
          Div4 = "R4")

#Create graphics for each time points

Bar.1 <- ggplot(subset(Df, Time.point == "1"), aes(x =Diversity , y = Px, fill = Feed_guild)) +
  geom_bar(stat="identity", position = position_dodge(1), width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), position = position_dodge(1), width = 0.2)+
  scale_y_continuous(name = "Average number larvae", 
                     expand = c(0, 0), limits = c(0, 1.3)) +
  scale_x_discrete(name = "", 
                   labels = Rich)+
  scale_fill_manual(values = MyColor, name = "Feeding guilds:",
                    labels = c("Control", "Chewer",
                               "Sap feeder", "Mixture"))+
  ggtitle("Time point 1") +
  geom_text(aes(y = Px + Px.se, label = as.factor(Sig)),  position = position_dodge(1), vjust = -0.5) + #to add text on your graphics
  facet_grid(. ~ Feed_guild
             ,scales = "free_x"
             , space = "free_x") +
  theme_classic2() +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        strip.text.x = element_blank(),
        #axis.text.x  = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"), 
        text = element_text(size = 12))
Bar.1


Bar.2 <- ggplot(subset(Df, Time.point == "2"), aes(x =Diversity , y = Px, fill = Feed_guild)) +
  geom_bar(stat="identity", position = position_dodge(1), width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), position = position_dodge(1), width = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) +
  scale_x_discrete(name = "", 
                   labels = Rich)+
  scale_fill_manual(values = MyColor, name = "Feeding guilds:",
                    labels = c("Control", "Chewer",
                               "Sap feeder", "Mixture"))+
  ggtitle("Time point 2") +
  geom_text(aes(y = Px + Px.se, label = as.factor(Sig)),  position = position_dodge(1), vjust = -0.5) +
  facet_grid(. ~ Feed_guild
             ,scales = "free_x"
             , space = "free_x") +
  theme_classic2() +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.y = element_blank(),
        #axis.text.x =  element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        strip.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 12))
Bar.2


#Arrange the graohic on the same page
Px.Abundance <- ggarrange(Bar.1, Bar.2, 
                          labels = c("", ""),
                          common.legend = TRUE,
                          legend = "bottom",
                          ncol = 2, nrow = 1)
Px.Abundance

#Export as pdf
ggsave("Px.abundance.pdf", plot = Px.Abundance, device = "pdf")

###For the treatments levels
Df1 <- Abund0%>%
  group_by(Feed_guild, Treatments, Time.point)%>%
  dplyr::summarise(Px = mean(Plutella_n), Px.se = sem(Plutella_n, na.rm = TRUE))%>%
  as.data.frame()

str(Df1)
summary(Df1)

Df1$Feed_guild <- factor(Df1$Feed_guild, levels = c("Ctrl", "Chewer", "Sap_feeder", "Mixture"))
levels(Df1$Feed_guild)

Df1$Treatments <- factor(Df1$Treatments, levels = c("Ctrl", "Ar", "Mb", "Pb", "Pc",
                                                    "Bb", "Le", "Mp", "Mpn", "MIX2", "MbPbPcAr", "BbMpMpnLe", "MIX4"))
levels(Df1$Treatments)

Species <- c(Ctrl = "Control", Ar = "Athalia rosae", Mb = "Mamestra brassicae", Mp = "Myzus persicae", Mpn = "Myzus persicae nicotianae",
             Pb = "Pieris brassicae", Pc = "Phaedon cochleariae", Bb = "Brevicoryne brassicae", Le = "Lipaphis erysimi",
             BbMpMpnLe = "All Aphid", MbPbPcAr = "All chewer", MIX2 = "Mixture of 2", MIX4 = "Mixture of 4")


BarT1 <- ggplot(subset(Df1, Time.point=="1"), aes(x = Treatments , y = Px, fill = Feed_guild)) +
  geom_bar(stat="identity", position = position_dodge(1), width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), position = position_dodge(1), width = 0.2)+
  scale_y_continuous(name = "Average number larvae", 
                     expand = c(0, 0), limits = c(0, 1.4)) +
  scale_x_discrete(name = "", 
                   labels = Species)+
  scale_fill_manual(values = MyColor, name = "Feeding guilds:",
                    labels = c("Control", "Chewer",
                               "Sap feeder", "Mixture"))+
  ggtitle("Time point 1") +
  theme_classic2() +
  theme(legend.direction = "horizontal",
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"), 
        text = element_text(size = 12))
BarT1

BarT2 <- ggplot(subset(Df1, Time.point == "2"), aes(x =Treatments , y = Px, fill = Feed_guild)) +
  geom_bar(stat="identity", position = position_dodge(1), width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), position = position_dodge(1), width = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.4)) +
  scale_x_discrete(name = "", 
                   labels = Species)+
  scale_fill_manual(values = MyColor, name = "Feeding guilds:",
                    labels = c("Control", "Chewer",
                               "Sap feeder", "Mixture"))+
  ggtitle("Time point 2") +
  theme_classic2() +
  theme(legend.direction = "horizontal",
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 12))
BarT2

Px.Abundance_Sp <- ggarrange(BarT1, BarT2, 
                             labels = c("", ""),
                             common.legend = TRUE,
                             legend = "none",
                             ncol = 2, nrow = 1)

Px.Abundance_Sp

ggsave("Px.abundance-Sp.pdf", plot = Px.Abundance_Sp, device = "pdf")


#########################################################################################################################################
######################### P. xylostella parasitism levels data

##Data loading
Para0 <- read.csv("Parasitism_Diversity_TC.csv", header=T, sep=",", 
                  na.strings=c("","NA"), stringsAsFactors = FALSE)

str(Para0) #Check the data strucuture
summary(Para0) #Check data stat summary

##Set the variable characteristics correctly for Abund0
Para0$Block<- as.factor(Para0$Block)
Para0$Time.point <- as.factor(Para0$Time.point)
Para0$Feed_guild <- as.factor(Para0$Feed_guild )
Para0$Diversity <- as.factor(Para0$Diversity)
Para0$Treatments <- as.factor(Para0$Treatments)
Para0$Px.start <- as.integer(Para0$Px.start)
Para0$Px.end <- as.integer(Para0$Px.end)
Para0$Px.para <- as.integer(Para0$Px.para )

str(Para0)
summary(Para0) #44 and 33 missing values for Px recovered and Px parasitised respectively?

subset(Para0, is.na(Px.end)) #4 observations have 1 parasitised Px recorded despite the Na values.
#Some NA values correspond to zero Px parasitised.
subset(Para0, is.na(Px.para)) #A 2 of missing values for parasitised Px had caterpillars recovered

#For now, I will remove the 44 NA values in Px.end, untill clarifications.
Para0 <- filter(Para0, !(is.na(Px.end)))

#9 observations have missing values for parasitised Px. Note only two of them have Px.end of 3 and 1 respectively.
Para0<-filter(Para0, !(is.na(Px.para)))

#Remove block 5 (due to unfortunate data loss)
Para0 <- subset(Para0, !(Block == "5"))

#Relevel factors
Para0$Block <- factor(Para0$Block)
levels(Para0$Block)

#Remove observations were no Px.end were recovered
Para0 <- subset(Para0, Px.end >0)#Remove observations were no Px.end were recovered
Para0 <- Para0[,-c(3)]
summary(Para0)

#Prepare the data for plotting
Dfp1 <- Para0 %>%
  group_by(Feed_guild)%>%
  dplyr::summarise(Px = mean(Px.para, na.rm = TRUE), Px.se = sem(Px.para, na.rm = TRUE))%>%
  as.data.frame()

str(Dfp1)
summary(Dfp1)

Dfp1$Feed_guild <- factor(Dfp1$Feed_guild, levels = c("Ctrl", "Chewer", "Sap_feeder", "Mixture"))
levels(Dfp1$Feed_guild)

#Put the compact letter based from the post hoc test (Done manually!!!!!!)
Dfp1$Sig <- 0

for (i in 1:4) {
  if (Dfp1$Feed_guild[[i]] == "Chewer" | Dfp1$Feed_guild[[i]] == "Sap_feeder" |  Dfp1$Feed_guild[[i]] == "Ctrl"){
    Dfp1$Sig[[i]] <- as.character("a")
  } else if (Dfp1$Feed_guild[[i]] == "Mixture") {
    Dfp1$Sig[[i]] <- as.character("b")
  } else { 
    Dfp1$Sig[[i]] <- as.character("")
  }
  
}

Dfp1

MyColor <- c("#66A61E",  "#666666", "#A6761D",  "#E6AB02")
names(MyColor) <- levels(factor(levels(Dfp1$Feed_guild)))
MyColor


Feed <- c(Ctrl = "Control",  Chewer  = " Chewer", Sap_feeder = "Sap feeder",
          Mixture = "Mixture")

BarF <- ggplot(Dfp1, aes(x = Feed_guild , y = Px, fill = Feed_guild)) +
  geom_bar(stat="identity",  alpha = 0.9) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), width = 0.2)+
  scale_y_continuous(name = "Average number parasitised larvae", 
                     expand = c(0, 0), limits = c(0, 2)) +
  scale_x_discrete(name = "", 
                   labels = Feed)+
  scale_fill_manual(values = MyColor, name = "Feeding guilds:",
                    labels = Feed)+
  geom_text(aes(y = Px + Px.se, label = as.factor(Sig)), vjust = -0.8) +
  theme_classic2() +
  theme(
        legend.direction = "horizontal",
        legend.position = "none",
        plot.margin = unit(c(1,-0.5,1,0.1), "cm"), #t, r, b, l  Dimensions of each margin. (To remember order, think trouble).
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"), 
        text = element_text(size = 12))
BarF

Dfp2 <- Para0 %>%
  group_by(Time.point)%>%
  dplyr::summarise(Px = mean(Px.para, na.rm = TRUE), Px.se = sem(Px.para, na.rm = TRUE))%>%
  as.data.frame()

str(Dfp2)
summary(Dfp2)

post_hocT <- emmeans(MP6, specs = pairwise ~ Time.point, type = "response")
summary(post_hocT)

complementary("steelblue")

MyColor2 <- c("#4682B4", "#B47846")
names(MyColor2) <- levels(factor(levels(Dfp2$Time.point)))
MyColor2

BarT <- ggplot(Dfp2, aes(x = Time.point , y = Px, fill = Time.point)) +
  geom_bar(stat="identity", position = "identity", alpha = 0.9, width =  0.5) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), width = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) +
  scale_x_discrete(name = "", 
                   labels = c("Time point 1", "Time point 2"))+
  scale_fill_manual(values = MyColor2)+
  geom_segment(aes(x = 1.1, xend = 1.9, y = 1.35, yend = 1.35)) +
   annotate("text",
           x = 1.5,
           y = 1.4,
           label = c("p < 0.001"),
           family = "", fontface = 3, size=4) + 
  theme_classic2() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,0.1,1,-0.5), "cm"), 
        axis.text.y = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 12))
BarT


PxF <- ggarrange(BarF, BarT, 
                 labels = c("A", "B"),
                 common.legend = FALSE,
                 legend = "none",
                 ncol = 2, nrow = 1)

PxF

ggsave("Px_para.pdf", plot = PxF, device = "pdf")

##Figures for the treatments levels
DfpS <- Para0%>%
  group_by(Feed_guild, Treatments)%>%
  dplyr::summarise(Px = mean(Px.para, na.rm = TRUE), Px.se = sem(Px.para, na.rm = TRUE))%>%
  as.data.frame()

str(DfpS)
summary(DfpS)

DfpS$Feed_guild <- factor(DfpS$Feed_guild, levels = c("Ctrl", "Chewer", "Sap_feeder", "Mixture"))
levels(DfpS$Feed_guild)

levels(DfpS$Treatments)

DfpS$Treatments <- factor(DfpS$Treatments, levels = c("Ctrl", "Ar", "Mb", "Pb", "Pc", "MbPbPcAr",
                                                      "Bb", "Le", "Mp", "Mpn", "BbMpMpnLe", "MIX2", "MIX4"))
levels(DfpS$Treatments)

Species <- c(Ctrl = "Control", Ar = "Athalia rosae", Mb = "Mamestra brassicae", Mp = "Myzus persicae", Mpn = "Myzus persicae nicotianae",
             Pb = "Pieris brassicae", Pc = "Phaedon cochleariae", Bb = "Brevicoryne brassicae", Le = "Lipaphis erysimi",
             BbMpMpnLe = "All Aphid", MbPbPcAr = "All chewer", MIX2 = "Mixture of 2", MIX4 = "Mixture of 4")


#Put the compact letter based from the post hoc test (Done manually!!!!!!)
DfpS$Sig <- 0

for (i in 1:13) {
  if (DfpS$Treatments[[i]] == "Ctrl" ){
    DfpS$Sig[[i]] <- as.character("a")
  } else if (DfpS$Treatments[[i]] == "Bb" | DfpS$Treatments[[i]] == "Le" |  DfpS$Treatments[[i]] == "Mb" | DfpS$Treatments[[i]] == "Pc" | 
             DfpS$Treatments[[i]] == "BbMpMpnLe" | DfpS$Treatments[[i]] == "MbPbPcAr") {
    DfpS$Sig[[i]] <- as.character("ab")
  } else if (DfpS$Treatments[[i]] == "MIX2" | DfpS$Treatments[[i]] == "Ar") {
    DfpS$Sig[[i]] <- as.character("abc")
  } else if (DfpS$Treatments[[i]] == "MIX4") {
    DfpS$Sig[[i]] <- as.character("c")
   } else { 
    DfpS$Sig[[i]] <- as.character("bc")
  }
  
}

DfpS

#.With compact leters CLD (see analysis)
tuk.cld$mcletters$Letters

DfpS <- cbind(DfpS, tuk.cld$mcletters$Letters)
DfpS


BarS <- ggplot(DfpS, aes(x = Treatments , y = Px, fill = Feed_guild)) +
  geom_bar(stat="identity", position = position_dodge(1), width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin= Px - Px.se , ymax= Px + Px.se), position = position_dodge(1), width = 0.2)+
  scale_y_continuous(name = "Average number of parasitised larvae", 
                     expand = c(0, 0), limits = c(0, 1.8)) +
  scale_x_discrete(name = "", 
                   labels = Species)+
  scale_fill_manual(values = MyColor, name = "Feeding guilds:",
                    labels = c("Control", "Chewer",
                               "Sap feeder", "Mixture"))+
  geom_text(aes(y = Px + Px.se, label = as.factor(tuk.cld$mcletters$Letters)), vjust = -0.8) +
  theme_classic2() +
  theme(legend.direction = "horizontal",
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"), 
        text = element_text(size = 12))
BarS

ggsave("Px.para-Sp_Bis.pdf", plot = BarS, device = "pdf")




