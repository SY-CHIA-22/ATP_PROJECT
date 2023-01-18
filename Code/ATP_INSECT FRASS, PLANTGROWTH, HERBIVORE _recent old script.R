dir() 
setwd("~/OneDrive/R-WORKING DIRECTORY/ANALYSIS DATA_2020-2024")
GrowthData<-read.table("leafarea.txt",header=TRUE)
leavesA<-read.table("leaves.txt",header=TRUE)
GrowthData_B<-read.table("leafareaB.txt",header=TRUE)
leavesB<-read.table("leavesB.txt",header=TRUE)
Delia_no<-read.table("Delia_recovered.txt",header=TRUE)
Delia_pupae<-read.table("Pupalweight.txt",header=TRUE)
DBM_no<-read.table("DBM_collected.txt",header=TRUE)
DBM_damage<-read.table("DBM_damagescore.txt",header=TRUE)
DBM_tot<-read.table("DBM_total.txt", header = TRUE)
LeafC<-read.table("leaves_C.txt", header = TRUE)# incubated frass in soil, no time added
LeafD<-read.table("leaves_D.txt", header=TRUE) # composted frass in lunch boxes
LeafE<-read.table("leaves_E.txt", header = TRUE)
GrowthIncub<-read.table("growth_incub_frass.txt", header = TRUE)# incubated frass in soil
Growth_comp<-read.table("growth_com_frass.txt", header=TRUE) # composted frass in lunch boxes

### # Modify legend titles
P + labs(color = "Soil amendment") # if you change "color" to "fill", it won't work
# P = plot object (can be anything or any letter or group of letters of choice)


###NOTE THAT SET A,B = Trial 1, 2 (Same trial conducted twice due to shortage of greenhouse space)

############### Leaf AREA DATA SET A
# save data in appropriate format
str(GrowthData)# view data structure, factors,etc.
GrowthData$Treatment <- as.factor(GrowthData$Treatment) # save Treatment as a factor
GrowthData$Time <- as.factor(GrowthData$Time) # save Time as a factor
GrowthData$Plant_ID <- as.factor(GrowthData$Plant_ID) # save Plant_ID as a factor

# model selection
Ma <- lmer(Area ~ Treatment * Time + (1|Plant_ID), data = GrowthData)
Mb <- glmmTMB(Area ~ Treatment * Time + (1|Plant_ID), data = GrowthData)
Md <- glmer(Area ~ Treatment * Time + (1|Plant_ID), family = Gamma(link = "log"), data = GrowthData)# Fit a Gamma model
Me <- glmer.nb(Area ~ Treatment * Time + (1|Plant_ID), data = GrowthData)
Mf <- glmer(Area ~ Treatment * Time + (1|Plant_ID), family = poisson (link = "log"), data = GrowthData)
AIC(Ma, Mb, Md, Me, Mf) 

#fit negative binomial model. this is best based on AIC VALUES. the gamma model fails to converge though with a similar AIC VALUE
library(lme4)
Me <- glmer.nb(Area ~ Treatment * Time + (1|Plant_ID), data = GrowthData)

#check for dispersion
P<-length(coef(Me))
N<-nrow(GrowthData)
E<-resid(Me,type="pearson")
Dispersion<-sum(E^2)/(N-P)
Dispersion
#  0.9929559 This is less than 1.9. therefore there is no overdispersion
Anova(Me)

emmeans # estimated marginal means. This tells you the mean response for each factor, adjusted for any other variables in the model.
library (emmeans) # load the emmeans package

# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Me, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments * timepoints
HSD_test <- emmeans(ref_grid(Me, transform = "response"), pairwise ~ Treatment * Time)
emmeans(ref_grid(Me, transform = "response"), pairwise ~ Treatment * Time) # the asterisk (*) in: pairwise ~ Treatment*Time makes several,less meaningful comparisons 
emmeans(ref_grid(Me, transform = "response"), pairwise ~ Treatment | Time) # changing "*" to "|" allows meaningful comparisons (i.e at experimental time points)
emmeans(ref_grid(Me, transform = "response"), pairwise ~ Time|Treatment)
library(ggplot2) # to plot

estimates <- emmeans(ref_grid(Me, transform = "response"), pairwise ~ Treatment * Time)
plot_output <- as.data.frame(estimates$emmeans) # plot estimates using emmeams package
ggplot(plot_output, aes(x = Time, y = response, color = Treatment, fill = Treatment)) + geom_point(size = 2)+ geom_errorbar(aes(ymin = response - SE, ymax = response + SE), 
width = 0.2) + geom_line(aes(group = Treatment)) + theme_bw () # plot the effect (response) of Treatment on leaf area over time or at a given time point

### number of plant leaves (B. rapa) leaves data
# save data in appropriate format
str(leavesA)# view data structure, factors,etc.
leavesA$Treatment <- as.factor(leavesA$Treatment) # save Treatment as a factor
leavesA$Time <- as.factor(leavesA$Time) # save Time as a factor
leavesA$Plant_ID <- as.factor(leavesA$Plant_ID) # save Plant_ID as a factor

# model selection for leaves data set A
Ma1 <- lmer(No_leaves ~ Treatment * Time + (1|Plant_ID), data = leavesA)
Mb1 <- glmmTMB(No_leaves ~ Treatment * Time + (1|Plant_ID), data = leavesA)
Md1 <- glmer(No_leaves ~ Treatment * Time + (1|Plant_ID), family = Gamma(link = "log"), data = leavesA)# Fit a Gamma model
Me1 <- glmer.nb(No_leaves ~ Treatment * Time + (1|Plant_ID), data = leavesA)
Mf1 <- glmer(No_leaves ~ Treatment * Time + (1|Plant_ID), family = poisson (link = "log"), data = leavesA)
AIC(Ma1, Mb1, Md1, Me1, Mf1) # Mb1, glmmTMB is the best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(No_leaves~Treatment,data=leavesA)) # run this and the next command to get plot
plot(lm(No_leaves~Time,data=leavesA)) # homogeneity of residuals confirmed here

library(car)
## tests for homogeneity (equality) of variances in data
leveneTest(No_leaves ~ Treatment, data = leavesA)
fligner.test(No_leaves ~ Treatment, data = leavesA)
bartlett.test(No_leaves ~ Treatment, data = leavesA) 
# variances not homogeneous according to the tests but visual assessment shows homogeneity, then lets proceed with GLMM

#fit selected model
Mb1 <- glmmTMB(No_leaves ~ Treatment * Time + (1|Plant_ID), data = leavesA)

# Check normality of residuals
qqnorm(resid(Mb1,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Mb1,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(leavesA$Treatment, resid(Mb1, type = "pearson"), main = "Treatment")
plot(leavesA$Time, resid(Mb1, type = "pearson"), main = "Time")

Anova(Mb1)


emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Mb1, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments * timepoints
HSD_test <- emmeans(ref_grid(Mb1, transform = "response"), pairwise ~ Treatment * Time)
emmeans(ref_grid(Mb1, transform = "response"), pairwise ~ Treatment * Time) 
emmeans(ref_grid(Mb1, transform = "response"), pairwise ~ Treatment | Time) 

library(ggplot2) # to plot

estimates <- emmeans(ref_grid(Mb1, transform = "response"), pairwise ~ Treatment * Time)
plot_output <- as.data.frame(estimates$emmeans) # plot estimates using emmeams package
ggplot(plot_output, aes(x = Time, y = emmean, color = Treatment, fill = Treatment)) + 
geom_point(size = 2)+ geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
width = 0.2) + geom_line(aes(group = Treatment)) + theme_bw ()


#######################################################################################
### Leaf AREA DATA SET B

# save data in appropriate format

GrowthData_B$Treatment <- as.factor(GrowthData_B$Treatment) # save Treatment as a factor
GrowthData_B$Time <- as.factor(GrowthData_B$Time) # save Time as a factor
GrowthData_B$Plant_ID <- as.factor(GrowthData_B$Plant_ID) # save Plant_ID as a factor
str(GrowthData_B)# view data structure, factors,etc.

# model selection
F1 <- lmer(area ~ Treatment * Time + (1|Plant_ID), data = GrowthData_B)
F2 <- glmmTMB(area ~ Treatment * Time + (1|Plant_ID), data = GrowthData_B)
F3 <- glmer(area ~ Treatment * Time + (1|Plant_ID), family = Gamma(link = "log"), data = GrowthData_B)# Fit a Gamma model
Fit4 <- glmer.nb(area ~ Treatment * Time + (1|Plant_ID), data = GrowthData_B)
F5 <- glmer(area ~ Treatment * Time + (1|Plant_ID), family = poisson (link = "log"), data = GrowthData_B)
AIC(F1, F2, F3, F5) 

#fit negative binomial model. this is best based on AIC VALUES. the gamma model fails to converge though with a similar AIC VALUE
library(lme4)
F3 <- glmer(area ~ Treatment * Time + (1|Plant_ID), family = Gamma(link = "log"), data = GrowthData_B)# Fit a Gamma model

# Check normality of residuals
qqnorm(resid(F3,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(F3,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(GrowthData_B$Treatment, resid(F3, type = "pearson"), main = "Treatment")
plot(GrowthData_B$Time, resid(F3, type = "pearson"), main = "Time")

library(car)
Anova(F3) # in car package

emmeans # estimated marginal means. This tells you the mean response for each factor, adjusted for any other variables in the model.
library (emmeans) # load the emmeans package

# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(F3, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments * timepoints
HSD_test <- emmeans(ref_grid(F3, transform = "response"), pairwise ~ Treatment * Time)
emmeans(ref_grid(F3, transform = "response"), pairwise ~ Treatment * Time) # the asterisk (*) in: pairwise ~ Treatment*Time makes several,less meaningful comparisons 
emmeans(ref_grid(F3, transform = "response"), pairwise ~ Treatment | Time) # changing "*" to "|" allows meaningful comparisons (i.e at experimental time points)
emmeans(ref_grid(F3, transform = "response"), pairwise ~ Time|Treatment)
library(ggplot2) # to plot

estimates <- emmeans(ref_grid(F3, transform = "response"), pairwise ~ Treatment * Time)
plot_output <- as.data.frame(estimates$emmeans) # plot estimates using emmeams package
ggplot(plot_output, aes(x = Time, y = response, color = Treatment, fill = Treatment)) + geom_point(size = 2)+ geom_errorbar(aes(ymin = response - SE, ymax = response + SE), 
                                                                                                                            width = 0.2) + geom_line(aes(group = Treatment)) + theme_bw () # plot the effect (response) of Treatment on leaf area over time or at a given time point

################################################
### number of plant leaves (B. rapa) leaves data, SET B
# save data in appropriate format

leavesB$Treatment <- as.factor(leavesB$Treatment) # save Treatment as a factor
leavesB$Time <- as.factor(leavesB$Time) # save Time as a factor
leavesB$Plant_ID <- as.factor(leavesB$Plant_ID) # save Plant_ID as a factor
str(leavesB)# view data structure, factors,etc.

# model selection for leaves data set A
F1a <- lmer(No_leavesb ~ Treatment * Time + (1|Plant_ID), data = leavesB)
F2a <- glmmTMB(No_leavesb ~ Treatment * Time + (1|Plant_ID), data = leavesB)
F3a <- glmer(No_leavesb ~ Treatment * Time + (1|Plant_ID), family = Gamma(link = "log"), data = leavesB)# Fit a Gamma model
F4 <- glmer.nb(No_leavesb ~ Treatment * Time + (1|Plant_ID), data = leavesB)
F5a <- glmer(No_leavesb ~ Treatment * Time + (1|Plant_ID), family = poisson (link = "log"), data = leavesB)
AIC(F1a, F2a, F3a, F4, F5a) # F3a, gamma is the best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(No_leavesb~Treatment,data=leavesB)) # run this and the next command to get plot
plot(lm(No_leavesb~Time,data=leavesB)) # homogeneity of residuals confirmed here

#fit selected model
F3a <- glmer(No_leavesb ~ Treatment * Time + (1|Plant_ID), family = Gamma(link = "log"), data = leavesB)# Fit a Gamma model

# Check normality of residuals
qqnorm(resid(F3a,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(F3a,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(leavesB$Treatment, resid(F3a, type = "pearson"), main = "Treatment")
plot(leavesB$Time, resid(F3a, type = "pearson"), main = "Time")

Anova(F3a)

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(F3a, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments * timepoints
HSD_test <- emmeans(ref_grid(F3a, transform = "response"), pairwise ~ Treatment * Time)
emmeans(ref_grid(F3a, transform = "response"), pairwise ~ Treatment * Time) 
emmeans(ref_grid(F3a, transform = "response"), pairwise ~ Treatment | Time) 

library(ggplot2) # to plot

estimates <- emmeans(ref_grid(Mb1, transform = "response"), pairwise ~ Treatment * Time)
plot_output <- as.data.frame(estimates$emmeans) # plot estimates using emmeams package
ggplot(plot_output, aes(x = Time, y = emmean, color = Treatment, fill = Treatment)) + 
  geom_point(size = 2)+ geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                                      width = 0.2) + geom_line(aes(group = Treatment)) + theme_bw ()


#######################################################################################

################# Number of Delia larvae plus pupae recovered, SET A DATA ########

# save data in appropriate format
str(Delia_no)# view data structure, factors,etc.
Delia_no$Treatment <- as.factor(Delia_no$Treatment) # save Treatment as a factor
Delia_no$Plant_ID <- as.factor(Delia_no$Plant_ID)

# model selection for number of Delia larvae and pupae data set A
str(DBM_tot)
Ma2b <- glm(Lar_Pup ~ Treatment, data = Delia_no)
Ma2c <- lm(Lar_Pup ~ Treatment, data = Delia_no)
Ma2d <- glm(Lar_Pup ~ Treatment, family = poisson (link = "log"), data = Delia_no)
Mb2 <- glmmTMB(Lar_Pup ~ Treatment + (1|Plant_ID), data = Delia_no)
Mc <- glmer.nb(Lar_Pup ~ Treatment + (1|Plant_ID), data = Delia_no)
Md2a <- glm.nb(Lar_Pup ~ Treatment, data = Delia_no)
Md2b <- glm.nb(Lar_Pup ~ Treatment+Plant_ID, data = Delia_no)
Me2 <- glmer(Lar_Pup ~ Treatment + (1|Plant_ID), family = poisson (link = "log"), data = Delia_no)
Mf2<-aov(Lar_Pup~Treatment*Plant_ID, data = Delia_no)
Mf2a <-aov(Lar_Pup~Treatment + Plant_ID, data = Delia_no)
Mf2b<-aov(Lar_Pup~Treatment, data = Delia_no)
AIC(Ma2b, Ma2c, Ma2d, Mb2, Mc, Md2a, Md2b, Me2, Mf2, Mf2a, Mf2b) # Ma2b, Mf2b, Ma2c are the best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(Lar_Pup~Treatment,data=Delia_no)) # run this and hit next command to get plot


#fit selected model: Ma2b
Ma2b <- glm(Lar_Pup ~ Treatment, data = Delia_no)# same results as the next
Mb2 <- glmmTMB(Lar_Pup ~ Treatment + (1|Plant_ID), data = Delia_no)
# Check normality of residuals
qqnorm(resid(Ma2b,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Ma2b,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(Delia_no$Treatment, resid(Ma2b, type = "pearson"), main = "Treatment")

Anova(Ma2b)

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Ma2b, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(Ma2b, transform = "response"), pairwise ~ Treatment)
emmeans(ref_grid(Ma2b, transform = "response"), pairwise ~ Treatment) 


################# barplots
## number of delia pupae plus larvae recovered, SET A DATA
library(plotrix)
summary.delia<-group_by(Delia_no, Treatment)  %>% 
  summarise(mean=mean(Lar_Pup),std_error=std.error(Lar_Pup)) %>% arrange(desc(mean))
summarised_Delia_no<-summarise_each(Delia_no,funs(mean=mean,std_error=std.error))
View(summary.delia)

ano1=aov(Lar_Pup ~ Treatment, data=Delia_no)
summary(ano1)

#tukey's test

T1<-TukeyHSD(ano1)
print(T1)

# compact letter display
library(multcompView)
T.cld1 = multcompLetters4(ano1,T1)
print(T.cld1)

cld1<-as.data.frame.list(T.cld1$'Treatment')
summary.delia$T1=cld1$Letters
View(summary.delia)
write_csv(summary.delia,"Delia_no_summary.csv") # save table as csv file

#Barplots
ggplot(summary.delia, aes(x=factor(Treatment), y = mean, fill=Treatment, color=Treatment))+
  geom_bar(stat = "identity", position = "dodge",alpha=0.5)+
  geom_errorbar(aes(ymin=mean-std_error, ymax=mean+std_error),
                position=position_dodge(0.9), width=0.25, show.legend = FALSE)+
  labs(x="Soil amendment",y="Number of CRF (larvae + pupae)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #leave legend outside
  geom_text(aes(label=T1), position=position_dodge(0.9),show.legend = FALSE, size=4,
            vjust=-1.4, hjust=-2.0, color="gray25")+
  ylim(0,8)+ # sets y axis to prevent cutting off error bars
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") 
ggsave("barplot.png", width = 4.5, height = 3, dpi = 1000)# SAVE PLOT



###### Number of Delia PUPAE recovered ##########

# model selection for number of Delia larvae and pupae data set A
Ma3b <- glm(Pupae ~ Treatment, data = Delia_no)
Ma3c <- lm(Pupae ~ Treatment, data = Delia_no)
Ma3d <- glm(Pupae ~ Treatment, family = poisson (link = "log"), data = Delia_no)
Mb3 <- glmmTMB(Pupae ~ Treatment + (1|Plant_ID), data = Delia_no)
MC3 <- glmer.nb(Pupae ~ Treatment + (1|Plant_ID), data = Delia_no)
Md3a <- glm.nb(Pupae ~ Treatment, data = Delia_no)
Md3b <- glm.nb(Pupae ~ Treatment+Plant_ID, data = Delia_no)
Me3 <- glmer(Pupae ~ Treatment + (1|Plant_ID), family = poisson (link = "log"), data = Delia_no)
Mf3<-aov(Pupae~Treatment*Plant_ID, data = Delia_no)
Mf3a <-aov(Pupae~Treatment + Plant_ID, data = Delia_no)
Mf3b<-aov(Pupae~Treatment, data = Delia_no)
AIC(Ma3b, Ma3c, Ma3d, Mb3, MC3, Md3a, Md3b, Me3, Mf3, Mf3a, Mf3b) # Ma3b best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(Pupae~Treatment,data=Delia_no)) # run this and hit next command to get plot

#fit selected model: Ma3b
Ma3b <- glm(Pupae ~ Treatment, data = Delia_no)# same results as the next

# Check normality of residuals
qqnorm(resid(Ma3b,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Ma3b,  type = "pearson"))

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Ma3b, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(Ma3b, transform = "response"), pairwise ~ Treatment)
emmeans(ref_grid(Ma3b, transform = "response"), pairwise ~ Treatment) 



###### Number of Delia LARVAE recovered, SET A DATA ##########

# model selection for number of Delia larvae and pupae data set A
M3b1 <- glm(tot_larvae ~ Treatment, data = Delia_no)
M3c <- lm(tot_larvae ~ Treatment, data = Delia_no)
M3d <- glm(tot_larvae ~ Treatment, family = poisson (link = "log"), data = Delia_no)
M3i <- glmmTMB(tot_larvae ~ Treatment, data = Delia_no)
M3 <- glm.nb(tot_larvae ~ Treatment, data = Delia_no)
M3a <- glm.nb(tot_larvae ~ Treatment, data = Delia_no)
M3b <- glm.nb(tot_larvae ~ Treatment, data = Delia_no)
M3c1 <- glm(tot_larvae ~ Treatment, family = poisson (link = "log"), data = Delia_no)
Mf3<-aov(tot_larvae~Treatment, data = Delia_no)
Mf3a <-aov(tot_larvae~Treatment, data = Delia_no)
Mf3i<-aov(tot_larvae~Treatment, data = Delia_no)
AIC(M3b1, M3c, M3d, M3i, M3, M3a, M3b, M3c1, Mf3, Mf3a, Mf3i) # Ma3d is the best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(tot_larvae~Treatment,data=Delia_no)) # run this and hit next command to get plot

#fit selected model: Ma3d
M3d <- glm(tot_larvae ~ Treatment, family = poisson (link = "log"), data = Delia_no)

#check for dispersion
P<-length(coef(M3d))
N<-nrow(Delia_no)
E<-resid(M3d,type="pearson")
Dispersion<-sum(E^2)/(N-P)
Dispersion # 0.7572016 looks fine


# Check normality of residuals
qqnorm(resid(M3d,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(M3d,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(Delia_no$Treatment, resid(M3d, type = "pearson"), main = "Treatment")

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(M3d, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(M3d, transform = "response"), pairwise ~ Treatment)
emmeans(ref_grid(M3d, transform = "response"), pairwise ~ Treatment)



###### Delia PUPAE WEIGHT, SET A DATA ##########

# model selection for number of Delia larvae and pupae data set A
M4a <- glm(Delia ~ Treatment, data = Delia_pupae)
M4b <- lm(Delia ~ Treatment, data = Delia_pupae)
M4c <- glm(Delia ~ Treatment, family = Gamma(link = "log"), data = Delia_pupae)
M4d <- glmmTMB(Delia ~ Treatment, data = Delia_pupae)
M4e <- glm.nb(Delia ~ Treatment, data = Delia_pupae)
M4f <- glm.nb(Delia ~ Treatment, data = Delia_pupae)
M4g <- glm.nb(Delia ~ Treatment, data = Delia_pupae)
M4h <- glm(Delia ~ Treatment, family = poisson (link = "log"), data = Delia_pupae)
M4i<-aov(Delia~Treatment, data = Delia_pupae)
M4j <-aov(Delia~Treatment, data = Delia_pupae)
M4k<-aov(Delia~Treatment, data = Delia_pupae)
AIC(M4a, M4b, M4c, M4d, M4e, M4f, M4g, M4h, M4i, M4j , M4k) # M4a is the best, least AIC, also similar to others

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(Delia~Treatment,data=Delia_pupae)) # run this and hit next command to get plot

#fit selected model: Ma3d
M4a <- glm(Delia ~ Treatment, data = Delia_pupae)
Anova(M4a) # anova function in "car" package, provides overall p value to indicate significance

# Check normality of residuals
qqnorm(resid(M4a,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(M4a,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(Delia_pupae$Treatment, resid(M4a, type = "pearson"), main = "Treatment")

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(M4a, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(M4a, transform = "response"), pairwise ~ Treatment)
emmeans(ref_grid(M4a, transform = "response"), pairwise ~ Treatment)

# plots, WEIGHT OF DELIA PUPAE, SET A DATA
# save data in appropriate format
str(Delia_pupae)# view data structure, factors,etc.
Delia_pupae$Treatment <- as.factor(Delia_pupae$Treatment) # save Treatment as a factor

## weight of delia pupae recovered

ano2=aov(Delia ~ Treatment, data=Delia_pupae)
summary(ano2)

#tukey's test

T2<-TukeyHSD(ano2)
print(T2)

# compact letter display
library(multcompView)
T.cld2 = multcompLetters4(ano2,T2)
print(T.cld2)

# Table with factors and third quantile(tk2 - Delia)
tk2<-group_by(Delia_pupae, Treatment)  %>% summarise(mean=mean(Delia), 
quant = quantile(Delia, probs = 0.75)) %>% arrange(desc(mean))

# Extract the compact letter display and add to the T table (tk2 - Delia)
cld2<-as.data.frame.list(T.cld2$Treatment)
tk2$cld2<-cld2$Letters
print(tk2)

## Boxplot - weight of Delia pupae, SET A DATA
pupawt <- ggplot(Delia_pupae, aes(Treatment, Delia))+ 
  geom_boxplot(aes(fill = Treatment), show.legend = FALSE) +
  theme_minimal() + labs(x="Soil amendment", y = "Mean weight (mg) of CRF pupae")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())+
  geom_text(data=tk2, aes(label=cld2, x=Treatment, y=quant),# brings in the letters
            vjust=-1, hjust=-0.9, size=4.5) + geom_jitter(width = 0.03)
pupawt + ylim(0,13)

#Save  the final figure
ggsave("boxplot.png", width = 4, height = 3, dpi = 1000)


############# Number of Plutella (DBM) LARVAE + PUPA COLLECTED, SET A

# save data in appropriate format
str(DBM_no)# view data structure, factors,etc.
DBM_no$Treatment <- as.factor(DBM_no$Treatment) # save Treatment as a factor
DBM_no$Plant_ID <- as.factor(DBM_no$Plant_ID)
DBM_no$Time <- as.factor(DBM_no$Time)
DBM_no$Tray <- as.factor(DBM_no$Tray)

# model selection for NUMBER OF DBM LARVAE data, set A
M <- glm(Larvae  ~ Treatment + Time, data = DBM_no)
N <- glmmTMB(Larvae  ~ Treatment * Time, data = DBM_no)
P<- glm.nb(Larvae  ~ Treatment * Time, data = DBM_no)
Q <- glm(Larvae  ~ Treatment * Time, family = poisson, data = DBM_no)
N1 <- glm(Larvae  ~ Treatment*Time, family = quasipoisson, data = DBM_no)
AIC(M, N, P, Q, N1) # Q is the best, NO OVERDISPERSION , AIC ALSO GOOD

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(Larvae~Treatment,data=DBM_no)) # run this and hit next command to get plot

#fit selected model: Ma3d
Q <- glm(Larvae ~ Treatment * Time, family = poisson, data = DBM_no)

#check for dispersion
# METHOD 1
P<-length(coef(Q))
N<-nrow(DBM_no)
E<-resid(Q,type="pearson")
Dispersion<-sum(E^2)/(N-P)
Dispersion # 0.9361205, which is fine.

#METHOD 2
Overdisp_fun <- function (model) {rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp * rp)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(Q) # ratio = 0.9361205, which is fine. then we continue with poisson model

# Check normality of residuals
qqnorm(resid(Q,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Q,  type = "pearson")) # RESIDUALS ARE NORMAL. most data points on the plot fall between -1 and 1.

# Assess homogeneity of variance (visual), model 
plot(DBM_no$Time, resid(Q, type = "pearson"), main = "Treatment"). # homoenous. -1 and 1

Anova(Q)# Provides anova table with chi square and p values of the model
# "Time" is significant
HSD.test(Q, "Time", console = T)

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Q, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(Q, transform = "response"), pairwise ~ Treatment*Time)
emmeans(ref_grid(Q, transform = "response"), pairwise ~ Treatment|Time)


##################### NUMBER OF DBM PUPAE SET A


# model selection for NUMBER OF DBM PUPAE data, set A
M1 <- glm(Pupae  ~ Treatment * Time, data = DBM_no)
N2 <- glmmTMB(Pupae  ~ Treatment * Time, data = DBM_no)
P1<- glm.nb(Pupae  ~ Treatment + Time, data = DBM_no)
Q1 <- glm(Pupae  ~ Treatment * Time, family = poisson, data = DBM_no)
N3 <- glm(Pupae  ~ Treatment*Time, family = quasipoisson, data = DBM_no)
AIC(M1, N2, P1, Q1, N3) # Q1 is the best, NO OVERDISPERSION , AIC ALSO GOOD

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(Pupae~Treatment,data=DBM_no)) # run this and hit next command to get plot

#fit selected model: Ma3d
Q1 <- glm(Pupae ~ Treatment * Time, family = poisson, data = DBM_no)

#check for dispersion
# METHOD 1
P<-length(coef(Q1))
N<-nrow(DBM_no)
E<-resid(Q1,type="pearson")
Dispersion<-sum(E^2)/(N-P)
Dispersion # 0.9601945, which is fine.

#METHOD 2
Overdisp_fun <- function (model) {rdf <- df.residual(model)
rp <- residuals(model, type = "pearson")
Pearson.chisq <- sum(rp * rp)
prat <- Pearson.chisq / rdf
pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(Q1) # ratio = 0.9601945, which is fine. then we continue with poisson model

# Check normality of residuals
qqnorm(resid(Q1,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Q1,  type = "pearson")) # RESIDUALS ARE NORMAL: Most data oints are found between -1 and 1 in the plot, which is fine

# Assess homogeneity of variance (visual), model 
plot(DBM_no$Time, resid(Q1, type = "pearson"), main = "Treatment") # most data oints are found between -1 and 1 in the plot, which is fine

Anova(Q1)# Provides analysis of deviance table with chi square and p values of the model

# check if data is stored as factors or not

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Q1, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(Q1, transform = "response"), pairwise ~ Treatment*Time)
emmeans(ref_grid(Q1, transform = "response"), pairwise ~ Treatment|Time)


######## DBM DAMGAE SCORES

# model selection for NUMBER OF DBM PUPAE data, set A
M2 <- glm(score ~ Treatment, data = DBM_damage)
N4 <- glmmTMB(score ~ Treatment, data = DBM_damage)
P2<- glm.nb(score ~ Treatment+Tray, data = DBM_damage)
Q2 <- glm(score ~ Treatment, family = poisson, data = DBM_damage)
N5 <- glm(score ~ Treatment, family=quasipoisson, data = DBM_damage)
N6 <- glm(score ~ Treatment, family = Gamma, data = DBM_damage)
AIC(M2, N4, P2, Q2, N5, N6) # N6 is the best, NO OVERDISPERSION , AIC ALSO GOOD

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. This shows homogeneity: most points fall between -2 and 2 
plot(lm(score~Treatment,data=DBM_damage)) # run this and hit next command to get plot

#fit selected model: Ma3d
N6 <- glm(score ~ Treatment, family = Gamma, data = DBM_damage)

# Check normality of residuals
qqnorm(resid(N6,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(N6,  type = "pearson")) # RESIDUALS ARE NORMAL: Most data oints are found between -2 and 2 in the plot, which is fine

Anova(N6)# Provides analysis of deviance table with chi square and p values of the model
# nO SIGNIFICANTEFFECT OF TREATMENT. Therefore, no need for posthoc test


## BARPLOT, FEEDING DAMAGE, DBM

library(plotrix)
summary.DBMDAMAGE<-group_by(DBM_damage, Treatment)  %>% 
  summarise(mean=mean(score),std_error=std.error(score)) %>% arrange(desc(mean))
summarised_DBM_damage<-summarise_each(DBM_damage,funs(mean=mean,std_error=std.error))
View(summary.DBMDAMAGE)

ano4=aov(score ~ Treatment, data=DBM_damage)
summary(ano4)

#tukey's test

T4<-TukeyHSD(ano4)
print(T4)

# compact letter display
library(multcompView)
T.cld4 = multcompLetters4(ano4,T4)
print(T.cld4)

cld4<-as.data.frame.list(T.cld4$'Treatment')
summary.DBMDAMAGE$T4=cld4$Letters
View(summary.DBMDAMAGE)

#Barplots
ggplot(summary.DBMDAMAGE, aes(Treatment, mean, fill=Treatment, color=Treatment))+
  geom_bar(stat = "identity", position = "dodge",alpha=0.5)+
  geom_errorbar(aes(ymin=mean-std_error, ymax=mean+std_error),
                position=position_dodge(0.9), width=0.25, show.legend = FALSE)+
  labs(x="Soil amendment",y="Feeding damage score (1-7)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #leave legend outside
    ylim(0,7)+ # sets y axis to prevent cutting off error bars
  scale_color_brewer(palette = "Dark2", "red", "blue")+
  scale_fill_brewer(palette = "Dark2") 
ggsave("barplot.png", width = 4.5, height = 3, dpi = 1000)# SAVE PLOT


######## Total DBM (larvae+pupae collected on day 7 after inoculation) collected when more than 50% of surviving larvae had pupated.

# save data in appropriate format
str(DBM_tot)# view data structure, factors,etc.
DBM_tot$Treatment <- as.factor(DBM_tot$Treatment) # save Treatment as a factor
DBM_tot$rep <- as.factor(DBM_tot$rep)

# model selection for NUMBER OF DBM LARVAE + PUPAE data, set A
M2a <- glm(tot_DBM ~ Treatment, data = DBM_tot)
N4a <- glmmTMB(tot_DBM ~ Treatment, data = DBM_tot)
P2a<- glm.nb(tot_DBM ~ Treatment, data = DBM_tot)
Q2a <- glm(tot_DBM ~ Treatment, family = poisson, data = DBM_tot)
N5a <- glm(tot_DBM ~ Treatment, family=quasipoisson, data = DBM_tot)
AIC(M2a, N4a, P2a, Q2a, N5a) # Q2a is the best, NO OVERDISPERSION , AIC ALSO GOOD

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. This shows homogeneity: most points fall between -1 and 1 
plot(lm(tot_DBM~Treatment,data=DBM_tot)) # run this and hit next command to get plot

#fit selected model: Ma3d
Q2a <- glm(tot_DBM ~ Treatment, family = poisson, data = DBM_tot)

#check for dispersion

P<-length(coef(Q2a))
N<-nrow(DBM_tot)
E<-resid(Q2a,type="pearson")
Dispersion<-sum(E^2)/(N-P)
Dispersion # 1.446527, < 1.9  which is fine.

# Check normality of residuals
qqnorm(resid(Q2a,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Q2a,  type = "pearson")) # RESIDUALS ARE NORMAL: Most data points are found between -2 and 2 in the plot, which is fine

# Get mean effects via anova function in car package 
Anova(Q2a)

# posthoc using HSD 
HSD.test(Q2a, "Treatment", console = T)


################## SET C DATA, INCUBATED FRASS 
#(Frass incubated in plastic lunch boxes under greenhouse conditons)

growth_incub<-read.table("growth_incub_frass.txt", header = TRUE)# incubated frass in soil

# STORE DATA AS FACTOR (Treatment)
growth_incub$Treatment<-as.factor(growth_incub$Treatment)
growth_incub$Plant_ID <- as.factor(growth_incub$Plant_ID)
growth_incub$Time <- as.factor(growth_incub$Time)
growth_incub$leaves_c <- as.numeric(growth_incub$leaves_c)# "c" means data set c = incubated frass expt
growth_incub$area_c <- as.numeric(growth_incub$area_c)
str(growth_incub)# see data structure
summary(growth_incub)

## Number of leaves_C: Frass incubated in soil before seed planting

# select model
afit<-glmer(leaves_c~Treatment*Time + (1|Plant_ID), family = poisson, data=growth_incub)
mod1 <- glmer(leaves_c ~ Treatment* Time + (1|Plant_ID), family = Gamma, data = growth_incub)# Fit a Gamma model
cfit3 <- glmer.nb(leaves_c ~ Treatment*Time + (1|Plant_ID), data = growth_incub)
dfit4<- glmer(leaves_c ~ Treatment*Time + (1|Plant_ID), data = growth_incub)
AIC(afit, mod1, cfit3, dfit4) # mod1, gamma model is best

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. 
plot(lm(leaves_c~Treatment,data=growth_incub))

# Check normality of residuals
qqnorm(resid(mod1,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(mod1,  type = "pearson"))

# Get mean effects via anova function in car package 
Anova(mod1) # no significant difference

# try emmeans method
Back_Trans <- ref_grid(mod1, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(mod1, regrid = "response"), pairwise ~ Treatment * Time)
emmeans(ref_grid(mod1, regrid = "response"), pairwise ~ Treatment|Time)
emmeans(ref_grid(mod1, regrid = "response"), pairwise ~ Time|Treatment)
# No significant effect of treatment, but there is an effect of Time confirmed.

# Plot leaves_c
# INSTALL and load package: Hmisc
L<- position_dodge(width=0.2)
Lf <- ggplot(data=growth_incub,mapping=aes(x=Time,y=leaves_c,color=Treatment)) +  
stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.2,position=L) + 
stat_summary(fun=mean,geom="line",aes(group=Treatment),position=L) +  
stat_summary(fun=mean,geom="point",position=L) + theme_bw ()+  
labs(x="Time (days)",y="Number of leaves") + 
theme(legend.position = c(0.2, 0.75), 
legend.background = element_rect(fill = "white"))

# Modify legend titles
Lf + labs(color = "Soil amendment") # if you change "color" to "fill", it won't work       
#save
ggsave("ggplot_Lf.jpeg", width = 4, height = 3, dpi = 1000)


####leaf Area_C: Frass incubated in soil before seed planting
# select model
fit5<-glmmTMB(area_c~Treatment*Time + (1|Plant_ID), data=growth_incub)
fit6 <- glmer(area_c ~ Treatment*Time + (1|Plant_ID), family = Gamma, data = growth_incub)# Fit a Gamma model
fit7 <- glmer.nb(area_c ~ Treatment*Time + (1|Plant_ID), data = growth_incub)
mod2 <- glmer(area_c ~ Treatment*Time + (1|Plant_ID), data = growth_incub)
AIC(fit5, fit6, fit7, mod2) # fit6, gamma model is best

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. 
plot(lm(Area_c~Treatment,data=growth_incub))

# Check normality of residuals
qqnorm(resid(mod2,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(mod2,  type = "pearson")) # Normal, most points on the plot within -1 and 1

# Get mean effects via anova function in car package 
Anova(mod2) # significant difference

# emmeans method
Back_Trans <- ref_grid(mod2, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(mod2, regrid = "response"), pairwise ~ Treatment*Time)
emmeans(ref_grid(mod2, regrid = "response"), pairwise ~ Treatment|Time)
emmeans(ref_grid(mod2, regrid = "response"), pairwise ~ Time|Treatment)

# Plot area_c
# INSTALL and load package: Hmisc
A <- position_dodge(width=0.2)
Ar <- ggplot(data=growth_incub,mapping=aes(x=Time,y=area_c,color=Treatment)) +  
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.2,position=A) + 
  stat_summary(fun=mean,geom="line",aes(group=Treatment),position=A) +  
  stat_summary(fun=mean,geom="point",position=A) + theme_bw ()+  
  labs(x="Time (days)",y="Leaf area (sq cm)") + 
  theme(legend.position = c(0.2, 0.75), 
  legend.background = element_rect(fill = "white"))

# Modify legend titles
Ar + labs(color = "Soil amendment") # if you change "color" to "fill", it won't work       
#save
ggsave("ggplot_Ar.jpeg", width = 4, height = 3, dpi = 1000)


############# SET D DATA, COMPOSTED FRASS, in lunch boxes

Growth_comp<-read.table("growth_com_frass.txt", header=TRUE) # composted frass in lunch boxes
Leave_comp<-read.table("LeavesD_com_frass.txt", header=TRUE)

# Save data in appropriate format
Growth_comp$Treatment<-as.factor(Growth_comp$Treatment)
Growth_comp$Plant_ID<-as.factor(Growth_comp$Plant_ID)
Growth_comp$Time<-as.factor(Growth_comp$Time)
Growth_comp$Length_d <-as.numeric(Growth_comp$Length_d)
Growth_comp$Area_d<-as.numeric(Growth_comp$Area_d)
summary(Growth_comp)
str(Growth_comp)# see data structure

#LeavesD
# Save data in appropriate format
Leave_comp$Treatment<-as.factor(Leave_comp$Treatment)
Leave_comp$Plant_ID<-as.factor(Leave_comp$Plant_ID)
Leave_comp$Time<-as.factor(Leave_comp$Time)
Leave_comp$Leaves_d <-as.numeric(Leave_comp$Leaves_d)

summary(Leave_comp)
str(Leave_comp)# see data structure

######### Number of leaves_D, frass composted in lunch boxes

# select model
fita<-glmer(Leaves_d~Treatment*Time + (1|Plant_ID), family = poisson, data=Leave_comp)
fitb <- glmmTMB(Leaves_d ~ Treatment*Time + (1|Plant_ID), data = Leave_comp)
Mod3<- lmer(Leaves_d ~ Treatment*Time + (1|Plant_ID), data = Leave_comp)
AIC(fita, fitb, Mod3) # Mod3 and fitb, best

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. 
plot(lm(Leaves_d~Treatment*Time,data=Leave_comp)) # homogeneous, most points between -1 and 1

# Check normality of residuals
qqnorm(resid(Mod3,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Mod3,  type = "pearson")) # normal

# Get mean effects via anova function in car package 
Anova(Mod3)

# emmeans method
Back_Trans <- ref_grid(Mod3, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(Mod3, regrid = "response"), pairwise ~ Treatment*Time)
emmeans(ref_grid(Mod3, regrid = "response"), pairwise ~ Treatment|Time)

#Plot LeavesD: Composted frass
# INSTALL and load package: Hmisc
Lcom <- position_dodge(width=0.2)
LL <- ggplot(data=Leave_comp,mapping=aes(x=Time,y=Leaves_d,color=Treatment)) +  
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.2,position=Lcom) + 
  stat_summary(fun=mean,geom="line",aes(group=Treatment),position=Lcom) +  
  stat_summary(fun=mean,geom="point",position=Lcom) + theme_bw ()+  
  labs(x="Time (days)",y="Number of leaves") + 
  theme(legend.position = c(0.2, 0.75), 
  legend.background = element_rect(fill = "white"))

# Modify legend titles
LL + labs(color = "Soil amendment") # if you change "color" to "fill", it won't work       
#save
ggsave("ggplot_LL.jpeg", width = 4, height = 3, dpi = 1000)


### ##leaf Area_D: composted frass in lunch boxes: Growth_comp

# select model
fite<-glmmTMB(Area_d~Treatment*Time + (1|Plant_ID), data=Growth_comp)
Mod4<- glmer.nb(Area_d ~ Treatment*Time + (1|Plant_ID), data = Growth_comp)
Mod5<- glmer(Area_d ~ Treatment*Time + (1|Plant_ID), data = Growth_comp)
fit61 <- glmer(Area_d ~ Treatment*Time + (1|Plant_ID), family = Gamma, data = Growth_comp)
AIC(fite, Mod4, Mod5, fit61) # fith and fite, have low AICs but fitg fits better and also with a similar AIC
###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. 
plot(lm(Area_d~Treatment*Time,data=Growth_comp))

# Check normality of residuals
qqnorm(resid(Mod5 ,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Mod5,  type = "pearson")) # Normal, most points on the plot within -1 and 1

# Get mean effects via anova function in car package 
Anova(Mod5) 

# emmeans method
Back_Trans <- ref_grid(Mod5, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(Mod5, regrid = "response"), pairwise ~ Treatment*Time)
emmeans(ref_grid(Mod5, regrid = "response"), pairwise ~ Treatment|Time)

#Plot Leafarea D: Composted frass
# INSTALL and load package: Hmisc
La <- position_dodge(width=0.2)
LA <- ggplot(data=Growth_comp,mapping=aes(x=Time,y=Area_d,color=Treatment)) +  
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.2,position=La) + 
  stat_summary(fun=mean,geom="line",aes(group=Treatment),position=La) +  
  stat_summary(fun=mean,geom="point",position=La) + theme_bw ()+  
  labs(x="Time (days)",y="Leaf area (sq cm)") + 
  theme(legend.position = c(0.2, 0.75), 
        legend.background = element_rect(fill = "white"))

# Modify legend titles
LA + labs(color = "Soil amendment") # if you change "color" to "fill", it won't work       
#save
ggsave("ggplot_LA.jpeg", width = 4, height = 3, dpi = 1000)


############# SET E DATA, BRASSICA OLERACEAE

# STORE DATA AS FACTOR (Treatment)
growth_olera$Treatment<-as.factor(growth_olera$Treatment)
growth_olera$Time<-as.factor(growth_olera$Time)
growth_olera$Plant_ID<-as.factor(growth_olera$Plant_ID)
# see data structure
str(growth_olera)


######### Number of leaves_E
# select model
fita1<-glmer(leaves_e~Treatment*Time + (1|Plant_ID), family = poisson, data=growth_olera)
Mod5 <-glmmTMB(leaves_e ~ Treatment*Time + (1|Plant_ID), data = growth_olera)
fitc1 <-glmer.nb(leaves_e ~ Treatment*Time + (1|Plant_ID), data = growth_olera)
fitd1<- glmer(leaves_e ~ Treatment*Time + (1|Plant_ID), data = growth_olera)
fite1<- glmer(leaves_e ~ Treatment*Time + (1|Plant_ID), family = Gamma(link = "log"), data = growth_olera)
AIC(fita1, Mod5, fitc1, fitd1, fite1) # fitd1 and fitb1, best

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. 
plot(lm(leavese~Treatment*Time,data=growth_olera)) # homogeneous, most points between -1 and 1

# Check normality of residuals
qqnorm(resid(Mod5,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Mod5,  type = "pearson")) # normal, most points within -2 and 2

# Get mean effects via anova function in car package 
Anova(Mod5) # NO significant difference 


## ##leaf Area, SET E Data

# select model
fite1<-glmmTMB(Area_e~Treatment*Time + (1|Plant_ID), data=growth_olera)
Mod6  <- glmer.nb(Area_e ~ Treatment*Time + (1|Plant_ID), data =growth_olera)
xx <- glmer(Area_e ~ Treatment*Time + (1|Plant_ID), data = growth_olera)
fiti1 <- glmer(Area_e ~ Treatment*Time + (1|Plant_ID), family = Gamma(link = "log"), data = growth_olera)
AIC(fite1, Mod6, xx, fiti1) # fiti1, has a low AIC but fitg fits better and also with a similar AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. 
plot(lm(Area_e~Treatment*Time,data=growth_olera))

# Check normality of residuals
qqnorm(resid(Mod6, type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Mod6, type = "pearson")) # Normal, most points on the plot within -1 and 1

# Get mean effects via anova function in car package 
Anova(Mod6) # no significant difference for treatment


# Try emmeans method
Back_Trans <- ref_grid(Mod6, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(Mod6, transform = "response"), pairwise ~ Treatment*Time)
emmeans(ref_grid(Mod6, transform = "response"), pairwise ~ Treatment|Time)
# no significant difference confirmed, CON - MWF close to a sig. diff. but P>0.05 (0.0529)


###### PLOTS OF THE DATA( A, B, C, D, E)

### Boxplots
#Leaf Damage scores, DBM (plutella)
boxplot(score~Treatment,data = DBM_damage,ylab= "Feeding Score (1-7)",xlab="Waste stream",main="",ylim=c(0, 7.5),las=1,outpch=NA)
means <- tapply(DBM_damage$score, DBM_damage$Treatment,function(x) mean(x,na.rm=T))
# If your data contains missing values, you might want to replace the last argument of the tapply function (i.e mean) with "function(x) mean(x,na.rm=T)"
points(means,col="blue",pch="*", cex=2)
stripchart(score~Treatment,data = DBM_damage, vertical=TRUE, method="jitter",pch=21,col="blue",bg="yellow",add=TRUE)

#Leaf Mean number of DBM PUPAE(plutella)
boxplot(Pupae~Treatment,data = DBM_no,ylab= "Number of DBM pupae",xlab="Waste stream",main="",ylim=c(0,3),las=1,outpch=NA)
means <- tapply(DBM_no$Pupae, DBM_no$Treatment,function(x) mean(x,na.rm=T))
# If your data contains missing values, you might want to replace the last argument of the tapply function (i.e mean) with "function(x) mean(x,na.rm=T)"
points(means,col="blue",pch="*", cex=2)
stripchart(Pupae~Treatment,data = DBM_no, vertical=TRUE, method="jitter",pch=21,col="blue",bg="yellow",add=TRUE)

boxplot(Pupae~Treatment,data = DBM_no,ylab= "Number of DBM pupae",xlab="Waste stream",main="",ylim=c(0, 7.5),las=1,outpch=NA)
means <- tapply(DBM_no$tot_DBM, DBM_no$Treatment,function(x) mean(x,na.rm=T))
# If your data contains missing values, you might want to replace the last argument of the tapply function (i.e mean) with "function(x) mean(x,na.rm=T)"
points(means,col="blue",pch="*", cex=2)
stripchart(DBM_no~Treatment,data = DBM_tot, vertical=TRUE, method="jitter",pch=21,col="blue",bg="yellow",add=TRUE)




###################################### 04 May 2021#############
dir() 
setwd("E:/ANALYSIS DATA_2020-2021")
growth_incub<-read.table("growth_incub_frass.txt",header=TRUE)
attach(growth_incub)
growth_com<-read.table("growth_com_frass.txt",header=TRUE)
leaves_com<-read.table("LeavesD_com_frass.txt",header=TRUE)
growth_olera<-read.table("growth_rawfrass_oler.txt",header=TRUE)
DBM_pup<-read.table("Pupation_DBM_B.txt",header=TRUE)
damagescore<-read.table("DBM_damagescoreB.txt",header=TRUE)
Delia<-read.table("Deliarecovered_B.txt",header=TRUE)
Pupweigt<-read.table("Pupalweight_B.txt",header=TRUE)
flowerA<-read.table("FloweringA.txt",header=TRUE)


############# Number of Plutella (DBM) LARVAE + PUPA COLLECTED, SET B

# save data in appropriate format
DBM_pup$Treatment <- as.factor(DBM_pup$Treatment) # save Treatment as a factor
str(DBM_pup)# view data structure, factors,etc.

# model selection for NUMBER OF DBM LARVAE data, set A
M <- glm(pupae  ~ Treatment, data = DBM_pup)
N <- glmmTMB(pupae  ~ Treatment, data = DBM_pup)
P<- glm.nb(pupae~ Treatment, data = DBM_pup)
Q <- glmmTMB(pupae  ~ Treatment, family = poisson, data = DBM_pup)
N1 <- glm(pupae  ~ Treatment, family = quasipoisson, data = DBM_pup)
AIC(M, N, P, Q, N1) # M and N are the best, NO OVERDISPERSION , AIC ALSO GOOD

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(pupae~Treatment,data=DBM_pup)) # run this and hit next command to get plot


# Check normality of residuals
qqnorm(resid(M,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(M,  type = "pearson")) # RESIDUALS ARE NORMAL. most data points on the plot fall between -1 and 1.

Anova(M)# Provides anova table with chi square and p values of the model
# "Time" is significant
HSD.test(M, "Treatment", console = T)
#no significant difference


######## DBM DAMGAE SCORES, SET B

# model selection for NUMBER OF DBM PUPAE data, set A
M2 <- glm(score ~ Treatment, data = damagescore)
N4 <- glmmTMB(score ~ Treatment, data = damagescore)
P2<- glm.nb(score ~ Treatment, data = damagescore)
Q2 <- glm(score ~ Treatment, family = poisson, data = damagescore)
N5 <- glm(score ~ Treatment, family=quasipoisson, data = damagescore)
N6 <- glm(score ~ Treatment, family = Gamma, data = damagescore)
AIC(M2, N4, P2, Q2, N5, N6) # M2 is the best, 

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values. This shows homogeneity: most points fall between -2 and 2 
plot(lm(score~Treatment,data=damagescore)) # run this and hit next command to get plot

# Check normality of residuals
qqnorm(resid(M2,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(M2,  type = "pearson")) # RESIDUALS ARE NORMAL: Most data oints are found between -2 and 2 in the plot, which is fine

Anova(M2)# Provides analysis of deviance table with chi square and p values of the model
# SIGNIFICANTEFFECT OF TREATMENT. Therefore, need for posthoc test
Anova(M2)
Anova(N4)#SHOWS A SIGNIFCANT DIFFRERENCE
Anova(P2)
Anova(Q2)
Anova(N5)
Anova(N6)
 ### NOTE, EMMEANS FAILS TO SHOW THE SIGNIFICANT DIFFRENCE OBTAINED ABOVE, THEREFORE, M2 IS FAVOURED

## BARPLOT, FEEDING DAMAGE, DBM, SET B DATA

library(plotrix)
summary.DBMDAMAGEb<-group_by(damagescore, Treatment)  %>% 
  summarise(mean=mean(score),std_error=std.error(score)) %>% arrange(desc(mean))
summarised_damagescore<-summarise_each(damagescore,funs(mean=mean,std_error=std.error))
View(summary.DBMDAMAGEb)

ano4b=aov(score ~ Treatment, data=damagescore)
summary(ano4b)

#tukey's test

T4b<-TukeyHSD(ano4b)
print(T4b)

# compact letter display
library(multcompView)
T.cld4b = multcompLetters4(ano4b,T4b)
print(T.cld4b)

cld4b<-as.data.frame.list(T.cld4b$'Treatment')
summary.DBMDAMAGEb$T4b=cld4b$Letters
View(summary.DBMDAMAGEb)

#Barplots
ggplot(summary.DBMDAMAGEb, aes(Treatment, mean, fill=Treatment, color=Treatment))+
  geom_bar(stat = "identity", position = "dodge",alpha=0.5)+
  geom_errorbar(aes(ymin=mean-std_error, ymax=mean+std_error),
                position=position_dodge(0.9), width=0.25, show.legend = FALSE)+
  labs(x="Soil amendment",y="Feeding damage score (1-7)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #leave legend outside
  ylim(0,6)+ # sets y axis to prevent cutting off error bars
  scale_color_brewer(palette = "Dark2", "red", "blue")+
  scale_fill_brewer(palette = "Dark2") 
ggsave("barplot.png", width = 4.5, height = 3, dpi = 1000)# SAVE PLOT



################# Number of Delia larvae and pupae recovered, SET B ########

# save data in appropriate format
str(Delia)# view data structure, factors,etc.
Delia$Treatment <- as.factor(Delia$Treatment) # save Treatment as a factor
Delia$Rep <- as.factor(Delia$Rep)
Delia$Larvae_b <- as.numeric(Delia$Larvae_b)
Delia$Pupae_b <- as.numeric(Delia$Pupae_b)
Delia$tot_LarvaPupa <- as.numeric(Delia$tot_LarvaPupa)
attach(Delia)

Ma2b <- glm(tot_LarvaPupa ~ Treatment, data = Delia)
Anova (Ma2b)

# model selection for number of Delia larvae and pupae data set B
Ma2b <- glm(tot_LarvaPupa ~ Treatment, data = Delia)
Ma2c <- lm(tot_LarvaPupa ~ Treatment, data = Delia)
Ma2d <- glm(tot_LarvaPupa ~ Treatment, family = poisson (link = "log"), data = Delia)
Mb2 <- glmmTMB(tot_LarvaPupa ~ Treatment + (1|Rep), data = Delia)
Mc <- glmer(tot_LarvaPupa ~ Treatment + (1|Rep), data = Delia)
Md2a <- glm.nb(tot_LarvaPupa ~ Treatment, data = Delia)
Md2b <- glm.nb(tot_LarvaPupa ~ Treatment+Rep, data = Delia)
Me2 <- glmer(tot_LarvaPupa ~ Treatment + (1|Rep), family = poisson (link = "log"), data = Delia)

AIC(Ma2b, Ma2c, Ma2d, Mb2, Mc, Md2a, Md2b, Me2) # Ma2b, Mf2b, Ma2c are the best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(tot_LarvaPupa~Treatment,data=Delia)) # run this and hit next command to get plot


# Check normality of residuals
qqnorm(resid(Ma2d,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Ma2d,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(Delia$Treatment, resid(Ma2d, type = "pearson"), main = "Treatment")

Anova(Ma2b)
Anova(Ma2b)
Anova(Ma2c)
Anova(Mb2) # SHOWS no SIGNIFICANT DIFFERENCE
Anova(Mc)
Anova(Md2a)
Anova(Md2b)
Anova(Me2)

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(Mb2, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(Mb2, regrid = "response"), pairwise ~ Treatment)
emmeans(ref_grid(Mb2, regrid= "response"), pairwise ~ Treatment) 

## BARPLOT, NUMBER OF DELIA LARVAE PLUS PUPAE RECOVERED, SET B DATA
library(plotrix)
summary.deliaB<-group_by(Delia, Treatment)  %>% 
  summarise(mean=mean(tot_LarvaPupa),std_error=std.error(tot_LarvaPupa)) %>% arrange(desc(mean))
summarised_Delia<-summarise_each(Delia,funs(mean=mean,std_error=std.error))
View(summary.deliaB)

ano1b=aov(tot_LarvaPupa ~ Treatment, data=Delia)
summary(ano1b)

#tukey's test

T1b<-TukeyHSD(ano1b)
print(T1b)

# compact letter display
library(multcompView)
T.cld1b = multcompLetters4(ano1b,T1b)
print(T.cld1b)

cld1b<-as.data.frame.list(T.cld1b$'Treatment')
summary.deliaB$T1b=cld1b$Letters
View(summary.deliaB)
write_csv(summary.deliaB,"Delia_summary.csv") # save table as csv file

#Barplots
ggplot(summary.deliaB, aes(x=factor(Treatment), y = mean, fill=Treatment, color=Treatment))+
  geom_bar(stat = "identity", position = "dodge",alpha=0.5)+
  geom_errorbar(aes(ymin=mean-std_error, ymax=mean+std_error),
                position=position_dodge(0.9), width=0.25, show.legend = FALSE)+
  labs(x="Soil amendment",y="Number of CRF (larvae + pupae)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #leave legend outside
  geom_text(aes(label=T1b), position=position_dodge(0.9),show.legend = FALSE, size=4,
            vjust=-1.4, hjust=-2.0, color="gray25")+
  ylim(0,7)+ # sets y axis to prevent cutting off error bars
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") 
ggsave("barplot.png", width = 4.5, height = 3, dpi = 1000)# SAVE PLOT



###### Number of Delia PUPAE WEIGHT, SET B ##########

str(Pupweigt)# view data structure, factors,etc.
Pupweigt$Treatment <- as.factor(Pupweigt$Treatment) # save Treatment as a factor
Pupweigt$Rep <- as.factor(Pupweigt$Rep)

# model selection for number of Delia larvae and pupae data set A
M4a <- glmmTMB(Delia_b ~ Treatment, data = Pupweigt)
M4b <- glm(Delia_b ~ Treatment, data = Pupweigt)
M4c <- glm(Delia_b ~ Treatment, family = Gamma, data = Pupweigt)
M4e <- glm.nb(Delia_b ~ Treatment, data = Pupweigt)
M4g <- glm.nb(Delia_b ~ Treatment, data = Pupweigt)
M4h <- glm(Delia_b ~ Treatment, family = poisson, data = Pupweigt)


AIC(M4a, M4b, M4c, M4e, M4g, M4h) # M4e and M4g are the same, best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(Delia_b~Treatment,data=Pupweigt)) # run this and hit next command to get plot

# Check normality of residuals
qqnorm(resid(M4a,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(M4a,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(Pupweigt$Treatment, resid(M4a, type = "pearson"), main = "Treatment")

Anova(M4a) 
Anova(M4b)
Anova(M4c)
Anova(M4d)
Anova(M4e) 
Anova(M4g) 
# SHOWS NO A SIGNIFICANT DIFFERENCE, m4a, m4b.best, AIC


################  plots, WEIGHT OF DELIA PUPAE, SET A DATA ###

# save data in appropriate format
str(Pupweigt)# view data structure, factors,etc.
Pupweigt$Treatment <- as.factor(Pupweigt$Treatment) # save Treatment as a factor
str(Pupweigt)

## weight of delia pupae recovered

ano2b=aov(Delia_b ~ Treatment, data=Pupweigt)
summary(ano2b)

#tukey's test

T2b<-TukeyHSD(ano2b)
print(T2b)

# compact letter display
library(multcompView)
T.cld2b = multcompLetters4(ano2b,T2b)
print(T.cld2b)

# Table with factors and third quantile(tk2 - Delia)
tk2b<-group_by(Pupweigt, Treatment)  %>% summarise(mean=mean(Delia_b), 
quant = quantile(Delia_b, probs = 0.75)) %>% arrange(desc(mean))

# Extract the compact letter display and add to the T table (tk2 - Delia)
cld2b<-as.data.frame.list(T.cld2b$Treatment)
tk2b$cld2b<-cld2b$Letters
print(tk2b)

## Boxplot - weight of Delia pupae, SET A DATA
pupawt <- ggplot(Pupweigt, aes(Treatment, Delia_b))+ 
  geom_boxplot(aes(fill = Treatment), show.legend = FALSE) +
  theme_minimal() + labs(x="Soil amendment", y = "Mean weight (mg) of CRF pupae")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())+
  geom_text(data=tk2b, aes(label=cld2b, x=Treatment, y=quant),# brings in the letters
            vjust=-1, hjust=-0.9, size=4.5) + geom_jitter(width = 0.03)
pupawt + ylim(0,13)

#Save  the final figure
ggsave("boxplot.png", width = 4, height = 3, dpi = 1000)


################# FLY EMERGENCE, DATA FOR SETS A and B (Trials 1 and 2)

setwd("E:/ANALYSIS DATA_2020-2021")
deliaflies<-read.table("Flyemergence_delia.txt",header=TRUE)

#see how many samples we have for each level of Trt
tapply(deliaflies$duration, INDEX=deliaflies$Trt, FUN=length) # design is not balanced (i.e. some groups have more samples)

# Check normality of residuals
qqnorm(resid(modfly,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(modfly,  type = "pearson")) # RESIDUALS ARE NORMAL

# further check for normality. skewness value below 0.5, we can consider the deviation from normality not big enough to transform the data
skewness(deliaflies$duration)

#MODELS
pois.mod <- glm(duration ~ Trt, data=deliaflies, family=c("poisson"))
pois.mod2 = glm(duration ~ Trt + Period, data=deliaflies, family=c("poisson"))
pois.mod3 = glm(duration ~ Trt + Period, data=deliaflies, family=c("quasipoisson"))
neg.mod <- glmmTMB(duration ~ Trt + (1|Period), data = deliaflies)# Fit a negative binomial model
neg.mod2 <- glm.nb(duration ~ Trt + Period, data = deliaflies)# Fit a negative binomial model

# select best model
AIC (pois.mod, pois.mod2, neg.mod, pois.mod3)# pois.mod2 is better, though underdispersed

# in poison, mean = variance to assume normal distribution
mean(deliaflies$duration)# CHECK MEAN OF DURATION
var(deliaflies$duration)  # CHECK VARIANCE OF DURATION
# Mean > variance

# check for dispersion
Overdisp_fun <- function (model) {rdf <- df.residual(model)
rp <- residuals(model, type = "pearson")
Pearson.chisq <- sum(rp * rp)
prat <- Pearson.chisq / rdf
pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(pois.mod2)  # data is underdispersed (ratio < 1)

# TEST FOR SIGNIFICANT DIFFERENCE among treatments
Anova(pois.mod2)
HSD.test(pois.mod2, "Trt", console = T)

# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(pois.mod2, transform = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(pois.mod2, transform = "response"), pairwise ~ Trt)
emmeans(ref_grid(pois.mod2, transform = "response"), pairwise ~ Trt)


################# FLower emergence, DATA FOR SETS A, B, C and D
setwd("E:/ANALYSIS DATA_2020-2021")


##FloweringA, SET A
flowerA <-read.table("FloweringA.txt",header=TRUE)
Mod7 <- glm(Days ~ Treatment, family = poisson, data = flowerA)
Mod8 <- glm(Days ~ Treatment, family = Gamma, data = flowerA)
Mod9 <- glm.nb(Days ~ Treatment, data = flowerA)
Mod10 <- glm(Days ~ Treatment, data = flowerA)
AIC (Mod7, Mod8, Mod9,Mod10 )

# check for dispersion
Overdisp_fun <- function (model) {rdf <- df.residual(model)
rp <- residuals(model, type = "pearson")
Pearson.chisq <- sum(rp * rp)
prat <- Pearson.chisq / rdf
pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(Mod7) #overdispersed, poisson

# Check normality of residuals
qqnorm(resid(Mod10,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Mod10,  type = "pearson")) # RESIDUALS ARE NORMAL

Anova(Mod10)# no sig. diff.

#PLOTS, create barplots with error bars in R using the ggpubr package
ggbarplot(flowerA, x = "Treatment", y = "Days", 
          add = c("mean_se"), fill = "#BF504D", ylab="Duration (days)", xlab="Soil amendment")

## FloweringB
flowerb <-read.table("FloweringB.txt",header=TRUE)
pois.mod5 <- glm(Days ~ Treatment, family = poisson, data = flowerb)
gamm.mod <- glm(Days ~ Treatment, family = Gamma, data = flowerb)
neg.mod4 <- glm.nb(Days ~ Treatment, data = flowerb)
glm.mod <- glm(Days ~ Treatment, data = flowerb)
AIC (pois.mod5, gamm.mod, neg.mod4,glm.mod )

# check for dispersion
Overdisp_fun <- function (model) {rdf <- df.residual(model)
rp <- residuals(model, type = "pearson")
Pearson.chisq <- sum(rp * rp)
prat <- Pearson.chisq / rdf
pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(pois.mod5) # check for overdisperse, poisson

# Check normality of residuals
qqnorm(resid(pois.mod5,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(pois.mod5,  type = "pearson")) # RESIDUALS ARE NORMAL

Anova(pois.mod5)

#PLOTS, create barplots with error bars in R using the ggpubr package
ggbarplot(flowerb, x = "Treatment", y = "Days", 
  add = c("mean_se"), fill = "#BF504D", ylab="Duration (days)", xlab="Soil amendment")

## FloweringC
flowerc <-read.table("FloweringC.txt",header=TRUE)
pois.mod6 <- glm(Days ~ Treatment, family = poisson, data = flowerc)
gamm.mod1 <- glm(Days ~ Treatment, family = Gamma, data = flowerc)
neg.mod5 <- glm.nb(Days ~ Treatment, data = flowerc)
glm.mod1 <- glm(Days ~ Treatment, data = flowerc)
AIC (pois.mod6, gamm.mod1, neg.mod5,glm.mod1 )

# check for dispersion
Overdisp_fun <- function (model) {rdf <- df.residual(model)
rp <- residuals(model, type = "pearson")
Pearson.chisq <- sum(rp * rp)
prat <- Pearson.chisq / rdf
pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(pois.mod6) # overdisperse, poisson

# Check normality of residuals
qqnorm(resid(pois.mod6,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(pois.mod6,  type = "pearson")) # RESIDUALS ARE NORMAL

Anova(pois.mod6)

#PLOTS, create barplots with error bars in R using the ggpubr package
ggbarplot(flowerc, x = "Treatment", y = "Days", 
          add = c("mean_se"), fill = "#BF504D", ylab="Duration (days)", xlab="Soil amendment")

## Floweringd
flowerd <-read.table("FloweringD.txt",header=TRUE)
pois.mod7 <- glm(Days ~ Treatment, family = poisson, data = flowerd)
gamm.mod2 <- glm(Days ~ Treatment, family = Gamma, data = flowerd)
neg.mod6 <- glm.nb(Days ~ Treatment, data = flowerd)
glm.mod2 <- glm(Days ~ Treatment, data = flowerd)
AIC (pois.mod7, gamm.mod2, neg.mod6, glm.mod2)

# Check normality of residuals
qqnorm(resid(glm.mod2,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(glm.mod2,  type = "pearson")) # RESIDUALS ARE NORMAL

Anova(glm.mod2)

#PLOTS, create barplots with error bars in R using the ggpubr package
ggbarplot(flowerd, x = "Treatment", y = "Days", position = position_dodge(0.8),
          add = c("mean_se"), fill = "#BF504D", space=1,  ylab="Duration (days)", xlab="Soil amendment")




#################################################################################
#################################################################################
################### COMBINED DATASETS_A & B   ###################################
              #####################################
## testing if trial (the one month interval between trials 1 & 2 or Set A & B) has an effect

setwd("~/OneDrive/R-WORKING DIRECTORY/ANALYSIS DATA_2020-2024")
setwd("C:/Users/chia002/OneDrive/R-WORKING DIRECTORY/ANALYSIS DATA_2020-2024")


#Leaf area (A & B)

areaAB<-read.table("Leafarea_trials_A_B.txt",header=TRUE)
str(areaAB)# view data structure, factors,etc.
areaAB$Treatt <- as.factor(areaAB$Treatt) # save Treatment as a factor
areaAB$Time <- as.factor(areaAB$Time) # save Time as a factor
areaAB$Plant_ID <- as.factor(areaAB$Plant_ID) # save Plant_ID as a factor
areaAB$Trial<- as.factor(areaAB$Trial) # save Trial as a factor
areaAB$Area<- as.numeric(areaAB$Area) # save Trial as a factor
areaAB$Length<- as.numeric (areaAB$Length) # save Trial as a factor
areaAB$Width<- as.numeric(areaAB$Width) # save Trial as a factor
str(areaAB)# view data structure, factors,etc.
summary(areaAB)
library(glmmTMB)

#fit model
f1 <- lmer(Area ~ Treatt * Time *Trial + (1|Plant_ID), data = areaAB)
f2<- glmmTMB(Area ~ Treatt * Time *Trial + (1|Plant_ID), data = areaAB)
f3 <- glmer(Area ~ Treatt * Time *Trial + (1|Plant_ID), family = Gamma(link = "log"), data = areaAB)# Fit a Gamma model
f4 <- glmmTMB(Area ~ Treatt * Time + Trial + (1|Plant_ID), data = areaAB)
AIC(f1, f2, f3, f4) # f3 is better

plot(lm(Area ~ Treatt,data=areaAB)) # run this and the next command to get plot
plot(lm(Area~Time,data=areaAB)) # homogeneity of residuals confirmed here
plot(lm(Area~Trial,data=areaAB))

library(car)
## tests for homogeneity (equality) of variances in data
leveneTest(Area ~ Treatt,data=areaAB) # P<0.05 = not homogenous
fligner.test(Area ~ Treatt,data=areaAB) # same
bartlett.test(Area ~ Treatt,data=areaAB) # same
# variances not homogeneous according to the tests but visual assessment shows 
 # homogeneity, then lets proceed with GLMM

# fit selected model
f2<- glmmTMB(Area ~ Treatt * Time *Trial + (1|Plant_ID), data = areaAB)
Anova(f2)

#emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(f2, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatt * trial * timepoints
HSD_test <- emmeans(ref_grid(f2, regrid = "response"), pairwise ~ Treatt * Time*Trial)
emmeans(ref_grid(f2, regrid = "response"), pairwise ~ Treatt * Time* Trial) 
emmeans(ref_grid(f2, regrid = "response"), pairwise ~ Treatt|Time|Trial) 
emmeans(ref_grid(f2, regrid = "response"), pairwise ~ Time|Treatt|Trial) 
emmeans(ref_grid(f2, regrid = "response"), pairwise ~ Trial|Treatt|Time) 

## Plot number of leafArea-combinedAB
library(ggplot2) # to plot

estimates <- emmeans(ref_grid(f2, regrid = "response"), pairwise ~ Treatt * Time*Trial)
plot_output <- as.data.frame(estimates$emmeans) # plot estimates using emmeans package
ggplot(plot_output, aes(x = Time, y = emmean, color = Treatt)) + facet_wrap(. ~ Trial)+
  geom_point(size = 0.8)+ geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
    width = 0.2) + geom_line(aes(group = Treatt)) + theme_bw ()+  
  labs(x="Time (days)",y="Leaf area (sq cm)") +
  theme(legend.position="bottom")
plot_output + labs(color = "Soil amendment")
ggsave("plot_output.jpeg", width = 5, height = 3, dpi = 1000)



## Number of leaves (A & B)

leavesAB<-read.table("Leaves_combinedAB.txt",header=TRUE) # READ IN DATA
str(leavesAB)# view data structure, factors,etc.
leavesAB$Treatt <- as.factor(leavesAB$Treatt) # save Treatment as a factor
leavesAB$Time <- as.factor(leavesAB$Time) # save Time as a factor
leavesAB$Plant_ID <- as.factor(leavesAB$Plant_ID) # save Plant_ID as a factor
leavesAB$Trial<- as.factor(leavesAB$Trial) # save Block as a factor
leavesAB$No_leaves<- as.numeric(leavesAB$No_leaves)

# model selection for leaves
lmer.mod1 <- lmer(No_leaves ~ Treatt * Time*Trial + (1|Plant_ID), data = leavesAB)
lmer.mod2 <- lmer(No_leaves ~ Treatt * Time+Trial + (1|Plant_ID), data = leavesAB)
glmmTMB.mod1 <- glmmTMB(No_leaves ~ Treatt * Time*Trial  + (1|Plant_ID), data = leavesAB)
glmmTMB.mod2 <- glmmTMB(No_leaves ~ Treatt * Time+Trial  + (1|Plant_ID), data = leavesAB)
gamma.mod1 <- glmer(No_leaves ~ Treatt * Time*Trial + (1|Plant_ID), family = Gamma(link = "log"), data = leavesAB)# Fit a Gamma model
gamma.mod2<- glmer(No_leaves ~ Treatt * Time+Trial + (1|Plant_ID), family = Gamma(link = "log"), data = leavesAB)# Fit a Gamma model
glmer.nb.mod1 <- glmer.nb(No_leaves ~ Treatt * Time*Trial + (1|Plant_ID), data = leavesAB)
glmer.nb.mod2 <- glmer.nb(No_leaves ~ Treatt * Time+Trial + (1|Plant_ID), data = leavesAB)
pois.mod1 <- glmer(No_leaves ~ Treatt * Time*Trial + (1|Plant_ID), family = poisson (link = "log"), data = leavesAB)
pois.mod2 <- glmer(No_leaves ~ Treatt * Time+Trial + (1|Plant_ID), family = poisson (link = "log"), data = leavesAB)

AIC(lmer.mod1, lmer.mod2, glmmTMB.mod1, glmmTMB.mod2, gamma.mod1, gamma.mod2, glmer.nb.mod1, glmer.nb.mod2, pois.mod1, pois.mod2) # glmmTMB.mod2 is the best, least AIC

###  checks for Homogeneity of variance (visual) in the data
#plotting residuals vs. fitted values
plot(lm(No_leaves~Treatment,data=leavesAB)) # run this and the next command to get plot
plot(lm(No_leaves~Time,data=leavesAB)) # homogeneity of residuals confirmed here

library(car)
## tests for homogeneity (equality) of variances in data
leveneTest(No_leaves ~ Treatment, data = leavesAB)
fligner.test(No_leaves ~ Treatment, data = leavesAB)
bartlett.test(No_leaves ~ Treatment, data = leavesAB) 
# variances not homogeneous according to the tests but visual assessment shows homogeneity, then lets proceed with GLMM

#fit selected model
glmmTMB.mod2 <- glmmTMB(No_leaves ~ Treatt * Time * Trial  + (1|Plant_ID), data = leavesAB)

# Check normality of residuals
qqnorm(resid(glmmTMB.mod2,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(glmmTMB.mod2,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(leavesAB$Treatment, resid(glmmTMB.mod2, type = "pearson"), main = "Treatment")
plot(leavesAB$Time, resid(glmmTMB.mod2, type = "pearson"), main = "Time")
plot(leavesAB$Trial, resid(glmmTMB.mod2, type = "pearson"), main = "Block")

Anova(glmmTMB.mod2)


emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(glmmTMB.mod2, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments * timepoints
HSD_test <- emmeans(ref_grid(glmmTMB.mod2, regrid = "response"), pairwise ~ Treatt*Time+Trial)
emmeans(ref_grid(glmmTMB.mod2, regrid = "response"), pairwise ~ Treatt|Time+Trial) 
emmeans(ref_grid(glmmTMB.mod2, regrid = "response"), pairwise ~   Time|Treatt|Trial) 
emmeans(ref_grid(glmmTMB.mod2, regrid = "response"), pairwise ~   Trial|Treatt|Time) 

## Plot number of leaves-combinedAB
library(ggplot2) # to plot

estimates <- emmeans(ref_grid(glmmTMB.mod2, regrid = "response"), pairwise ~ Treatt * Time+Trial)
plot_output_Lf <- as.data.frame(estimates$emmeans) # plot estimates using emmeams package
ggplot(plot_output_Lf, aes(x = Time, y = emmean, color = Treatt)) + facet_wrap(. ~ Trial)+
  geom_point(size = .8)+ geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
    width = 0.2) + geom_line(aes(group = Treatt)) + theme_bw ()+  
  labs(x="Time (days)",y="Number of leaves")+
theme(legend.position="bottom")

# Modify legend titles
plot_output_Lf + labs(color = "Soil amendment")
ggsave("plot_output_Lf.jpeg", width = 5, height = 3, dpi = 1000)


#### Larval development and survival: Number of Delia radicum retrieved from pot soil

setwd("C:/Users/chia002/OneDrive/R-WORKING DIRECTORY/ANALYSIS DATA_2020-2024")

delia<-read.table("Delia_recovered_trials_A_B.txt",header=TRUE)

# Safe data in appropriate forms
delia$Treatt <- as.factor(delia$Treatt) # save Treatt as a factor
delia$Trial <- as.factor(delia$Trial)
delia$Plant_ID <- as.factor(delia$Plant_ID)
delia$Larvae <- as.numeric(delia$Larvae)
delia$Pupae <- as.numeric(delia$Pupae)
delia$tot_Larv_Pupa <- as.numeric(delia$tot_Larv_Pupa)
str(delia)# view data structure, factors,etc.
summary(delia)

#fit model
d1 <- glm(tot_Larv_Pupa ~ Treatt+Trial, data = delia)
d1a <- glm(tot_Larv_Pupa ~ Treatt*Trial, data = delia)
AIC (d1, d1a)
Anova(d1)


# Check normality of residuals
qqnorm(resid(d1,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(d1,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(delia$Treatt, resid(d1, type = "pearson"), main = "Treatt") # # appear not to be homogeneous
plot(delia$Trial, resid(d1, type = "pearson"), main = "Trial") # homogeneous

# no sig difference, therefore, no post hoc

#Plot number of delia recovered

library(plotrix) # needed to view summary of data
library(tidyverse)
summary.delia<-group_by(delia, Treatt)  %>% 
  summarise(mean=mean(tot_Larv_Pupa),std_error=std.error(tot_Larv_Pupa)) %>% arrange(desc(mean))
summarised_delia<-summarise_each(delia,funs(mean=mean,std_error=std.error))
View(summary.delia)

ano_delia=aov(tot_Larv_Pupa ~ Treatt, data=delia)
summary(ano_delia)

#tukey's test

T<-TukeyHSD(ano_delia)
print(T)

# compact letter display
library(multcompView)
T.cld = multcompLetters4(ano_delia,T)
print(T.cld)

cld<-as.data.frame.list(T.cld$'Treatt')
summary.delia$T=cld$Letters
View(summary.delia)
write_csv(summary.delia,"delia_summary.csv") # save table as csv file

#Barplots
ggplot(summary.delia, aes(x=factor(Treatt), y = mean, fill=Treatt, color=Treatt))+
  geom_bar(stat = "identity", position = "dodge",alpha=0.5)+
  geom_errorbar(aes(ymin=mean-std_error, ymax=mean+std_error),
                position=position_dodge(0.9), width=0.25, show.legend = FALSE)+
  labs(x="Soil amendment",y="Number of CRF (larvae + pupae)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #leave legend outside
  geom_text(aes(label=T), position=position_dodge(0.9),show.legend = FALSE, size=4,
            vjust=-1.4, hjust=-2.0, color="gray25")+
  ylim(0,8)+ # sets y axis to prevent cutting off error bars
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") 
ggsave("barplot.jpeg", width = 4.5, height = 3, dpi = 1000)# SAVE PLOT

detach(delia) # to free the R environment of similar headers in next file/scipt


#### Delia biomass: fresh weight of pupae retrieved from pot soil

setwd("C:/Users/chia002/OneDrive/R-WORKING DIRECTORY/ANALYSIS DATA_2020-2024")
delia_wt<-read.table("Delia_Pupalweight_trials_A_B.txt",header=TRUE)

# save data in appropriate format
delia_wt$Treatt <- as.factor(delia_wt$Treatt) # save Treatt as a factor
delia_wt$Trial <- as.factor(delia_wt$Trial)
delia_wt$sample <- as.factor(delia_wt$sample)
delia_wt$pupalweight <- as.numeric(delia_wt$pupalweight)
str(delia_wt)# view data structure, factors,etc.
summary(delia_wt)

#fit model
d2 <- glm(pupalweight ~ Treatt * Trial, family = Gamma (link = "inverse"), data = delia_wt)
d5 <- glm(pupalweight ~ Treatt*Trial, data = delia_wt)
Anova(d5) # anova function in "car" package, provides overall p value to indicate significance

AIC(d2, d5)# d5 is better than d2 

# Check normality of residuals
qqnorm(resid(d5,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(d5,  type = "pearson")) # RESIDUALS ARE NORMAL

# Assess homogeneity of variance (visual), model 
plot(delia_wt$Treatt, resid(d5, type = "pearson"), main = "Treatt") # homogeneous
plot(delia_wt$Trial, resid(d5, type = "pearson"), main = "Trial") # homogeneous

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(d5, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments
HSD_test <- emmeans(ref_grid(d5, regrid = "response"), pairwise ~ Treatt+Trial)
emmeans(ref_grid(d5, regrid = "response"), pairwise ~ Treatt|Trial)
emmeans(ref_grid(d5, regrid = "response"), pairwise ~ Trial|Treatt)

# plots, WEIGHT OF DELIA PUPAE

## Boxplot - weight of Delia pupae
pupwt <- ggplot(delia_wt, aes(Treatt, pupalweight))+ 
  geom_boxplot(aes(fill = Trial), notch = TRUE, position = position_dodge(0.9), 
               show.legend = TRUE) +
  theme_minimal() + labs(x="Soil amendment", y = "Wet weight (mg) of CRF pupae")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())+ 
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))
pupwt + ylim(0,17)
#Save  the final figure
ggsave("boxplot_pupwt.jpeg", width = 4, height = 3, dpi = 1000)



#### Plutella: Larval development - Number of Plutella xylostella pupae


DBM<-read.table("DBM_pupaeRecovered_trials_A_B.txt",header=TRUE)
DBM1<-read.table("DBM_pupaeRecovered_trials_A_B_days_6_10_added_zeros.txt",header=TRUE)

# save data in appropriate format
DBM$Treatt <- as.factor(DBM$Treatt) # save Treatt as a factor
DBM$Trial <- as.factor(DBM$Trial)
DBM$Date <- as.factor(DBM$Date)
DBM$Time <- as.factor(DBM$Time)
DBM$Plant_ID <- as.factor(DBM$Plant_ID)
DBM$Pupae <- as.numeric(DBM$Pupae)
str(DBM)# view data structure, factors,etc.
summary(DBM)
detach(DBM1)

p1 <- glm(Pupae ~ Treatt * Time+Trial, family = poisson, data = DBM1)
p1a <- glm(Pupae ~ Treatt * Time*Trial, family = poisson, data = DBM)
p1b <- glm(Pupae ~ Treatt * Time, family = poisson, data = DBM1)
AIC (p1, p1a, p1b)
Anova(p1a)
      
#check for dispersion
# METHOD 1
P<-length(coef(p1a))
N<-nrow(DBM)
E<-resid(p1a,type="pearson")
Dispersion<-sum(E^2)/(N-P)
Dispersion # 0.9207607, which is fine.

#METHOD 2
Overdisp_fun <- function (model) {rdf <- df.residual(model)
rp <- residuals(model, type = "pearson")
Pearson.chisq <- sum(rp * rp)
prat <- Pearson.chisq / rdf
pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}
Overdisp_fun(p1) # ratio = 1.009416, which is fine. then we continue with poisson model

# Check normality of residuals
qqnorm(resid(p1a,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(p1a,  type = "pearson")) # RESIDUALS ARE NORMAL: 
#most data points are found between -1 and 1 in the plot, which is fine

# Assess homogeneity of variance (visual), model 
plot(DBM$Treatt, resid(p1a, type = "pearson"), main = "Treatt")
plot(DBM$Time, resid(p1a, type = "pearson"), main = "Time") 
# most data points are found between -1 and 1 in the plot, which is fine

Anova(p1a)# Provides analysis of deviance table with chi square and p values of the model

emmeans # estimated marginal means. 
library (emmeans) # load the emmeans package
# Recalculate estimations from log link function to the original response variable
Back_Trans <- ref_grid(p1a, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatments*Time points
HSD_test <- emmeans(ref_grid(p1a, regrid = "response"), pairwise ~ Treatt*Time)
emmeans(ref_grid(p1a, regrid = "response"), pairwise ~ Treatt|Time)
emmeans(ref_grid(p1a, regrid = "response"), pairwise ~ Time|Trial|Treatt)
emmeans(ref_grid(p1a, regrid = "response"), pairwise ~ Trial|Treatt|Time)



#Plot number of DMB recovered over time

library(dplyr) 
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

# Compute the summary statistics of "Pupae" grouped by "Time" and "Treatt":
DBM.summary <- DBM %>%
  group_by(Time, Treatt) %>%
  summarise(
    sd = sd(Pupae, na.rm = TRUE),
    Pupae = mean(Pupae))
DBM.summary

# Line plot + error bars
P <- ggplot(DBM.summary, aes(Time, Pupae)) +
  geom_line(aes(linetype = Treatt, group = Treatt))+
  geom_point()+
  geom_errorbar(
    aes(ymin = Pupae-sd, ymax = Pupae+sd, group = Treatt),
    width = 0.2)+ labs(x="Time (days)", y = "Number of pupae")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())

# Modify legend titles
P + labs(fill = "Soil amendment") # Did not work yet

ggsave("ggplot_P.jpeg", width = 4, height = 3, dpi = 1000)



#### Plant damage by larvae of P. xylostella

damagescore<-read.table("DBM_leafdamagescore_trials_A_B.txt",header=TRUE)

# save data in appropriate format
damagescore$Treatt <- as.factor(damagescore$Treatt) # save Treatt as a factor
damagescore$Trial <- as.factor(damagescore$Trial)
damagescore$score <- as.numeric(damagescore$score)
damagescore$Plant_ID <- as.factor(damagescore$Plant_ID)
str(damagescore)# view data structure, factors,etc.
summary(damagescore)

#fit model
Dscore <- glm(score ~ Treatt+Trial, data=damagescore) # Trial has no effect
Dscore1 <- glm(score ~ Treatt, data=damagescore) # better, lower AIC value
AIC(Dscore, Dscore1)
Anova(Dscore1)

#post hoc
Back_Trans <- ref_grid(Dscore1, regrid = "response") 

# Combine ref_grid with emmeans and a post-hoc test among all treatment points
HSD_test <- emmeans(ref_grid(Dscore1, regrid = "response"), pairwise ~ Treatt)
emmeans(ref_grid(Dscore1, regrid = "response"), pairwise ~ Treatt)

# PLOT, leave Damage, DBM larvae
ds <- ggplot(damagescore, aes(Treatt, score))+ 
  geom_boxplot(aes(fill = Treatt), notch = FALSE, position = position_dodge(0.9), 
               show.legend = TRUE) +
  stat_summary(fun=mean, geom="point", shape=19, size=5, color="red", fill="red")+
  theme_minimal() + labs(x="Soil amendment", y = "Feeding damage score (1-7)")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())
ds + ylim(2,8)+ theme(legend.position="bottom")

#Save  the final figure
ggsave("boxplot_ds.jpeg", width = 4, height = 3, dpi = 1000)



#### Flower development time for B. rapa plants in the greenhouse

Flower<-read.table("Flowering_trials_A_B.txt",header=TRUE)

# save data in appropriate format
Flower$ Treatment <- as.factor(Flower$ Treatment) # save  Treatment as a factor
Flower$Trial <- as.factor(Flower$Trial)
Flower$Days <- as.numeric(Flower$Days)
Flower$Plant_ID <- as.factor(Flower$Plant_ID)
str(Flower)# view data structure, factors,etc.

Flwer_ABCD$Frass <- as.factor(Flwer_ABCD$ Frass) # save  Treatment as a factor
Flwer_ABCD$tt <- as.factor(Flwer_ABCD$tt)
Flwer_ABCD$DaysFF <- as.numeric(Flwer_ABCD$DaysFF)
str(Flwer_ABCD)# view data structure, factors,etc.
summary(Flwer_ABCD)

#select model
Md7 <- glm(Days ~ Treatment+Trial, family = poisson, data = Flower) #shows a sig. diff.
Md8 <- glm(Days ~ Treatment+Trial, family = Gamma, data = Flower)
Md9 <- lm(Days ~ Treatment+Trial, data = Flower)# best, 
Md10 <- glm(Days ~ Treatment*Trial, data = Flower) # best, 
AIC (Md7, Md8, Md9,Md10) # All similar, but Md7 is used here, it produces asig. diff
AIC (Md7,Md10)
# Check normality of residuals
qqnorm(resid(Md10,  type = "pearson"), ylab="Std Residuals", xlab="Normal Scores")
qqline(resid(Md10,  type = "pearson")) # RESIDUALS ARE NORMAL ??

Anova(Md10)# sig. diff. for "Trial" but not for "Treatment"
Anova(Md7) # better with low AIC 

#Post hoc
HSD_test <- emmeans(ref_grid(Md7, regrid = "response"), pairwise ~ Treatment+Trial)
emmeans(ref_grid(Md7, regrid = "response"), pairwise ~ Treatment|Trial)
emmeans(ref_grid(Md7, regrid = "response"), pairwise ~ Trial|Treatment)

#PLOTS, Flower development time

flwer <- ggplot(Flower, aes(Treatment, Days))+ 
  geom_boxplot(aes(fill = Trial), notch = FALSE, position = position_dodge(0.9), 
               show.legend = TRUE) +
  theme_minimal() + labs(x="Soil amendment", y = "Flowering time (days)")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())+ 
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))
flwer + ylim (20, 100) + theme(legend.position="right")

#Save  the final figure
ggsave("boxplot_flwer.jpeg", width = 4, height = 3, dpi = 1000)


#### Flowering, A, B, C, D
Flwer_ABCD<-read.table("Flowering_trials_A_B_C_D.txt",header=TRUE)
Flw_ABCD<-read.table("Flowering_trials_A_B_C_D_no control.txt",header=TRUE) # no control

Flw_ABCD$Frass <- as.factor(Flw_ABCD$ Frass) # save Treatment as a factor
Flw_ABCD$tt <- as.factor(Flw_ABCD$tt)
Flw_ABCD$DaysFF <- as.numeric(Flw_ABCD$DaysFF)
str(Flw_ABCD)# view data structure, factors,etc.
summary(Flw_ABCD)


Mf <- glm(DaysFF ~ Frass+tt, family = poisson, data = Flw_ABCD)
Anova(Mf)
MF <- glm(DaysFF ~ Frass*tt, family = poisson, data = Flw_ABCD)
Anova(MF)
AIC(Mf, MF)# Mf better with low AIC, poor post hoc, hence MF is prepared

Anova(MF) 

#Post hoc
HSD_test <- emmeans(ref_grid(MF, regrid = "response"), pairwise ~ Frass + tt)
emmeans(ref_grid(MF, regrid = "response"), pairwise ~ Frass|tt)
emmeans(ref_grid(MF, regrid = "response"), pairwise ~ tt|Frass)

#Plot Flowering, A, B, C, D ### not that in plotting, the control was included, 
   #but left out in the analysis
flw <- ggplot(Flwer_ABCD, aes(Frass, DaysFF))+ 
  geom_boxplot(aes(fill = tt), notch = FALSE, position = position_dodge(0.9), 
               show.legend = TRUE) +
  theme_minimal() + labs(x="Soil amendment", y = "Flowering time (days)")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank()) + 
  scale_fill_discrete(labels=c('No frass (control)', 'Composted frass', 'Incubated frass', 'Raw frass'))
flw + ylim (0, 100) + theme(legend.position="right") 

# Modify legend titles
flw + labs(fill = "Status of frass used") + ylim (20, 100) + theme(legend.position="right") 


#Save  the final figure
ggsave("boxplot_flw.jpeg", width = 4, height = 2.5, dpi = 1000)


#### Fly Emergence from Delia radicum pupae retrieved from pot soil and incubated at 20 degrees, in KLIMA WUR
setwd("~/OneDrive/R-WORKING DIRECTORY/ANALYSIS DATA_2020-2024")
DFly<-read.table("FlyEmergence_Delia_trials_A_B.txt",header=TRUE)

# save data in appropriate format
DFly$treatt <- as.factor(DFly$treatt) # save Treatt as a factor
DFly$trial <- as.factor(DFly$trial)
DFly$time <- as.numeric(DFly$time)
DFly$pupae <- as.numeric(DFly$pupae)
DFly$flies <- as.numeric(DFly$flies)
DFly$percent <- as.numeric(DFly$percent)
str(DFly)# view data structure, factors,etc.
summary(DFly)

#see how many samples we have for each level of Treatt

tapply(DFly$percent, INDEX=DFly$treatt, FUN=length) # design is not balanced (i.e. some groups have more samples)???

#Data 
# Calculate percentage emergence of D. radicum

percent.em <- DFly$flies/DFly$pupae

DFly <- cbind(percent.em, DFly)

#look at New percent data
DFly

#fit model: logistic regression model using the binary function

#percent emergence
fly_glm <- glm(percent.em ~ treatt*trial, 
data = DFly, family = "binomial", weights = pupae)
Anova(fly_glm)

fly_glm_1 <- glm(percent.em ~ treatt+trial, 
data = DFly, family = "binomial", weights = pupae)
Anova(fly_glm_1)

AIC(fly_glm, fly_glm) # fly_glm is better, output agrees with visual observation of plots

#development time
fly_glm1 <- glm(time ~ treatt+trial, family = Gamma (link ="log"), data = DFly)
fly_glm2 <- glm(time ~ treatt*trial, family = Gamma (link ="log"), data = DFly)

AIC(fly_glm1, fly_glm2) # fly_glm1 better

Anova(fly_glm1) # no effect of treatt, but able to show sig. diff. of trial in posthoc
Anova(fly_glm2) # no effect of treatt, not able to show sig. diff. of trial in post hoc

#Post hoc test
library (emmeans)

# percent emergence
HSD_test <- emmeans(ref_grid(fly_glm, regrid = "response"), pairwise ~ treatt * trial)    
emmeans(ref_grid(fly_glm, regrid="response"), pairwise~treatt|trial)
emmeans(ref_grid(fly_glm, regrid="response"), pairwise~trial|treatt)

#development time
HSD_test <- emmeans(ref_grid(fly_glm1, regrid = "response"), pairwise ~ treatt *trial) 
emmeans(ref_grid(fly_glm1, regrid="response"), pairwise~trial|treatt) 

# plot for Percent Emergence of Delia adults- DFly

em <- ggplot(DFly, aes(treatt, percent))+ 
  geom_boxplot(aes(fill = trial), notch = FALSE, position = position_dodge(0.9), 
               show.legend = TRUE) +
  theme_minimal() + labs(x="Soil amendment", y = "Fly emergence (%)")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())+ 
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))
em + ylim(0,80)

#Adjust text and font size of legend, percent adult emergence

em + theme(legend.text=element_text(size = 10),
              legend.title=element_text(size=12, hjust=0.5),
              legend.key.height=unit(0.5,"cm"),
              legend.key.width=unit(0.5,"cm"),
              legend.position = "bottom") + ylim(0,80)
ggsave("boxplot_em.jpeg", width = 4, height = 3, dpi = 1000)

# plot for Time taken to Emerge (development time) of Delia adults- DFly

dev<- ggplot(DFly, aes(treatt, time))+ 
  geom_boxplot(aes(fill = trial), notch = FALSE, position = position_dodge(0.9), 
               show.legend = TRUE) +
  theme_minimal() + labs(x="Soil amendment", y = "Fly development time (days)")+
  theme_bw() + theme(panel.grid.major = element_blank(),# removes the grid background 
                     panel.grid.minor = element_blank())+ 
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))

dev + ylim(4,16)

#Adjust text and font size of legend, percent adult emergence

dev + theme(legend.text=element_text(size = 10),
           legend.title=element_text(size=12, hjust=0.5),
           legend.key.height=unit(0.5,"cm"),
           legend.key.width=unit(0.5,"cm"),
           legend.position = "bottom") + ylim(4,16)
ggsave("boxplot_dev.jpeg", width = 4, height = 3, dpi = 1000)

