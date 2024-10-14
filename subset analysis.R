library(knitr)
library(reporttools)
library(xtable)
library(ICC)
library(readxl)
library(forestplot)
library(dplyr)
library(glmtoolbox)
library(ggplot2)
library(nlme) 
library(lme4)
library(MASS)
library(splines)
library(multcomp)
library(car)
library(ggpubr)


##############################
##subset analysis
##############################

#read data
data.bt008.used<-read.csv("/Users/chixiang.chen/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/UMSOM_collaboration/woodworth/BT008tech/manuscript/revision/data.bt008.used.csv")
data.bt008.used$HD3l<-ifelse(data.bt008.used$HD<=0.5,"1",ifelse(
  data.bt008.used$HD<=1.5,"2","3"
))

#subset the data
table(data.bt008.used$controller)
data.bt008.used<-data.bt008.used[data.bt008.used$version==2,]
dim(data.bt008.used)

#do univariate analysis
y<-data.bt008.used$perTarEnhanced
id<-data.bt008.used$PatientIDNumber
Site<-data.bt008.used$subject.id
names(data.bt008.used)
data.univariate<-data.bt008.used[,c(11:13,18:21,25:26,
                                    28,32,37:41,59
)]
head(data.univariate)

p_all<-rep()
library(lme4)
library(nlme)
for(i in 1:ncol(data.univariate))
{
  data.univariate.i<-data.frame(y=y,x=as.numeric(data.univariate[,i]),id=id,Site=Site)
  data.univariate.i<-data.univariate.i[complete.cases(data.univariate.i),]
  fit<-lme(fixed=y ~ x,
           random=~1|Site/id,
           data=data.univariate.i)
  result<-summary(fit) 
  p_i<-coef(result)[2,c(1,5)]
  ci<-intervals(fit,which = "fixed")
  p_all<-cbind(p_all,c(p_i,ci$fixed[2,c(1,3)]))
}
colnames(p_all)<-names(data.univariate)
Table_uni<-t(round(p_all,3)) #no include NumberOfPulses and cycle duty
colnames(Table_uni)<-c("Effect size","P-value","95%LL","95%UL")
Table_uni
#generate results for supplementary table 3 
Table_uni[c("GridSize","AvgSDR",
            "ActElWithAngleBelow15",
            "NoElOn","SubspotPowerMedian",
            "HD","sum.sonit","Tumor_size"),] 


#do multivariable analysis
data.subset<-data.bt008.used[,c("perTarEnhanced","Site","subject.id","PatientIDNumber",
                                "HD","HD3l","SubspotPowerMedian","NoElOn","Tumor_size")]

Subs2<-data.subset[complete.cases(data.subset),]
length(unique(Subs2$subject.id))  #6 subjects in analysis


#with subspotmedpower
fit<-lme(perTarEnhanced ~ HD + SubspotPowerMedian + NoElOn 
         ,
         random=~1|subject.id/PatientIDNumber,
         data=Subs2)
result<-summary(fit)
ci<-intervals(fit,which = "fixed")
multi_reduced<-cbind(coef(result),ci$fixed)
multi_reduced<-multi_reduced[,c(1,2,4,5,6,8)]
multi_reduced<-data.frame(multi_reduced)
colnames(multi_reduced)<-c("Effect size",	"Std.Error",	"t-statistics",
                           "P-value",	"95%LL",	"95%UL")

round(multi_reduced,3)

multi_reduced<-multi_reduced[-1,c(1,4,5,6)]
#generat results in Table 4
round(multi_reduced,3) 


#stratefied analysis based on 3-level HD
Subs2$HD2l<-ifelse(Subs2$HD<=0.5,0,1)
Subs2$HD3l<-as.factor(Subs2$HD3l)
fit<-lme(perTarEnhanced ~ relevel(HD3l, ref = "1"), 
         # (1 | id:Site) +
         #  (1 | Site),
         # random = list(
         #   Site = ~1,
         #   id = ~HD),
         random=~1|subject.id/PatientIDNumber,
         # random=~1|Site/id,
         
         # random=~HD|Site,
         #correlation = corAR1(),
         data=Subs2)
dim(Subs2)
result<-summary(fit)
ci<-intervals(fit,which = "fixed")
multi_reduced<-cbind(coef(result),ci$fixed)
multi_reduced<-multi_reduced[,c(1,2,4,5,6,8)]
multi_reduced<-data.frame(multi_reduced)
colnames(multi_reduced)<-c("Effect size",	"Std.Error",	"t-statistics",
                           "P-value",	"95%LL",	"95%UL")
#get significant level for supplementary figure 2
round(multi_reduced,3)


#generate supplementary figure 2
stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
  #"1",     "2", "ns",
  #"1",     "3", "ns",
  "1",     "2", "**",
  "1",     "3", "**",
  "2",     "3", "ns"
)

library(ggpubr)
ggboxplot(Subs2, x = "HD3l", y = "perTarEnhanced",fill="HD3l",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.4,lwd=0.7) +
  stat_pvalue_manual(
    stat.test,
    size=9,
    y.position = 105, step.increase = 0.1,
    label = "p.adj"
  ) +
  ylim(0,130) +
  xlab("SD (A.U.)") + ylab("T1c (% enhance)") +
  theme(  axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 18),
          axis.title=element_text(size=20,face="bold"))+
  scale_fill_brewer(type = "seq", palette = "Blues")+
  scale_x_discrete(labels=c('(0,0.5]', '(0.5,1.5]','(1.5,2]')) +
  geom_smooth(method = "lm", se=FALSE, color="purple",aes(group="HD3l")
  )

















