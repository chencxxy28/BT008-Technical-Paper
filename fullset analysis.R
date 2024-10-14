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
##full analysis
##############################

#read data file
data.bt008.used<-read.csv("/Users/chixiang.chen/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/UMSOM_collaboration/woodworth/BT008tech/manuscript/revision/data.bt008.used.csv")

#generate extra variables for exploration
table(data.bt008.used$WSVer)
data.bt008.used$controller<-ifelse(data.bt008.used$WSVer<7.4,"0",
                                   ifelse(data.bt008.used$WSVer<7.41,"1","1"))
data.bt008.used$controller<-ifelse(data.bt008.used$WSVer<7.4,"0",
                                   ifelse(data.bt008.used$WSVer<7.41,"1","1"))
table(data.bt008.used$controller)
tapply(data.bt008.used$perTarEnhanced,data.bt008.used$controller,mean)
tapply(data.bt008.used$perTarEnhanced,data.bt008.used$version,mean)
summary(data.bt008.used$version)

data.bt008.used$HD3l<-ifelse(data.bt008.used$HD<=0.5,"1",ifelse(
  data.bt008.used$HD<=1.5,"2","3"
))
data.bt008.used$HD3l <- factor(data.bt008.used$HD3l,levels = c("1", "2", "3"))

data.bt008.used$HD4l<-ifelse(data.bt008.used$HD<=0.5,"1",ifelse(
  data.bt008.used$HD<=1,"2", ifelse(data.bt008.used$HD<=1.5,"3","4")
))
data.bt008.used$HD4l <- factor(data.bt008.used$HD4l,levels = c("1", "2", "3","4"))

summary(data.bt008.used$sum.sonit)
summary(data.bt008.used$T.sum.sonit)

tapply(data.bt008.used$perTarEnhanced,data.bt008.used$T.sum.sonit,mean)
tapply(data.bt008.used$perTarEnhanced,data.bt008.used$version,mean)


#conduct univariate analysis for the full set analysis
names(data.bt008.used)
y<-data.bt008.used$perTarEnhanced
id<-data.bt008.used$PatientIDNumber
Site<-data.bt008.used$subject.id

unique(data.bt008.used$subject.id)
unique(data.bt008.used$PatientIDNumber)
names(data.bt008.used)
table(data.bt008.used$version)
data.univariate<-data.bt008.used[,c(11:13,18:21,25:26,
                                    28,32,37:43,48,58,3,
                                    59
                                    )]
head(data.univariate)
p_all<-rep()
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

#generate results for supplementary table 2 
Table_uni[c("GridSize","AvgSDR",
            "ActElWithAngleBelow15",
            "NoElOn","MedPower",
            "HD","version","sum.sonit","infusion","Tumor_size"),] 




#conduct 3-level multivariable analysis
names(data.bt008.used)
data.subset<-data.bt008.used[,c("perTarEnhanced","Site","id", "MedPower","subject.id","PatientIDNumber",
                                "HD3l","HD4l", "sum.sonit", "GridSize", "AvgSDR","ActElWithAngleBelow15",
                                "NoElOn","sum.sonit",
                                "controller","version",
                                "HD","infusion","Tumor_size")]
Subs1<-data.subset[complete.cases(data.subset),]
length(unique(Subs1$subject.id))  #28 subjects in analysis

##final analysis
fit<-lme(perTarEnhanced ~  HD + MedPower + GridSize + ActElWithAngleBelow15 +
           sum.sonit + AvgSDR +
           NoElOn + 
            infusion +
        Tumor_size,
         random=~1|subject.id/PatientIDNumber,
         method="ML",
         data=Subs1)
result<-summary(fit)
ci<-intervals(fit,which = "fixed")
multi_full<-cbind(coef(result),ci$fixed)
multi_full<-multi_full[,c(1,2,4,5,6,8)]
multi_full<-data.frame(multi_full)
colnames(multi_full)<-c("Effect size",	"Std.Error",	"t-statistics",
                     "P-value",	"95%LL",	"95%UL")
round(multi_full,3)

multi_full<-multi_full[-1,c(1,4,5,6)]
used.names<-rownames(multi_full)
multi_full<-multi_full[c(1,5,2,7,4,6,3,8,9),]
#generate the results for Table 2
round(multi_full,3) 

#exploratory analyses 
##variable selection
bestfit <- stepAIC(fit)
fit<-lme(perTarEnhanced ~  HD  +
           sum.sonit +
           NoElOn + 
           Tumor_size,
         random=~1|subject.id/PatientIDNumber,
         method="ML",
         data=Subs1)
result<-summary(fit)
ci<-intervals(fit,which = "fixed")
multi_full<-cbind(coef(result),ci$fixed)
multi_full<-multi_full[,c(1,2,4,5,6,8)]
multi_full<-data.frame(multi_full)
colnames(multi_full)<-c("Effect size",	"Std.Error",	"t-statistics",
                        "P-value",	"95%LL",	"95%UL")
round(multi_full,3)

#analyses with B-spline expansion 
Bbase<-bs(Subs1$HD,degree=2)
Subs1$base1<-Bbase[,1]
Subs1$base2<-Bbase[,2]
fit<-lme(perTarEnhanced ~  base1 + base2 +
           MedPower + GridSize + ActElWithAngleBelow15 +
           sum.sonit + AvgSDR +
           NoElOn + 
           infusion +
           Tumor_size,
         random=~1|subject.id/PatientIDNumber,
         method="ML",
         data=Subs1)
result<-summary(fit)
ci<-intervals(fit,which = "fixed")
multi_full<-cbind(coef(result),ci$fixed)
multi_full<-multi_full[,c(1,2,4,5,6,8)]
multi_full<-data.frame(multi_full)
colnames(multi_full)<-c("Effect size",	"Std.Error",	"t-statistics",
                        "P-value",	"95%LL",	"95%UL")
round(multi_full,3)

hypothesis_matrix<-rbind(
  c(0,1,0,rep(0,8)),
  c(0,0,1,rep(0,8)))
wald_test <- linearHypothesis(fit, hypothesis_matrix, contrast = contrast_vector)
#test results for figure 5(A)
wald_test 

##stratefied analysis based on 2-level HD
Subs1$HD2l<-as.factor(ifelse(Subs1$HD<=1, 0,1))
fit<-lme(perTarEnhanced ~ relevel(HD2l, ref = "0") +
         MedPower + AvgSDR  + GridSize + ActElWithAngleBelow15 +
           NoElOn 
         + sum.sonit + 
           + infusion +
           Tumor_size,
         random=~1|subject.id/PatientIDNumber,
         data=Subs1)
result<-summary(fit)
ci<-intervals(fit,which = "fixed")
multi_full_2l<-cbind(coef(result),ci$fixed)
multi_full_2l<-multi_full_2l[,c(1,2,4,5,6,8)]
multi_full_2l<-data.frame(multi_full_2l)
colnames(multi_full_2l)<-c("Effect size",	"Std.Error",	"t-statistics",
                        "P-value",	"95%LL",	"95%UL")
# get p-value for figure 5 (B)
round(multi_full_2l,3) 



#generate figure 5(B)
stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
  #"1",     "2", "ns",
  #"1",     "3", "ns",
  "0",     "1", "**",
)

library(ggpubr)
p1<-ggboxplot(Subs1, x = "HD2l", y = "perTarEnhanced",fill="HD2l",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.4,lwd=0.7) +
  stat_pvalue_manual(
    stat.test,
    size=8,
    y.position = 102, step.increase = 0.05,
    label = "p.adj"
  ) +
  ylim(0,110) +
  xlab("SD (A.U.)") + ylab("BBB Opening (% New T1c)") +
  theme(  axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=12,face="bold")
          # axis.text = element_text(size = 18),
          # axis.title=element_text(size=20,face="bold")
          )+
  scale_fill_brewer(type = "seq", palette = "Blues")+
  scale_x_discrete(labels=c('(0,1]', '(1,2]')) +
  scale_y_continuous(
    breaks = c(0, 50, 100),
  ) 
p1


# #for scatter plot combining both data
# names(Subs1)
# Subs1$Group<-ifelse(Subs1$version==2,"Subsonication Level","Target Level")
# 
# ggplot(Subs1, aes(x=HD, y=perTarEnhanced, group=Group,col=Group)) +
#   geom_point(alpha=.3,size=1,aes(group = Group)) +
#   geom_smooth(size=1.5,method = "lm",fill=NA) +
#   scale_color_manual(values=c( "#E69F00", "skyblue")) +
#   xlab("SD (A.U.)")+
#   ylab("BBB Opening (% New T1c)") +
#   #theme(legend.position="")
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title = element_text(size=14))
# 
# ggplot(Subs1, aes(x=HD, y=perTarEnhanced)) +
#   # geom_point(data = Subs1[Subs1$perTarEnhanced>=50 |Subs1$perTarEnhanced<=10,],alpha=.3,size=1, aes(color = (subject.id))) +
#   geom_point(alpha=.3,size=1, aes(color = (subject.id))) +
#   # geom_smooth(size=1.5,method = "lm",fill=NA, color="skyblue",fullrange = TRUE) +
#   geom_smooth(size=1.5,method = lm,formula = y ~ splines::bs(x, degree=3),fill=NA, color="skyblue",fullrange = TRUE) +
#   #scale_color_manual(values=c( "skyblue")) +
#   xlab("SD (A.U.)")+
#   ylab("BBB Opening (% New T1c)") +
#   #ylim(c(50,100)) +
#   # scale_y_continuous(
#   #   breaks = c(0, 80, 100),
#   # ) +
#   #theme(legend.position="")
#   theme_classic() +
#   theme(axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=12,face="bold"),
#     legend.position = "none")



# generate figure 5(A)
library(mgcv)
library(broom)
model1 <- gam(perTarEnhanced ~ s(HD,k=4), data=Subs1)
broom::tidy(model1)
plot(model1)
p2<-ggplot(Subs1, aes(x=HD, y=perTarEnhanced)) +
  # geom_point(data = Subs1[Subs1$perTarEnhanced>=50 |Subs1$perTarEnhanced<=10,],alpha=.3,size=1, aes(color = (subject.id))) +
  # geom_point(alpha=.3,size=1, aes(color = (PatientIDNumber))) +
  geom_point(alpha=.4,size=1,color="lightblue") +
  # geom_smooth(size=1.5,method = "lm",fill=NA, color="skyblue",fullrange = TRUE) +
  geom_smooth(size=1.5,method = lm,formula = y ~ splines::bs(x, degree=3),se=T,color="grey46",fullrange = TRUE,alpha = 0.2) +
  #scale_color_manual(values=c( "skyblue")) +
  ylim(0,100) +
  xlab("SD (A.U.)")+
  ylab("BBB Opening (% New T1c)") +
  #ylim(c(50,100)) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = c(0, 50, 100),
  ) +
  #theme(legend.position="")
  theme_classic() +
  theme(axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "none")
p2


#combine Figure 5(A) and Figure 5(B)
ggarrange(p2, p1, ncol = 2, nrow = 1, widths=c(1.4,1))










