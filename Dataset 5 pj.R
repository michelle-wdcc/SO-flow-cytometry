my_data<-read.csv("C:/Uni/BioSci220/Sempj/Southern_Ocean_flow_cytometry.csv")
#Clean data
my_data<-my_data[1:263,c(1:12,18)]
my_data$Start.Longitude[37]<-NA
my_data$Depth..m. <- as.numeric(my_data$Depth..m.)
my_data$Total.bacteria.mL<-as.numeric(my_data$Total.bacteria.mL)
my_data$Active.bacteria.mL<-as.numeric(my_data$Active.bacteria.mL)

my_data <- my_data[,-c(13),] #Remove Integrated Phytoplankton as there are very few measurements for this variable
my_data <- na.omit(my_data) #removes remaining samples with missing values 

summary(my_data)
str(my_data)
install.packages("s20x")
install.packages("tidyverse")
library("s20x")
library("tidyverse")
#Investigate the relationship between variables
pairs20x(pj_data)

##Choose variables to reesearch about
pj_data=cbind(my_data[,7:12],Transect=my_data$Transect.number)
pj_data$Transect=as.factor(pj_data$Transect)

##Exploring plots
#Phytoplankton concentration distribution across transects
ggplot(pj_data, aes(x=Phytoplankton.mL)) + geom_histogram()+facet_wrap(~Transect)+
  ggtitle("Phytoplankton distribution")+theme_bw()
#Protozoa concentration distribution across transects
ggplot(pj_data, aes(x=Protozoa.mL)) + geom_histogram()+facet_wrap(~Transect)+
  ggtitle("Protozoa distribution")+theme_bw()
#Phytoplankton concentration and protozoa concentration at different depth
ggplot(pj_data, aes(x=Depth..m., y=log(Phytoplankton.mL), col=Transect)) + geom_point()
ggplot(pj_data, aes(x=Depth..m., y=log(Protozoa.mL),col=Transect)) + geom_point()
ggplot(pj_data, aes(x=log(Protozoa.mL), y=log(Phytoplankton.mL),col=Transect)) + geom_point()

ggplot(pj_data, aes(x=log(Protozoa.mL), y=log(Phytoplankton.mL),col=Depth..m.)) + geom_point()+geom_smooth(method = 'lm')+facet_wrap(~Transect)+
  xlab("Log(Protozoa/mL)")+ylab("Log(Phytoplankton/mL)")+
  labs(title="Relationship between log(Protozoa/mL) and log(Phytoplankton/mL) across 6 transects in Southern Ocean", subtitle="The chosen transects were where they put CTD stations ", col="Depth (m)")+
  theme_bw() + theme(axis.title = element_text(face="bold"))



#Interactive model with all 3 explanatory variable
phyto<-lm(log(Phytoplankton.mL)~Depth..m.*Transect*log(Protozoa.mL),data = pj_data)
anova(phyto)

#Interactive model with Depth and Transect only
#Chosen model
phyto.lm<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)+Depth..m.*Transect,data = pj_data)

#Other model
phyto.lm1<-lm(log(Phytoplankton.mL)~log(Protozoa.mL),data = pj_data) #simple linear regression
phyto.lm2<-lm(log(Phytoplankton.mL)~Transect*Depth..m.+(Depth..m.>70),data = pj_data)
phyto.lm3<-lm(log(Phytoplankton.mL)~Transect*Depth..m.*(Depth..m.>70),data = pj_data)
phyto.lm4<-lm(log(Phytoplankton.mL)~Transect+Depth..m.+(Depth..m.>70),data = pj_data)
phyto.lm5<-lm(log(Phytoplankton.mL)~(Transect+Depth..m.+(Depth..m.>70))^2,data = pj_data)
phyto.lm6<-lm(log(Phytoplankton.mL)~Depth..m.*Transect+(Depth..m.>70)+Depth..m.:(Depth..m.>70),data = pj_data)             
phyto.lm7<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)+(Depth..m.>70), data = pj_data)
phyto.lm8<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)+(Depth..m.>70)+Transect, data = pj_data)
phyto.lm9<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)+Depth..m.*(Depth..m.>70)+Transect+Transect:Depth..m., data = pj_data)
phyto.lm10<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)*Transect+(Depth..m.>70), data=pj_data)
phyto.lm11<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)*Transect*g70*Depth..m., data=pj_data)
g70=pj_data$Depth..m.
phyto.lm12<-lm(log(Protozoa.mL)~log(Phytoplankton.mL)*Transect*g70*Depth..m., data=pj_data)

summary(phyto.lm2)
summary(phyto.lm3)
anova(phyto.lm2,phyto.lm3)
anova(phyto.lm3)
anova(phyto.lm5)
anova(phyto.lm4)
anova(phyto.lm6)
summary(phyto.lm6)
summary(phyto.lm7)
summary(phyto.lm8)
summary(phyto.lm9)
anova(phyto.lm9)
summary(phyto.lm10)

phyto.lm
#anova & summary
library(MuMIn)
options(na.action=na.fail)
all.fit=dredge(phyto.lm11)
head(all.fit)
mod1=get.models(all.fit,1)[[1]]
summary(mod1)

all.fit2=dredge(phyto.lm12)
mod2=get.models(all.fit2,1)[[1]]
summary(mod2)

ggplot(pj_data, aes(x=Depth..m.,y=log(Total.bacteria.mL)))+
  geom_point(color="turquoise")+facet_wrap(~Transect)+theme_bw()

#ANCOVA
phyto.lm3<-lm(log(Phytoplankton.mL)~log(Protozoa.mL)+Depth..m.+Transect,data = pj_data)#Additive model
anova(phyto.lm,phyto.lm3) #Investigate if interaction term is significant

AIC(phyto.lm1,phyto.lm,phyto,phyto.lm2,phyto.lm3)



#Choose the smallest AIC model and look at the summary
summary(phyto.lm)

#Check assumption
plot(phyto.lm,which=1) #relatively constant scatter
normcheck(phyto.lm,shapiro.wilk=T) #normality checked

#Quantify findings
confint(phyto.lm)
100*((2^confint(phyto.lm)[2,])-1)


#Prediction 
transect1<- data.frame(Depth..m. = 100, Transect="1", Protozoa.mL=150) ## create new data frame with data we want to predict to
exp(predict(phyto.lm, newdata = transect1, interval = "prediction"))
transect1_x2<- data.frame(Depth..m. = 100, Transect="1", Protozoa.mL=300) ## create new data frame with double Protozoa concentration

transect3<- data.frame(Depth..m. = 100, Transect="3", Protozoa.mL=150) ## create new data frame with different transect
exp(predict(phyto.lm, newdata = transect3, interval = "prediction"))

transect7<- data.frame(Depth..m. = 100, Transect="7", Protozoa.mL=150) ## create new data frame with different transect
exp(predict(phyto.lm, newdata = transect7, interval = "prediction"))

compare=data.frame(Depth..m.=c(1,10), Transect=c("1","1","3","3"), Protozoa.mL=150)## create new data frame with different depth and transect
a=exp(predict(phyto.lm, newdata = compare, interval = "conf"))
compare
a[c(1,2),] ## Phytoplankton concentration when depth =1,10m with Protozoa.ml=150 in transect 1
a[c(3,4),] ## Phytoplankton concentration when depth =1,10m with Protozoa.ml=150 in transect 1