
# packages
library(AICcmodavg)
library(corrplot)
library(ggplot2)
library(interactions)
library(lme4)
library(MASS)
library(performance)
library(tidyverse)
library(visreg) 

# data
NOCAspatial <- read.csv("~/NOCAspatial_short.csv",header=TRUE)

# independent variable selection
colnames(NOCAspatial)
corNOCAspatial <- NOCAspatial[, c(24:36)] # select independent variables

# change column names to match manuscript text
colnames(corNOCAspatial) <- c("ShrubDen", "MSden","OSden","GapWeight", "NumGaps", 
                              "DisEdge", "UrbanOpen", "UrbanDev", "UrbanMed", 
                              "UrbanHi", "DecForest","PatchArea", "NestHeight")

# correlation plot
cor <- cor(corNOCAspatial)
corrplot(cor, type="upper", order="hclust", tl.cex = 0.9, cl.cex = 1, mar=c(0,0,0,0), tl.col='black') 

# scale 10 retained indepedent variables
NOCAspatial$scShrubDtc <- c(scale(NOCAspatial$ShrubDtc, scale=TRUE))
NOCAspatial$scMidstoryDtc <- c(scale(NOCAspatial$MidstoryDtc, scale=TRUE))
NOCAspatial$scOverstoryDtc <- c(scale(NOCAspatial$OverstoryDtc, scale=TRUE)); 
NOCAspatial$scTotalGapDomin <- c(scale(NOCAspatial$totalGapDomin, scale=TRUE))
NOCAspatial$scNumGaps <- c(scale(NOCAspatial$numGaps, scale=TRUE))
NOCAspatial$scDisToEdge <- c(scale((NOCAspatial$DisToEdge), scale=TRUE))
NOCAspatial$scDeveloped_Low <- c(scale((NOCAspatial$Developed_Low), scale=TRUE))
NOCAspatial$scDeciduous_Forest <- c(scale(NOCAspatial$Deciduous_Forest, scale=TRUE))
NOCAspatial$scAreaStudySite <- c(scale(NOCAspatial$AreaStudySite, scale=TRUE))
NOCAspatial$NestHeight[is.na(NOCAspatial$NestHeight)]<-mean(NOCAspatial$NestHeight,na.rm=TRUE) #interpolate nest height NA's
NOCAspatial$scNestHeight <- c(scale(as.numeric(NOCAspatial$NestHeight), scale=TRUE))
colnames(NOCAspatial)

# exposure link function
# library(MASS)
remove(..exposure) # remove mean exposure after plotting interaction if code is run again

logexp <- function(exposure = 1) {
  get_exposure <- function() {
    if (exists("..exposure", env=.GlobalEnv))
      return(get("..exposure", envir=.GlobalEnv))
    exposure
  }
  linkfun <- function(mu) qlogis(mu^(1/get_exposure()))
  linkinv <- function(eta) plogis(eta)^get_exposure()
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta)>30,.Machine$double.eps,
           exp(eta)/(1+exp(eta))^2)
  }
  mu.eta <- function(eta) {       
    get_exposure() * plogis(eta)^(get_exposure()-1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

# candidate models
Cand <- list()

# null model
Cand[[1]] <- lme4::glmer(Success ~ + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[1]])

## local-scale models

# overstory density (Dtc = detections)
Cand[[2]] <- lme4::glmer(Success ~ scOverstoryDtc + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"), nAGQ = 4); summary(Cand[[2]]) 

# shrub density 
Cand[[3]] <- lme4::glmer(Success ~ scShrubDtc + (1|Year:Site),data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[3]])

# midstory density
Cand[[4]] <- lme4::glmer(Success ~ scMidstoryDtc + (1|Year:Site),data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[4]])

# shrub and oveerstory density
Cand[[5]] <- lme4::glmer(Success ~ scOverstoryDtc + scShrubDtc + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"), nAGQ = 4); summary(Cand[[5]]) 

# overstory and midstory density
Cand[[6]] <- lme4::glmer(Success ~ scOverstoryDtc + scMidstoryDtc + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"), nAGQ = 4); summary(Cand[[6]]) 

# gap number
Cand[[7]] <- lme4::glmer(Success ~ scNumGaps + (1|Year:Site),data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[7]])

# gap weight
Cand[[8]] <- lme4::glmer(Success ~ scTotalGapDomin + (1|Year:Site),data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[8]])

# nest height
Cand[[9]] <- lme4::glmer(Success ~ scNestHeight + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                         control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[9]])

# overstory density and gap number
Cand[[10]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[10]]) 
# top AIC local-model passed to patch-and landscape-scale models 
performance::check_collinearity(Cand[[10]]) # low VIF 

# overstory density and gap weight
Cand[[11]] <- lme4::glmer(Success ~ scOverstoryDtc + scTotalGapDomin + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[11]]) 

# overstory density and nest height
Cand[[12]] <- lme4::glmer(Success ~ scOverstoryDtc + scNestHeight + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[12]])

# gap number and nest height
Cand[[13]] <- lme4::glmer(Success ~ scNumGaps + scNestHeight + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)),
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[13]])

## patch- and landscape-scale models added to best best local-scale model (cand10)

# overstory density, gap number and patch size
Cand[[14]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + scAreaStudySite + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[14]]) 

# overstory density, gap number and edge distance
Cand[[15]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + scDisToEdge + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[15]]) 

# overstory density, gap number and urban development
Cand[[16]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps +  scDeveloped_Low + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[16]]) 
# overstory density, gap number and forest cover
Cand[[17]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps +  scDeciduous_Forest + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[17]]) 

# overstory density, gap number,urban development and edge distance
Cand[[18]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + scDeveloped_Low + scDisToEdge + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)),
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[18]]) 
# overstory density, gap number,forest cover and edge distance
Cand[[19]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + scDeciduous_Forest + scDisToEdge + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[19]]) 

## ecologically relevant interactions

# overstory density, gap number,edge distance, urban development and edge distance * urban development (tests if edge effects on nest success vary with urban development level)
Cand[[20]] <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + scDisToEdge * scDeveloped_Low + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)),
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary(Cand[[20]]) 

# overstory density, edge distance, gap number, urban develoement and  gap number * urban development (tests if number of gaps effect on nest success varies with urban development level)
Cand[[21]] <- lme4::glmer(Success ~ scOverstoryDtc + scDisToEdge + scNumGaps * scDeveloped_Low + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary( Cand[[21]]) 

# overstory density, urban development, gap number, edge distance and gap number * edge distance (tests if number of gaps effect on nest success varies with edge distance)
Cand[[22]] <- lme4::glmer(Success ~ scOverstoryDtc + scDeveloped_Low + scNumGaps * scDisToEdge + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary( Cand[[22]]) 

# gap number, edge distance, overstory density, urban development and  overstory density * urban development (tests if overstory density effect on nest success varies with urban development level)
Cand[[23]] <- lme4::glmer(Success ~ scNumGaps + scDisToEdge + scOverstoryDtc * scDeveloped_Low + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary( Cand[[23]]) 

# gap number, urban development, overstory density, edge distance and overstory density * edge distance (test if overstory density effect effect on nest success varies with edge distance)  
Cand[[24]] <- lme4::glmer(Success ~  scNumGaps + scDeveloped_Low + scOverstoryDtc*scDisToEdge + (1|Year:Site), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                          control = glmerControl(optimizer = "bobyqa"),nAGQ = 4); summary( Cand[[23]]) 

# AIC table
AICtable <-   AICcmodavg::aictab(Cand, sort = TRUE); AICtable
AICt<-as.data.frame(AICtable)
write_csv(AICt, "~/Downloads/AICtable.csv")

# top AIC model variance inflation factors (VIF)
performance::check_collinearity(Cand[[20]]) # low VIF

# top model AIC convergence 
performance::check_convergence(Cand[[20]]) # true

# set exposure variable to mean exposure for plotting
..exposure <- mean(NOCAspatial$Exposure)

# main effects
a<-visreg(Cand[[20]], "scOverstoryDtc", scale="response" ,rug=2, band=F); #a$fit
b<-visreg(Cand[[20]], "scNumGaps", scale="response" ,rug=2, band=F);# b$fit
c<-visreg(Cand[[20]], "scDeveloped_Low", scale="response" ,rug=2); c$fit
d<-visreg(Cand[[20]], "scDisToEdge", scale="response" ,rug=2); d$fit

# interaction plot
e<-visreg(Cand[[20]], "scDisToEdge", by="scDeveloped_Low", overlay=T, type="contrast", breaks=2); e$fit # doesn't average across sites, not ideal

# interaction plot using interact_plot (visreg convergence issues avoided)

# partion 8 sites into 4 most urban development and 4 least (easier to interpret; used in publication)
NOCAspatial$Urban <- cut(NOCAspatial$scDeveloped_Low,
                         breaks=c(-5, 0.23, 5),
                         labels=c('Less', 'More'))

fit4 <- lme4::glmer(Success ~ scOverstoryDtc + scNumGaps + Urban * scDisToEdge + (1|Year), data = NOCAspatial, family = binomial(link=logexp(NOCAspatial$Exposure)), 
                    control = glmerControl(optimizer = "bobyqa"),nAGQ = 4) # 1|Year random effect used here because 1|Year:Site does not converge with urban development as categories

interactions::interact_plot(fit4, pred=scDisToEdge, modx=Urban, colors=c("cornflowerblue","black"), int.width = .95,
                            lty=1,  plot.points=T, point.shape=T,#modx.values="plus-minus",
                            x.label="\n Distance to Edge", y.label = "Probability of Success \n",
                            legend.label = "Development\n Level", modx.labels = c("Less Developed", "More Developed"),
                            point.size = 3, line.thickness = 1.5, vary.lty = T, jitter=0.02)+
  theme_bw()+
  theme(panel.grid=element_blank()) 


# Based on mean and sd's (regression line intersection more accurately depicted; not in publication)
interactions::interact_plot(Cand[[20]], pred=scDisToEdge, modx=scDeveloped_Low, colors=c("cornflowerblue","black"),#interval = TRUE,int.type = "prediction", int.width = .95,
                            lty=1,  plot.points=T, point.shape=T, modx.values = c(-1, 1), #Set at +-1 sd or modx.values="plus-minus",
                            x.label="\n Distance to Edge", y.label = "Probability of Success \n",
                            legend.label = "Development\n Level", modx.labels = c("Less Developed", "More Developed"),
                            point.size = 4, line.thickness = 1.5, vary.lty = T, jitter=0.01)+ 
  scale_color_gradient(low = "cornflowerblue", high = "black")+
  theme_bw()+
  theme(panel.grid=element_blank()) #+

# distances among sites
AmongSite <- c(9500	,1390	,3970	,6980,	6120	,5210	,5630,
               10050	,13540,	16320,	15610	,10770,	15090,
               4140,	6250	,5640,	6710,	5100,
               4360,	2690	,6140	,1650,
               1800,	10690,	4460,
               8890,	2530,
               6780)


