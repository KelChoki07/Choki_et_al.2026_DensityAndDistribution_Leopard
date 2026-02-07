
# remove.packages("RPresence")
# 
# if ("RPresence" %in% installed.packages()) {
#   remove.packages("RPresence")
# } else {
#   message("Package not installed.")
# }



# install.packages("unmarked")
# install.packages('RPresence',repo='https://www.mbr-pwrc.usgs.gov/mbrCRAN')	
library("RPresence"); library(dplyr); library(ggplot2); 
library(kableExtra); library(tidyverse)


csv1 <- read.csv("overall_combined_DetHist.csv", header=T)

#impute blanks with data- cannot run AIC when number of datapoints differs between models
csv1$elevation<-ifelse(is.na(csv1$elevation)==T,mean(csv1$elevation,na.rm=T),csv1$elevation)
csv1$treecover<-ifelse(is.na(csv1$treecover)==T,mean(csv1$treecover,na.rm=T),csv1$treecover)
csv1$strmden<-ifelse(is.na(csv1$strmden)==T,mean(csv1$strmden,na.rm=T),csv1$strmden)
csv1$housden<-ifelse(is.na(csv1$housden)==T,mean(csv1$housden,na.rm=T),csv1$housden)
csv1$prey_count<-ifelse(is.na(csv1$prey_count)==T,mean(csv1$prey_count,na.rm=T),csv1$prey_count)

sum(is.na(csv1$elevation))
summary(csv1)


nsites=nrow(csv1); nsrvys=ncol(csv1[, c(2:17)])  #  set number of sites,surveys from det. data
dethist=matrix(as.integer(unlist(csv1[,c(2:17)])),nrow=nsites) # replace missing values(-) with NA


# create site covariate data frame, scaling the covariates to help optimization...
cov1= data.frame(apply(csv1[,18:22], 2, scale))

##create input "pao" object, for use with occMod function
data=createPao(dethist,
               unitcov=cov1,
               nsurveyseason = c(8,8),
               title="2-species multi-season occupancy data")


str(data)

#  _____            _   _ _  _______ _   _  _____ 
# |  __ \     /\   | \ | | |/ /_   _| \ | |/ ____|
# | |__) |   /  \  |  \| | ' /  | | |  \| | |  __ 
# |  _  /   / /\ \ | . ` |  <   | | | . ` | | |_ |
# | | \ \  / ____ \| |\  | . \ _| |_| |\  | |__| |
# |_|  \_\/_/    \_\_| \_|_|\_\_____|_| \_|\_____|
                                                 
# quick tips:
#look at the unique probabilities for psi estimates
# lapply(dot$real, FUN = unique)


#  _____       _            _   _             
# |  __ \     | |          | | (_)            
# | |  | | ___| |_ ___  ___| |_ _  ___  _ __  
# | |  | |/ _ \ __/ _ \/ __| __| |/ _ \| '_ \ 
# | |__| |  __/ ||  __/ (__| |_| | (_) | | | |
# |_____/ \___|\__\___|\___|\__|_|\___/|_| |_|
                                             
#Model 1
NULL.mod.p <- occMod(model=list( 
  psi~SP,                
  gamma~SP,    
  epsilon~SP,            
  p~SP),   
  data=data,
  type="do.2sp.1")                                             

INT.mod.p1 <- occMod(model=list( 
  psi~SP,                
  gamma~SP,    
  epsilon~SP,            
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

aic.p <- createAicTable(
  list(NULL.mod.p, #1,
       INT.mod.p1),
  use.aicc = F
)

aic_table.p <- aic.p$table
aic_table.p #INT.mod.p is the only top model!!

coef(INT.mod.p1, param="p")
lapply(INT.mod.p1$real, FUN = unique)


#  _____   _____ _____ 
# |  __ \ / ____|_   _|
# | |__) | (___   | |  
# |  ___/ \___ \  | |  
# | |     ____) |_| |_ 
# |_|    |_____/|_____|

#  defining a function to print the real estimates of a model in
#   a table.  Only the 1st site of each parameter is printed.
prnt <- function(mod) {
  x=matrix(unlist(lapply(mod$real,function(a) a[1,])),ncol=4,byrow=F)
  rownames(x)=rep(names(mod$real))
  colnames(x)=colnames(mod$real[[1]])
  cat('model:',mod$modname,'\n'); print(unique(x))
}

#Model 1
NULL.mod.psi <- occMod(model=list(  #with no covariates used
  psi~SP,                 #  PsiA,psiB,psiBa
  gamma~SP,    #  species effect on colonization
  epsilon~SP,             #  interaction among all effects on extinction
  p~SP+INT_o),   #  pA,pB,rA,rBA,rBa
  data=data,
  type="do.2sp.1")

prnt(NULL.mod.psi) 
NULL.mod.psi$beta$psi #only parameter psiA and psiBA (bc psiBA=psiBa)

INT.mod.psi <- occMod(model=list(  #with no covariates used
  psi~SP+INT,                 #  PsiA,psiB,psiBa
  gamma~SP,    #  species effect on colonization
  epsilon~SP,             #  interaction among all effects on extinction
  p~SP+INT_o),   #  pA,pB,rA,rBA,rBa
  data=data,
  type="do.2sp.1")

INT.mod.psi$beta$psi  #Three parameters (psiA and psiBA, psiBa)

aic.psi <- createAicTable(
  list(NULL.mod.psi, #1
       INT.mod.psi),
  use.aicc = T
)

aic_table.psi <- aic.psi$table
aic_table.psi  #NULL.mod.psi is top but followed up by INT.mod.psi 
#but we used both these models as per our hypotheses testing,

coef(NULL.mod.psi, param ="psi", prob = 0.95)
coef(INT.mod.psi, param ="psi", prob = 0.95)

#Actual model based on the hypothesis 
psi~SP+INT+elevation #Tiger model
psi~SP+INT+housden+elevation #topdown model

psi~SP+INT*housden+elevation #Human Shield
psi~SP+INT*prey_count+elevation #Human Shield

psi~SP+prey+elevation #prey model
psi~SP+houden+elevation #human model
psi~SP+treecover+stremden+elevation #habitat model
psi~SP+prey_count+treecover+stremden+elevation #bottomup model

Psi~SP+INT*housden+INT*prey_count+prey_count+treecover+stremden+elevation #global model

  #   _____              __  __ 
  # / ____|      /\     |  \/  |
  # | |  __     /  \    | \  / |
  # | | |_ |   / /\ \   | |\/| |
  # | |__| |  / ____ \  | |  | |
  # \_____|  /_/    \_\ |_|  |_|
  

gam.null <- occMod(model=list( 
  psi~SP,            
  gamma~SP,    
  epsilon~SP,             
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

coef(gam.null, param="gamma")


gam.INT_A_INT_B<- occMod(model=list( 
  psi~SP,                 
  gamma~SP+INT_A+INT_B,    
  epsilon~SP,             
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

gam.INT_A_INT_B$beta$gamma
coef(gam.INT_A_INT_B, param="gamma")

coef(gam.null, param ="gamma", prob = 0.95)
coef(gam.INT_A_INT_B, param ="gamma", prob = 0.95)
coef(gam.preyNtree, param ="gamma", prob = 0.95)

#gam with cov 
#You can just use only INT_A but not just INT_B bc it doesn't really make much sense to use it alone


gam.preyNtree <- occMod(model=list( 
  psi~SP,                 
  gamma~SP+prey_count+treecover,    
  epsilon~SP,             
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

coef(gam.preyNtree, param="gamma")

aic.gam.cov <- createAicTable(
  list(gam.null,
       gam.INT_A_INT_B, #tiger (topdown)
       gam.preyNtree    #prey+tree(bottomup)
      ),
  use.aicc = T
)

aic.gam.cov.ALL <- aic.gam.cov$table
aic.gam.cov.ALL #all 3 models are top models


#build 3*each parameter in top model for epsilon
#add gamma on psi (botttom and top down) model

#Its EFFECT PARAMETERIZATION
#so you just add up the parameter to get the est for second para bc first is intercept


#gam has six parameters
lapply(gam.null$real, FUN = unique)
lapply(NULL.mod.psi$real, FUN = unique)


gam.null$dmat$gamma.dm
gam.null$beta$gamma #only two para gamAB and gamBAA bc rest 4 para are same to it
gam.preyNtree$beta$gamma

gam.INT_A_INT_B$beta$gamma
lapply(gam.null$real, FUN = unique)

plogis(gam.null$beta$gamma$est)[1] #gamAB #tiger 
plogis(gam.null$beta$gamma$est)[2] #
plogis(gam.null$beta$gamma$est[1]+gam.null$beta$gamma$est[2]) # gamBAA=gamBAa=gamBaA=gamBaa #leopard

#4 unique para= gamAB, gamAb,gamBAA,gamBaA

gam.INT_A_INT_B$dmat$gamma.dm
gam.INT_A_INT_B$beta$gamma
unique(gam.INT_A_INT_B$real$gamAb)
unique(gam.INT_A_INT_B$real$gamAB)


lapply(gam.INT_A_INT_B$real, FUN = unique)

plogis(gam.INT_A_INT_B$beta$gamma$est)[1] #0.1288 #GamAb
plogis(gam.INT_A_INT_B$beta$gamma$est)[2] #0.741632 # this is wrong way of doing since it not mean parameterization
plogis(gam.INT_A_INT_B$beta$gamma$est[1]+gam.INT_A_INT_B$beta$gamma$est[2]) #0.29 gamBAA
plogis(gam.INT_A_INT_B$beta$gamma$est[1]+gam.INT_A_INT_B$beta$gamma$est[3]) #GamAB
plogis(gam.INT_A_INT_B$beta$gamma$est[1]+gam.INT_A_INT_B$beta$gamma$est[2]+gam.INT_A_INT_B$beta$gamma$est[4]) # 0.235 GamBaA
# you are adding [2] here maybe bc of relatedness of gamBAA to gamBaA



gam.preyNtree$dmat$gamma.dm
gam.preyNtree$beta$gamma #2 for para # 2 for covs
View(lapply(gam.preyNtree$real, FUN = unique))
lapply(gam.preyNtree$real, FUN = unique)[14] #gamAB=gamAb, gamBAA=gamBAa=gamBaA=gamBaa
#its going be unique for each site but all para is going be same

plogis(gam.preyNtree$beta$gamma$est)[1]

#95% CI
0.09*1.96
0.45*1.96
0.45*1.44 #85% CI #1.49


  #  ______   _____     _____ 
  # |  ____| |  __ \   / ____|
  # | |__    | |__) | | (___  
  # |  __|   |  ___/   \___ \ 
  # | |____  | |       ____) |
  # |______| |_|      |_____/ 
  # 
  
#gam~SP
esp.NULL <- occMod(model=list( 
  psi~SP,                 
  gamma~SP,    
  epsilon~SP,             
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

esp.INT_A_INT_B <- occMod(model=list( 
  psi~SP,                 
  gamma~SP,   
  epsilon~SP+INT_A+INT_B,             
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")


esp.preyNtree <- occMod(model=list( #doesnt work as same cov has to be used in gamma
  psi~SP,                 
  gamma~SP,   
  epsilon~SP+prey_count+treecover,             
  p~SP+INT_o),   # 
  data=data,
  type="do.2sp.1")


#gam~SP+INT_A+INT_B

esp.INT_A_INT_B2 <- occMod(model=list( 
  psi~SP,                 
  gamma~SP+INT_A+INT_B,   
  epsilon~SP+INT_A+INT_B,             
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

esp.preyNtree2 <- occMod(model=list( #doesnt work
  psi~SP,                 
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,             
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

#gam~SP+prey_count+treecover
esp.preyNtree3 <- occMod(model=list( 
  psi~SP,                 
  gamma~SP+prey_count+treecover,    
  epsilon~SP+prey_count+treecover,             
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

aic.esp <- createAicTable(
  list(esp.NULL,
       esp.INT_A_INT_B,
       esp.preyNtree, #1
       esp.INT_A_INT_B2,
       esp.preyNtree2, #2
       esp.preyNtree3),
  use.aicc = T
)

aic_table.esp <- aic.esp$table
aic_table.esp 



#  _    _          ____ _____ _______    _______    _____ ______      _______ 
# | |  | |   /\   |  _ \_   _|__   __|/\|__   __|  / ____/ __ \ \    / / ____|
# | |__| |  /  \  | |_) || |    | |  /  \  | |    | |   | |  | \ \  / / (___  
# |  __  | / /\ \ |  _ < | |    | | / /\ \ | |    | |   | |  | |\ \/ / \___ \ 
# | |  | |/ ____ \| |_) || |_   | |/ ____ \| |    | |___| |__| | \  /  ____) |
# |_|  |_/_/    \_\____/_____|  |_/_/    \_\_|     \_____\____/   \/  |_____/ 
# The intercept is A1, 
# A2 is the effect of species (difference in occupancy between species A and species B on the logit scale), and 
# A3 is the effect of species A being absent on the occupancy of species B.  
# A3 and A4 are the covariate effects.  
# Since A1 is the intercept, qlogis(A1) is the estimate of occupancy for species A when housden and elev are zero.
# Since A2 is the difference in occupancy between species B and species A, a positive value means that occupancy of species B > occupancy of species A.  If A3 is positive, then the occupancy of species B is higher when species A is absent.  The confidence intervals will indicate if the effect is "significant".


NULL.mod <- occMod(model=list( 
  psi~SP,                
  gamma~SP,    
  epsilon~SP,            
  p~SP),   
  data=data,
  type="do.2sp.1")

# BaseMod <- occMod(model=list(  
#   psi~SP,               
#   epsilon~SP, 
#   p~SP+INT_o),   
#   data=data,
#   type="do.2sp.1")

BaseMod1 <- occMod(model=list(  
  psi~SP+INT, 
  gamma~SP,#you forgot gamma here which is why the mod avg function is not working
  epsilon~SP, 
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")
# 
# summary(BaseMod1)
# BaseMod1$beta
# BaseMod1$derived
# printModResults(BaseMod1)
# print(BaseMod1)
# lapply(BaseMod1$real, FUN = unique)
# lapply(BaseMod1$derived, FUN = unique)

#Using the top model from epsilon  psi(SP)gamma(SP+INT_A+INT_B)epsilon(SP+INT_A+INT_B)p(SP+INT_o)

# gamma~SP+INT_A+INT_B, 
# gamma~SP, 
# epsilon~SP+prey_count+treecover


#Prey model
Prey.mod <- occMod(model=list( 
  psi~SP+prey_count+elevation,            
  gamma~SP+INT_A+INT_B,   
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

Prey.mod1 <- occMod(model=list( 
  psi~SP+prey_count+elevation,            
  gamma~SP,   
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

#Bottom Model
Bottom.mod <- occMod(model=list( 
  psi~SP+prey_count+treecover+strmden+elevation,               
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

Bottom.mod1 <- occMod(model=list( 
  psi~SP+prey_count+treecover+strmden+elevation,               
  gamma~SP,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")


#Habitat Model
Habitat.mod <- occMod(model=list( 
  psi~SP+treecover+strmden+elevation,                
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")

Habitat.mod1 <- occMod(model=list( 
  psi~SP+treecover+strmden+elevation,                
  gamma~SP,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")


#Human model
Human.mod <- occMod(model=list( 
  psi~SP+housden +elevation,            
  gamma~SP+INT_A+INT_B,   
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")

Human.mod1 <- occMod(model=list( 
  psi~SP+housden +elevation,            
  gamma~SP,   
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")


#Top down model
Topdown.mod <- occMod(model=list( 
  psi~SP+INT+housden+elevation,             
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,           
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

Topdown.mod1 <- occMod(model=list( 
  psi~SP+INT+housden+elevation,             
  gamma~SP,    
  epsilon~SP+prey_count+treecover,           
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")


#Tiger top down model
Tiger.mod <- occMod(model=list( 
  psi~SP+INT+elevation,                
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

Tiger.mod1 <- occMod(model=list( 
  psi~SP+INT+elevation,                
  gamma~SP,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

#Human shield model
HumanShield.mod <- occMod(model=list( 
  psi~SP+INT*housden+elevation,                
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,  
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")

HumanShield.mod1 <- occMod(model=list( 
  psi~SP+INT*housden+elevation,                
  gamma~SP,    
  epsilon~SP+prey_count+treecover,  
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")


#resource competition model
ResourceComp.mod <- occMod(model=list(  
  psi~SP+INT*prey_count+elevation,                 
  gamma~SP+INT_A+INT_B,  
  epsilon~SP+prey_count+treecover,           
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

ResourceComp.mod1 <- occMod(model=list(  
  psi~SP+INT*prey_count+elevation,                 
  gamma~SP,  
  epsilon~SP+prey_count+treecover,           
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

#global model
Global.mod <- occMod(model=list( 
  psi~SP+INT*prey_count+INT*housden+treecover+strmden+elevation,               
  gamma~SP+INT_A+INT_B,  
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")


aic.ALL <- createAicTable(
  list(NULL.mod, #
       BaseMod1, #
        Prey.mod, #
       Prey.mod1, #
       Bottom.mod, #
       Bottom.mod1, #6
       Habitat.mod,#7
       Habitat.mod1,#8
       Human.mod, #
       Human.mod1, #10
       Topdown.mod, #
       Topdown.mod1,
       Tiger.mod, #
       Tiger.mod1, #14
       HumanShield.mod, #
       HumanShield.mod1, #16
       ResourceComp.mod, #
       ResourceComp.mod1,
       Global.mod),#
  use.aicc = T
)


 aic_table.ALL <- aic.ALL$table
aic_table.ALL #8,14,16,10,7,6

#Relative likelihood
Rel.like<-exp(-0.5*aic_table.ALL$DAICc)

aic_table.ALL$RelLik<-round(Rel.like, digits=3)
aic_table.ALL



#elevation
elevation_seq <- seq(min(cov1$elevation),
                     max(cov1$elevation),length.out = 100)

SP_seq <- c(rep(0,length(elevation_seq)),rep(1,length(elevation_seq)))
INT_seq<-c(rep(0,length(SP_seq)),rep(1,length(SP_seq)))

# Create a new data frame for prediction, holding other covariates at their mean
newdata.elevation <- data.frame(elevation = elevation_seq,
                                SP=SP_seq,
                                INT=INT_seq,
                                treecover = mean(cov1$treecover, na.rm =T ),
                                strmden = mean(cov1$strmden, na.rm = T),
                                housden= mean(cov1$housden, na.rm=T),
                                prey_count = mean(cov1$prey_count, na.rm= T))

#head(newdata.elevation)


#forest cover
treecover_seq <- seq(min(cov1$treecover),
                     max(cov1$treecover),length.out = 100)

SP_seq <- c(rep(0,length(elevation_seq)),rep(1,length(elevation_seq)))
INT_seq<-c(rep(0,length(SP_seq)),rep(1,length(SP_seq)))

# Create a new data frame for prediction, holding other covariates at their mean
newdata.treecover <- data.frame(treecover = treecover_seq,
                                SP=SP_seq,
                                INT=INT_seq,
                                elevation = mean(cov1$elevation, na.rm =T ),
                                strmden = mean(cov1$strmden, na.rm = T),
                                housden= mean(cov1$housden, na.rm=T),
                                prey_count = mean(cov1$prey_count, na.rm= T))

#head(newdata.treecover)

#stream density
streamden_seq <- seq(min(cov1$strmden),
                     max(cov1$strmden),length.out = 100)

SP_seq <- c(rep(0,length(elevation_seq)),rep(1,length(elevation_seq)))
INT_seq<-c(rep(0,length(SP_seq)),rep(1,length(SP_seq)))

# Create a new data frame for prediction, holding other covariates at their mean
newdata.strmden <- data.frame(strmden= streamden_seq,
                                SP=SP_seq,
                                INT=INT_seq,
                                elevation = mean(cov1$elevation, na.rm =T ),
                                treecover = mean(cov1$treecover, na.rm = T),
                                housden= mean(cov1$housden, na.rm=T),
                                prey_count = mean(cov1$prey_count, na.rm= T))

##head(newdata.strmden)

#housingdensity
housden_seq <- seq(min(cov1$housden),
                     max(cov1$housden),length.out = 100)

SP_seq <- c(rep(0,length(elevation_seq)),rep(1,length(elevation_seq)))
INT_seq<-c(rep(0,length(SP_seq)),rep(1,length(SP_seq)))

# Create a new data frame for prediction, holding other covariates at their mean
newdata.housden <- data.frame(housden= housden_seq,
                              SP=SP_seq,
                              INT=INT_seq,
                              elevation = mean(cov1$elevation, na.rm =T ),
                              treecover = mean(cov1$treecover, na.rm = T),
                              strmden= mean(cov1$strmden, na.rm=T),
                              prey_count = mean(cov1$prey_count, na.rm= T))

#head(newdata.housden)

#preycount
prey_seq <- seq(min(cov1$prey_count),
                   max(cov1$prey_count),length.out = 100)

SP_seq <- c(rep(0,length(elevation_seq)),rep(1,length(elevation_seq)))
INT_seq<-c(rep(0,length(SP_seq)),rep(1,length(SP_seq)))

# Create a new data frame for prediction, holding other covariates at their mean
newdata.preycount<- data.frame(prey_count= prey_seq,
                              SP=SP_seq,
                              INT=INT_seq,
                              elevation = mean(cov1$elevation, na.rm =T ),
                              treecover = mean(cov1$treecover, na.rm = T),
                              strmden= mean(cov1$strmden, na.rm=T),
                            housden = mean(cov1$housden, na.rm= T))

#  __  __    ____    _____    ______   _                    __      __  ______   _____                _____   ______     ______   _    _   _   _    _____   _______   _____    ____    _   _   
# |  \/  |  / __ \  |  __ \  |  ____| | |            /\     \ \    / / |  ____| |  __ \      /\      / ____| |  ____|   |  ____| | |  | | | \ | |  / ____| |__   __| |_   _|  / __ \  | \ | |  
# | \  / | | |  | | | |  | | | |__    | |           /  \     \ \  / /  | |__    | |__) |    /  \    | |  __  | |__      | |__    | |  | | |  \| | | |         | |      | |   | |  | | |  \| |  
# | |\/| | | |  | | | |  | | |  __|   | |          / /\ \     \ \/ /   |  __|   |  _  /    / /\ \   | | |_ | |  __|     |  __|   | |  | | | . ` | | |         | |      | |   | |  | | | . ` |  
# | |  | | | |__| | | |__| | | |____  | |____     / ____ \     \  /    | |____  | | \ \   / ____ \  | |__| | | |____    | |      | |__| | | |\  | | |____     | |     _| |_  | |__| | | |\  |  
# |_|  |_|  \____/  |_____/  |______| |______|   /_/    \_\     \/     |______| |_|  \_\ /_/    \_\  \_____| |______|   |_|       \____/  |_| \_|  \_____|    |_|    |_____|  \____/  |_| \_|  


dynocc2sp.Modavg<-function (aic.tab, param = "psi", param2 = "psiA",parmgrp = "real", index = 1:nrow(aic.tab$table), 
                            conf = 0.95, predict = FALSE, newdata, replaceNaN = TRUE) 
{
  linkfn = "logit"
  if (length(grep("phi|nu", param)) > 0) 
    linkfn = "exp"
  aic.table = list(table = aic.tab$table[index, ], models = aic.tab$models[index])
  class(aic.table) = "aic.table"
  jaic = grep("AIC", colnames(aic.table$table))[1]
  jdaic = grep("AIC", colnames(aic.table$table))[2]
  aic.table$table$DAIC = aic.table$table[, jaic] - min(aic.table$table[, 
                                                                       jaic])
  aic.table$table$modlike = exp(-aic.table$table[, jdaic]/2)
  aic.table$table$wgt = aic.table$table$modlike/sum(aic.table$table$modlike)
  #added this to rename p columns to est and se from V1 and V2
  for (i in seq_along(aic.table$models)) {
    names(aic.table$models[[i]]$beta$p) <- c("est","se")
  }
  
  #added this to rename p columns to est and se from V1 and V2
  for (i in seq_along(aic.table$models)) {
    names(aic.table$models[[i]]$beta$epsilon) <- c("est","se")
    print(i)
  }
  
  #not needed below
  # #added this to rename gamma intercept to Ab, which is the same as gamAB where not parameterized
  # for (i in seq_along(aic.table$models)) {
  #   gamma_block <- aic.table$models[[i]]$beta$gamma
  #   
  #   # only proceed if gamma_block is a 2D object (matrix or data.frame)
  #   if (!is.null(gamma_block) && !is.null(dim(gamma_block))) {
  #     rn <- rownames(gamma_block)
  #     rn[rn == "B1_gamAB[1]"] <- "B1_gamAb[1]"
  #     rownames(aic.table$models[[i]]$beta$gamma) <- rn
  #   }
  # }
  
  
  
  # if (!predict) {
  #   i = which(names(aic.table$models[[1]]) == parmgrp)
  #   nr = nrow(get(param, aic.table$models[[1]][[i]]))
  #   est = array(unlist(lapply(aic.table$models, function(xx) get(param, 
  #                                                                xx[[i]])[, 1])), dim = c(nr, length(aic.table$models)))
  #   se = array(unlist(lapply(aic.table$models, function(xx) get(param, 
  #                                                               xx[[i]])[, 2])), dim = c(nr, length(aic.table$models)))
  # }
  # else {
  pred = lapply(aic.table$models, function(xx) predict(xx, 
                                                       newdata, param = param))
  est = array(unlist(lapply(pred, function(xx) xx[, 1])), 
              dim = c(nrow(newdata), length(aic.table$models)))
  se = array(unlist(lapply(pred, function(xx) xx[, 2])), 
             dim = c(nrow(newdata), length(aic.table$models)))
  #}
  if (replaceNaN) 
    se[is.na(se)] = 0
  ma = est %*% aic.table$table$wgt
  ma.se = sqrt((se^2 + (est - as.vector(ma))^2) %*% aic.table$table$wgt)
  if (linkfn == "logit") {
    logit.est = log(ma/(1 - ma))
    logit.se = ma.se
    on.bounds = ma %in% 1 | ma %in% 0
    logit.se[on.bounds] = 0
    logit.se[!on.bounds] = ma.se[!on.bounds]/(ma[!on.bounds] * 
                                                (1 - ma[!on.bounds]))
  }
  else {
    logit.est = log(ma)
    logit.se = ma.se
    logit.se[ma == 0] = 0
  }
  alpha = (1 - conf)/2
  z = -qnorm(alpha)
  lower = as.matrix(sapply(z, function(zz) logit.est - zz * 
                             logit.se))
  colnames(lower) = paste("lower", conf, sep = "_")
  upper = as.matrix(sapply(z, function(zz) logit.est + zz * 
                             logit.se))
  colnames(upper) = paste("upper", conf, sep = "_")
  if (linkfn == "logit") {
    lower = plogis(lower)
    upper = plogis(upper)
  }
  else {
    lower = exp(lower)
    upper = exp(upper)
  }
  result = data.frame(est = ma, se = ma.se, lower, upper)
  if (!predict) {
    rownames(result) = rownames(get(param, aic.table$models[[1]][[i]]))
  }
  else {
    rownames(result) = rownames(newdata)
  }
  return(result)
}




#Model average for all psi (the design matrix is included in the newdata)
#we only have to look at the specific combo of SP and Int to get the psiA, psiBA, psiBa
#See below
psi_avg <- modAvg(aic.tab = aic.ALL, 
                             param = "psi",
                             #parmgrp = "derived",
                             predict=T,
                             newdata=newdata.elevation,
                             index = 1:nrow(aic.ALL$table),
                             conf = 0.95)
# nrow(psi_avg) #400
Global.mod$dmat$psi#inspect the design matrix to fill out below

psi.avg<-cbind(psi_avg,newdata.elevation)#combind newdata with psi_avg
# head(psi.avg)
# nrow(psi.avg)

psi.avg$Psi<-ifelse(psi.avg$SP==0&psi.avg$INT==0,"psiA",
                    ifelse(psi.avg$SP==1&psi.avg$INT==1,"psiBa",
                    ifelse(psi.avg$SP==1&psi.avg$INT==0,"psiBA",NA)))
# remove the SP=0 INT=1 case that I do not think makes sense
psi.avg<-psi.avg[!is.na(psi.avg$Psi)==T,]

# head(psi.avg)

#and plot the model averaged predictions :)
mean_elev <- mean(csv1$elevation,na.rm=T) #2302.298
sd_elev <- sd(csv1$elevation,na.rm=T) #1049.267
mean_treecov<-mean(csv1$treecover,na.rm=T) #0.8993651
sd_treecov<-sd(csv1$treecover,na.rm=T) #0.1319367
mean_strmden<-mean(csv1$strmden,na.rm=T) #3.748604
sd_strmden<-sd(csv1$strmden,na.rm=T) #0.564268
mean_housden<-mean(csv1$housden,na.rm=T) #10.32206
sd_housden<-sd(csv1$housden,na.rm=T) #20.55243
mean_prey<-mean(csv1$prey_count,na.rm=T) # 17.21864
sd_prey<-sd(csv1$prey_count,na.rm=T) #24.37


#Plotting
library(ggplot2)
# Set the factor level order so they plot in desired sequence
psi.avg$Psi <- factor(psi.avg$Psi, levels = c("psiA", "psiBA", "psiBa"))


ggplot(data = psi.avg, aes(x = elevation * sd_elev + mean_elev, y = est)) + 
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95, fill = Psi), alpha = 0.6, lwd = 1) +
  geom_line(aes(colour = Psi), lwd = 1) + 
  scale_fill_manual(
    values = c("psiA" = "#c86b4a", 
               "psiBa" = "#eed7a1", 
               "psiBA" = "#D55E00"),
    labels = c(
      "psiA" = expression(psi[Tiger]), 
      "psiBa" = expression(psi[Leopard*"|"*Tiger_Absence]), 
      "psiBA" = expression(psi[Leopard*"|"*Tiger_Presence])
    )
  ) +
  scale_colour_manual(
    values = c("psiA" = "#c86b4a", 
               "psiBa" = "#eed7a1", 
               "psiBA" = "#D55E00"),
    labels = c(
      "psiA" = expression(psi[Tiger]), 
      "psiBa" = expression(psi[Leopard*"|"*Tiger_Absence]), 
      "psiBA" = expression(psi[Leopard*"|"*Tiger_Presence])
    )
  ) +
  labs(
    x = "Elevation (m)", 
    y = expression("Habitat-use (" * psi * ")") #Removed "Initial"
  ) +
  coord_sf(ylim = c(0.05, 0.7)) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.line = element_line(color = "grey50"),
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)  #
  )


#avg code for gamma 
#need to create newdata for gamma including all the vars from global gamma model
#SP+INT_A+INT_B
# aic.table$models[[5]]$dmat$gamma.dm
# aic.table$models[[1]]$beta$gamma

gam.newdata<-data.frame(SP=c(0,1,0,1,0,1,0,1),
                        INT_A=c(0,0,0,0,1,1,1,1),
                        INT_B=c(0,0,1,1,0,0,1,1))
# aic.table$table
# Global.mod$beta$p

gamma_avg <- modAvg(aic.tab = aic.ALL, 
param = "gamma",
#parmgrp = "derived",
predict=T,
newdata=gam.newdata,
index = 1:nrow(aic.ALL$table),
conf = 0.95)


#yay works! Just need to check it
# unique(Global.mod$dmat$gamma)
# 
# Global.mod$beta$gamma
# Habitat.mod1$beta$gamma

gamma_avg<-cbind(gamma_avg,gam.newdata)
gamma_avg$Gam<-ifelse(gamma_avg$SP==0&gamma_avg$INT_A==1&gamma_avg$INT_B==0,"Tiger",#gamAB
                  ifelse(gamma_avg$SP==0&gamma_avg$INT_A==0&gamma_avg$INT_B==0,"gamAb",
                         ifelse(gamma_avg$SP==1&gamma_avg$INT_A==1&gamma_avg$INT_B==0,"Leopard",#gamBAA
                                ifelse(gamma_avg$SP==1&gamma_avg$INT_A==1&gamma_avg$INT_B==1,"gamBaA",NA))))

gamma_avg<-gamma_avg[!is.na(gamma_avg$Gam)==T,]
gamma_avg<-gamma_avg[!gamma_avg$Gam=="gamAb",]
gamma_avg<-gamma_avg[!gamma_avg$Gam=="gamBaA",]


ggplot(gamma_avg, aes(x = Gam, y = est, fill = Gam)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_0.95, ymax = upper_0.95),
                width = 0.2, size = 0.7) +
  labs(
    x     = "",
    y     = "Colonization Probability(± 95% CI)",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x    = element_text(angle = 25, hjust = 1),
  )




##compute model-averaged estimates for parameters appearing in top
##models
Topmods <- aic_table.ALL[c(1:6),] #selecting the top six mods
Topmods

Topmodlist <- list()
Topmodlist[[1]] <- occMod(model=list( 
  psi~SP+treecover+strmden+elevation,                
  gamma~SP,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")

Topmodlist[[2]] <- occMod(model=list( 
  psi~SP+INT+elevation,                
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

Topmodlist[[3]] <-occMod(model=list( 
  psi~SP+INT*housden+elevation,                
  gamma~SP,    
  epsilon~SP+prey_count+treecover,  
  p~SP+INT_o),   
  data=data,
  type="do.2sp.1")
  
Topmodlist[[4]] <- occMod(model=list( 
  psi~SP+housden +elevation,            
  gamma~SP,   
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")
  
  
Topmodlist[[5]]  <- occMod(model=list( 
  psi~SP+treecover+strmden+elevation,                
  gamma~SP+INT_A+INT_B,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o), 
  data=data,
  type="do.2sp.1")

Topmodlist[[6]]  <- occMod(model=list( 
  psi~SP+prey_count+treecover+strmden+elevation,               
  gamma~SP,    
  epsilon~SP+prey_count+treecover,            
  p~SP+INT_o),  
  data=data,
  type="do.2sp.1")

Modnames <- c("Habitat.mod1",#8
              "Tiger.mod1", #14
              "HumanShield.mod1",#16
              "Human.mod1", #10
              "Habitat.mod",#7
              "Bottom.mod1")#6
sapply(Topmodlist, is.null)
sapply(Topmodlist, class)
str(Topmodlist[[1]])

aic.ALL <- createAicTable(list(Habitat.mod1,#8
                                 Tiger.mod1, #14
                                 HumanShield.mod1,#16
                                 Human.mod1, #10
                                 Habitat.mod,#7
                                 Bottom.mod1))

aic.table <- aic.ALL$table
aic.table
aic.ALL$chat
??modAvg
names(Habitat.mod1$real)

str(Habitat.mod1)


aic.test <- createAicTable(list(Habitat.mod1, Tiger.mod1))
modAvg(aic.tab = aic.test, param = "psiBA", parmgrp = "real")

get("psi",Tiger.mod1$beta)[1,]
View(modAvg)


nrow(newdata.elevation)
psiA_avg <- dynocc2sp.Modavg(aic.tab = aic.ALL, 
                    param = "psi",
                    param2="psiA",
                    #parmgrp = "derived",
                   predict=T,
                   newdata=newdata.elevation,
                    index = 1:nrow(aic.ALL$table),
                    conf = 0.95)
head(psiA_avg)
psiBA_avg <- modAvg(aic.tab = aic.ALL, 
                   param = "psiBA",
                   parmgrp = "real",
                   index = 1:nrow(aic.ALL$table),
                   conf = 0.95)

psiBa_avg <- modAvg(aic.tab = aic.ALL, 
                    param = "psiBa",
                   # parmgrp = "real",
                    #index = 1:nrow(aic.ALL$table),
                    predict=T,
                    newdata=newdata.elevation,
                    conf = 0.95)

rBA_avg <- modAvg(aic.tab = aic.ALL, 
                    param = "rBA",
                    parmgrp = "real",
                    index = 1:nrow(aic.ALL$table),
                    conf = 0.95)

rBa_avg <- modAvg(aic.tab = aic.ALL, 
                  param = "rBa",
                  parmgrp = "real",
                  index = 1:nrow(aic.ALL$table),
                  conf = 0.95)

Human.mod$beta
plot(cov1$elevation,Human.mod$real$psiA[,1])
plot(cov1$housden,Human.mod$real$psiA[,1])

plot(cov1$elevation,psiBA_avg$est)
plot(cov1$housden,psiBA_avg$est)

#ModelAveragePrediction

#elevation
elevation_seq <- seq(min(cov1$elevation),
                     max(cov1$elevation),length.out = 100)

SP_seq <- c(rep(0,length(elevation_seq)),rep(1,length(elevation_seq)))
INT_seq<-c(rep(0,length(SP_seq)),rep(1,length(SP_seq)))

# Create a new data frame for prediction, holding other covariates at their mean
newdata.elevation <- data.frame(elevation = elevation_seq,
                                SP=SP_seq,
                                INT=INT_seq,
                                treecover = mean(cov1$treecover, na.rm =T ),
                                strmden = mean(cov1$strmden, na.rm = T),
                                housden= mean(cov1$housden, na.rm=T),
                                prey_count = mean(cov1$prey_count, na.rm= T))

head(newdata.elevation)

Habitat.mod1$beta
# Predict using the top model
elevation_pred <- predict(Habitat.mod1, newdata = newdata.elevation, param = "psi", appendData=T)
elevation_pred2 <- predict(Tiger.mod1, newdata = newdata.elevation, param = "psi", appendData=T)
elevation_pred3 <- predict(HumanShield.mod1, newdata = newdata.elevation, param = "psi", appendData=T)
elevation_pred4 <- predict(Human.mod1, newdata = newdata.elevation, param = "psi", appendData=T)
elevation_pred5 <- predict(Habitat.mod, newdata = newdata.elevation, param = "psi", appendData=T)
elevation_pred6 <- predict(Bottom.mod1, newdata = newdata.elevation, param = "psi", appendData=T)


??rPresence
?modAvg

# Create an empty data frame to store results
results.elevation <- data.frame(elevation = rep(elevation_seq,4), 
                                est = numeric(length(elevation_seq)),
                                se  = numeric(length(elevation_seq)),
                                lower_0.95 = numeric(length(elevation_seq)), 
                                upper_0.95 = numeric(length(elevation_seq))
)
results.elevation 

# Fill in the results.elevation 
for (i in 1:length(elevation_pred6)) {
  results.elevation $est <- elevation_pred6$est
  results.elevation $se <- elevation_pred6$se
  results.elevation $lower_0.95<- elevation_pred6$lower_0.95
  results.elevation $upper_0.95<- elevation_pred6$upper_0.95
}

results.elevation


# Extract predictions and standard errors
elevation_pred.est <- results.elevation$est
elevation_pred.est2 <- results.elevation$est
elevation_pred.est3 <- results.elevation$est
elevation_pred.est4 <- results.elevation$est
elevation_pred.est5 <- results.elevation$est
elevation_pred.est6 <- results.elevation$est

elevation_pred.se <- results.elevation$se
elevation_pred.se2 <- results.elevation$se
elevation_pred.se3 <- results.elevation$se
elevation_pred.se4 <- results.elevation$se
elevation_pred.se5 <- results.elevation$se
elevation_pred.se6 <- results.elevation$se

# Compute the model-averaged prediction
AICc.wt <- aic.table[,9]

Mod.Avg.Pred.elevation <- elevation_pred.est*AICc.wt[1]+
                          elevation_pred.est2*AICc.wt[2]+
                          elevation_pred.est3*AICc.wt[3]+
                          elevation_pred.est4*AICc.wt[4]+
                          elevation_pred.est5*AICc.wt[5]+
                          elevation_pred.est6*AICc.wt[6]

# Step 4: Compute the unconditional variance and SE
elev_var_avg <- AICc.wt[1] * (elevation_pred.se^2 + (elevation_pred.est- Mod.Avg.Pred.elevation)^2) +
  AICc.wt[2] * (elevation_pred.se2^2 + (elevation_pred.est2- Mod.Avg.Pred.elevation)^2)+
  AICc.wt[3]* (elevation_pred.se3^2 + (elevation_pred.est3- Mod.Avg.Pred.elevation)^2)+
  AICc.wt[4] *(elevation_pred.se4^2 + (elevation_pred.est4- Mod.Avg.Pred.elevation)^2)+
  AICc.wt[5]* (elevation_pred.se5^2 + (elevation_pred.est5- Mod.Avg.Pred.elevation)^2)+
  AICc.wt[6]* (elevation_pred.se6^2 + (elevation_pred.est6- Mod.Avg.Pred.elevation)^2)
  
  
elev_se_avg <- sqrt(elev_var_avg)

#CI
# Compute 95% Confidence Intervals
elev_LCL <- Mod.Avg.Pred.elevation - (1.96 * elev_se_avg)
elev_UCL <- Mod.Avg.Pred.elevation + (1.96 * elev_se_avg)


# model-averaged prediction and its SE
results.elevation.all <- data.frame(elevation = rep(elevation_seq,4), 
                                    est = Mod.Avg.Pred.elevation,
                                    se = elev_se_avg,
                                    lcl = elev_LCL, 
                                    ucl = elev_UCL
)
head(results.elevation.all)


Habitat.mod1$beta
Tiger.mod1$beta
HumanShield.mod1$beta
Human.mod1$beta
Habitat.mod$beta
Bottom.mod1$beta


??modAvg
#  defining a function to print the real estimates of a model in
#   a table.  Only the 1st site of each parameter is printed.

prnt <- function(mod) {
  x=matrix(unlist(lapply(mod$real,function(a) a[1,])),ncol=4,byrow=F)
  rownames(x)=rep(names(mod$real))
  colnames(x)=colnames(mod$real[[1]])
  cat('model:',mod$modname,'\n'); print(unique(x))
}

prnt(Habitat.mod1)

View(lapply(Habitat.mod1$real, FUN = unique))
lapply(Habitat.mod1$real, FUN = unique)

#1
Habitat.mod1$beta$psi
Habitat.mod1$beta$epsilon

dplyr::distinct(fitted(Habitat.mod1, param = "epsilon"))
coef(Habitat.mod1, param='psi', prob = 0.95)
coef(Habitat.mod1, param='gamma', prob = 0.95)
coef(Habitat.mod1, param='epsilon', prob = 0.95)
printModResults(Habitat.mod1)


#2
prnt(Tiger.mod1)
Tiger.mod1$beta$psi
coef(Tiger.mod1, param='psi', prob = 0.95)
coef(Tiger.mod1, param='gamma', prob = 0.95)
coef(Tiger.mod1, param='epsilon', prob = 0.95)

#3
prnt(HumanShield.mod1)
coef(HumanShield.mod1, param='psi', prob = 0.95)
coef(HumanShield.mod1, param='gamma', prob = 0.95)
coef(HumanShield.mod1, param='epsilon', prob = 0.95)

#4
prnt(Human.mod1)
coef(Human.mod1, param='psi', prob = 0.95)
coef(Human.mod1, param='gamma', prob = 0.95)
coef(Human.mod1, param='epsilon', prob = 0.95)

#5
prnt(Habitat.mod)
coef(Habitat.mod, param='psi', prob = 0.95)
coef(Habitat.mod, param='gamma', prob = 0.95)
coef(Habitat.mod, param='epsilon', prob = 0.95)


#6
prnt(Bottom.mod1)
coef(Bottom.mod1, param='psi', prob = 0.95)
coef(Bottom.mod1, param='gamma', prob = 0.95)
coef(Bottom.mod1, param='epsilon', prob = 0.95)
Bottom.mod1$beta


#model averaging
names(aic_table.ALL)
ncol(aic_table.ALL)

AICc.wt <- aic_table.ALL["wgt"]




























































































































