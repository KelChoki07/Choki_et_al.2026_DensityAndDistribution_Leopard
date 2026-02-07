
library(secr);library(lubridate);library(raster)
library(sp);library(sf);library(terra); library(ggplot2)



# Load traps and capture history built from the initial codes
cap_hist <- read.capthist(
  "capture_history_120.txt",
  "traps.txt",
  detector = "proximity",
  sep = "\t"
)

summary(cap_hist)
initialsigma <- RPSV(cap_hist, CC = TRUE) #4332.366 [initialsigma * 3=12997.1]
initialsigma

# Plot the habitat mask
plot(hab_mask, col = "lightgray", main = "Habitat Mask with Polygon Overlay")

# Overlay the habitat polygon
plot(spc_mask, add = TRUE, border = "black")  # Use `add = TRUE` to overlay
points(traps(cap_hist), col = "darkgreen", pch = 16, cex=0.5) # Plot trap locations
plot(cap_hist, add=TRUE, tracks = TRUE, title = "", main = "leopard capture history", cex=0.8)  


m <- unlist(moves(cap_hist))
par(mar = c(3.2,4,1,1), mgp = c(2.1,0.6,0)) # reduce margins
hist(m,  xlab = "Movement m", main = "")
# ______       _               _    _                ___  ___            _        _ 
# |  _  \     | |             | |  (_)               |  \/  |           | |      | |
# | | | | ___ | |_  ___   ___ | |_  _   ___   _ __   | .  . |  ___    __| |  ___ | |
# | | | |/ _ \| __|/ _ \ / __|| __|| | / _ \ | '_ \  | |\/| | / _ \  / _` | / _ \| |
# | |/ /|  __/| |_|  __/| (__ | |_ | || (_) || | | | | |  | || (_) || (_| ||  __/| |
# |___/  \___| \__|\___| \___| \__||_| \___/ |_| |_| \_|  |_/ \___/  \__,_| \___||_|


#12/5/2024   12:36 am
args <- list("cap_hist", buffer = 15000, mask = hab_mask,
             detectfn = 'HN', ncores=25)

models <-  list(
  list(D ~ 1, g0 ~ 1, sigma ~ 1),   
  list(D ~ 1, g0 ~ b, sigma ~ 1)            
)
detection_fit <- list.secr.fit(model = models, constant = args)

# Model selection
det.modelsel <- AIC(detection_fit, criterion = 'AICc')
det.modelsel


Rel.like1<-exp(-0.5*AIC(detection_fit)$dAIC)

det.modelsel$RelLik1<-round(Rel.like1, digits=3)

det.modelsel

#plotting detection probability
plot(detection_fit$fit1, limits = T,  xv = 0:15000, ylim = c(0, 0.006), main= "Detection Probability") 

top_mod_det$realnames


models <-  list(
  list(D ~ housden * tigerden + prey_count * tigerden + treecover + strmden + elevation, g0 ~ 1, sigma ~ 1),   
  list(D ~ housden * tigerden + prey_count * tigerden + treecover + strmden + elevation, g0 ~ b, sigma ~ 1)            
)
detection_fit <- list.secr.fit(model = models, constant = args) #1;40pm 5/20/25

# Model selection
det.modelsel <- AIC(detection_fit, criterion = 'AICc')
det.modelsel

#Detection probability estimates
summary(detection_fit$fit1)
plogis(-5.805697) #0.003 Just to double check
plogis(0.07614186 ) #0.519
plogis(-5.954932 ) #0.00258
plogis(-5.656461 ) #0.00348

initialsigma <- RPSV(cap_hist, CC = TRUE) #4332.366 [initialsigma * 3=12997.1]


out0 <- secr.fit(
  cap_hist,
  model = list(D ~ 1, g0 ~ 1, sigma ~ 1),
  mask = NULL, # can upload a shapefile or CSV mask (this needs workaround)
  buffer = 15000, # make it sufficiently large to encompass all activity centres
  # buffer = initialsigma * 3, # normally for half-normal models
  ncores = 7 # depends on your machine
)
# 1st-Warning message:
# In bias.D(buffer, temptrps, detectfn = output$detectfn, detectpar = dpar,  :
# bias.D() does not allow for variable effort (detector usage)
# out0 #D=0.91

# To check if the density estimate plateaus with the chosen buffer width
esa.plot(out0)
# abline(v = 3 * initialsigma, lty = 2, col = 'red')
abline(v = 15000, lty = 2, col = 'red')



# ______                   _  _            ___  ___            _        _      
# |  _  \                 (_)| |           |  \/  |           | |      | |     
# | | | | ___  _ __   ___  _ | |_  _   _   | .  . |  ___    __| |  ___ | | ___ 
# | | | |/ _ \| '_ \ / __|| || __|| | | |  | |\/| | / _ \  / _` | / _ \| |/ __|
# | |/ /|  __/| | | |\__ \| || |_ | |_| |  | |  | || (_) || (_| ||  __/| |\__ \
# |___/  \___||_| |_||___/|_| \__| \__, |  \_|  |_/ \___/  \__,_| \___||_||___/
#                                   __/ |                                      
#                                  |___/                                       
#

tigerden  <- raster("TigerDensity_1km.tif")
tigerden<-raster(aggregate(rast(tigerden),fact=5,fun="mean",ncores=10))
spc_mask_o<-as.polygons(rast(tigerden))
spc_mask_o<-st_as_sf(aggregate(spc_mask_o),1,dissolve=T)

#spc_mask_o <- st_read("tiger_boundary_AC.shp")
plot(spc_mask_o)

crs(spc_mask_o)

mask_spdf <- as(spc_mask_o, "Spatial") #class: SpatialPolygonsDataFrame

spc_mask <- spChFIDs(mask_spdf, "SA")

# Make SECR mask
hab_mask <- make.mask(
  traps(cap_hist),
  buffer = 15000,
  type = "trapbuffer",
  spacing = 5000, # 5km
  poly = spc_mask,
  poly.habitat = TRUE,
  keep.poly = TRUE
)
str(hab_mask)
elevation <- raster("elevation5km.tif")
strmden <- raster("treecover5km.tif")
housden <- raster("housden5km.tif")
strmden   <- raster("strmden5km.tif")
tigerden  <- raster("tigerden5km.tiff")

crs(elevation, describe=T)


mshp<-st_read("prey_2025_feb.shp") #there are 58 stations no det muntjac

mshp$prey_count
mshp$prey_count<-scale(as.numeric(mshp$prey_count))
mshp

# mshp$prey_count<-as.numeric(mshp$prey_count)
# mean(mshp$prey_count, na.rm=T) #17.93
# sd(mshp$prey_count, na.rm=T) #26.63

hab_mask_with_covariates <- st_join(st_as_sf(hab_mask,coords = c("x","y"),
                                             crs=crs(mshp)), mshp["prey_count"],
                                    join = st_nearest_feature) #34470 obs


str(hab_mask_with_covariates)

hab_mask <- addCovariates(hab_mask, hab_mask_with_covariates,columns="prey_count")

#34470 obs
sum(is.na(covariates(hab_mask))) #2170 (127 + 2006)

head(attr(hab_mask, "covariates"))
sum(is.na(attr(hab_mask, "covariates")$prey_count)) #2006

length(attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$elevation), "elevation"] )
length(attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$treecover), "treecover"])
length(attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$housden), "housden"] )
length(attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$strmden), "strmden"] )
length(attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$tigerden), "tigerden"] )
length(attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$prey_count), "prey_count"])

attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$elevation), "elevation"] <- mean(attr(hab_mask, "covariates")$elevation,na.rm=T)
attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$treecover), "treecover"] <- mean(attr(hab_mask, "covariates")$treecover,na.rm=T)
attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$housden), "housden"] <- mean(attr(hab_mask, "covariates")$housden,na.rm=T)
attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$strmden), "strmden"] <- mean(attr(hab_mask, "covariates")$strmden,na.rm=T)
attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$tigerden), "tigerden"] <- mean(attr(hab_mask, "covariates")$tigerden,na.rm=T)
# attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$prey_count), "prey_count"] <- 0
attr(hab_mask, "covariates")[is.na(attr(hab_mask, "covariates")$prey_count), "prey_count"] <-   mean(attr(hab_mask, "covariates")$prey_count,na.rm=T)



# ______                   _  _            ___  ___            _        _      
# |  _  \                 (_)| |           |  \/  |           | |      | |     
# | | | | ___  _ __   ___  _ | |_  _   _   | .  . |  ___    __| |  ___ | | ___ 
# | | | |/ _ \| '_ \ / __|| || __|| | | |  | |\/| | / _ \  / _` | / _ \| |/ __|
# | |/ /|  __/| | | |\__ \| || |_ | |_| |  | |  | || (_) || (_| ||  __/| |\__ \
# |___/  \___||_| |_||___/|_| \__| \__, |  \_|  |_/ \___/  \__,_| \___||_||___/
#                                   __/ |                                      
#                                  |___/                                       
#args <- list("cap_hist", mask = hab_mask,ncores=25) 

models2 <-  list(
  list(D~1, g0~1, sigma~1), 
  list(D~prey_count+elevation, g0~1, sigma~1),#prey model
  list(D~prey_count*tigerden+elevation, g0~1, sigma~1),#resouce competition
  list(D~housden*tigerden+elevation, g0~1, sigma~1),#human shield
  list(D~housden+tigerden+elevation, g0~1, sigma~1),#top down
  list(D~housden+elevation, g0~1, sigma~1),#human impacts
  list(D~tigerden+elevation, g0~1, sigma~1),#tiger impacts
  list(D~treecover+strmden+elevation, g0~1, sigma~1), #habitat
  list(D~treecover+strmden+prey_count+elevation, g0~1, sigma~1), #Bottom up
  list(D~housden*tigerden+prey_count*tigerden+
         treecover+strmden+elevation, g0~1, sigma~1) #Global
)

density_fit2 <- list.secr.fit(model = models2, constant = args)

AIC(density_fit2, criterion = 'AICc') #model 10, 8,9 #criterion = "AIC"

#Model ranking and selection 

modelsel <- AIC(density_fit2, criterion = 'AICc')
modelsel

names(modelsel)
aic_table <- modelsel[c(1,3,5,7,8, 4)]
aic_table$AIC <- round(aic_table$AIC, digit=2)
aic_table$AICwt <- round(aic_table$AICwt, digit=2)
aic_table$dAIC <- round(aic_table$dAIC, digit=2)
aic_table$logLik <- round(aic_table$logLik, digit=2)
colnames(aic_table) <- c( "Model", 
                          "K",
                          "AICc",
                          "?? AICc", "AICWt", "logLik")


library(kableExtra)

aic_table %>%
  kbl(caption = "_______ Table 2. Summary of the density model selection by AIC for the 10 competing models based on different
covariate combinations on the density _______(D=individuals/100 km2), null detection probability (g0) and the scale parameter (sigma) using half
normal (HN) detection function") %>%
  kable_classic(full_width =T , html_font = "Cambria") %>%
  footnote(general = "______ K is the number of parameters in the model.
______ AIC is the Akaike information criterion.
______ ?? AIC is the difference in the AIC between given model and the top model.
______ AICWt is the AICc weights, or the probability that of the models tested the given model fits the data best.")


summary(density_fit2$fit10)
summary(density_fit2$fit8)
summary(density_fit2$fit9)

#Model Average Density Estimate
fxsurface8 <- fx.total(density_fit2$fit8, mask = hab_mask,ncores=25)
class(fxsurface8)
plot(fxsurface8, covariate = "D.sum", poly = FALSE) #here it shows

model2.rast<-rast(fxsurface8,values=attr(fxsurface8,"covariates")$D.sum)

fxsurface9 <- fx.total(density_fit2$fit9, mask = hab_mask,ncores=25)
class(fxsurface9)
plot(fxsurface9, covariate = "D.sum", poly = FALSE) #here it shows
model3.rast<-rast(fxsurface9,values=attr(fxsurface9,"covariates")$D.sum)


fxsurface10 <- fx.total(density_fit2$fit10, mask = hab_mask,ncores=25)
class(fxsurface10)
plot(fxsurface10, covariate = "D.sum", poly = FALSE) #here it shows
model1.rast<-rast(fxsurface10,values=attr(fxsurface10,"covariates")$D.sum)
plot(model1.rast)

plot(model1.rast)

aic.out<-AIC(density_fit2, criterion = 'AICc') #model 10, 8,9 #criterion = "AIC"
names(aic.out)
mod.avg.ras<-model1.rast
values(mod.avg.ras)<-(values(model1.rast)*aic.out$AICcwt[1])+
  (values(model2.rast)*aic.out$AICcwt[2])+
  (values(model3.rast)*aic.out$AICcwt[3])

plot(mod.avg.ras)

dsurf_rast <- raster(fxsurface10, covariate = "D.sum")
class(dsurf_rast)
plot(dsurf_rast)


plot(
  mod.avg.ras * 5E3, # spacing in 1km, convert density to per 100 km2
  axes = FALSE,
  box = FALSE,
  col = rev(terrain.colors(200)),
  horizontal = TRUE,
  legend.args = list(
    text = as.expression(bquote(paste("Leopard AC density (100 km"^{-2}, ")")))
  )
)

plot(spc_mask_sf, add=T)

# Save the raster
writeRaster(dsurf_rast * 5E3, 
            "Leopard_density_5km_modelaverageprediction.tif", 
            overwrite = TRUE)


leop <- raster("Leopard_density_5km_modelaverageprediction.tif")
range(leop)


# Convert from sp to sf (as spc_mask is sp object so you need sf library)

library(sf)
spc_mask_sf <- st_as_sf(spc_mask)
plot(spc_mask_sf)
# Save as a shapefile
st_write(spc_mask_sf, "Leopard_density_5km_boundary.shp", delete_layer = TRUE)

boundary <- vect("Leopard_density_5km_boundary.shp")
plot(boundary, axes = FALSE,
     box = FALSE)

leopard <- raster("Leopard_density_5km.tif")
plot(leopard,
     axes = FALSE,
     box = FALSE,)
plot(spc_mask_sf, add=T)

# install.packages("viridis")  # Install the package if not already installed
library(viridis)  # Load the package


plot(
  dsurf_rast * 5E3, # Adjusting density scale
  axes = FALSE,
  box = FALSE,
  col = viridis(200), # Use viridis colormap for better contrast
  horizontal = TRUE,
  legend.args = list(
    text = as.expression(bquote(paste("Leopard AC density (100 km"^{-2}, ")")))
  )
)

# library(ggplot2)
# ggplot(dsurf_rast * 5E3) +
#   geom_hex() + coord_fixed() +
#   scale_fill_viridis() + theme_bw()

# plot(spc_mask, add = TRUE)

top_mod <- density_fit2$fit10
top_mod2 <- density_fit2$fit8
top_mod3 <- density_fit2$fit9

# Estimate abundance
plot(top_mod, limits = T,  xv = 0:15000, ylim = c(0, 0.006)) 


# Summary of the top model
topmod <- summary(density_fit2$fit10) #topmod <- summary(density_fits$fit.secr_elev)$predicted
topmod2 <- summary(density_fit2$fit8)
topmod3 <- summary(density_fit2$fit9)


est <- as.data.frame(topmod$coef)
est2 <- as.data.frame(topmod2$coef)
est3 <- as.data.frame(topmod3$coef)
class(est)

est
est2
est3
# topmod <- est[c(2,3,4,5)]

est$beta      <- round(est$beta, digits = 2)
est$SE.beta <- round(est$SE.beta, digits = 2)
est$lcl <- round(est$lcl, digits = 2)
est$ucl <- round(est$ucl, digits = 2)
# options(digits=3)
row.names(est) <- c( "Intercept of D(log)", "housden", "tigerden", 
                     "prey_count", "treecover", "strmden", "elevation", 
                     "housden:tigerden", "tigerden:prey_count", "g0(log)", "sigma(logit)")
colnames(est) <- c("Estimates", "SE", "lcl", "ucl")
est

est %>%
  kbl(caption = "Table 3. Maximum likelihood estimates for leopards detected by camera trapping in the Bhutan, 2022") %>%
  kable_classic(full_width =T , html_font = "Cambria")





#population estimates
popest_all1 <- region.N(
  top_mod,
  hab_mask,
  spacing = 5000,
  se.N = TRUE
)

popest_all1

polybound <- spc_mask 

popest_all2 <- region.N(
  top_mod2,
  hab_mask,
  spacing = 5000,
  se.N = TRUE
)

popest_all2

popest_all3 <- region.N(
  top_mod3,
  hab_mask,
  spacing = 5000,
  se.N = TRUE
)

popest_all3

#Population model prediction estimates
popest_all <- popest_all1*aic.out$AICcwt[1]+
              popest_all2*aic.out$AICcwt[2]+
              popest_all3*aic.out$AICcwt[3]


popest_all

polybound <- spc_mask 

# Convert to animals per 1000 km^2 (and calculate 95% CIs)
a <- data.frame(
  Study = names(hab_mask), popest_all[1, ],
  area_km2 = expanse(vect(polybound))/(1000^2)
)
a

# per 100km2 #Dhat=n/esa
a$Dest_100km2 <- a$estimate/(a$area_km2/(100))
a$Dse_100km2 <- a$SE.estimate/(a$area_km2/(100)) 
a$Dlcl_100km2 <- a$lcl/(a$area_km2/(100))
a$Ducl_100km2 <- a$ucl/(a$area_km2/(100))
a

ncol(a)
names(a)
b <- a[1,2:11] #E.N
b


b$estimate    <- round(b$estimate, digits = 2)
b$SE.estimate <- round(b$SE.estimate, digits = 2)
b$lcl <- round(b$lcl, digits = 2)
b$ucl <- round(b$ucl , digits = 2)
b$Dest_100km2 <- round(b$Dest_100km2, digits = 2)
b$Dse_100km2 <- round(b$Dse_100km2, digits = 2)
b$Dlcl_100km2 <- round(b$Dlcl_100km2, digits = 2)
b$Ducl_100km2 <- round(b$Ducl_100km2, digits=2)

colnames(b) <- c("Abundance", "SE", "Lower 95% CI", "Upper 95% CI", "n", "Area(km2)", "Density/100km2",
                 "SE", "Lower 95% CI", "Upper 95% CI")
b
library(kableExtra)
b %>%
  kbl(caption = "     Population and density estimates of leopard in Bhutan, 2022") %>%
  kable_classic(full_width =T , html_font = "Cambria")


save.image(paste("Leopard_Density_models_AEC_5km",Sys.Date(),".RData"))


#model averaging detection prob
#model AIC weights (aic_weights) that sum to 1,
# Extract g0 from each model
coef(density_fit2[1])

g0_estimates <- sapply(density_fit2, function(mod) coef(mod)["g0", "beta"])
g0_SEs  <- sapply(density_fit2, function(mod) coef(mod)["g0", "SE.beta"])
g0_lcls <- sapply(density_fit2, function(mod) coef(mod)["g0", "lcl"])
g0_ucls <- sapply(density_fit2, function(mod) coef(mod)["g0", "ucl"])

# Weighted average of g0
g0_avg <- sum(g0_estimates * modelsel$AICcwt)
g0_SE_avg <- sum(g0_SEs * modelsel$AICcwt)
g0_lcl_avg <- sum(g0_lcls * modelsel$AICcwt)
g0_ucl_avg <- sum(g0_ucls * modelsel$AICcwt)


g0_avg_prob <- plogis(g0_avg)  # plogis = inverse logit
g0_SE_avg_prob<- plogis(g0_SE_avg)
g0_lcl_avg_prob <- plogis(g0_lcl_avg)
g0_ucl_avg_prob<- plogis(g0_ucl_avg)

g0_mod_avg <- cbind(g0_avg_prob,g0_SE_avg_prob,g0_lcl_avg_prob,g0_ucl_avg_prob)
g0_mod_avg
#[1,] 0.003027946      0.5189484     0.002610858      0.00351143

#model average the predicted detection probability 
#you can just predict detection at the mean value of all covariates.

det.modelsel #AIC 
detection_fit #model list
det.AICc.wt <- det.modelsel[8]
det.top <- detection_fit$fit1


newdata.det <- data.frame(housden= mean(covariate_data$housden, na.rm =T),
                               tigerden = mean(covariate_data$tigerden, na.rm =T),
                               treecover = mean(covariate_data$treecover, na.rm =T),
                               elevation = mean(covariate_data$elevation, na.rm =T),
                               strmden = mean(covariate_data$strmden, na.rm = T),
                               prey_count = mean(covariate_data$prey_count, na.rm=T))

newdata.det.pred <- predict(det.top, newdata = newdata.det, appendData=T)
head(newdata.det.pred)

results.g0  <- data.frame(      estimate = numeric(length(housden)),
                                SE.estimate = numeric(length(housden)),
                                lcl = numeric(length(housden)), 
                                ucl = numeric(length(housden))
)
head(results.g0)

# Fill in the results.treecover 

  results.g0$estimate <- newdata.det.pred$estimate[2]
  results.g0$SE.estimate <- newdata.det.pred$SE.estimate[2]
  results.g0$lcl<- newdata.det.pred$lcl[2]
  results.g0$ucl <- newdata.det.pred$ucl[2]


head(results.g0)[1,]


 #  __  __           _      _                       _____              _ _      _   _             
 # |  \/  |         | |    | |     /\              |  __ \            | (_)    | | (_)            
 # | \  / | ___   __| | ___| |    /  \__   ____ _  | |__) | __ ___  __| |_  ___| |_ _  ___  _ __  
 # | |\/| |/ _ \ / _` |/ _ \ |   / /\ \ \ / / _` | |  ___/ '__/ _ \/ _` | |/ __| __| |/ _ \| '_ \ 
 # | |  | | (_) | (_| |  __/ |  / ____ \ V / (_| | | |   | | |  __/ (_| | | (__| |_| | (_) | | | |
 # |_|  |_|\___/ \__,_|\___|_| /_/    \_\_/ \__, | |_|   |_|  \___|\__,_|_|\___|\__|_|\___/|_| |_|
 #                                           __/ |                                                
 #                                          |___/                                                 

AICc.wt <- modelsel[8]

covariate_data <- attr(hab_mask, "covariates")
head(covariate_data)

range(covariate_data$housden) #-0.3517039 18.5438257
hist(covariate_data$housden)

range(covariate_data$tigerden) #-0.6943142  6.6498452
mean(covariate_data$tigerden) #0
hist(covariate_data$tigerden)
median(covariate_data$tigerden) #-0.29
  
#tiger
# To answer the below - you need to input the actual numbers
# (as numeric) here -make sure it is scaled or unscaled as it
# was input into the model (i.e.- the numbers selected are on 
# the same scale as the variable that went 
# into the model): 
#   
  tigerden.s <- factor(
  rep(c("tigerden.low", "tigerden.med", "tigerden.high"), each = 100),
  levels = levels(covariate_data$tigerden_cat)
)
# For example:
#   Tigerden.s <- rep(c(0,1,5),100)

# covariate_data$tigerden ranges from -0.69 to 6.65
covariate_data$tigerden_cat <- cut(
  covariate_data$tigerden,
  breaks = c(-Inf, 0, 3.3, Inf),  # adjust 3.3 to a value that makes sense as the "medium" cutoff
  labels = c("tigerden.low", "tigerden.med", "tigerden.high"),
  right = TRUE
)


library(classInt)

# Get the natural breaks
breaks <- classIntervals(covariate_data$tigerden, n = 3, style = "jenks")

# Show the breakpoints
breaks$brks #-0.6943142  0.2913168  2.6117459  6.6498452

# Ensure factor levels match those used in the model
tigerden.s <- c(rep(0.29,100),rep(2.61,100),rep(5.55,100))

housden.s<-seq(min(covariate_data$housden),max(covariate_data$housden),
             (max(covariate_data$housden)-min(covariate_data$housden))/99.99)


# Create a new data frame for prediction, holding other covariates at their mean
newdata.hous.tig <- data.frame(housden= housden.s,
                               tigerden = tigerden.s,
                               treecover = mean(covariate_data$treecover, na.rm =T ),
                               elevation = mean(covariate_data$elevation, na.rm =T ),
                               strmden = mean(covariate_data$strmden, na.rm = T),
                              prey_count = mean(covariate_data$prey_count, na.rm= T))

head(newdata.hous.tig)
nrow(newdata.hous.tig)

# Predict using the top model
hous.tig_pred1 <- predict(top_mod, newdata = newdata.hous.tig, appendData=T)
hous.tig_pred2 <- predict(top_mod2, newdata = newdata.hous.tig, appendData=T)
hous.tig_pred3 <- predict(top_mod3, newdata = newdata.hous.tig, appendData=T)

# Create an empty data frame to store results
results.hous.tig <- data.frame(housden= rep(housden.s,3), 
                               tigerden=tigerden.s,
                                estimate = numeric(length(housden.s)*3),
                                SE.estimate = numeric(length(housden.s)*3),
                                lcl = numeric(length(housden.s)*3), 
                                ucl = numeric(length(housden.s)*3)
)
head(results.hous.tig) 
nrow(results.hous.tig)
 ##if there is err0r in plotting, just rerun- Yes you are using just the hous.tig_pred1 model
# Fill in the results.treecover 
for (i in 1:length(hous.tig_pred1)) {
  results.hous.tig $estimate[i] <- hous.tig_pred1[[i]]$estimate[1]
  results.hous.tig $SE.estimate[i] <- hous.tig_pred1[[i]]$SE.estimate[1]
  results.hous.tig $lcl[i] <- hous.tig_pred1[[i]]$lcl[1]
  results.hous.tig $ucl[i] <- hous.tig_pred1[[i]]$ucl[1]
}

head(results.hous.tig)
# Error in `$<-.data.frame`(`*tmp*`, "estimate", value = c(8.56714382884571e-05,  : 
#                                                            replacement has 201 rows, data has 200
#work from here
#Housing:tiger were plot is from top

# Extract predictions and standard errors
# hous.tig.est <- results.hous.tig$estimate
# hous.tig.est2 <- results.hous.tig$estimate
# hous.tig.est3 <- results.hous.tig$estimate
# 
# hous.tig.se <- results.hous.tig$SE.estimate
# hous.tig.se2 <- results.hous.tig$SE.estimate
# hous.tig.se3 <- results.hous.tig$SE.estimate
# 
# # Compute the model-averaged prediction
# Mod.Avg.Pred.treecover <- hous.tig.est*AICc.wt[1,]+
#   hous.tig.est2*AICc.wt[2,]+
#   hous.tig.est3*AICc.wt[3,]
# 
# 
# # Step 4: Compute the unconditional variance and SE
# tree_var_avg <- AICc.wt[1,] * (hous.tig.se^2 + (hous.tig.est- Mod.Avg.Pred.treecover)^2) +
#   AICc.wt[2,] * (hous.tig.se2^2 + (hous.tig.est2- Mod.Avg.Pred.treecover)^2) +
#   AICc.wt[3,]* (hous.tig.se3^2 + (hous.tig.est3- Mod.Avg.Pred.treecover)^2)
# 
# tree_se_avg <- sqrt(tree_var_avg)
# 
# #CI
# # Compute 95% Confidence Intervals
# tree_LCL <- Mod.Avg.Pred.treecover - (1.96 * tree_se_avg)
# tree_UCL <- Mod.Avg.Pred.treecover + (1.96 * tree_se_avg)
# 
# # Add CIs to results dataframe
# results.hous.tig$model_avg_estimate <- Mod.Avg.Pred.treecover
# results.hous.tig$model_avg_se <- tree_se_avg
# results.hous.tig$model_avg_lcl <- tree_LCL
# results.hous.tig$model_avg_ucl <- tree_UCL
# 
# # View final results
# head(results.hous.tig)
# 
# 
# # model-averaged prediction and its SE
# Mod.Avg.Pred.hous.tig
# tree_se_avg
# Mod.Avg.Pred.hous.tig.all <- data.frame(cbind(Mod.Avg.Pred.hous.tig,tree_se_avg))
# Mod.Avg.Pred.hous.tig.all
# 
# results.hous.tig.all <- data.frame(treecover = treecover_seq, 
#                                     estimate = Mod.Avg.Pred.treecover,
#                                     SE.estimate = tree_se_avg,
#                                     lcl = tree_LCL, 
#                                     ucl = tree_UCL
# )


results.hous.tig$tigerden.cat<-ifelse(results.hous.tig$tigerden==min(results.hous.tig$tigerden),"Low Tiger Density",
                                      ifelse(results.hous.tig$tigerden==max(results.hous.tig$tigerden),"High Tiger Density",
                                                                            "Mod Tiger Density"))

# Filter out "Mod Tiger Density"
head(results.hous.tig)
str(results.hous.tig)

results.hous.tig <- results.hous.tig[results.hous.tig$tigerden.cat != "Mod Tiger Density", ]

# Define new colors
line_colors <- c("High Tiger Density" = "#d95f02",   # burnt orange
                 "Low Tiger Density" = "brown")    # goldenrod

fill_colors <- c("High Tiger Density" = "#f4b183",    # soft orange fill
                 "Low Tiger Density" = "brown")     # soft yellow fill

ggplot(results.hous.tig, aes(x = housden*sd_housden+housden_mean, y = estimate*1000))+
  geom_line(aes(color = tigerden.cat), linewidth = 1) +
  geom_ribbon(aes(min=lcl*1000, max=ucl*1000,fill = tigerden.cat), alpha = 0.4)+
  #geom_ribbon(aes(min=(estimate-SE.estimate)*1000, max=(estimate+SE.estimate)*1000,fill = tigerden.cat), alpha = 0.4)+
  geom_line(aes(y = lcl*1000, color = tigerden.cat), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = ucl*1000, color = tigerden.cat), linetype = "dashed", linewidth = 0.6) +
  labs(
    x = expression("Housing Density (km"^2*")"),
    y = expression("Leopard Density (per 100 km"^2*")"))+
  #ylim(0, 10)+ #tho range is 0.07188756 0.99998481
  coord_sf(ylim=c(0,1.5),)+
  scale_color_manual(values = line_colors) +
  scale_fill_manual(values = fill_colors) +
  theme_classic()+
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 14),  # Increase axis number size
    axis.title = element_text(size = 15),  # Increase axis label size
    axis.line = element_line(color = "grey50"),  # Default axis color
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.8),  # move legend inside plot
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA), # semi-transparent background
    legend.key.size = unit(0.8, "cm")
)

    

#prey:tiger
prey.s<-seq(min(covariate_data$prey_count),max(covariate_data$prey_count),
               (max(covariate_data$prey_count)-min(covariate_data$prey_count))/99.99)


# Create a new data frame for prediction, holding other covariates at their mean
newdata.prey.s <- data.frame(prey_count= prey.s,
                               tigerden = tigerden.s,
                               treecover = mean(covariate_data$treecover, na.rm =T ),
                               elevation = mean(covariate_data$elevation, na.rm =T ),
                               strmden = mean(covariate_data$strmden, na.rm = T),
                               housden = mean(covariate_data$housden, na.rm= T))

head(newdata.prey.s)
nrow(newdata.prey.s) #300

# Predict using the top model
prey.pred <- predict(top_mod, newdata = newdata.prey.s, appendData=T)

# Create an empty data frame to store results
results.prey <- data.frame(prey_count = rep(prey.s,3), 
                                tigerden=tigerden.s,
                                estimate = numeric(length(prey.s)*3),
                                SE.estimate = numeric(length(prey.s)*3),
                                lcl = numeric(length(prey.s)*3), 
                                ucl = numeric(length(prey.s)*3)
)
head(results.prey) 
nrow(results.prey) #300

# Fill in the results.treecover 
for (i in 1:length(prey.pred)) {
  results.prey $estimate[i] <- prey.pred[[i]]$estimate[1]
  results.prey $SE.estimate[i] <- prey.pred[[i]]$SE.estimate[1]
  results.prey $lcl[i] <- prey.pred[[i]]$lcl[1]
  results.prey $ucl[i] <- prey.pred[[i]]$ucl[1]
}

head(results.prey)
nrow(results.prey)

results.prey$tigerden.cat1<-ifelse(results.prey$tigerden==min(results.prey$tigerden),"Low Tiger Density",
                                      ifelse(results.prey$tigerden==max(results.prey$tigerden),"High Tiger Density",
                                             "Mod Tiger Density"))

# Filter out "Mod Tiger Density"
results.prey <- results.prey[results.prey$tigerden.cat != "Mod Tiger Density", ]

head(results.prey)
nrow(results.prey)

# Define new colors
line_colors <- c("High Tiger Density" = "#634081",   # deep purple
                 "Low Tiger Density" = "#a690b4")    # soft lavender


line_colors <- c("High Tiger Density" = "#a678de",   # Lavender
                 "Low Tiger Density"  = "#e573b3")   # Rosy pink

fill_colors <- c("High Tiger Density" = "#d1b3ff",   # Soft lavender fill
                 "Low Tiger Density"  = "#f8c1d9")

#used this
fill_colors <- c("High Tiger Density" = "#d1b3ff",   # soft lavender
                 "Low Tiger Density"  = "#e573b3")   # soft pink

# Darker shades for the lines:
line_colors <- c("High Tiger Density" = "#6a00a8",   # dark purple
                 "Low Tiger Density"  = "#b30059")   # dark magenta/pink

ggplot(results.prey, aes(x = prey_count*sd_prey_count+mean_prey_count, y = estimate*1000))+
  geom_line(aes(color = tigerden.cat1), linewidth = 1) +
  geom_ribbon(aes(min=lcl*1000, max=ucl*1000,fill = tigerden.cat1), alpha = 0.4)+
  #geom_ribbon(aes(min=(estimate-SE.estimate)*1000, max=(estimate+SE.estimate)*1000,fill = tigerden.cat), alpha = 0.4)+
  geom_line(aes(y = lcl*1000, color = tigerden.cat1), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = ucl*1000, color = tigerden.cat1), linetype = "dashed", linewidth = 0.6) +
  labs(
    x = expression("Prey Count"),
    y = expression("Leopard Density (per 100 km"^2*")"))+
  #ylim(0, 10)+ #tho range is 0.07188756 0.99998481
  coord_sf(ylim=c(0,1.5),)+
  scale_color_manual(values = line_colors) +
  scale_fill_manual(values = fill_colors) +
  theme_classic()+
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 14),  # Increase axis number size
    axis.title = element_text(size = 15),  # Increase axis label size
    axis.line = element_line(color = "grey50"),  # Default axis color
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_blank(),
    # legend.text = element_text(size = 12),
    legend.position = c(0.25, 0.8),  # move legend inside plot
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA), # semi-transparent background
    legend.key.size = unit(0.8, "cm")
    )


    
#treecover
treecover_seq <- seq(min(covariate_data$treecover),max(covariate_data$treecover),length.out = 100)


# Create a new data frame for prediction, holding other covariates at their mean
newdata.treecover <- data.frame(treecover = treecover_seq,
                      elevation = mean(covariate_data$elevation, na.rm =T ),
                      strmden = mean(covariate_data$strmden, na.rm = T),
                      housden= mean(covariate_data$housden, na.rm=T),
                      tigerden= mean(covariate_data$tigerden, na.rm=T),
                      prey_count = mean(covariate_data$prey_count, na.rm= T))

head(newdata.treecover)

# Predict using the top model
tree_pred <- predict(top_mod, newdata = newdata.treecover, appendData=T)
tree_pred2 <- predict(top_mod2, newdata = newdata.treecover, appendData=T)
tree_pred3 <- predict(top_mod3, newdata = newdata.treecover, appendData=T)

# Create an empty data frame to store results
results.treecover <- data.frame(treecover = treecover_seq, 
                      estimate = numeric(length(treecover_seq)),
                      SE.estimate = numeric(length(treecover_seq)),
                      lcl = numeric(length(treecover_seq)), 
                      ucl = numeric(length(treecover_seq))
)
results.treecover 

# Fill in the results.treecover 
for (i in 1:length(tree_pred)) {
  results.treecover $estimate[i] <- tree_pred[[i]]$estimate[1]
  results.treecover $SE.estimate[i] <- tree_pred[[i]]$SE.estimate[1]
  results.treecover $lcl[i] <- tree_pred[[i]]$lcl[1]
  results.treecover $ucl[i] <- tree_pred[[i]]$ucl[1]
}

head(results.treecover)

# Extract predictions and standard errors
tree.pred.est <- results.treecover$estimate
tree.pred.est2 <- results.treecover$estimate
tree.pred.est3 <- results.treecover$estimate

tree.pred.se <- results.treecover$SE.estimate
tree.pred.se2 <- results.treecover$SE.estimate
tree.pred.se3 <- results.treecover$SE.estimate

# Compute the model-averaged prediction
Mod.Avg.Pred.treecover <- tree.pred.est*AICc.wt[1,]+
                      tree.pred.est2*AICc.wt[2,]+
                    tree.pred.est3*AICc.wt[3,]


# Step 4: Compute the unconditional variance and SE
tree_var_avg <- AICc.wt[1,] * (tree.pred.se^2 + (tree.pred.est- Mod.Avg.Pred.treecover)^2) +
  AICc.wt[2,] * (tree.pred.se2^2 + (tree.pred.est2- Mod.Avg.Pred.treecover)^2) +
  AICc.wt[3,]* (tree.pred.se3^2 + (tree.pred.est3- Mod.Avg.Pred.treecover)^2)

tree_se_avg <- sqrt(tree_var_avg)

#CI
# Compute 95% Confidence Intervals
tree_LCL <- Mod.Avg.Pred.treecover - (1.96 * tree_se_avg)
tree_UCL <- Mod.Avg.Pred.treecover + (1.96 * tree_se_avg)

# Add CIs to results dataframe
results.treecover$model_avg_estimate <- Mod.Avg.Pred.treecover
results.treecover$model_avg_se <- tree_se_avg
results.treecover$model_avg_lcl <- tree_LCL
results.treecover$model_avg_ucl <- tree_UCL

# View final results
head(results.treecover)


# model-averaged prediction and its SE
Mod.Avg.Pred.treecover
tree_se_avg
Mod.Avg.Pred.treecover.all <- data.frame(cbind(Mod.Avg.Pred.treecover,tree_se_avg))
Mod.Avg.Pred.treecover.all

results.treecover.all <- data.frame(treecover = treecover_seq, 
                                estimate = Mod.Avg.Pred.treecover,
                                SE.estimate = tree_se_avg,
                                lcl = tree_LCL, 
                                ucl = tree_UCL
)
head(results.treecover.all)


#   _____  _       _
 # |  __ \| |     | |  
 # | |__) | | ___ | |_ 
 # |  ___/| |/ _ \| __|
 # | |    | | (_) | |_ 
 # |_|    |_|\___/ \__|


elevation_mean <- cellStats(elevation, stat = 'mean', na.rm = TRUE)
elevation_mean #3042.964
sd_elevation <- cellStats(elevation, stat = 'sd', na.rm = TRUE)
sd_elevation #1596.076
treecover_mean <- cellStats(treecover, stat = 'mean', na.rm = TRUE) #86.38
treecover_mean
sd_treecover <-  cellStats(treecover, stat = 'sd', na.rm = TRUE) #17.76
housden_mean <- cellStats(housden, stat = 'mean', na.rm = TRUE) #9.791
housden_mean
sd_housden <- cellStats(housden, stat = 'sd', na.rm = TRUE) #27.84
sd_housden
strmden_mean <- cellStats(strmden, stat = 'mean', na.rm = TRUE) #2.68
strmden_mean
sd_strmden <- cellStats(strmden, stat = 'sd', na.rm = TRUE) #1.81
sd_strmden
tigerden_mean <- cellStats(tigerden, stat = 'mean', na.rm = TRUE) #0.26
tigerden_mean
sd_tigerden <- cellStats(tigerden, stat = 'sd', na.rm = TRUE) #0.37
sd_tigerden
mshp<-st_read("prey_2025_feb.shp")
mshp$prey_count <- as.numeric(mshp$prey_count)
class(mshp$prey_count)
prey_count_mean <- mean(mshp$prey_count, na.rm = TRUE) #18
sd_prey_count<- sd(mshp$prey_count, na.rm = TRUE) #26.63


#ggplot
#treecover
tree <- (results.treecover.all$treecover *sd_treecover +treecover_mean)
range(tree*100)
range(results.treecover.all$estimate*1000)
range(results.treecover.all$lcl*1000)

#tree in percentage
ggplot(results.treecover.all, aes(x = tree*100, y = estimate*1000))+
  geom_line(color = '#4c633c', linewidth = 1) +
  geom_ribbon(aes(ymin=lcl*1000, ymax=ucl*1000),fill = '#8ba974', alpha = 0.4)+
  geom_line(aes(y = lcl*1000), color = "#8ba974", linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = ucl*1000), color = "#8ba974", linetype = "dashed", linewidth = 0.6) +
    labs(
    x = expression("Tree cover (%)"),
    y = expression("Leopard Density (per 100 km"^2*")"))+
  #ylim(0, 1.5) +
  coord_sf(xlim=c(min(tree*100),100),ylim=c(0,1.5))+
  theme_classic()+
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 14),  # Increase axis number size
    axis.title = element_text(size = 15),  # Increase axis label size
    axis.line = element_line(color = "grey50"),  # Default axis color
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1)  # Match default axis color
  )                   

#Tree from TOP MODEL Only 

tree <- (results.treecover$treecover *sd_treecover +treecover_mean)
#tree in percentage
range(tree)
ggplot(results.treecover, aes(x = tree*100, y = estimate*1000))+
  geom_line(color = '#4c633c', linewidth = 1) +
  geom_ribbon(aes(min=lcl*1000, max=ucl*1000),fill = '#8ba974', alpha = 0.4)+
  labs(
    x = expression("Tree Cover (%)"),
    y = expression("Leopard Density (per 100 km"^2*")"))+
  ylim(0, 1.5)+  xlim(20, 100)+ #tho range is 0.07188756 0.99998481
  theme_classic()+
  theme(
    axis.text = element_text(size = 21),  # Increase axis number size
    axis.title = element_text(size = 21),  # Increase axis label size
    axis.line = element_line(color = "grey50"),  # Default axis color
    panel.border = element_rect(color = "grey50", fill = NA, size = 1)  # Match default axis color
  )


