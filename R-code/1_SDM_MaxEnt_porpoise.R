

# Date: 19 Dec 2025
# Author: Floris M. van Beest (flbe[at]ecos.au.dk)
# Why: The below code was generated during the OBAMA-NEXT project to fulfil IP 18. 
# Reference: The code follows the procedure described in the article: van Beest FM, Carstensen J, Dietz R, Nabe-Nielsen J, Sveegaard S, Teilmann J (2025) Shifts in habitat suitability for harbour porpoises leads to reduced importance of marine protected areas. Biological Conservation 302:111009

# Brief overview: The code below reads the harbour porpoise tracking data (based on ARGOS tags), cleans the Argos data, and performs a MaxEnt analyses to generate habitat suitability layers for one season and period (out of 6). To generate habitat suitability layers for all season and periods described in the article, adjust the code to the season and period of interest and rerun the procedure.  


#----------------------------------------------------#
1_load_format_ARGOS_data <- function()
#----------------------------------------------------#

library("readxl")
argos_porp <- read_excel("~/SDM_porpoise_IP18/data/AllposHP_1997_2022.xlsx")
str(argos_porp)
argos_porp <- argos_porp[,c("animal", "pttid","date","time","longitud","latitude","lc94","sex","deployed","deploc", "month", "day", "year", "hour", "minute", "second")]

# generate datetime column
argos_porp$date <- paste(argos_porp$year, argos_porp$month, argos_porp$day, sep="-")
argos_porp$time <- paste(argos_porp$hour, argos_porp$minute, argos_porp$second, sep=":")
argos_porp$datetime <- paste(argos_porp$date, argos_porp$time, sep=" ")
argos_porp $datetime <- as.POSIXct(argos_porp$datetime,format="%Y-%m-%d %H:%M:%S", tz="UTC")
table(argos_porp$pttid)

# location quality index
table(argos_porp$lc94)
argos_porp$lc94 <- as.factor(argos_porp$lc94)

## get year
argos_porp$Year <- as.numeric(format(argos_porp$datetime,'%Y'))
table(argos_porp$pttid,argos_porp$Year)
table(argos_porp$Year)

# get month
argos_porp$Month <- as.numeric(format(argos_porp$datetime,'%m'))
table(argos_porp$Month)


##---------------------------------------------------------------##
# remove location 24 hour after tagging to reduce behavioural bias
##---------------------------------------------------------------##
argos_porp$hours.since.deployed <- difftime(argos_porp $datetime, argos_porp $deployed, units="hours")
argos_porp <- subset(argos_porp, hours.since.deployed>24)


##-----------------------------------------------##
# remove_locations_on_land 
##-----------------------------------------------##
library(terra)
bathy <- rast("~/SDM_porpoise_IP18/rasters/bathymetry.grd")
bathy <- project(bathy, "+proj=longlat +ellps=WGS84 +datum=WGS84")

Bathymetry<-extract(bathy, argos_porp[,c('longitud','latitude')],ID=F)
argos_porp <- cbind(argos_porp, Bathymetry)
argos_porp <- subset(argos_porp, !is.na(Bathymetry))

plot(bathy)
points(argos_porp $longitud, argos_porp $latitude)

names(argos_porp)
argos_porp <- argos_porp[,c("animal", "pttid","datetime","longitud","latitude","lc94","sex","deployed","deploc", "Year", "Month")]
names(argos_porp)

##-----------------------------------------------##
# create Seasons and Periods
##-----------------------------------------------##

argos_porp$Season <- ifelse(argos_porp$Month %in%c(5,6,7,8,9,10), "Summer", "Winter")

argos_porp1 <- subset(argos_porp, Year<2005)
argos_porp1$Period <- "period_1997_2004"

argos_porp2 <- subset(argos_porp, Year>2004&Year<2013)
argos_porp2$Period <- "period_2005_2012"

argos_porp3 <- subset(argos_porp, Year>2012)
argos_porp3$Period <- "period_2013_2022"

argos_porp <- rbind(argos_porp1, argos_porp2, argos_porp3)
 
 
##-----------------------------------------------##
# change projection as used in OBAMA-NEXT
##-----------------------------------------------##

library(sp)
latslongs <- SpatialPointsDataFrame(coords= argos_porp[,c(4,5)],data= argos_porp,proj4string =CRS("+proj=longlat +datum=WGS84")) 
latslongs <- spTransform(latslongs, CRS("+init=epsg:3035"))

argos_porp <- as.data.frame(latslongs)
names(argos_porp)[14] <- "X"
names(argos_porp)[15] <- "Y"


##---------------------------------------------------------------##
# make a new column date to be used for sampling with replacement
##---------------------------------------------------------------##

library(lubridate)
argos_porp$date <- as.Date(ymd_hms(argos_porp$datetime))

#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################

#-------------------------------------------#
2_setup_for_MaxEnt_summer_97_04 <- function()
#-------------------------------------------#


## subset the period and season to model (1 out of 6)
dat1 <- subset(argos_porp, Season=="Summer" & Period=="period_1997_2004")


##-----------------------------------------------##
# subsample the Argos data
##-----------------------------------------------##

dat1 $ddat_id <- paste(dat1 $date, dat1 $pttid, sep="_")

## keep only 1 location per day with replacement
library(dplyr)
dat2 <- dat1 %>%
		    group_by(ddat_id) %>%
		    slice_sample(n=1, replace=T) 

## keep only IDs with at least 10 tracking days
dat3 <- dat2 %>%
		    group_by(pttid) %>%
		    filter(n() >=10) 
summary(dat3)
head(dat3)

# dataframe of locations for model pruning
dat0 <- dat3[,c(14:15)]


##----------------------------------------------------------------------##
# load environmental layers for season period to be used in MaxEnt model
##----------------------------------------------------------------------##
library(raster)
predictors <- stack("~/SDM_porpoise_IP18/rasters/stack_summer_97_04.grd")

## reclass sediment type
sed <- raster::subset(predictors, grep('Sediment', names(predictors), value = T))
s <- setValues(raster(sed), sed[])
s
cls <- data.frame(ID=0:3, Sediment=c("Bedrock", "Hard bottom complex", "Mud", "Sand"))
levels(s) <- cls
is.factor(s)
s
sediments <-s
names(sediments) <- "Sediments"

# create classification matrix
reclass_df <- c(0, 0.5, 0,
              0.5, 1.5, 1,
             1.5, 2.5, 2,
             2.5, Inf, 3)
reclass_m <- matrix(reclass_df,
                ncol = 3,
                byrow = TRUE)
# reclassify the raster using the reclass object
sediments <- reclassify(sediments,
                     reclass_m)

cls <- data.frame(ID=0:3, Sediment=c("Bedrock", "Hard bottom complex", "Mud", "Sand"))
levels(sediments) <- cls
is.factor(sediments)
sediments
sediments <-sediments
names(sediments) <- "Sediments"

##-----------------------------------------------##
# check for collinearity among environmental layers
##-----------------------------------------------##

## remove sediment for collinearity test
predictors <- dropLayer(predictors, c(13))

#predictors <- scale(predictors)

predictors2 <- terra::rast(predictors[[1:14]])
library(usdm)
v1 <-  vif(predictors2)
v1
mean(v1$VIF, na.rm=T)

v2 <- vifstep(predictors2, th=3)
v2

library(raster)    
jnk=layerStats(predictors, 'pearson', na.rm=T)
jnk

## remove correlated variables
names(predictors)

library(ENMTools)

raster.cor.plot(predictors)
predictors <- dropLayer(predictors, c(1,2,5,7,9,12))

## add sediment back in
predictors <- stack(predictors, sediments)
predictors <- terra::rast(predictors)


#-----------------------------------------------------#
2_1_generate_background_points_using_sampling_bias <- function()
#-----------------------------------------------------#

## generate sampling bias file
dat1 <- subset(argos_porp, season=="summer" & period=="period_1997_2004")

pres.locs <- coordinates(dat1[,c(17,18)])
 
library(MASS) # for 2D kernel density function
h_max  <-cbind(bandwidth.nrd(pres.locs[,1]), bandwidth.nrd(pres.locs[,2]))
h_max

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(680,260),h= h_max , lims = c(range(4280000, 4540000), range(3420000, 4100000)))
dens.ras <- raster(dens)

# crop to raster of the study area extent
refRaster <- raster("~/SDM_porpoise_IP18/rasters/stack_summer_97_04.grd")
dens.ras <- crop(dens.ras, extent(refRaster), snap="near")
# mask
dens.ras <- resample(dens.ras, refRaster[[1]],method="bilinear")
dens.ras <- mask(dens.ras, refRaster[[1]])

plot(dens.ras)

# replacing NA's by zero
dens.ras[is.na(dens.ras[])] <- 0 
set.seed(25)

library(SDMtune)
bg <- xyFromCell(dens.ras, sample(ncell(dens.ras), 10000, replace=T, prob=values(dens.ras))) 
bg <- prepareSWD(species = "Porpoise", a = bg, env = predictors)


#-----------------------------------------------------#
2_2_MaxEnt_model_pruning <- function()
#-----------------------------------------------------#

library(ENMeval)
library(ecospat)

bg <- as.matrix(bg@coords)
colnames(bg) <- c("X", "Y")
ps <- list(orientation = "lon_lat")
tune.args <- list(fc = c("LQ"), rm = seq(0.5,5,0.5))
enmeval_results <- ENMevaluate(dat0, predictors, partitions ="block", bg=bg,
                                 algorithm='maxent.jar',overlap=F, partition.settings = ps, categoricals = "Sediments", tune.args = tune.args, doClamp=TRUE )

enmeval_results@results

# find optimal model settings (tobe used in loop below)
library(tibble)
library(dplyr)
res <- eval.results(enmeval_results)
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
mod.seq <- eval.models(enmeval_results)[[opt.seq$tune.args]]

# use terra instead of raster
predictors <- terra::rast(predictors)


#-----------------------------------------------------#
2_3_loop_best_model_response_VI_predict <- function()
#-----------------------------------------------------#


rast.now<- rast()
l.model <- list()
null_all <- data.frame()
threshold_all <- data.frame()
AUC <- data.frame()
OR10 <- data.frame()
VI_dat <- data.frame()
response_dat <- data.frame()
lat.change.kappa <- data.frame()
lat.change.spec_sens <- data.frame()
lat.change.no_omission <- data.frame()
size.change.kappa <- data.frame()
size.change.spec_sens <- data.frame()
size.change.no_omission <- data.frame()
nni.change.kappa <- data.frame()
nni.change.spec_sens <- data.frame()
nni.change.no_omission <- data.frame()

#------------------------------------#
# run the loop 100 times to quantify uncertainty
#------------------------------------#
run = 100
for(i in 1:run){

	# generate dataset for each run
	## keep only 1 location per day with replacement
	library(dplyr)
	dat2 <- dat1 %>%
		    group_by(ddat_id) %>%
		    slice_sample(n=1, replace=T) 

	## keep only IDs with at least 10 tracking days
	dat3 <- dat2 %>%
		    group_by(pttid) %>%
		    filter(n() >=10) 

	dat4 <- dat3[,c(14:15)]
	
	# generate training testing dataset for each run
	smp_size <- floor(0.8 * nrow(dat4))
	train_ind <- sample(seq_len(nrow(dat4)), size = smp_size)
	train <- dat4[train_ind, ]
	test <- dat4[-train_ind, ]
	
	#------------------------------------#
	# run the best ENMeval model
	#------------------------------------#
	#ps <- list(aggregation.factor =2)
	ps <- list(kfolds=5,orientation = "lon_lat")
	tune.args <- list(fc = c("LQ"), rm = 3.5)

	best.mod <- ENMevaluate(train[,c(1:2)], predictors, partitions ="block", bg =bg, occs.testing=test[,c(1:2)], algorithm='maxent.jar',overlap=F, partition.settings = ps, tune.args = tune.args, doClamp=FALSE, categoricals = "Sediments")
	l.model <- list(l.model, best.mod)
	

	#------------------------------------#
	# null model with optmal settings
	#------------------------------------#
	 mod.null <- ENMnulls(best.mod, mod.settings = list(fc = "LQ", rm = 3.5), no.iter = 1)
	 null <- data.frame()
	 null <- null.emp.results(mod.null)
	 null$Period <- "1997-2004"
	 null$Season <- "Summer"
	 null $model <- i
	 null <- null[c(1,3),c(1,4,6,8,9,10,11)]
     null_all <- rbind(null_all,null)
	
	#------------------------------------#
	# Model Evaluation
	#------------------------------------#
	mod.evT1 <- evaluate(best.mod@models[[1]], p = test[,c(1:2)], a = bg, x = predictors)
	
	#------------------------------------#
	# Capturing the AUC values etc
	#------------------------------------#
	auc.T1 <- as.data.frame(best.mod@results$auc.train)
	names(auc.T1)[1] <- "AUC"
	auc.T1$model <- i
	auc.T1 $Period <- "1997-2004"
	auc.T1 $Season <- "Summer"
	AUC <- rbind(AUC, auc.T1)

	#------------------------------------#
	# Capturing the OR10 values etc
	#------------------------------------#
	auc.T1 <- as.data.frame(best.mod@results$or.10p.avg)
	names(auc.T1)[1] <- "OR10"
	auc.T1$model <- i
	auc.T1 $Period <- "1997-2004"
	auc.T1 $Season <- "Summer"
	OR10 <- rbind(OR10, auc.T1)
	
	
	#------------------------------------#
	# extract variable importance
	#------------------------------------#
	
	VI <- as.data.frame(best.mod@variable.importance)
	names(VI)[1] <- "Variable"
	names(VI)[2] <- "Percent.contribution"
	names(VI)[3] <- "Permutation importance"
	VI $Period <- "1997-2004"
	VI $Season <- "Summer"
	VI $model <- i
	VI_dat <- rbind(VI_dat, VI)

	#------------------------------------#
	# predict habitat suitability and stack
	#------------------------------------#
	predict.now <- predict(predictors,best.mod@models[[1]], type = 'cloglog', na.rm=T)
	rast.now <- c(rast.now, predict.now)

}

## generate the mean and sd and cv habitat suitability maps
mean_map <- mean(rast.now)
sd_map <- stdev(rast.now)
cv_map <- cv(raster::stack(rast.now))


#-----------------------------------------------------#
2_4_save_output_files <- function()
#-----------------------------------------------------#


# save NULL data
 write.table(null_all, file = "~/SDM_porpoise_IP18/results/NULL data/dat97_04_summer.txt", row.names = FALSE, append = FALSE, col.names = TRUE, sep="\t")
 
# save AUC data
 write.table(AUC, file = "~/SDM_porpoise_IP18/results/AUC data/dat97_04_summer.txt", row.names = FALSE, append = FALSE, col.names = TRUE, sep="\t")

# save OR10 data
 write.table(OR10, file = "~/SDM_porpoise_IP18/results/OR10 data/dat97_04_summer.txt", row.names = FALSE, append = FALSE, col.names = TRUE, sep="\t")

# save VI data
 write.table(VI_dat, file = "~/SDM_porpoise_IP18/results/Variable importance data/dat97_04_summer.txt", row.names = FALSE, append = FALSE, col.names = TRUE, sep="\t")

 # save rasters
 mean_map = setNames(mean_map[[1]], "summer1997-2004")
writeRaster(mean_map, filename = "~/SDM_porpoise_IP18/results/SDM maps/mean_SDM_stack_97_04_summer.grd", overwrite=TRUE)

 sd_map = setNames(sd_map[[1]], "summer1997-2004")
writeRaster(sd_map, filename = "~/SDM_porpoise_IP18/results/SDM maps/sd_SDM_stack_97_04_summer.grd", overwrite=TRUE)

 cv_map = setNames(cv_map[[1]], "summer1997-2004")
writeRaster(cv_map, filename = "~/SDM_porpoise_IP18/results/SDM maps/cv_SDM_stack_97_04_summer.grd", overwrite=TRUE)


#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################

## repeat the process above for any other season period combination of interest
## data (Argos and environmental layers) are available for summer and winter for the "period_1997_2004" & "period_2005_2012" & "period_2013_2022"



