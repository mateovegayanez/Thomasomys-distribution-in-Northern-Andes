### Thomasomys distribution in northern Andes ###
#SDM
library(readr)
library(raster)
library(sf)
library(usdm)
library(sdm)
library(ecospat)
xy <- read.csv("Thomasomys_species.csv")
bios <- stack(list.files(path = "D:/Desktop/Thomasomys/SDM//Bios", pattern = "tif", full.names = T))
plot(bios$bio_4)
bios_df <- as.data.frame(bios)
v1 <- vifcor(bios_df, th = 0.7)
bios.red <- exclude(bios, v1)
v1
abs_obs <- 500
pseudo_abs_train <- ecospat.rand.pseudoabsences(nbabsences=abs_obs[1],
                                                glob=bios.all,colxyglob=1:2, colvar=1:2, presence=presences,
                                                colxypresence=1:2, mindist=0.04166667)
DataSpecies <- SpatialPointsDataFrame(coords = xy, data = DataSpeciesTrain.spp,
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
DataSpecies <- DataSpecies[,c(-1,-2)]

d1 <- sdmData(formula=Rep1~., train=DataSpecies, predictors=v1)
m1 <- sdm(Rep1~, data=d1, methods= c("glm", "maxent"), test.percent=30, n=10)
m1
bios_ps <- stack(list.files(path = "D:/Desktop/Thomasomys/SDM//Bios", pattern = "tif", full.names = T))
names(bios_ps) <- names(bios)
plot(bios_ps[[2]])
e2a <- ensemble(m1,newdata=bios_ps,filename='ps.tif', setting=list(method='weighted',stat='TSS'))
#Binary threshold of 0.9 occurrences
mpa.e1 <- ecospat.mpa(Model, Coordinates, perc = .9)
mpa.e1
bin.e1 <- ecospat.binary.model(Model, mpa.e1)
library(letsR)
PAM <- lets.presab(
  bin.e1,
  xmn = -82, xmx = -75,
  ymn = -7,  ymx = 3,
  resol = 0.1,
  show.matrix = TRUE,
  crs = CRS("+proj=longlat +datum=WGS84"),
  cover = 0
)
PAM
#PD
Phy <- read.nexus("Phy_Thomasomys")
pd <- pd(PAM, Phy, include.root = TRUE)
head(pd)
modelPD <- lm(PD ~ SR, data = pd)
summary(modelPD)
str(modelPD)
pd_res$resPD <- residuals(modelPD)
head(pd_res)
#Distribution centroids
library(raster)
library(sp)
library(sf)
species_df <- read.csv("Thomasomys_species.csv", stringsAsFactors = FALSE)

crs_utm <- CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

calc_centroid_shift <- function(species_name) {
  
  past_file    <- paste0("Past/", species_name, "_ps.tif")
  present_file <- paste0("Present/", species_name, "_p.tif")
  
  if (!file.exists(past_file) | !file.exists(present_file)) {
    return(NULL)
  }
  
  sdm_past    <- raster(past_file)
  sdm_present <- raster(present_file)
  
  sdm_past_utm <- projectRaster(sdm_past, crs = crs_utm, method = "ngb")
  sdm_pres_utm <- projectRaster(sdm_present, crs = crs_utm, method = "ngb")
  
  coords_past <- xyFromCell(
    sdm_past_utm,
    which(values(sdm_past_utm) == 1)
  )
  
  coords_pres <- xyFromCell(
    sdm_pres_utm,
    which(values(sdm_pres_utm) == 1)
  )
  
  if (nrow(coords_past) == 0 | nrow(coords_pres) == 0) {
    return(NULL)
  }
  
  centroid_past <- colMeans(coords_past)
  centroid_pres <- colMeans(coords_pres)
  
  dx <- centroid_pres[1] - centroid_past[1]
  dy <- centroid_pres[2] - centroid_past[2]
  
  distance_km <- sqrt(dx^2 + dy^2) / 1000
  angle_deg   <- (atan2(dy, dx) * 180 / pi) %% 360
  
  data.frame(
    species = species_name,
    x_past = centroid_past[1],
    y_past = centroid_past[2],
    x_present = centroid_pres[1],
    y_present = centroid_pres[2],
    distance_km = distance_km,
    angle_deg = angle_deg
  )
}

results <- bind_rows(
  lapply(species_df$species, calc_centroid_shift)
)
#PGLS
library(ape)
library(caper) 

#Explore variables
Clim <- read.csv("Values_Thomasomys.csv")

Phy <- read.nexus("Phy_Thomasomys")

comp_data <- comparative.data(phy = Phy, data = clima, vcv = TRUE)
pgls_t <- pgls(temp_p ~ temp_ps, data = comp_data, lambda = "ML")
summary(pgls_t)
pgls_p <- pgls(pre_p ~ pre_ps, data = comp_data, lambda = "ML")
summary(pgls_p)
#Residuals
res_t <- residuals(pgls_t)
res_t
res_p <- residuals(pgls_p)
res_p
#Scale variables 
comp_data$data$clim_pas  <- scale(comp_data$data$temp_ps) + scale(comp_data$data$pre_ps)
comp_data$data$clim_pres <- scale(comp_data$data$temp_p) + scale(comp_data$data$pre_p)
pgls_model <- pgls(clim_pres ~ clim_pas, data = comp_data, lambda = "ML")
summary(pgls_model)
#Residuals
res_scale <- residuals(pgls_model)
res_scale




