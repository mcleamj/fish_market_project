

if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(funrar)){install.packages("funrar"); library(funrar)}

###########################
## IMPORT SPECIES TRAITS ##
###########################

traits <- read.table("data/clean traits.txt")

####################################################
## IMPORT AND WRANGLE REEF AND LANDINGS DATA SETS ##
####################################################

reef_data <- read.csv("data/clean reef data.csv")
reef_meta <- reef_data %>%
  select(Site_SPC_year : geographic)
reef_fish <- reef_data %>%
  select(Acanthurus.blochii : Siganus.punctatus)
reef_log <- log10(reef_fish+1)
colnames(reef_log) <- gsub("\\.", " ", colnames(reef_log))
reef_log <- as.matrix(reef_log)
rownames(reef_log) <- paste("com",seq(1,nrow(reef_log),1))

market_data <- read.csv("data/clean market data.csv")
market_data$num_fishers_lines <- as.numeric(market_data$num_fishers_lines)
market_meta <- market_data %>%
  dplyr::select('Monitoring.Code':'num_fishers_lines')
market_fish <-market_data %>%
  dplyr::select('Acanthurus.blochii' : 'Siganus.punctatus')
market_log <- log10(market_fish+1)
colnames(market_log) <- gsub("\\.", " ", colnames(market_log))
market_log <- as.matrix(market_log)
rownames(market_log) <- paste("com",seq(1,nrow(market_log),1))

#######################
## BUILD TRAIT SPACE ##
#######################

trait_space <- prcomp(traits, scale. = TRUE)
graphics.off()
biplot(trait_space)
summary(trait_space)

PC1 <- trait_space$x[,1]
PC2 <- trait_space$x[,2]
PC3 <- trait_space$x[,3]
PC4 <- trait_space$x[,4]

axes <- data.frame(PC1, PC2, PC3, PC4)

###############################################################
## CALCULATE MEAN BIOMASS PER SPECIES PER GEOGRAPHIC SECTION
## FOR PLOTTING TRAIT SPACES
###############################################################

# FIRST BY SITE YEAR, THEN SITE, THEN GEOGAPRAPHIC SECTION

mean_biomass_reefs <- 
  data.frame( select(reef_meta, geographic, site, site_year), reef_log) %>%
  group_by(geographic, site, site_year) %>%
  summarise_all(mean) %>%
  select(-site_year) %>%
  group_by(geographic, site) %>%
  summarise_all(mean) %>%
  select(-site) %>%
  group_by(geographic) %>%
  summarise_all(mean)

# FIRST BY FISHER THEN GEOGRAPHIC SECTION

mean_biomass_landings <- 
  data.frame( select(market_meta, geographic, Fisher.Name), market_log) %>%
  group_by(geographic, Fisher.Name) %>%
  summarise_all(mean) %>%
  select(-Fisher.Name) %>%
  group_by(geographic) %>%
  summarise_all(mean) %>%
  filter(geographic %in% mean_biomass_reefs$geographic)



###############################################################
## CALCULATE MEAN BIOMASS PER SPECIES PER GEOGRAPHIC SECTION
## FOR PLOTTING TRAIT SPACES
###############################################################

# JUST CALCULATE MEAN BY GEOGRAPHIC

mean_biomass_reefs <- 
  data.frame( select(reef_meta, geographic), reef_log) %>%
  group_by(geographic) %>%
  summarise_all(mean)


# FIRST BY FISHER THEN GEOGRAPHIC SECTION

mean_biomass_landings <- 
  data.frame( select(market_meta, geographic), market_log) %>%
  group_by(geographic) %>%
  summarise_all(mean) %>%
  filter(geographic %in% mean_biomass_reefs$geographic)


##################################################
## PLOT TRAIT SPACE FOR EACH GEOGRAPHIC SECTION ##
##################################################

# REEFS 

#par(mfrow=c(2,10))

par(mfrow=c(1,5))

for(i in unique(mean_biomass_reefs$geographic)){
#for(i in "Northwest"){

  geo_biomass <- mean_biomass_reefs %>%
    filter(geographic == i) %>%
    select(-geographic) %>%
    as.matrix() %>%
    make_relative() %>%
    t() %>%
    as.data.frame() %>%
    rename(Biomass=V1) %>% 
    cbind( axes) %>%
    arrange(desc(Biomass))

  rownames(geo_biomass) <- gsub("\\."," ", rownames(geo_biomass))
  identical(rownames(geo_biomass), rownames(axes))
  
  plot(PC1, PC2, cex=0)
  
  geo_biomass <- geo_biomass %>%
    filter(Biomass != 0)
  
  points(geo_biomass$PC1, geo_biomass$PC2, pch=21, bg="grey", col="black", cex=geo_biomass$Biomass*20)
  title(paste(i,"- Reef Observations"))
  
  hpts <- chull(cbind(geo_biomass$PC1,geo_biomass$PC2))
  hpts <- c(hpts, hpts[1])
  polygon(cbind(geo_biomass$PC1,geo_biomass$PC2)[hpts, ], col = adjustcolor("grey",alpha.f=0.3), border="grey")
  

}


# LANDINGS

#par(mfrow=c(2,5))
for(i in unique(mean_biomass_reefs$geographic)){
#for(i in "Northwest"){
  
  geo_biomass <- mean_biomass_landings %>%
    filter(geographic == i) %>%
    select(-geographic) %>%
    as.matrix() %>%
    make_relative() %>%
    t() %>%
    as.data.frame() %>%
    rename(Biomass=V1) %>% 
    cbind( axes) %>%
    arrange(desc(Biomass))
    
  rownames(geo_biomass) <- gsub("\\."," ", rownames(geo_biomass))
  identical(rownames(geo_biomass), rownames(axes))
  
  plot(PC1, PC2, cex=0)
  
  geo_biomass <- geo_biomass %>%
    filter(Biomass != 0)
  
  points(geo_biomass$PC1, geo_biomass$PC2, pch=21, bg="grey", col="black", cex=geo_biomass$Biomass*20)
  title(paste(i,"- Fisheries Landings"))
  
  hpts <- chull(cbind(geo_biomass$PC1,geo_biomass$PC2))
  hpts <- c(hpts, hpts[1])
  polygon(cbind(geo_biomass$PC1,geo_biomass$PC2)[hpts, ], col = adjustcolor("grey",alpha.f=0.3), border="grey")
  
}



