

#################################################################
#################################################################
## CODE FOR TOUCHIE ET AL. TRAIT SELECTIVITY IN REEF FISHERIES ##
#################################################################
#################################################################

###########################################
# MUST INSTALL ELBOW PACKAGE FROM GITHUB ##
###########################################

#devtools::install_github("ahasverus/elbow", build_vignettes = TRUE)

##################################
## INSTALL AND LIBRARY PACKAGES ##
##################################

if(!require(reshape2)){install.packages("reshape2"); library(reshape2)}
if(!require(sjPlot)){install.packages("sjPlot"); library(sjPlot)}
if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(cluster)){install.packages("cluster"); library(cluster)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(ggpubr)){install.packages("ggpubr"); library(ggpubr)}
if(!require(FD)){install.packages("FD"); library(FD)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(raster)){install.packages("raster"); library(raster)}
if(!require(sf)){install.packages("sf"); library(sf)}
if(!require(Matrix)){install.packages("Matrix"); library(Matrix)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(funrar)){install.packages("funrar"); library(funrar)}
if(!require(performance)){install.packages("performance"); library(performance)}
if(!require(rfishbase)){install.packages("rfishbase"); library(rfishbase)}
if(!require(mFD)){install.packages("mFD"); library(mFD)}
if(!require(forcats)){install.packages("forcats"); library(forcats)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}

library("elbow")

#############################################################
## DEFINE FUNCTION FOR STANDARDIZATION USING 2 SD'S 
## GELMAN AND HILL RECOMMEND SETTING PREDICTOR VARIABLES
## TO MEAN 0 AND SD 0.5 SO THAT EFFECT SIZES BETWEEN
## CONTINUOUS AND CATEGORICAL VARIABLES ARE MORE COMPARABLE
#############################################################

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

############################
## LOAD KOSRAE SHAPEFILES ##
############################

FSM <- readRDS("data/gadm36_FSM_0_sp.rds")
st_crs(FSM)

kos_shoreline <- st_read("data/1-01_Shoreline_and_base_layer/kos_shoreline.shp")
kos_shoreline <- st_transform(kos_shoreline, st_crs(FSM))

kos_reefs <- st_read("data/2-01_Reef_base_layer/kos_coral_reefs.shp")
kos_reefs <-  st_transform(kos_reefs, st_crs(FSM))

kos_merged <- st_union(kos_shoreline, kos_reefs)

kos_polys <- st_read("data/geo_polygons_layer/kos_geo_polygons.shp")
kos_polys <- st_transform(kos_polys, st_crs(FSM))

kos_utm <- st_transform(kos_merged, 32655)

kos_buffer <- st_buffer(st_make_valid(kos_utm), dist = 1000)  # 1000 m

kos_buffer <- st_transform(kos_buffer, 4326)

marina_coords <- read.csv("data/marina coords.csv")

###################################################
## IMPORT SPECIES TRAITS AND TROPHIC INFORMATION ##
###################################################

traits <- read.table("data/clean traits 2025.txt") %>%
  tibble::rownames_to_column("species") %>%
  arrange(species) %>%
  column_to_rownames("species")

species_info <- read.table("data/clean species info 2025.txt") %>%
  arrange(species)
diet_cols <- kovesi.rainbow(length(unique(species_info$Diet)))
names(diet_cols) <- c(
  "Grazer",
  "Microphage",
  "Planktivore",
  "Omnivore",
  "Invertivore",
  "Piscivore"
)

####################################################
## IMPORT AND WRANGLE REEF AND LANDINGS DATA SETS ##
####################################################

# REEF DATA 
reef_data <- read.csv("data/clean reef biomass data 2025.csv")

reef_meta <- reef_data %>%
  select(c(Site_SPC_year, site, SPC_Rep, year, 
           Observer, Lat, Lon, site_year, geographic))

reef_fish <- reef_data %>%
  select(-c(Site_SPC_year, site, SPC_Rep, year, 
            Observer, Lat, Lon, site_year, geographic)) %>%
  select(sort(names(.)))

reef_log <- log10(reef_fish+1)

colnames(reef_log) <- gsub("\\.", " ", colnames(reef_log))

reef_log <- as.matrix(reef_log)

rownames(reef_log) <- paste("com",seq(1,nrow(reef_log),1))

# MARKET DATA
market_data <- read.csv("data/clean market data.csv")
market_data$num_fishers_lines <- as.numeric(market_data$num_fishers_lines)

market_meta <- market_data %>%
  dplyr::select(c(Monitoring.Code, geographic, Fisher.Name, CPUE, gear, 
                  season, month, moon_phase, moon_light, wind, num_fishers_lines))

market_fish <-market_data %>%
  dplyr::select(-c(Monitoring.Code, geographic, Fisher.Name, CPUE, gear, 
                   season, month, moon_phase, moon_light, wind, num_fishers_lines))  %>%
  select(sort(names(.)))

market_log <- log10(market_fish+1)

colnames(market_log) <- gsub("\\.", " ", colnames(market_log))

market_log <- as.matrix(market_log)

rownames(market_log) <- paste("com",seq(1,nrow(market_log),1))

##################################
## LOAD GEOGRRAPHIC COORDINATES ##
##################################

geo_coords <-read.csv("data/geographic coordinates.csv")

#######################
## BUILD TRAIT SPACE ##
#######################

trait_space <- prcomp(traits, scale. = TRUE)
graphics.off()
biplot(trait_space)
summary(trait_space)

######################################################
## CALCULATE OPTIMAL NUMBER OF AXES FOR TRAIT SPACE ##
######################################################

# auc_tab <- NULL
# for(i in 1:ncol(trait_space$x)){
#   auc <- AUC_ln_K(R_NX(coranking(traits, trait_space$x[,1:i])))
#   print(1:i)
#   auc_tab <- rbind(auc_tab, auc)
# }
# 
# auc_tab <- data.frame(axes=seq(1,10), auc=auc_tab)
# 
# elbow(as.data.frame(auc_tab))

PC1 <- trait_space$x[,1]
PC2 <- trait_space$x[,2]
PC3 <- trait_space$x[,3]
PC4 <- trait_space$x[,4]

axes <- data.frame(PC1, PC2, PC3, PC4)

vectors_1_2 <- envfit(axes[,1:2], traits)
vectors_3_4 <- envfit(axes[,3:4], traits)

######################################################
## CALCUALTE FUNCTIONAL DIVERSITY FOR BOTH DATASETS ##
######################################################

###########
## REEFS ##
###########


reef_fide <- alpha.fd.multidim(
  as.matrix(axes),
  reef_log,
  ind_vect = c("fide"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE,
  verbose = TRUE
)

reef_fide <- reef_fide$functional_diversity_indices
reef_fide$com <- rownames(reef_fide)
reef_FD <- reef_fide

sp_dist <- as.matrix(vegdist(scale(traits), method = "euclidean"))
reef_entropy <- alpha.fd.hill(reef_log, sp_dist, tau="mean", q=1)
reef_entropy <- reef_entropy$asb_FD_Hill

reef_FD <- data.frame(reef_FD, reef_entropy = reef_entropy)

# CWM TRAITS 
cwm_reef <- functcomp(traits, reef_log)
mean(cwm_reef$max_length)
sd(cwm_reef$max_length)

##############
## LANDINGS ##
##############

market_fide <- alpha.fd.multidim(
  as.matrix(axes),
  market_log,
  ind_vect = c("fide"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE,
  verbose = TRUE
)

market_fide <- market_fide$functional_diversity_indices
market_fide$com <- rownames(market_fide)
market_FD <- market_fide

market_entropy <- alpha.fd.hill(market_log, sp_dist, tau="mean", q=1)
market_entropy <- market_entropy$asb_FD_Hill

market_FD <- data.frame(market_FD, market_entropy=market_entropy)

# CWM
market_cwm <- functcomp(traits, market_log)
mean(market_cwm$max_length)
sd(market_cwm$max_length)

######################################
## SPECIES VULNERABILITY TO FISHING ##
######################################

# NOTE - Carangoides orthogrammus WAS RECENTLY CHANGED TO Ferdauia orthogrammus
# NOTE - Carangoides ferdau WAS RECENTLY CHANGED TO Ferdauia ferdau
# AS OF JULY 2026 FISH BASE HAS NOT UPDATED THESE
# NOTE FISHBASE VULNERABILITY VALUES ARE NOW BINNED TO 0, 25, 50, ETC. FOR MANY SPECIES
# SO HERE WE USE HISTORICAL ARCHIVE OF VALUES AT ORIGINAL TIME OF STUDY (2021)
# USING VERSION 21.06
vuln <- rfishbase::species(species_info$species, fields = c("Species", "Vulnerability"), version = "21.06")
vuln <- vuln %>% dplyr::rename(species = Species)
vuln <- merge(species_info[,c("species","family")], vuln, by="species",
              all=TRUE)
vuln$Vulnerability[vuln$species=="Ferdauia orthogrammus"] <- species("Carangoides orthogrammus", fields="Vulnerability", version = "21.06")$Vulnerability
vuln$Vulnerability[vuln$species=="Ferdauia ferdau"] <- species("Carangoides ferdau", fields="Vulnerability", version = "21.06")$Vulnerability
rownames(vuln) <- vuln$species
vuln$family <- NULL

############################################################################
## FILL MYRPRISTIS SP. WITH MEAN OF CONGENERIC SPECIES OBSERVED IN KOSRAE ##
############################################################################

myripristis_vuln <- species(c("Myripristis adusta","Myripristis berndti","Myripristis murdjan"),
                            fields="Vulnerability", version = "21.06")

vuln$Vulnerability[vuln$species=="Myripristis sp"] <- mean(myripristis_vuln$Vulnerability)
sum(is.na(vuln$Vulnerability))


################################
## FUNCTIONAL DISTINCTIVENESS ##
################################

dist_euc <- vegdist(scale(traits), method = "euclidean")
pres <- matrix(ncol=nrow(traits),nrow=1,rep(1))
colnames(pres) <- rownames(traits)
di <- as.data.frame(t(distinctiveness(pres, as.matrix(dist_euc))));colnames(di)<-"di"
di$di <-round(di$di, 2)
di <- di %>%
  rownames_to_column("species")

# WHAT TRAITS DRIVE DISTINCTIVENESS?

di_glm <- glm(di$di ~ scale(traits), family = Gamma(link="log") )
summary(di_glm)

#######################################################
## CLASSIFY TOP DISTINCT AND TOP VULNERABLE SPECIES  ##
#######################################################

# DISTINCT 
quars <- quantile(di$di,probs=c(0.25,0.5,0.75))
di$quartile <- ifelse(di$di <= quars[1],"Q1",
                      ifelse(di$di > quars[1] & 
                               di$di <= quars[2],"Q2",
                             ifelse(di$di > quars[2] & 
                                      di$di <= quars[3],"Q3",
                                    ifelse(di$di > quars[3],"Q4",NA))))

boxplot(di$di ~ di$quartile)
table(di$quartile)
hist(di$di, xlab="Distinctiveness", main="Histogram of Distinctiveness",col="darkgrey",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
# abline(v=quars[1],lty=1,col=2,lwd=3)
# abline(v=quars[2],lty=1,col=2,lwd=3)
abline(v=quars[3],lty=1,col=2,lwd=3)
polygon(x=c(quars[3],quars[3],par("usr")[2],par("usr")[2]),
        y=c(par("usr")[3],par("usr")[4],par("usr")[4],par("usr")[3]),
        border=NA, col=adjustcolor("red",alpha.f = 0.2))
text(5.75,18,"Top 25% Most Distinct",cex=1.5)

# VULNERABLE
quars <- quantile(vuln$Vulnerability,probs=c(0.25,0.5,0.75))
vuln$quartile <- ifelse(vuln$Vulnerability <= quars[1],"Q1",
                        ifelse(vuln$Vulnerability > quars[1] & 
                                 vuln$Vulnerability <= quars[2],"Q2",
                               ifelse(vuln$Vulnerability > quars[2] & 
                                        vuln$Vulnerability <= quars[3],"Q3",
                                      ifelse(vuln$Vulnerability > quars[3],"Q4",NA))))

boxplot(vuln$Vulnerability ~ vuln$quartile)
table(vuln$quartile)
hist(vuln$Vulnerability, xlab="Vulnerability", main="Histogram of Vulnerability",
     col="darkgrey",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
# abline(v=quars[1],lty=1,col=2,lwd=3)
# abline(v=quars[2],lty=1,col=2,lwd=3)
abline(v=quars[3],lty=1,col=2,lwd=3)
polygon(x=c(quars[3],quars[3],par("usr")[2],par("usr")[2]),
        y=c(par("usr")[3],par("usr")[4],par("usr")[4],par("usr")[3]),
        border=NA, col=adjustcolor("red",alpha.f = 0.2))
text(65,20,"Top 25% Most Vulnerable",cex=1.5)


######################################################################
## CALCULATE PROPORTION OF DI SPECIES PER FISHING EVENT AND PER SPC ##
######################################################################

market_di <- market_log[,colnames(market_log) %in% rownames(di)[di$quartile=="Q4"]]
market_proportion_di <- rowSums(market_di)/rowSums(market_log)

reef_di <- reef_log[,colnames(reef_log) %in% rownames(di)[di$quartile=="Q4"]]
reef_proportion_di <- rowSums(reef_di)/rowSums(reef_log)

##############################################################################
## CALCULATE PROPORTION OF VULNERABLE SPECIES PER FISHING EVENT AND PER SPC ##
##############################################################################

market_vuln <- market_log[,colnames(market_log) %in% rownames(vuln)[vuln$quartile=="Q4"]]
market_proportion_vuln <- rowSums(market_vuln)/rowSums(market_log)

reef_vuln <- reef_log[,colnames(reef_log) %in% rownames(vuln)[vuln$quartile=="Q4"]]
reef_proportion_vuln <- rowSums(reef_vuln)/rowSums(reef_log)

######################################################################
## CALCULATE TROPHIC GROUP PROPORTION PER FISHING EVENT AND PER SPC ##
######################################################################

# Trait prep
diet_trait <- species_info %>%
  select(species, Diet) %>%
  column_to_rownames("species")

identical(rownames(diet_trait), colnames(market_log))

guild_names <- c("grazer", "invertivore", "microphage", "omnivore", "piscivore", "planktivore")

# --- MARKET TROPHIC ---
market_trophic <- functcomp(diet_trait, market_log, CWM.type = "all")
names(market_trophic) <- guild_names

# Squeeze all columns simultaneously & re-normalize
n_market <- nrow(market_trophic)
market_trophic <- (market_trophic * (n_market - 1) + 0.5) / n_market
market_trophic <- as.data.frame(make_relative(as.matrix(market_trophic)))

rowSums(market_trophic)
min(market_trophic)

# --- REEF TROPHIC ---
reef_trophic <- functcomp(diet_trait, reef_log, CWM.type = "all")
names(reef_trophic) <- guild_names

# Squeeze all columns simultaneously & re-normalize
n_reef <- nrow(reef_trophic)
reef_trophic <- (reef_trophic * (n_reef - 1) + 0.5) / n_reef
reef_trophic <- as.data.frame(make_relative(as.matrix(reef_trophic)))

rowSums(reef_trophic)
min(reef_trophic)


######################################################################
## INTERCEPT-ONLY HIERARCHICAL MODELS TO CALCULATE GEOGRAPHIC MEANS ##---------------------------------------------------------------------------
######################################################################

###############
## FOR REEFS ##----------------------------------------------------------------------------------------------------------------------------------
###############

reef_raw_avg <- aggregate(cbind(reef_FD,reef_proportion_di,reef_proportion_vuln, reef_trophic), 
                          by=list(reef_meta$geographic),FUN=mean,na.rm=TRUE)
colnames(reef_raw_avg)[1] <- "geographic"

# TRANSFORMATION TO CONSTRAIN 0 AND 1 VALUES FOR BETA DISTRIBUTION MODEL
reef_models <- data.frame(reef_meta, reef_FD, reef_proportion_di, reef_proportion_vuln, reef_trophic)
reef_models$reef_proportion_di <- (reef_models$reef_proportion_di*(nrow(reef_models) - 1) + 0.5) / nrow(reef_models)
reef_models$reef_proportion_vuln <- (reef_models$reef_proportion_vuln*(nrow(reef_models) - 1) + 0.5) / nrow(reef_models)


#################
# PC 1 CENTROID #
#################

reef_PC1_model <- brm(log(fide_PC1+3) ~ (1 | geographic/site/site_year), 
                      set_prior(class="Intercept", "normal(0,1)"),
                      data=reef_models,
                      family=gaussian, chains=4, iter=2000)

reef_PC1_draws <- data.frame(as.matrix(reef_PC1_model))
reef_PC1 <- reef_PC1_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_PC1 <- exp(reef_PC1 + reef_PC1_draws$b_Intercept) - 3
colnames(reef_PC1) <- geo_coords$geographic
PC1_centroid_reefs <- apply(reef_PC1, 2, median)
plot(reef_raw_avg$fide_PC1, PC1_centroid_reefs, pch=19)

#################
# PC 2 CENTROID #
#################

reef_PC2_model <- brm(fide_PC2 ~ (1 | geographic/site/site_year), 
                      set_prior(class="Intercept", "normal(0,1)"),
                      data=reef_models,
                      family=gaussian, chains=4, iter=2000)

reef_PC2_draws <- data.frame(as.matrix(reef_PC2_model))
reef_PC2 <- reef_PC2_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_PC2 <- reef_PC2 + reef_PC2_draws$b_Intercept
colnames(reef_PC2) <- geo_coords$geographic
PC2_centroid_reefs <- apply(reef_PC2, 2, median)
plot(reef_raw_avg$fide_PC2, PC2_centroid_reefs,pch=19)

######################
# FUNCTIONAL ENTROPY #
######################

reef_entropy_model <- brm((FD_q1+1) ~ (1 | geographic/site/site_year), 
                          set_prior(class="Intercept", "normal(0,1)"),
                          data=reef_models,
                          family=Gamma(link="log"),chains=4, iter=2000)

reef_entropy_draws <- data.frame(as.matrix(reef_entropy_model))
reef_entropy <- reef_entropy_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_entropy <- (exp(reef_entropy + reef_entropy_draws$b_Intercept)) - 1
colnames(reef_entropy) <- geo_coords$geographic
reef_entropy <- apply(reef_entropy, 2, median)
plot(reef_raw_avg$FD_q1, reef_entropy,pch=19)


#######################
# PROPORTION DISTINCT #
#######################

reef_di_model <- brm(reef_proportion_di ~ (1 | geographic/site/site_year), 
                     set_prior(class="Intercept", "normal(0,1)"),
                     data=reef_models,
                     family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)

reef_di_draws <- data.frame(as.matrix(reef_di_model))
reef_di <- reef_di_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_di <- inv_logit_scaled(reef_di + reef_di_draws$b_Intercept)
colnames(reef_di) <- geo_coords$geographic
reef_di <- apply(reef_di, 2, median)
plot(reef_raw_avg$reef_proportion_di, reef_di, pch=19)


#########################
# PROPORTION VULNERABLE #
#########################

reef_vuln_model <- brm(reef_proportion_vuln ~ (1 | geographic/site/site_year), 
                       set_prior(class="Intercept", "normal(0,1)"),
                       data=reef_models,
                       family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)

reef_vuln_draws <- data.frame(as.matrix(reef_vuln_model))
reef_vuln <- reef_vuln_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_vuln <- inv_logit_scaled(reef_vuln + reef_vuln_draws$b_Intercept)
colnames(reef_vuln) <- geo_coords$geographic
reef_vuln <- apply(reef_vuln, 2, median)
plot(reef_raw_avg$reef_proportion_vuln, reef_vuln, pch=19)


################################
# PROPORTION OF TROPHIC GROUPS 
# USING DIRICHLET MODEL 
################################

# Dirichlet model
reef_troph_model <- brm(
  cbind(grazer, invertivore, microphage,
        omnivore, piscivore, planktivore) ~ 
    (1 | geographic/site/site_year),
  data = reef_models,
  family = dirichlet(),
  chains = 4, iter = 2000
)

reef_geo_grid <- data.frame(geographic = unique(reef_models$geographic))

reef_ep <- posterior_epred(
  reef_troph_model,
  newdata = reef_geo_grid,
  re_formula = ~(1 | geographic) # Retain geographic, set site & site_year to 0
)

dimnames(reef_ep)[[2]] <- reef_geo_grid$geographic
dimnames(reef_ep)[[3]] <- c(
  "grazer", "invertivore", "microphage",
  "omnivore", "piscivore", "planktivore"
)

# Averages across all 4000 draws -> 11 rows (regions) x 6 columns (categories)
reef_troph <- colMeans(reef_ep)

reef_troph <- as.data.frame(reef_troph) %>%
  rownames_to_column("geographic") %>%
  filter(!is.na(geographic) & geographic != "NA.") %>%
  arrange(geographic) %>%
  mutate(geographic = gsub("\\.", " ", geographic)) %>%
  column_to_rownames("geographic")

rowSums(reef_troph)



################## ##---------------------------------------------------------------------------------------------------------------------
##################
## FOR LANDINGS ##-------------------------------------------------------------------------------------------------------------------------
##################
################## ##----------------------------------------------------------------------------------------------------------------------

landings_raw_avg <- aggregate(cbind(market_FD,market_proportion_di,market_proportion_vuln, market_trophic),
                              by=list(market_meta$geographic),FUN=mean,na.rm=TRUE)
colnames(landings_raw_avg)[1] <- "geographic"

landings_models <- data.frame(market_meta, market_FD, market_proportion_di, market_proportion_vuln, market_trophic)
landings_models$market_proportion_di <- (landings_models$market_proportion_di*(nrow(landings_models) - 1) + 0.5) / nrow(landings_models)
landings_models$market_proportion_vuln <- (landings_models$market_proportion_vuln*(nrow(landings_models) - 1) + 0.5) / nrow(landings_models)


#################
# PC 1 CENTROID #
#################

landings_PC1_model <- brm(fide_PC1 ~ (1 | geographic) + (1 | Fisher.Name),
                          set_prior(class="Intercept", "normal(0,1)"),
                          data=landings_models,
                          family=gaussian, chains=4, iter=2000)

landings_PC1_draws <- data.frame(as.matrix(landings_PC1_model))
landings_PC1 <- landings_PC1_draws %>%
  select(starts_with("r_geographic"))
landings_PC1 <- landings_PC1 + landings_PC1_draws$b_Intercept
colnames(landings_PC1) <- geo_coords$geographic
PC1_centroid_landings <- apply(landings_PC1, 2, median)
plot(landings_raw_avg$fide_PC1, PC1_centroid_landings,pch=19)
plot(landings_raw_avg$fide_PC1, PC1_centroid_landings,pch=19,
     xlim=range(c(landings_raw_avg$fide_PC1, PC1_centroid_landings)),
     ylim=range(c(landings_raw_avg$fide_PC1, PC1_centroid_landings)))
abline(0,1)


#################
# PC 2 CENTROID #
#################

landings_PC2_model <- brm(fide_PC2 ~ (1 | geographic) + (1 | Fisher.Name),
                          set_prior(class="Intercept", "normal(0,1)"),
                          data=landings_models,
                          family=gaussian, chains=4, iter=2000)

landings_PC2_draws <- data.frame(as.matrix(landings_PC2_model))
landings_PC2 <- landings_PC2_draws %>%
  select(starts_with("r_geographic"))
landings_PC2 <- landings_PC2 + landings_PC2_draws$b_Intercept
colnames(landings_PC2) <- geo_coords$geographic
PC2_centroid_landings <- apply(landings_PC2, 2, median)
plot(landings_raw_avg$fide_PC2, PC2_centroid_landings,pch=19)


######################
# FUNCTIONAL ENTROPY #
######################

landings_entropy_model <- brm((FD_q1+1) ~ (1 | geographic) + (1 | Fisher.Name),
                              set_prior(class="Intercept", "normal(0,1)"),
                              data=landings_models,
                              family=Gamma(link="log"),chains=4, iter=2000 )

landings_entropy_draws <- data.frame(as.matrix(landings_entropy_model))

landings_entropy <- landings_entropy_draws %>%
  select(starts_with("r_geographic"))
landings_entropy <- (exp(landings_entropy + landings_entropy_draws$b_Intercept)) - 1
colnames(landings_entropy) <- geo_coords$geographic
landings_entropy <- apply(landings_entropy, 2, median)
plot(landings_raw_avg$FD_q1, landings_entropy,pch=19)


#################
# PROPORTION DI #
#################

landings_di_model <- brm(market_proportion_di~ (1 | geographic) + (1 | Fisher.Name), 
                         set_prior(class="Intercept", "normal(0,1)"),
                         data=landings_models,
                         family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )
landings_di_draws <- data.frame(as.matrix(landings_di_model))

landings_di <- landings_di_draws %>%
  select(starts_with("r_geographic"))
landings_di <- inv_logit_scaled(landings_di + landings_di_draws$b_Intercept)
colnames(landings_di) <- geo_coords$geographic
landings_di <- apply(landings_di, 2, median)
plot(landings_raw_avg$market_proportion_di, landings_di, pch=19)


#########################
# PROPORTION VULNERABLE #
#########################

landings_vuln_model <- brm(market_proportion_vuln ~ (1 | geographic) + (1 | Fisher.Name),
                           set_prior(class="Intercept", "normal(0,1)"),
                           data=landings_models,
                           family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )

landings_vuln_draws <- data.frame(as.matrix(landings_vuln_model))
landings_vuln <- landings_vuln_draws %>%
  select(starts_with("r_geographic"))
landings_vuln <- inv_logit_scaled(landings_vuln + landings_vuln_draws$b_Intercept)
colnames(landings_vuln) <- geo_coords$geographic
landings_vuln <- apply(landings_vuln, 2, median)
plot(landings_raw_avg$market_proportion_vuln, landings_vuln, pch=19)


#############################
# PROPORTION TROPHIC GROUPS #
#############################

# Dirichlet model for fisheries landings
landings_troph_model <- brm(
  cbind(grazer, invertivore, microphage, 
        omnivore, piscivore, planktivore) ~ (1 | geographic) + (1 | Fisher.Name),
  data = landings_models,
  family = dirichlet(),
  chains = 4, iter = 2000
)

geo_grid <- data.frame(geographic = unique(landings_models$geographic))

landings_ep <- posterior_epred(
  landings_troph_model,
  newdata = geo_grid,
  re_formula = ~(1 | geographic)
)

dimnames(landings_ep)[[2]] <- geo_grid$geographic
dimnames(landings_ep)[[3]] <- c(
  "grazer", "invertivore", "microphage",
  "omnivore", "piscivore", "planktivore"
)

# Averages across all 4000 draws -> 11 rows (regions) x 6 columns (categories)
landings_troph <- colMeans(landings_ep)

landings_troph <- as.data.frame(landings_troph) %>%
  rownames_to_column("geographic") %>%
  filter(!is.na(geographic) & geographic != "NA.") %>%
  arrange(geographic) %>%
  mutate(geographic = gsub("\\.", " ", geographic)) %>%
  column_to_rownames("geographic")

rowSums(landings_troph)


#############################################################
## PREPARE DATA FORD BAR CHARTS PER ZONE BY TROPHIC GROUPS ##--------------------------------------------------------------------------------------
#############################################################

##########
## REEF ##
##########

reef_diet <- melt(reef_troph)
colnames(reef_diet)[1] <- "diet"
reef_diet$location <- as.factor(geo_coords$geographic)
reef_diet <- reef_diet[ order(reef_diet$diet, reef_diet$location),]

# SUMMARY 
reef_diet %>%
  group_by(diet) %>%
  summarise(mean(value))

reef_diet %>%
  group_by(diet) %>%
  summarise(max(value))

##############
## LANDINGS ##
##############

landings_diet <- melt(landings_troph)
colnames(landings_diet)[1] <- "diet"
landings_diet$location <- as.factor(geo_coords$geographic)
landings_diet <- landings_diet[ order(landings_diet$diet, landings_diet$location),]

# SUMMARY 
landings_diet %>%
  group_by(diet) %>%
  summarise(mean(value))

landings_diet %>%
  group_by(diet) %>%
  summarise(max(value))

# ##############################
# TWO PANEL STACKED BAR CHART ##
################################

source("scripts/figures/Figure stacked bars.R")

##################
## SAVE METRICS ##
##################

reef_metrics <- data.frame(PC1_centroid_reefs, PC2_centroid_reefs, reef_entropy, reef_di, reef_vuln)
saveRDS(reef_metrics, "reef_metrics_original.rds")

landings_metrics <- data.frame(PC1_centroid_landings, PC2_centroid_landings, landings_entropy, landings_di, landings_vuln)
saveRDS(landings_metrics, "landings_metrics_original.rds")


###########################################################
## MODELS FOR DRIVERS OF REEF AND MARKET TRAIT DIVERSITY ##---------------------------------------------------------------------------
###########################################################

reef_drivers <- read.table("data/site_covariates.txt")
reef_drivers$Lat <- NULL; reef_drivers$Lon <- NULL; reef_drivers$geographic <- NULL
length(unique(reef_drivers$site))

#########################################
## IMPORT COVARIATES FROM GOOGLE EARTH ##
#########################################

google_earth <- read.csv("data/KOS google earth covariates.csv")
reef_drivers <- merge(reef_drivers, google_earth, by="site")

##########################################
## CALCULATE FISHING EVENTS PER SECTION ##
##########################################

events_per_section <- market_meta %>%
  group_by(geographic) %>%
  summarise(n_distinct(Monitoring.Code))
names(events_per_section) <- c("geographic","fishing_events")

#################
## REEF MODELs ##
#################

reef_drivers <- merge(reef_models, reef_drivers, by="site")
length(unique(reef_drivers$site))
reef_drivers <- merge(reef_drivers, events_per_section, by="geographic")

#########################################
## ADD BENTHIC DATA AT SITE-YEAR LEVEL ##
#########################################

benthic_site_year <- read.table("data/all_benthic_site_year_estimates.txt")
reef_drivers <- merge(reef_drivers, benthic_site_year, by="site_year")

#################
# PC1 CENTROID ##
#################

drivers <- reef_drivers %>%
  select(Lat, PC1_all, PC2_all, reef_area_5km,
         dist_MPA, fishing_events,
         tot_grav_pop, wave_energy)
corrplot::corrplot(cor(drivers))

reef_PC1_drivers_model <- brm(log(fide_PC1+3) ~ 
                                z_score_2sd(fishing_events) + 
                                z_score_2sd(tot_grav_pop) +
                                z_score_2sd(reef_area_5km) + 
                                z_score_2sd(PC1_all) +
                                z_score_2sd(PC2_all) +
                                (1 | geographic/site/site_year),
                              c(set_prior(class="Intercept", "normal(0,1)"),
                                set_prior(class="b", "normal(0,1)")),
                              data=reef_drivers,
                              family=gaussian, chains=4, iter=3000,
                              control = list(adapt_delta=0.95))

reef_PC1_drivers <- data.frame(as.matrix(reef_PC1_drivers_model))
summary(reef_PC1_drivers_model, prob=0.90)
reef_PC1_drivers <- reef_PC1_drivers %>%
  dplyr::select('b_z_score_2sdtot_grav_pop', 'b_z_score_2sdfishing_events', 
                'b_z_score_2sdreef_area_5km', 'b_z_score_2sdPC1_all', 'b_z_score_2sdPC2_all')
color_scheme_set("darkgray")
mcmc_intervals(reef_PC1_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)


#########################
# FUNCTIONAL ENTROPY #
#########################

reef_entropy_drivers_model <- brm((FD_q1 +1)  ~ 
                                    z_score_2sd(PC1_all) +
                                    z_score_2sd(PC2_all) +
                                    z_score_2sd(reef_area_5km) + 
                                    z_score_2sd(fishing_events) + 
                                    z_score_2sd(tot_grav_pop) +
                                    (1 | geographic/site/site_year), 
                                  c(set_prior(class="Intercept", "normal(0,1)"),
                                    set_prior(class="b", "normal(0,1)")),
                                  data=reef_drivers,
                                  family=Gamma(link="log"),chains=4, iter=3000,
                                  control = list(adapt_delta=0.95))

summary(reef_entropy_drivers_model, prob=0.90)
reef_entropy_drivers <- data.frame(as.matrix(reef_entropy_drivers_model))
reef_entropy_drivers <-  reef_entropy_drivers %>%
  dplyr::select('b_z_score_2sdtot_grav_pop', 'b_z_score_2sdfishing_events', 
                'b_z_score_2sdreef_area_5km', 'b_z_score_2sdPC1_all', 'b_z_score_2sdPC2_all')
mcmc_intervals(reef_entropy_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)


##########################
## MODEL SUMMARY TABLES ##
##########################

tab_model(reef_PC1_drivers_model, show.ci=0.90, title="Reef Functional Identity PC1 Drivers",
          transform = NULL,
          file="reef_PC1_drivers_model.html")
webshot("reef_PC1_drivers_model.html",
        "reef_PC1_drivers_model.png",
        zoom = 3)

tab_model(reef_entropy_drivers_model, show.ci=0.90, title="Reef Functional Entropy Drivers",
          transform = NULL,
          file="reef_entropy_drivers_model.html")
webshot("reef_entropy_drivers_model.html",
        "reef_entropy_drivers_model.png",
        zoom = 3)

saveRDS(reef_PC1_drivers_model, "reef_PC1_drivers_model.rds")
saveRDS(reef_entropy_drivers_model, "reef_entropy_drivers_model.rds")

#####################
## LANDINGS MODELs ##--------------------------------------------------------------------------------
#####################

reef_intercepts <- data.frame(PC1_centroid_reefs, reef_entropy)
reef_intercepts$geographic <- rownames(reef_intercepts)

landings_drivers <- merge(landings_models, reef_intercepts, by="geographic")

landings_drivers$fishing_events <- with(events_per_section,
                                        fishing_events[match(landings_drivers$geographic, geographic)])

#################################################
## ADD WAVE ENERGY AT THE GEO SECTION MIDPOINT ##
#################################################

geo_wave <- read.table("data/Geographic wave energy.txt")
geo_wave <- geo_wave %>%
  dplyr::select(geographic, energy) %>%
  dplyr::rename(wave_energy = energy)
landings_drivers <- merge(landings_drivers, geo_wave, by="geographic")

####################################################
## ADD MARKET GRAVITY AT THE GEO SECTION MIDPOINT ##
####################################################

geo_gravity <- read.csv("data/Geographic section gravity.csv")
geo_gravity <- geo_gravity %>%
  select(geographic, tot_grav_pop)
landings_drivers <- merge(landings_drivers, geo_gravity, by="geographic")

#################################
## SET FACTOR VARIABLE LEVELES ##
#################################

landings_drivers$moon_phase <- as.factor(landings_drivers$moon_phase)
landings_drivers$moon_phase <- forcats::fct_relevel(landings_drivers$moon_phase, "low moon","medium moon","big moon")
levels(landings_drivers$moon_phase)

landings_drivers$gear <- as.factor(landings_drivers$gear)
landings_drivers$gear <- relevel(landings_drivers$gear, "Spear")
levels(landings_drivers$gear)

################
# PC1 CENTROID #
################

landings_PC1_drivers_model <- brm(fide_PC1 ~ 
                                    z_score_2sd(PC1_centroid_reefs) + 
                                    z_score_2sd(fishing_events) + 
                                    z_score_2sd(wave_energy) + 
                                    z_score_2sd(tot_grav_pop) +
                                    gear + 
                                    z_score_2sd(num_fishers_lines) + 
                                    (1 | geographic) + (1 | Fisher.Name),
                                  c(set_prior(class="Intercept", "normal(0,1)"),
                                    set_prior(class="b", "normal(0,1)")),
                                  family=gaussian, data=landings_drivers, chains=4, iter=4000,
                                  control = list(adapt_delta=0.95))

summary(landings_PC1_drivers_model, prob=0.90)
landings_PC1_drivers <- data.frame(as.matrix(landings_PC1_drivers_model))
landings_PC1_drivers$spearguns <- rep(0)
landings_PC1_drivers <- landings_PC1_drivers %>%
  select('b_z_score_2sdPC1_centroid_reefs','b_z_score_2sdtot_grav_pop', 'b_z_score_2sdfishing_events', 
         'b_z_score_2sdwave_energy', 
         'b_gearShallowBottomFishing','spearguns', 'b_z_score_2sdnum_fishers_lines')
mcmc_intervals(landings_PC1_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)


#########################
# FUNCTIONAL ENTROPY  #
#########################

landings_entropy_drivers_model <- brm((FD_q1+1) ~ 
                                        z_score_2sd(reef_entropy+1) + 
                                        z_score_2sd(fishing_events) + 
                                        z_score_2sd(wave_energy) + 
                                        z_score_2sd(tot_grav_pop) +
                                        gear + 
                                        z_score_2sd(num_fishers_lines) + 
                                        (1 | geographic) + (1 | Fisher.Name),
                                      c(set_prior(class="Intercept", "normal(0,1)"),
                                        set_prior(class="b", "normal(0,1)")),
                                      family=Gamma(link="log"), data=landings_drivers, chains=4, iter=2000,
                                      control = list(adapt_delta=0.95))

summary(landings_entropy_drivers_model, prob=0.90)
landings_entropy_drivers <- data.frame(as.matrix(landings_entropy_drivers_model))
landings_entropy_drivers$spearguns <- rep(0)
landings_entropy_drivers <- landings_entropy_drivers %>%
  select('b_z_score_2sdreef_entropyP1','b_z_score_2sdtot_grav_pop', 'b_z_score_2sdfishing_events', 
         'b_z_score_2sdwave_energy', 
         'b_gearShallowBottomFishing','spearguns', 'b_z_score_2sdnum_fishers_lines')
mcmc_intervals(landings_entropy_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

