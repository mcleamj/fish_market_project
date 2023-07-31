

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

if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(coRanking)){install.packages("coRanking"); library(coRanking)}
if(!require(cluster)){install.packages("cluster"); library(cluster)}
if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(ggpubr)){install.packages("ggpubr"); library(ggpubr)}
if(!require(FD)){install.packages("FD"); library(FD)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(rgdal)){install.packages("rgdal"); library(rgdal)}
if(!require(raster)){install.packages("raster"); library(raster)}
if(!require(sp)){install.packages("sp"); library(sp)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(Matrix)){install.packages("Matrix"); library(Matrix)}
if(!require(rstanarm)){install.packages("rstanarm"); library(rstanarm)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(funrar)){install.packages("funrar"); library(funrar)}
if(!require(performance)){install.packages("performance"); library(performance)}
if(!require(rfishbase)){install.packages("rfishbase"); library(rfishbase)}
if(!require(mFD)){install.packages("mFD"); library(mFD)}
if(!require(forcats)){install.packages("forcats"); library(forcats)}
if(!require(sjPlot)){install.packages("sjPlot"); library(sjPlot)}

library("elbow")

####################################
## DEFINE VARIABLE COLOR FUNCTION ##
####################################

variablecol <- function(colvar, col, clim) {
  ncol <- length(col)
  colvar[colvar < min(clim)] <- NA
  colvar[colvar > max(clim)] <- NA
  rn <- clim[2] - clim[1]
  ifelse (rn != 0, Col <- col[1 + trunc((colvar - clim[1])/rn *
                                          (ncol - 1)+1e-15)], Col <- rep(col[1], ncol))               
  return(Col)
}

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
proj4string(FSM)

kos_shoreline <- raster::shapefile("data/1-01_Shoreline_and_base_layer/kos_shoreline.shp")
kos_shoreline <- spTransform(kos_shoreline, FSM@proj4string)

kos_reefs <- raster::shapefile("data/2-01_Reef_base_layer/kos_coral_reefs.shp")
kos_reefs <- spTransform(kos_reefs, FSM@proj4string)

kos_merged <- raster::union(kos_shoreline, kos_reefs)

kos_buffer <- raster::buffer(kos_merged, width=0.01)

###################################################
## IMPORT SPECIES TRAITS AND TROPHIC INFORMATION ##
###################################################

traits <- read.table("data/clean traits.txt")

species_info <- read.table("data/clean species info.txt")
diet_cols <- kovesi.rainbow(length(unique(species_info$Diet)))

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

auc_tab <- NULL
for(i in 1:ncol(trait_space$x)){
  auc <- AUC_ln_K(R_NX(coranking(traits, trait_space$x[,1:i])))
  print(1:i)
  auc_tab <- rbind(auc_tab, auc)
}

auc_tab <- data.frame(axes=seq(1,10), auc=auc_tab)

elbow(as.data.frame(auc_tab))

PC1 <- trait_space$x[,1]
PC2 <- trait_space$x[,2]
PC3 <- trait_space$x[,3]
PC4 <- trait_space$x[,4]

axes <- data.frame(PC1, PC2, PC3, PC4)

vectors_1_2 <- envfit(trait_space$x[,1:2], traits)
vectors_3_4 <- envfit(trait_space$x[,3:4], traits)

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

########################################################################
## PROXY FOR SPECIES DESIRABILITY USING THE PREDATOR PREFERENCE INDEX
## FROM CHESSON 1978 'Measuring Preference in Selective Predation'
## alphai = (ri/ni) * [1/sum(rj/nj)] 
## where prey type ri or rj is the proportion in the diet and
## ni or nj the proportion in the environment
########################################################################

mean_reef_biomass <- colMeans(reef_fish)
mean_reef_biomass <- mean_reef_biomass/sum(mean_reef_biomass)
sum(mean_reef_biomass)
hist(mean_reef_biomass)
hist(log(mean_reef_biomass))
mean_reef_biomass <- log(mean_reef_biomass)

mean_market_biomass <- colMeans(market_fish)
mean_market_biomass <- mean_market_biomass/sum(mean_market_biomass)
hist(mean_market_biomass)
hist(log(mean_market_biomass))
mean_market_biomass <- log(mean_market_biomass)

# ADD CONSTANT TO MAKE VALUES POSITIVE
constant <- ceiling(abs(min(c(mean_reef_biomass, mean_market_biomass))))

mean_reef_biomass <- mean_reef_biomass + constant
mean_market_biomass <- mean_market_biomass + constant

alpha_index <- data.frame(alpha = (mean_market_biomass / mean_reef_biomass)/sum(mean_market_biomass / mean_reef_biomass))
rownames(alpha_index) <- rownames(traits)

#######################################################
## CALCULATE MEAN ALPHA INDEX FOR REEFS AND LANDINGS ##
#######################################################

reef_pref <- functcomp(as.matrix(alpha_index), as.matrix(reef_log))
reef_FD$reef_pref <- reef_pref$alpha

landings_pref <- functcomp(as.matrix(alpha_index), as.matrix(market_log))
market_FD$market_pref <- landings_pref$alpha

plot(density(reef_pref$alpha),xlim=range(c(reef_pref$alpha,landings_pref$alpha)),
     #ylim=c(0,max(density(landings_pref$alpha)$x)),
     lty="blank",
     xlab=NA, ylab=NA,main=NA,bty="n")
polygon(density(reef_pref$alpha),col=adjustcolor("blue",alpha.f = 0.4),border="blue")
polygon(density(landings_pref$alpha),col=adjustcolor("red",alpha.f = 0.4),border="red")
title("Mean Preference")
legend("topright", legend=c("Reefs","Landings"),
       pch=19,col=c("blue","red"))

######################################
## SPECIES VULNERABILITY TO FISHING ##
######################################

# NOTE - Carangoides orthogrammus was just changed to Ferdauia orthogrammus 
vuln <- species(species_info$species, fields = c("Species", "Vulnerability"))
vuln <- vuln %>% dplyr::rename(species = Species)
vuln <- merge(species_info[,c("species","family")], vuln, by="species",
              all=TRUE)
vuln$Vulnerability[vuln$species=="Carangoides orthogrammus"] <- species("Ferdauia orthogrammus", fields="Vulnerability")$Vulnerability
rownames(vuln) <- vuln$species
vuln$family <- NULL

#################################################################
## FILL MYRPRISTIS SP. WITH MEAN OF OBSERVED SPECIES IN GENERA ##
#################################################################

myripristis_vuln <- species(c("Myripristis adusta","Myripristis berndti","Myripristis murdjan"),
                            fields="Vulnerability")

vuln$Vulnerability[vuln$species=="Myripristis sp"] <- mean(myripristis_vuln$Vulnerability)
sum(is.na(vuln$Vulnerability))

################################
## FUNCTIONAL DISTINCTIVENESS ##
################################

dist_euc <- vegdist(scale(traits), method = "euclidean")
pres <- matrix(ncol=nrow(traits),nrow=1,rep(1))
colnames(pres) <- rownames(traits)
di <- as.data.frame(t(distinctiveness(pres, as.matrix(dist_euc))));colnames(di)<-"di"


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

diet_trait <- as.data.frame(species_info$Diet)
rownames(diet_trait) <- species_info$species

market_trophic <- functcomp(diet_trait, market_log, CWM.type = "all")
names(market_trophic) <- c("grazer","invertivore","microphage","omnivore","piscivore","planktivore")

# TRANSFORMATION TO CONSTRAIN 0 AND 1 VALUES FOR BETA DISTRIBUTION MODEL
market_trophic$grazer <- (market_trophic$grazer*(nrow(market_trophic) - 1) + 0.5) / nrow(market_trophic)
market_trophic$invertivore <- (market_trophic$invertivore*(nrow(market_trophic) - 1) + 0.5) / nrow(market_trophic)
market_trophic$microphage <- (market_trophic$microphage*(nrow(market_trophic) - 1) + 0.5) / nrow(market_trophic)
market_trophic$omnivore <- (market_trophic$omnivore*(nrow(market_trophic) - 1) + 0.5) / nrow(market_trophic)
market_trophic$piscivore <- (market_trophic$piscivore*(nrow(market_trophic) - 1) + 0.5) / nrow(market_trophic)
market_trophic$planktivore <- (market_trophic$planktivore*(nrow(market_trophic) - 1) + 0.5) / nrow(market_trophic)

market_trophic <- as.data.frame(make_relative(as.matrix(market_trophic)))
rowSums(market_trophic)
min(market_trophic)

reef_trophic <- functcomp(diet_trait, reef_log, CWM.type = "all")
names(reef_trophic) <- c("grazer","invertivore","microphage","omnivore","piscivore","planktivore")

# TRANSFORMATION TO CONSTRAIN 0 AND 1 VALUES FOR BETA DISTRIBUTION MODEL
reef_trophic$grazer <- (reef_trophic$grazer*(nrow(reef_trophic) - 1) + 0.5) / nrow(reef_trophic)
reef_trophic$invertivore <- (reef_trophic$invertivore*(nrow(reef_trophic) - 1) + 0.5) / nrow(reef_trophic)
reef_trophic$microphage <- (reef_trophic$microphage*(nrow(reef_trophic) - 1) + 0.5) / nrow(reef_trophic)
reef_trophic$omnivore <- (reef_trophic$omnivore*(nrow(reef_trophic) - 1) + 0.5) / nrow(reef_trophic)
reef_trophic$piscivore <- (reef_trophic$piscivore*(nrow(reef_trophic) - 1) + 0.5) / nrow(reef_trophic)
reef_trophic$planktivore <- (reef_trophic$planktivore*(nrow(reef_trophic) - 1) + 0.5) / nrow(reef_trophic)

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

y <- log(reef_models$fide_PC1+3)
yrep <- posterior_predict(reef_PC1_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Reef Functional Identity (PC1)"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

performance::check_model(reef_PC1_model)

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

y <- reef_models$fide_PC2
yrep <- posterior_predict(reef_PC2_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Reef Functional Identity (PC2)"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

check_model(reef_PC2_model)


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

y <- reef_models$FD_q1+1
yrep <- posterior_predict(reef_entropy_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Reef Functional Entropy"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

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

y <- reef_models$reef_proportion_di
yrep <- posterior_predict(reef_di_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Reef Proportion Distinct"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)


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

y <- reef_models$reef_proportion_vuln
yrep <- posterior_predict(reef_vuln_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Reef Proportion Vulnerable"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)


################################
# PROPORTION OF TROPHIC GROUPS #
################################

reef_troph_model <-  brm(mvbind(grazer, invertivore, microphage,
                                omnivore, piscivore, planktivore) ~ (1  | geographic/site/site_year), 
                         set_prior(class="Intercept", "normal(0,1)"),
                         data=reef_models,
                         family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)
reef_troph_draws <- data.frame(as.matrix(reef_troph_model))

reef_grazer <- reef_troph_draws %>%
  select(starts_with("r_geographic__grazer"))
reef_grazer <- inv_logit_scaled(reef_grazer + reef_troph_draws$b_grazer_Intercept)

reef_invertivore <- reef_troph_draws %>%
  select(starts_with("r_geographic__invertivore"))
reef_invertivore <- inv_logit_scaled(reef_invertivore + reef_troph_draws$b_invertivore_Intercept)

reef_microphage <- reef_troph_draws %>%
  select(starts_with("r_geographic__microphage"))
reef_microphage <- inv_logit_scaled(reef_microphage + reef_troph_draws$b_microphage_Intercept)

reef_omnivore <- reef_troph_draws %>%
  select(starts_with("r_geographic__omnivore"))
reef_omnivore <- inv_logit_scaled(reef_omnivore + reef_troph_draws$b_omnivore_Intercept)

reef_piscivore <- reef_troph_draws %>%
  select(starts_with("r_geographic__piscivore"))
reef_piscivore <- inv_logit_scaled(reef_piscivore + reef_troph_draws$b_piscivore_Intercept)

reef_planktivore <- reef_troph_draws %>%
  select(starts_with("r_geographic__planktivore"))
reef_planktivore <- inv_logit_scaled(reef_planktivore + reef_troph_draws$b_planktivore_Intercept)

# ppc_dens_overlay(reef_models$grazer,  posterior_predict(reef_troph_model)[1:100,,1])
# ppc_dens_overlay(reef_models$invertivore,  posterior_predict(reef_troph_model)[1:100,,2])
# ppc_dens_overlay(reef_models$microphage,  posterior_predict(reef_troph_model)[1:100,,3])
# ppc_dens_overlay(reef_models$omnivore,  posterior_predict(reef_troph_model)[1:100,,4])
# ppc_dens_overlay(reef_models$piscivore,  posterior_predict(reef_troph_model)[1:100,,5])
# ppc_dens_overlay(reef_models$planktivore,  posterior_predict(reef_troph_model)[1:100,,6])

reef_troph <- cbind(apply(reef_grazer,2,median), apply(reef_microphage,2,median), 
                    apply(reef_planktivore,2,median),  apply(reef_omnivore,2,median),
                    apply(reef_invertivore,2,median),apply(reef_piscivore,2,median))
colnames(reef_troph) <- c("1_grazer", "2_microphage", "3_planktivore",
                          "4_omnivore", "5_invertivore", "6_piscivore")
rowSums(reef_troph)


##########################
## MODEL SUMMARY TABLES ##
##########################

tab_model(reef_PC1_model,show.ci=0.90, title="Reef Functional Identity PC1")
tab_model(reef_PC2_model,show.ci=0.90, title="Reef Functional Identity PC2")
tab_model(reef_entropy_model,show.ci=0.90, title="Reef Functional Entropy")
tab_model(reef_di_model,show.ci=0.90, title="Reef Functional Distinctiveness")
tab_model(reef_vuln_model,show.ci=0.90, title="Reef Fishing Vulnerability")
tab_model(reef_troph_model, show.ci=0.90, title="Reef Trophic Guilds")

##################
## FOR LANDINGS ##-------------------------------------------------------------------------------------------------------------------------
##################

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

model_data <- na.omit(landings_models[,c("fide_PC1","geographic","Fisher.Name")])
y <- model_data$fide_PC1
yrep <- posterior_predict(landings_PC1_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Landings Functional Identity (PC1)"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)


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

model_data <- na.omit(landings_models[,c("fide_PC2","geographic","Fisher.Name")])
y <- model_data$fide_PC2
yrep <- posterior_predict(landings_PC2_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Landings Functional Identity (PC2)"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

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

model_data <- na.omit(landings_models[,c("FD_q1","geographic","Fisher.Name")])
y <- model_data$FD_q1+1
yrep <- posterior_predict(landings_entropy_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Landings Functional Entropy"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

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

model_data <- na.omit(landings_models[,c("market_proportion_di","geographic","Fisher.Name")])
y <- model_data$market_proportion_di
yrep <- posterior_predict(landings_di_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Landings Proportion Distinct"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

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

model_data <- na.omit(landings_models[,c("market_proportion_vuln","geographic","Fisher.Name")])
y <- model_data$market_proportion_vuln
yrep <- posterior_predict(landings_vuln_model, draws=100)
ggarrange( ppc_dens_overlay(y, yrep) + ggtitle("Landings Proportion Vulnerable"), 
           ppc_stat(y, yrep, stat="mean"),
           ncol=2, nrow=1)

#############################
# PROPORTION TROPHIC GROUPS #
#############################

landings_troph_model  <- brm(mvbind(grazer, invertivore, microphage,
                                    omnivore, piscivore, planktivore) ~ (1 | geographic) + (1 | Fisher.Name),
                             set_prior(class="Intercept", "normal(0,1)"),
                             data=landings_models,
                             family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )

landings_troph_draws <- data.frame(as.matrix(landings_troph_model))

landings_grazer <- landings_troph_draws %>%
  select(starts_with("r_geographic__grazer"))
landings_grazer <- inv_logit_scaled(landings_grazer + landings_troph_draws$b_grazer_Intercept)

landings_invertivore <- landings_troph_draws %>%
  select(starts_with("r_geographic__invertivore"))
landings_invertivore <- inv_logit_scaled(landings_invertivore + landings_troph_draws$b_invertivore_Intercept)

landings_microphage <- landings_troph_draws %>%
  select(starts_with("r_geographic__microphage"))
landings_microphage <- inv_logit_scaled(landings_microphage + landings_troph_draws$b_microphage_Intercept)

landings_omnivore <- landings_troph_draws %>%
  select(starts_with("r_geographic__omnivore"))
landings_omnivore <- inv_logit_scaled(landings_omnivore + landings_troph_draws$b_omnivore_Intercept)

landings_piscivore <- landings_troph_draws %>%
  select(starts_with("r_geographic__piscivore"))
landings_piscivore <- inv_logit_scaled(landings_piscivore + landings_troph_draws$b_piscivore_Intercept)

landings_planktivore <- landings_troph_draws %>%
  select(starts_with("r_geographic__planktivore"))
landings_planktivore <- inv_logit_scaled(landings_planktivore + landings_troph_draws$b_planktivore_Intercept)

# model_data <- na.omit(landings_models[,c("grazer","invertivore","microphage",
#                                  "omnivore","piscivore","planktivore",
#                                  "geographic","Fisher.Name")])
# ppc_dens_overlay(model_data$grazer,  posterior_predict(landings_troph_model)[1:100,,1])
# ppc_dens_overlay(model_data$invertivore,  posterior_predict(landings_troph_model)[1:100,,2])
# ppc_dens_overlay(model_data$microphage,  posterior_predict(landings_troph_model)[1:100,,3])
# ppc_dens_overlay(model_data$omnivore,  posterior_predict(landings_troph_model)[1:100,,4])
# ppc_dens_overlay(model_data$piscivore,  posterior_predict(landings_troph_model)[1:100,,5])
# ppc_dens_overlay(model_data$planktivore,  posterior_predict(landings_troph_model)[1:100,,6])

landings_troph <- cbind(apply(landings_grazer,2,median), apply(landings_microphage,2,median),
                        apply(landings_planktivore,2,median), apply(landings_omnivore,2,median), 
                        apply(landings_invertivore,2,median), apply(landings_piscivore,2,median))
colnames(landings_troph) <- c("1_grazer", "2_microphage", "3_planktivore",
                              "4_omnivore", "5_invertivore", "6_piscivore")
rowSums(landings_troph)

troph_test <- as.data.frame(landings_troph - reef_troph)
troph_test <- stack(troph_test)
troph_test$geographic <- rep(landings_raw_avg$geographic)

##########################
## MODEL SUMMARY TABLES ##
##########################

tab_model(landings_PC1_model,show.ci=0.90, title="Landings Functional Identity PC1")
tab_model(landings_PC2_model,show.ci=0.90, title="Landings Functional Identity PC2")
tab_model(landings_entropy_model,show.ci=0.90, title="Landings Functional Entropy")
tab_model(landings_di_model,show.ci=0.90, title="Landings Functional Distinctiveness")
tab_model(landings_vuln_model,show.ci=0.90, title="Landings Fishing Vulnerability")
tab_model(landings_troph_model, show.ci=0.90, title="Landings Trophic Guilds")

#################################
## SCALE TO 1 FOR STACKED BARS ##
#################################

reef_troph_scaled <- make_relative(reef_troph)
rowSums(reef_troph_scaled)

landings_troph_scaled <- make_relative(landings_troph)
rowSums(landings_troph_scaled)

######################################################
## DO STACKED BAR CHARTS PER ZONE BY TROPHIC GROUPS ##--------------------------------------------------------------------------------------
######################################################

##########
## REEF ##
##########

reef_diet <- melt(reef_troph_scaled)
colnames(reef_diet)[2] <- "diet"
reef_diet$location <- as.factor(geo_coords$geographic)
reef_diet <- reef_diet[ order(reef_diet$diet, reef_diet$location),]

ggplot(reef_diet, aes(fill=diet,x=location,y=value)) +
  geom_bar(position = "fill",stat="identity") +
  scale_fill_manual(values=diet_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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

landings_diet <- melt(landings_troph_scaled)
colnames(landings_diet)[2] <- "diet"
landings_diet$location <- as.factor(geo_coords$geographic)
landings_diet <- landings_diet[ order(landings_diet$diet, landings_diet$location),]

ggplot(landings_diet, aes(fill=diet,x=location,y=value)) +
  geom_bar(position = "fill",stat="identity") +
  scale_fill_manual(values=diet_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SUMMARY 
landings_diet %>%
  group_by(diet) %>%
  summarise(mean(value))

landings_diet %>%
  group_by(diet) %>%
  summarise(max(value))

#######################################
## TRAIT SPACE FIGURE WITH CENTROIDS ##---------------------------------------------------------------------------------------------------------
#######################################

###########################
## PANEL A - TRAIT SPACE ##
###########################

graphics.off()
par(mfrow=c(2,2))

plot(PC1, PC2, cex=0, cex.axis=1.2, cex.lab=1.2)

hpts <- chull(cbind(PC1,PC2))
hpts <- c(hpts, hpts[1])
polygon(cbind(PC1,PC2)[hpts, ], col = adjustcolor("grey",alpha.f=0.3), border="grey")

points(PC1,PC2,pch=21,col=1,bg="grey",cex=1.5)

mtext("A", font=2, cex=1.5, adj=-0.1, line=0.5)

#######################
## PANEL B - VECTORS ##
#######################

plot(PC1, PC2, cex=0, cex.axis=1.2, cex.lab=1.2)
plot(vectors_1_2,col=1)

mtext("B", font=2, cex=1.5, adj=-0.1, line=0.5)

#################################
## PANEL C - DIET CONVEX HULLS ##
#################################

species_info$diet_col <- ifelse(species_info$Diet=="Grazer",diet_cols[1],
                                ifelse(species_info$Diet=="Microphage",diet_cols[2],
                                       ifelse(species_info$Diet=="Planktivore",diet_cols[3],
                                              ifelse(species_info$Diet=="Omnivore",diet_cols[4],
                                                     ifelse(species_info$Diet=="Invertivore",diet_cols[5],
                                                            ifelse(species_info$Diet=="Piscivore",diet_cols[6],NA))))))


plot(PC1, PC2, cex=0, cex.axis=1.2, cex.lab=1.2)

legend("bottomleft",legend=c("Grazers","Microphages","Planktivores","Omnivores","Invertivores","Piscivores"),
       col=diet_cols,
       pch=19,
       cex=1)

trophic_scores <- data.frame(species_info$Diet, PC1, PC2)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Grazer")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[1],alpha.f=0.1), border=diet_cols[1])

trophic_scores <- data.frame(species_info$Diet, PC1, PC2)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Microphage")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[2],alpha.f=0.1), border=diet_cols[2])

trophic_scores <- data.frame(species_info$Diet, PC1, PC2)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Planktivore")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[3],alpha.f=0.1), border=diet_cols[3])

trophic_scores <- data.frame(species_info$Diet, PC1, PC2)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Omnivores")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[4],alpha.f=0.1), border=diet_cols[4])

trophic_scores <- data.frame(species_info$Diet, PC1, PC2)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Invertivore")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[5],alpha.f=0.1), border=diet_cols[5])

trophic_scores <- data.frame(species_info$Diet, PC1, PC2)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Piscivore")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[6],alpha.f=0.1), border=diet_cols[6])

points(PC1, PC2,
       pch=19,col=species_info$diet_col, cex=1.5)

mtext("C", font=2, cex=1.5, adj=-0.1, line=0.5)


#########################
## PANEL D - CENTROIDS ##
#########################

# NOTE - WE CAN SET ZOOM TO ANY SCALE ON PANEL D
# BY PLAYING WITH X AND Y LIMITS

plot(PC1, PC2, cex=0, cex.lab=1.2, cex.axis=1.2,
     xlim=range(1.5*c(PC1_centroid_reefs, PC1_centroid_landings)),
     ylim=range(2*c(PC2_centroid_reefs, PC2_centroid_landings)))

points(PC1_centroid_reefs, PC2_centroid_reefs,
       pch=21, col=1, bg="blue", cex=2)

points(PC1_centroid_landings, PC2_centroid_landings,
       pch=21, col=1, bg="red", cex=2)

legend("topright", legend=c("Reef Observations", "Market Landings"),
       pch=19, col=c("blue","red"), pt.cex=1.5)

mtext("D", font=2, cex=1.5, adj=-0.1, line=0.5)



##############################
## MAP FUNCTIONAL DIVERSITY ##---------------------------------------------------------------------------------------------------------
##############################


################################
# FUNCTIONAL IDENTITY/CENTROID #
################################

graphics.off()
par(mfrow=c(1,2))

## REEF 
map_color <- variablecol(PC1_centroid_reefs, col=brewer.blues(length(PC1_centroid_reefs)),
                         clim=range(c(PC1_centroid_reefs, PC1_centroid_landings)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = PC1_centroid_reefs,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.blues(length(PC1_centroid_reefs)),
          clim=range(c(PC1_centroid_reefs, PC1_centroid_landings)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"PC1 Centroid")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = PC1_centroid_reefs,
          colkey = FALSE, xlab=NA,ylab=NA,
          cex=4,
          col = 1, bg=map_color,
          add=TRUE)
plot(kos_buffer, add=TRUE)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(PC1_centroid_landings, col=brewer.blues(length(PC1_centroid_landings)),
                         clim=range(c(PC1_centroid_reefs, PC1_centroid_landings)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = PC1_centroid_landings,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.blues(length(PC1_centroid_landings)),
          clim=range(c(PC1_centroid_reefs, PC1_centroid_landings)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"PC1 Centroid")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = PC1_centroid_landings,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Fisheries Landings")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(PC1_centroid_landings),
                      diff = PC1_centroid_reefs - PC1_centroid_landings)
di_diff$sign <- ifelse(di_diff$diff < 0, "negative","positive")

ggplot(di_diff,aes(geographic,diff)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green"))


#########################
# FUNCTIONAL ENTROPY  #
#########################

graphics.off()
par(mfrow=c(1,2))

## REEF 
map_color <- variablecol(reef_entropy, col=brewer.oranges(length(reef_entropy)),
                         clim=range(c(reef_entropy, landings_entropy)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = reef_entropy,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.oranges(length(reef_entropy)),
          clim=range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Functional Entropy")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = reef_entropy,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(reef_entropy/max(c(reef_entropy,landings_entropy)))*4,
          #cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_entropy, col=brewer.oranges(length(landings_entropy)),
                         clim=range(c(reef_entropy, landings_entropy)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = landings_entropy,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.oranges(length(landings_entropy)),
          clim=range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Functional Entropy")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = landings_entropy,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_entropy/max(c(reef_entropy,landings_entropy)))*4,
          #cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Fisheries Landings")
mtext("(b)",line=1,font=2,adj=-0.1,cex=1.75)



## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(landings_entropy),
                      diff = reef_entropy - landings_entropy)
di_diff$sign <- ifelse(di_diff$diff < 0, "negative","positive")

ggplot(di_diff,aes(geographic,diff)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green"))


#######################
# PROPORTION DISTINCT #
#######################

graphics.off()
par(mfrow=c(1,2))

## REEF 
map_color <- variablecol(reef_di, col=brewer.purples(length(reef_di)),
                         clim=range(c(reef_di, landings_di)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = reef_di,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.purples(length(reef_di)),
          clim=range(c(reef_di, landings_di)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1.5,"Proportion Distinct Species")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = reef_di,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(reef_di/max(c(reef_di,landings_di)))*4,
          #cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_di, col=brewer.purples(length(landings_di)),
                         clim=range(c(reef_di, landings_di)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = landings_di,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.purples(length(landings_di)),
          clim=range(c(reef_di, landings_di)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Proportion Distinct Species")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = landings_di,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_di/max(c(reef_di,landings_di)))*4,
          #cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Fisheries Landings")
mtext("(b)",line=1,font=2,adj=-0.1,cex=1.75)


## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(landings_di),
                      diff = reef_di - landings_di)
di_diff$sign <- ifelse(di_diff$diff < 0, "negative","positive")

ggplot(di_diff,aes(geographic,diff)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green"))


#########################
# PROPORTION VULNERABLE #
#########################

graphics.off()
par(mfrow=c(1,2))

## REEF 
map_color <- variablecol(reef_vuln, col=brewer.reds(length(reef_vuln)),
                         clim=range(c(reef_vuln, landings_vuln)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = reef_vuln,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.reds(length(reef_vuln)),
          clim=range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Proportion Vulnerable Species")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = reef_vuln,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(reef_vuln/max(c(reef_vuln,landings_vuln)))*4,
          #cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_vuln, col=brewer.reds(length(landings_vuln)),
                         clim=range(c(reef_vuln, landings_vuln)))
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = landings_vuln,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.reds(length(landings_vuln)),
          clim=range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1.5,"Proportion Vulnerable Species")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=21, colvar = landings_vuln,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_vuln/max(c(reef_vuln,landings_vuln)))*4,
          #cex=4,
          col = 1, bg=map_color)
plot(kos_buffer, add=TRUE)
title("Fisheries Landings")
mtext("(b)",line=1,font=2,adj=-0.1,cex=1.75)


## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(landings_vuln),
                      diff = reef_vuln - landings_vuln)
di_diff$sign <- ifelse(di_diff$diff < 0, "negative","positive")

ggplot(di_diff,aes(geographic,diff)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green"))


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

model_data <- reef_drivers[,c("fide_PC1","PC1_all","PC2_all","reef_area_5km","fishing_events",
                              "tot_grav_pop","geographic","site","site_year")]

ppc_dens_overlay(log(reef_PC1_drivers_model$data$fide_PC1+3), posterior_predict(reef_PC1_drivers_model,draws=100)[1:100,])
ppc_error_scatter_avg(log(reef_PC1_drivers_model$data$fide_PC1+3), posterior_predict(reef_PC1_drivers_model,draws=100))

y <- log(reef_PC1_drivers_model$data$fide_PC1+3)
yrep <- posterior_predict(reef_PC1_drivers_model,draws=100)
plot(colMeans(yrep), y)
ppc_stat(y, yrep, stat="mean")

#########################
# FUNCTIONAL DISPERSION #
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

ppc_dens_overlay(reef_entropy_drivers_model$data$FD_q1+1, posterior_predict(reef_entropy_drivers_model,draws=100)[1:100,])

y <- (reef_entropy_drivers_model$data$FD_q1 + 1)
yrep <- posterior_predict(reef_entropy_drivers_model,draws=100)
plot(colMeans(yrep), y)
ppc_stat(y, yrep, stat="mean")

##########################
## MODEL SUMMARY TABLES ##
##########################

tab_model(reef_PC1_drivers_model, show.ci=0.90, title="Reef Functional Identity PC1 Drivers")
tab_model(reef_entropy_drivers_model, show.ci=0.90, title="Reef Functional Identity PC1 Drivers")

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
                                    #z_score_2sd(wave_energy) + 
                                    z_score_2sd(tot_grav_pop) +
                                    #moon_phase + 
                                    z_score_2sd(wind) +
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
         #'b_z_score_2sdwave_energy', 
         'b_z_score_2sdwind',
         #'b_moon_phasemediummoon', 'b_moon_phasebigmoon',
         'b_gearShallowBottomFishing','spearguns', 'b_z_score_2sdnum_fishers_lines')
mcmc_intervals(landings_PC1_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

model_data <- landings_drivers[,c("fide_PC1","PC1_centroid_reefs","fishing_events",
                                  "wave_energy","tot_grav_pop","moon_phase","wind",
                                  "gear","num_fishers_lines","geographic","Fisher.Name")]
model_data <- na.omit(model_data)
y <- model_data$fide_PC1
yrep <- posterior_predict(landings_PC1_drivers_model,draws=100)
ppc_dens_overlay(y, yrep[1:100,])
ppc_stat(y, yrep, stat="mean")
plot(colMeans(yrep), y)

#########################
# FUNCTIONAL ENTROPY  #
#########################

landings_entropy_drivers_model <- brm((FD_q1+1) ~ 
                                           z_score_2sd(reef_entropy+1) + 
                                           z_score_2sd(fishing_events) + 
                                           #z_score_2sd(wave_energy) + 
                                           z_score_2sd(tot_grav_pop) +
                                           #moon_phase + 
                                           z_score_2sd(wind) +
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
         #'b_z_score_2sdwave_energy', 
         'b_z_score_2sdwind',
         #'b_moon_phasemediummoon', 'b_moon_phasebigmoon',
         'b_gearShallowBottomFishing','spearguns', 'b_z_score_2sdnum_fishers_lines')
mcmc_intervals(landings_entropy_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

model_data <- landings_drivers[,c("FD_q1","reef_entropy","fishing_events",
                                  "wave_energy","tot_grav_pop","moon_phase","wind",
                                  "gear","num_fishers_lines","geographic","Fisher.Name")]
model_data <- na.omit(model_data)
y <- model_data$FD_q1+1
yrep <- posterior_predict(landings_entropy_drivers_model,draws=100)
ppc_dens_overlay(y, yrep[1:100,])
ppc_stat(y, yrep, stat="mean")
plot(colMeans(yrep), y)

##########################
## MODEL SUMMARY TABLES ##
##########################

tab_model(landings_PC1_drivers_model, show.ci=0.90, title="Landings Functional Identity PC1 Drivers")
tab_model(landings_entropy_drivers_model, show.ci=0.90, title="Landings Functional Identity PC1 Drivers")

######################################################
## RELATION OF MEAN PREFERENCE TO DIVERSITY METRICS ##
######################################################

pref_model <- brm(market_pref ~ z_score_2sd(fide_PC1) + z_score_2sd(FD_q1) +
                    (1 | geographic) + (1 | Fisher.Name),
                  family=Gamma(link="log"), data=landings_drivers, chains=4, iter=2000)

summary(pref_model, prob=0.90)
posterior_summary(pref_model, pars=c("fide_PC1","FD_q1"), prob=c(0.1,0.9))
pref_drivers <- data.frame(as.matrix(pref_model)) %>%
  dplyr::select(b_z_score_2sdfide_PC1, b_z_score_2sdFD_q1 )
mcmc_intervals(pref_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + 
  geom_vline(xintercept = 0) +
  scale_y_discrete(
    labels = c("b_z_score_2sdFD_q1" = "Functional Entropy",
               "b_z_score_2sdfide_PC1" = "Functional Identity (PC1)")) +
  ggtitle("Mean Fisher Preference")
tab_model(pref_model, show.ci=0.90, title="Preference and Trait Diversity Model")


############################################
# CPUE HIGHER OR LOWER FOR DESIRABLE FISH? #
############################################

CPUE <- brm(log(CPUE) ~ z_score_2sd(market_pref) +
              (1 | geographic) + (1 | Fisher.Name),
            family=gaussian, data=landings_drivers, chains=4, iter=2000)

summary(CPUE, prob=0.9)
marginal_effects(CPUE)
CPUE_drivers <- as.matrix(CPUE)
tab_model(CPUE, show.ci=0.90, title="CPUE and Preferece Model")

# CPUE_drivers <- (CPUE_drivers[,1:2]) 
# mcmc_intervals(CPUE_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
#                outer_size = 1,
#                inner_size = 4,
#                point_size = 6) + geom_vline(xintercept = 0)






######################################################
## CODE TO MAKE TRAIT SPACE FIGURE WITH PC3 AND PC4 ##
######################################################


#################
# PC 3 CENTROID #
#################

reef_PC3_model <- brm(log(fide_PC3+1) ~ (1 | geographic/site/site_year), 
                      set_prior(class="Intercept", "normal(0,1)"),
                      data=reef_models,
                      family=gaussian, chains=4, iter=2000)

reef_PC3_draws <- data.frame(as.matrix(reef_PC3_model))
reef_PC3 <- reef_PC3_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_PC3 <- exp(reef_PC3 + reef_PC3_draws$b_Intercept) - 1
colnames(reef_PC3) <- geo_coords$geographic
PC3_centroid_reefs <- apply(reef_PC3, 2, median)
plot(reef_raw_avg$fide_PC3, PC3_centroid_reefs, pch=19)

#################
# PC 4 CENTROID #
#################

reef_PC4_model <- brm(fide_PC4 ~ (1 | geographic/site/site_year), 
                      set_prior(class="Intercept", "normal(0,1)"),
                      data=reef_models,
                      family=gaussian, chains=4, iter=2000)

reef_PC4_draws <- data.frame(as.matrix(reef_PC4_model))
reef_PC4 <- reef_PC4_draws %>%
  dplyr::select("r_geographic.East.Intercept.":"r_geographic.West.Intercept.")
reef_PC4 <- reef_PC4 + reef_PC4_draws$b_Intercept
colnames(reef_PC4) <- geo_coords$geographic
PC4_centroid_reefs <- apply(reef_PC4, 2, median)
plot(reef_raw_avg$fide_PC4, PC4_centroid_reefs,pch=19)





#################
# PC 3 CENTROID #
#################

landings_PC3_model <- brm(fide_PC3 ~ (1 | geographic) + (1 | Fisher.Name),
                          set_prior(class="Intercept", "normal(0,1)"),
                          data=landings_models,
                          family=gaussian, chains=4, iter=2000)

landings_PC3_draws <- data.frame(as.matrix(landings_PC3_model))
landings_PC3 <- landings_PC3_draws %>%
  select(starts_with("r_geographic"))
landings_PC3 <- landings_PC3 + landings_PC3_draws$b_Intercept
colnames(landings_PC3) <- geo_coords$geographic
PC3_centroid_landings <- apply(landings_PC3, 2, median)
plot(landings_raw_avg$fide_PC3, PC3_centroid_landings,pch=19)
plot(landings_raw_avg$fide_PC3, PC3_centroid_landings,pch=19,
     xlim=range(c(landings_raw_avg$fide_PC3, PC3_centroid_landings)),
     ylim=range(c(landings_raw_avg$fide_PC3, PC3_centroid_landings)))
abline(0,1)


#################
# PC 4 CENTROID #
#################

landings_PC4_model <- brm(fide_PC4 ~ (1 | geographic) + (1 | Fisher.Name),
                          set_prior(class="Intercept", "normal(0,1)"),
                          data=landings_models,
                          family=gaussian, chains=4, iter=2000)

landings_PC4_draws <- data.frame(as.matrix(landings_PC4_model))
landings_PC4 <- landings_PC4_draws %>%
  select(starts_with("r_geographic"))
landings_PC4 <- landings_PC4 + landings_PC4_draws$b_Intercept
colnames(landings_PC4) <- geo_coords$geographic
PC4_centroid_landings <- apply(landings_PC4, 2, median)
plot(landings_raw_avg$fide_PC4, PC4_centroid_landings,pch=19)


###########################
## PANEL A - TRAIT SPACE ##
###########################

graphics.off()
par(mfrow=c(2,2))

plot(PC3, PC4, cex=0, cex.axis=1.2, cex.lab=1.2)

hpts <- chull(cbind(PC3,PC4))
hpts <- c(hpts, hpts[1])
polygon(cbind(PC3,PC4)[hpts, ], col = adjustcolor("grey",alpha.f=0.3), border="grey")

points(PC3,PC4,pch=21,col=1,bg="grey",cex=1.5)

mtext("A", font=2, cex=1.5, adj=-0.1, line=0.5)

#######################
## PANEL B - VECTORS ##
#######################

plot(PC3, PC4, cex=0, cex.axis=1.2, cex.lab=1.2)
plot(vectors_3_4,col=1)

mtext("B", font=2, cex=1.5, adj=-0.1, line=0.5)

#################################
## PANEL C - DIET CONVEX HULLS ##
#################################

species_info$diet_col <- ifelse(species_info$Diet=="Grazer",diet_cols[1],
                                ifelse(species_info$Diet=="Microphage",diet_cols[2],
                                       ifelse(species_info$Diet=="Planktivore",diet_cols[3],
                                              ifelse(species_info$Diet=="Omnivore",diet_cols[4],
                                                     ifelse(species_info$Diet=="Invertivore",diet_cols[5],
                                                            ifelse(species_info$Diet=="Piscivore",diet_cols[6],NA))))))


plot(PC3, PC4, cex=0, cex.axis=1.2, cex.lab=1.2)

legend("bottomleft",legend=c("Grazers","Microphages","Planktivores","Omnivores","Invertivores","Piscivores"),
       col=diet_cols,
       pch=19,
       cex=1)

trophic_scores <- data.frame(species_info$Diet, PC3, PC4)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Grazer")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[1],alpha.f=0.1), border=diet_cols[1])

trophic_scores <- data.frame(species_info$Diet, PC3, PC4)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Microphage")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[2],alpha.f=0.1), border=diet_cols[2])

trophic_scores <- data.frame(species_info$Diet, PC3, PC4)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Planktivore")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[3],alpha.f=0.1), border=diet_cols[3])

trophic_scores <- data.frame(species_info$Diet, PC3, PC4)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Omnivores")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[4],alpha.f=0.1), border=diet_cols[4])

trophic_scores <- data.frame(species_info$Diet, PC3, PC4)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Invertivore")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[5],alpha.f=0.1), border=diet_cols[5])

trophic_scores <- data.frame(species_info$Diet, PC3, PC4)
diet_hull <- subset(trophic_scores, trophic_scores$species_info.Diet=="Piscivore")[,2:3]
hpts <- chull(diet_hull)
hpts <- c(hpts, hpts[1])
polygon(diet_hull[hpts, ], col = adjustcolor(diet_cols[6],alpha.f=0.1), border=diet_cols[6])

points(PC3, PC4,
       pch=19,col=species_info$diet_col, cex=1.5)

mtext("C", font=2, cex=1.5, adj=-0.1, line=0.5)


#########################
## PANEL D - CENTROIDS ##
#########################

# NOTE - WE CAN SET ZOOM TO ANY SCALE ON PANEL D
# BY PLAYING WITH X AND Y LIMITS

plot(PC3, PC4, cex=0, cex.lab=1.2, cex.axis=1.2,
     xlim=range(1.5*c(PC3_centroid_reefs, PC3_centroid_landings)),
     ylim=range(2*c(PC4_centroid_reefs, PC4_centroid_landings)))

points(PC3_centroid_reefs, PC4_centroid_reefs,
       pch=21, col=1, bg="blue", cex=2)

points(PC3_centroid_landings, PC4_centroid_landings,
       pch=21, col=1, bg="red", cex=2)

legend("topright", legend=c("Reef Observations", "Market Landings"),
       pch=19, col=c("blue","red"), pt.cex=1.5)

mtext("D", font=2, cex=1.5, adj=-0.1, line=0.5)
