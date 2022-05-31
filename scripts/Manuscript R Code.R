

#################################################################
#################################################################
## CODE FOR TOUCHIE ET AL. TRAIT SELECTIVITY IN REEF FISHERIES ##
#################################################################
#################################################################

###########################################
# MUST INSTALL THESE PACKAGES FROM GITHUB #
###########################################

remotes::install_github("ropensci/rfishbase")
remotes::install_github("CmlMagneville/mFD")
devtools::install_github("ahasverus/elbow", build_vignettes = TRUE)

###############################################
## NOTE - multidimFD FUNCTION MUST BE LOADED 
## THE CODE IS AVAILABLE IN THE
## SCRIPTS FOLDER ON THE GITHUB REPOSITORY
###############################################

##################################
## INSTALL AND LIBRARY PACKAGES ##
##################################

if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(cluster)){install.packages("cluster"); library(cluster)}
if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(FD)){install.packages("FD"); library(FD)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(rethinking)){install.packages("rethinking"); library(rethinking)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(rstanarm)){install.packages("rstanarm"); library(rstanarm)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(funrar)){install.packages("funrar"); library(funrar)}
if(!require(coRanking)){install.packages("coRanking"); library(coRanking)}
if(!require(phylolm)){install.packages("phylolm"); library(phylolm)}
if(!require(fishtree)){install.packages("fishtree"); library(fishtree)}
#if(!require(ggridges)){install.packages("ggridges"); library(ggridges)}
if(!require(performance)){install.packages("performance"); library(performance)}
library("rfishbase")
library("mFD")
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

###############################
## LOAD MICRONESIA SHAPEFILE ##
###############################

FSM <- readRDS("data/gadm36_FSM_0_sp.rds")

###################################################
## IMPORT SPECIES TRAITS AND TROPHIC INFORMATION ##
###################################################

traits <- read.table("data/clean traits.txt")

species_info <- read.table("data/clean species info.txt")

####################################################
## IMPORT AND WRANGLE REEF AND LANDINGS DATA SETS ##
####################################################

reef_data <- read.csv("data/clean reef data.csv")
reef_meta <- reef_data[,1:6]
reef_fish <- reef_data[,7:ncol(reef_data)]
reef_log <- log10(reef_fish+1)
colnames(reef_log) <- gsub("\\.", " ", colnames(reef_log))
reef_log <- as.matrix(reef_log)
rownames(reef_log) <- paste("com",seq(1,nrow(reef_log),1))

market_data <- read.csv("data/clean market data.csv")
market_data$num_fishers_lines <- as.numeric(market_data$num_fishers_lines)
market_meta <- market_data[,1:10]
market_fish <- market_data[,11:ncol(market_data)]
market_log <- log10(market_fish+1)
colnames(market_log) <- gsub("\\.", " ", colnames(market_log))
market_log <- as.matrix(market_log)
rownames(market_log) <- paste("com",seq(1,nrow(market_log),1))

##################################
## LOAD GEOGRRAPHIC COORDINATES ##
##################################

geo_coords <-read.table("data/geographic coordinates.txt")

#######################
## BUILD TRAIT SPACE ##
#######################

trait_space <- prcomp(traits, scale. = TRUE)
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

reef_FD <- multidimFD(as.matrix(axes),weight=reef_log, check_species_pool=TRUE, verb=TRUE,
                      folder_plot=NULL, nm_asb_plot=NULL, Faxes_plot=NULL, Faxes_nm_plot=NULL, 
                      plot_pool=FALSE, col_bg="grey90", col_sp_pool="grey30", pch_sp_pool="+", cex_sp_pool=1,
                      pch_sp=21, col_sp="#1145F0", transp=50 ) 

sp_dist <- as.matrix(vegdist(scale(traits), method = "euclidean"))

reef_entropy <- alpha.fd.hill(reef_log, sp_dist, tau="mean", q=1)
reef_entropy <- reef_entropy$asb_FD_Hill

reef_FD <- data.frame(reef_FD, reef_entropy = reef_entropy)

# LANDINGS #
market_FD <- multidimFD(as.matrix(axes),weight=market_log, check_species_pool=TRUE, verb=TRUE,
                        folder_plot=NULL, nm_asb_plot=NULL, Faxes_plot=NULL, Faxes_nm_plot=NULL, 
                        plot_pool=FALSE, col_bg="grey90", col_sp_pool="grey30", pch_sp_pool="+", cex_sp_pool=1,
                        pch_sp=21, col_sp="#1145F0", transp=50 ) 

market_entropy <- alpha.fd.hill(market_log, sp_dist, tau="mean", q=1)
market_entropy <- market_entropy$asb_FD_Hill

market_FD <- data.frame(market_FD, market_entropy=market_entropy)


########################################################################
## PROXY FOR SPECIES DESIRABILITY AS DEVIATION FROM 1:1 CATCH PATTERN ##
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

deviation <- data.frame( deviation = mean_market_biomass - mean_reef_biomass)
rownames(deviation) <- rownames(traits)

graphics.off()
plot(mean_reef_biomass, mean_market_biomass,cex=0,
     pch=19, ylim=range(mean_reef_biomass, mean_market_biomass),
     xlim=range(mean_reef_biomass, mean_market_biomass),
     xlab="Reef Biomass Proportion (log scale)",
     ylab="Market Biomass Proportion (log scale)")
abline(0,1)
segments(x0=mean_reef_biomass, y0=mean_reef_biomass,
         x1=mean_reef_biomass, y1=mean_market_biomass)
points(mean_reef_biomass, mean_market_biomass,pch=21,
       col=1,bg=ifelse(deviation$deviation>0,"red","green"))
legend("topleft", legend=c("Over-targeted","Under-targeted"),
       pch=19,col=c("red","green"))

alpha_index <- data.frame(alpha = (mean_market_biomass / mean_reef_biomass)/sum(mean_market_biomass / mean_reef_biomass))
rownames(alpha_index) <- rownames(traits)

plot(deviation$deviation, alpha_index$alpha)

#####################################################
## CALCULATE MEAN DEVIATION FOR REEFS AND LANDINGS ##
#####################################################

reef_dev <- functcomp(as.matrix(deviation), as.matrix(reef_log))
reef_FD$reef_deviation <- reef_dev$deviation

landings_dev <- functcomp(as.matrix(deviation), as.matrix(market_log))
market_FD$market_deviation <- landings_dev$deviation

plot(density(reef_dev$deviation),xlim=range(c(reef_dev-1,landings_dev+1)),
     ylim=c(0,max(density(landings_dev$deviation)$y)),lty="blank",
     xlab=NA, ylab=NA,main=NA,bty="n")
polygon(density(reef_dev$deviation),col=adjustcolor("blue",alpha.f = 0.4),border="blue")
polygon(density(landings_dev$deviation),col=adjustcolor("red",alpha.f = 0.4),border="red")
title("Mean Deviation/Preference")
legend("topleft", legend=c("Reefs","Landings"),
       pch=19,col=c("blue","red"))

######################################
## SPECIES VULNERABILITY TO FISHING ##
######################################

vuln <- data.frame(Species=colnames(reef_fish),
                   Vulnerability= species(species_info$species, fields = "Vulnerability"))
rownames(vuln) <- rownames(traits)

#################################################################
## FILL MYRPRISTIS SP. WITH MEAN OF OBSERVED SPECIES IN GENERA ##
#################################################################

myripristis_vuln <- species(c("Myripristis adusta","Myripristis berndti","Myripristis murdjan"),
                            fields="Vulnerability")

vuln$Vulnerability[vuln$Species=="Myripristis.sp"] <- mean(myripristis_vuln$Vulnerability)
vuln$Vulnerability[vuln$Species=="Parupeneus.bifasciatus"] <- 31

################################
## FUNCTIONAL DISTINCTIVENESS ##
################################

dist_euc <- vegdist(scale(traits), method = "euclidean")
pres <- matrix(ncol=nrow(traits),nrow=1,rep(1))
colnames(pres) <- rownames(traits)
di <- as.data.frame(t(distinctiveness(pres, as.matrix(dist_euc))));colnames(di)<-"di"

#######################################
##  WHICH TRAITS ARE MOST DISTINCT ? ##
#######################################

di_trait_data <- data.frame(di=di$di, scale(traits))

di_traits_model <- brm(di ~ ., family = lognormal(), data=di_trait_data, chains=4)
di_traits <- as.matrix(di_traits_model)
di_traits <- (di_traits[,2:11]) 
di_traits <- di_traits[,order(abs(colMeans(di_traits)), decreasing = TRUE)]
color_scheme_set("darkgray")
mcmc_intervals(di_traits, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)


#####################################################
## ARE DISTINCTIVENESS AND VULNERABILITY REALTED ? ##
#####################################################

graphics.off()
plot(vuln$Vulnerability ~ di$di, pch=21,col=1,bg="grey",
     xlab="Distinctiveness", ylab="Vulnerability",cex=2,
     cex.axis=1.5,cex.lab=1.5)
abline(lm(vuln$Vulnerability ~ di$di),lty=2,lwd=3,col="skyblue")
cor.test(vuln$Vulnerability, di$di)
legend("topleft",legend="r = 0.44",cex=1.5,pt.cex=0)

##########################################
## PHYLOGENTICALLY-CORRECTED REGRESSION ##
##########################################

phy_lm_data <- data.frame(di, vuln)
sp_list <- rownames(phy_lm_data)

tree <- fishtree_complete_phylogeny(sp_list)
tree <- tree[[1]]
tree$tip.label
tree$tip.label <- gsub("_",".",tree$tip.label)

phy_lm_data <- phy_lm_data[phy_lm_data$Species %in% tree$tip.label,]
rownames(phy_lm_data) <- phy_lm_data$Species

fit = phylolm(Vulnerability ~ di, data=phy_lm_data,phy=tree,
              model="BM")
summary(fit) # r2 = 0.2, r=0.44

############################################################################
## IS REEF BIOMASS GREATER OR LESSER FOR DISTINCT OR VULNERABLE SPECIES ? ##
############################################################################

par(mfrow=c(1,1))

plot(log(colMeans(reef_fish)) ~ vuln$Vulnerability, pch=19)
vuln_reef_model <- lm(log(colMeans(reef_fish)) ~ vuln$Vulnerability)
summary(vuln_reef_model)
hist(resid(vuln_reef_model))
cor.test(log(colMeans(reef_fish)), vuln$Vulnerability)

plot(log(colMeans(reef_fish)) ~ di$di, pch=19)
di_reef_model <- lm(log(colMeans(reef_fish)) ~ di$di)
summary(di_reef_model)
hist(resid(di_reef_model))
cor.test(log(colMeans(reef_fish)), di$di)

#########################################################################
## IS BIOMASS OF LANDINGS GREATER FOR DISTINCT OR VULNERABLE SPECIES ? ##
#########################################################################

par(mfrow=c(1,1))

plot(log(colMeans(market_fish)) ~ vuln$Vulnerability, pch=21,col=1,bg="grey",cex=2,
     xlab="Vulnerability", ylab="Market Biomass (log scale)", cex.axis=1.5, cex.lab=1.5)
vuln_catch_model <- lm(log(colMeans(market_fish)) ~ vuln$Vulnerability)
abline(vuln_catch_model, lty=2, col="skyblue", lwd=3)
summary(vuln_catch_model)
#hist(resid(vuln_catch_model))
cor.test(log(colMeans(market_fish)), vuln$Vulnerability)
legend("bottomright",legend="r = 0.39",cex=1.5,pt.cex=0)

plot(log(colMeans(market_fish)) ~ di$di, pch=19)
di_catch_model <- lm(log(colMeans(market_fish)) ~ di$di)
summary(di_catch_model)
hist(resid(di_catch_model))
cor.test(log(colMeans(market_fish)), di$di)

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
market_trophic[market_trophic==0] <- 0.001
market_trophic[market_trophic==1] <- 0.999
market_trophic <- as.data.frame(make_relative(as.matrix(market_trophic)))
rowSums(market_trophic)
min(market_trophic)

reef_trophic <- functcomp(diet_trait, reef_log, CWM.type = "all")
names(reef_trophic) <- c("grazer","invertivore","microphage","omnivore","piscivore","planktivore")
reef_trophic[reef_trophic==0] <- 0.001
reef_trophic[reef_trophic==1] <- 0.999
reef_trophic <- as.data.frame(make_relative(as.matrix(reef_trophic)))
rowSums(reef_trophic)
min(reef_trophic)

######################################################################
## INTERCEPT-ONLY HIERARCHICAL MODELS TO CALCULATE GEOGRAPHIC MEANS ##
######################################################################

###############
## FOR REEFS ##
###############

reef_raw_avg <- aggregate(cbind(reef_FD,reef_proportion_di,reef_proportion_vuln, reef_trophic), 
                          by=list(reef_meta$geographic),FUN=mean,na.rm=TRUE)
colnames(reef_raw_avg)[1] <- "geographic"

reef_models <- data.frame(reef_meta, reef_FD, reef_proportion_di, reef_proportion_vuln, reef_trophic)
reef_models$reef_proportion_di <- (reef_models$reef_proportion_di*(nrow(reef_models) - 1) + 0.5) / nrow(reef_models)
reef_models$reef_proportion_vuln <- (reef_models$reef_proportion_vuln*(nrow(reef_models) - 1) + 0.5) / nrow(reef_models)


#################
# PC 1 CENTROID #
#################

reef_PC1_model <- stan_glmer(log(FIde_PC1+3) ~ (1 | geographic/Site), data=reef_models,
                             prior_intercept = normal(0,10),
                       family=gaussian, chains=4, iter=2000)
reef_PC1 <- as.matrix(reef_PC1_model)
reef_PC1 <- exp(reef_PC1[,18:27] + reef_PC1[,1]) - 3
colnames(reef_PC1) <- geo_coords$geographic
mcmc_areas(reef_PC1)
PC1_centroid_reefs <- apply(reef_PC1, 2, median)
plot(reef_raw_avg$FIde_PC1, PC1_centroid_reefs, pch=19)
y <- log(reef_models$FIde_PC1+3)
yrep <- posterior_predict(reef_PC1_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#################
# PC 2 CENTROID #
#################

reef_PC2_model <- stan_glmer(FIde_PC2 ~ (1 | geographic/Site), data=reef_models,
                             family=gaussian, chains=4, iter=2000)
reef_PC2 <- as.matrix(reef_PC2_model)
reef_PC2 <- reef_PC2[,18:27] + reef_PC2[,1]
colnames(reef_PC2) <- geo_coords$geographic
mcmc_areas(reef_PC2)
PC2_centroid_reefs <- apply(reef_PC2, 2, median)
plot(reef_raw_avg$FIde_PC2, PC2_centroid_reefs,pch=19)
ppc_dens_overlay(reef_models$FIde_PC2,  posterior_predict(reef_PC2_model, draws=100)[1:100,])
y <- reef_models$FIde_PC2
yrep <- posterior_predict(reef_PC2_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

##################
# entropy DIVERSITY #
##################

reef_entropy_model <- stan_glmer((FD_q1+1) ~ (1 | geographic/Site), data=reef_models,
                             family=Gamma(link="log"),chains=4, iter=2000)
reef_entropy <- as.matrix(reef_entropy_model)
reef_entropy <- (exp((reef_entropy[,18:27]) + reef_entropy[,1])) - 1
colnames(reef_entropy) <- geo_coords$geographic
mcmc_areas(reef_entropy)
reef_entropy <- apply(reef_entropy, 2, median)
plot(reef_raw_avg$FD_q1, reef_entropy,pch=19,
     xlim=range(c(reef_raw_avg$FD_q1, reef_entropy)),
     ylim=range(c(reef_raw_avg$FD_q1, reef_entropy)))
abline(0,1)
y <- reef_models$FD_q1+1
yrep <- posterior_predict(reef_entropy_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

############
# EVENNESS #
############

reef_feve_model <- brm(FEve ~ (1 | geographic/Site), data=reef_models,
                       family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)
reef_feve <- as.matrix(reef_feve_model)
reef_feve <- invlogit((reef_feve[,5:14]) + reef_feve[,1])
colnames(reef_feve) <- geo_coords$geographic
mcmc_areas(reef_feve)
reef_feve <- apply(reef_feve, 2, median)
plot(reef_raw_avg$FEve, reef_feve,pch=19)
model_data <- na.omit(reef_models[,c("FEve","geographic","Site")])
y <- model_data$FEve
yrep <- posterior_predict(reef_feve_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#######################
# PROPORTION DISTINCT #
#######################

reef_di_model <- brm(reef_proportion_di ~ (1 | geographic/Site), data=reef_models,
                       family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)

reef_di <- as.matrix(reef_di_model)
reef_di <- invlogit((reef_di[,5:14]) + reef_di[,1])
colnames(reef_di) <- geo_coords$geographic
mcmc_areas(reef_di)
reef_di <- apply(reef_di, 2, median)
plot(reef_raw_avg$reef_proportion_di, reef_di, pch=19)
y <- reef_models$reef_proportion_di
yrep <- posterior_predict(reef_di_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#########################
# PROPORTION VULNERABLE #
#########################

reef_vuln_model <- brm(reef_proportion_vuln ~ (1 | geographic/Site), data=reef_models,
                     family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)

reef_vuln <- as.matrix(reef_vuln_model)
reef_vuln <- invlogit((reef_vuln[,5:14]) + reef_vuln[,1])
colnames(reef_vuln) <- geo_coords$geographic
mcmc_areas(reef_vuln)
reef_vuln <- apply(reef_vuln, 2, median)
plot(reef_raw_avg$reef_proportion_vuln, reef_vuln, pch=19)
y <- reef_models$reef_proportion_vuln
yrep <- posterior_predict(reef_vuln_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

################################
# PROPORTION OF TROPHIC GROUPS #
################################

reef_troph_model <-  brm(mvbind(grazer, invertivore, microphage,
                                omnivore, piscivore, planktivore) ~ (1  | geographic/Site), data=reef_models,
                         family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000)
reef_troph <- as.matrix(reef_troph_model)
reef_grazer <- invlogit((reef_troph[,25:34]) + reef_troph[,1])
reef_invertivore <- invlogit((reef_troph[,51:60]) + reef_troph[,2])
reef_microphage <- invlogit((reef_troph[,77:86]) + reef_troph[,3])
reef_omnivore <- invlogit((reef_troph[,103:112]) + reef_troph[,4])
reef_piscivore <- invlogit((reef_troph[,129:138]) + reef_troph[,5])
reef_planktivore <- invlogit((reef_troph[,155:164]) + reef_troph[,6])

ppc_dens_overlay(reef_models$grazer,  posterior_predict(reef_troph_model)[1:100,,1])
ppc_dens_overlay(reef_models$invertivore,  posterior_predict(reef_troph_model)[1:100,,2])
ppc_dens_overlay(reef_models$microphage,  posterior_predict(reef_troph_model)[1:100,,3])
ppc_dens_overlay(reef_models$omnivore,  posterior_predict(reef_troph_model)[1:100,,4])
ppc_dens_overlay(reef_models$piscivore,  posterior_predict(reef_troph_model)[1:100,,5])
ppc_dens_overlay(reef_models$planktivore,  posterior_predict(reef_troph_model)[1:100,,6])

reef_troph <- cbind(apply(reef_grazer,2,median), apply(reef_microphage,2,median), 
                    apply(reef_planktivore,2,median),  apply(reef_omnivore,2,median),
                    apply(reef_invertivore,2,median),apply(reef_piscivore,2,median))
colnames(reef_troph) <- c("1_grazer", "2_microphage", "3_planktivore",
                                                      "4_omnivore", "5_invertivore", "6_piscivore")
rowSums(reef_troph)


##################
## FOR LANDINGS ##
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

landings_PC1_model <- stan_glmer(FIde_PC1 ~ (1 | geographic) + (1 | Fisher.Name),
                                  data=landings_models,
                             family=gaussian, chains=4, iter=2000)
landings_PC1 <- as.matrix(landings_PC1_model)
# GEOGRAPHIC INTERCEPTS
landings_PC1 <- landings_PC1[,117:126] + landings_PC1[,1]
colnames(landings_PC1) <- geo_coords$geographic
mcmc_areas(landings_PC1)
PC1_centroid_landings <- apply(landings_PC1, 2, median)
plot(landings_raw_avg$FIde_PC1, PC1_centroid_landings,pch=19)
plot(landings_raw_avg$FIde_PC1, PC1_centroid_landings,pch=19,
     xlim=range(c(landings_raw_avg$FIde_PC1, PC1_centroid_landings)),
     ylim=range(c(landings_raw_avg$FIde_PC1, PC1_centroid_landings)))
abline(0,1)
model_data <- na.omit(landings_models[,c("FIde_PC1","geographic","Fisher.Name")])
y <- model_data$FIde_PC1
yrep <- posterior_predict(landings_PC1_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#################
# PC 2 CENTROID #
#################

landings_PC2_model <- stan_glmer(FIde_PC2 ~ (1 | geographic) + (1 | Fisher.Name),
                                data=landings_models,
                             family=gaussian, chains=4, iter=2000)
landings_PC2 <- as.matrix(landings_PC2_model)
# GEOGRAPHIC INTERCEPTS
landings_PC2 <- landings_PC2[,117:126] + landings_PC2[,1]
colnames(landings_PC2) <- geo_coords$geographic
mcmc_areas(landings_PC2)
PC2_centroid_landings <- apply(landings_PC2, 2, median)
plot(landings_raw_avg$FIde_PC2, PC2_centroid_landings,pch=19)
model_data <- na.omit(landings_models[,c("FIde_PC2","geographic","Fisher.Name")])
y <- model_data$FIde_PC2
yrep <- posterior_predict(landings_PC2_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

##################
# entropy DIVERSITY #
##################

landings_entropy_model <- stan_glmer((FD_q1+1) ~ (1 | geographic) + (1 | Fisher.Name),
                                   data=landings_models,
                                 family=Gamma(link="log"),chains=4, iter=2000 )
landings_entropy <- as.matrix(landings_entropy_model)
landings_entropy <- (exp((landings_entropy[,117:126]) + landings_entropy[,1])) - 1
colnames(landings_entropy) <- geo_coords$geographic
mcmc_areas(landings_entropy)
landings_entropy <- apply(landings_entropy, 2, median)
plot(landings_raw_avg$FD_q1, landings_entropy,pch=19)
model_data <- na.omit(landings_models[,c("FD_q1","geographic","Fisher.Name")])
y <- model_data$FD_q1+1
yrep <- posterior_predict(landings_entropy_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

############
# EVENNESS #
############

landings_feve_model <- brm(FEve ~ (1 | geographic) + (1 | Fisher.Name),
                             data=landings_models,
                           family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )
landings_feve <- as.matrix(landings_feve_model)
landings_feve <- invlogit((landings_feve[,112:121]) + landings_feve[,1])
colnames(landings_feve) <- geo_coords$geographic
mcmc_areas(landings_feve)
landings_feve <- apply(landings_feve, 2, median)
plot(landings_raw_avg$FEve, landings_feve, pch=19)
model_data <- na.omit(landings_models[,c("FEve","geographic","Fisher.Name")])
y <- model_data$FEve
yrep <- posterior_predict(landings_feve_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#################
# PROPORTION DI #
#################

landings_di_model <- brm(market_proportion_di~ (1 | geographic) + (1 | Fisher.Name), 
                         data=landings_models,
                           family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )
landings_di <- as.matrix(landings_di_model)
landings_di <- invlogit((landings_di[,120:129]) + landings_di[,1])
colnames(landings_di) <- geo_coords$geographic
mcmc_areas(landings_di)
landings_di <- apply(landings_di, 2, median)
plot(landings_raw_avg$market_proportion_di, landings_di, pch=19)
model_data <- na.omit(landings_models[,c("market_proportion_di","geographic","Fisher.Name")])
y <- model_data$market_proportion_di
yrep <- posterior_predict(landings_di_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#########################
# PROPORTION VULNERABLE #
#########################

landings_vuln_model <- brm(market_proportion_vuln ~ (1 | geographic) + (1 | Fisher.Name),
                           data=landings_models,
                         family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )
landings_vuln <- as.matrix(landings_vuln_model)
landings_vuln <- invlogit((landings_vuln[,120:129]) + landings_vuln[,1])
colnames(landings_vuln) <- geo_coords$geographic
mcmc_areas(landings_vuln)
landings_vuln <- apply(landings_vuln, 2, median)
plot(landings_raw_avg$market_proportion_vuln, landings_vuln, pch=19)
model_data <- na.omit(landings_models[,c("market_proportion_vuln","geographic","Fisher.Name")])
y <- model_data$market_proportion_vuln
yrep <- posterior_predict(landings_vuln_model, draws=100)
ppc_dens_overlay(y, yrep)
ppc_stat(y, yrep, stat="mean")

#############################
# PROPORTION TROPHIC GROUPS #
#############################

landings_troph_model  <- brm(mvbind(grazer, invertivore, microphage,
                                  omnivore, piscivore, planktivore) ~ (1 | geographic) + (1 | Fisher.Name),
                           data=landings_models,
                           family=Beta(link = "logit", link_phi = "log"),chains=4, iter=2000 )
landings_troph <- as.matrix(landings_troph_model)
landings_grazer <- invlogit((landings_troph[,140:149]) + landings_troph[,1])
landings_invertivore <- invlogit((landings_troph[,265:274]) + landings_troph[,2])
landings_microphage <- invlogit((landings_troph[,390:399]) + landings_troph[,3])
landings_omnivore <- invlogit((landings_troph[,515:524]) + landings_troph[,4])
landings_piscivore <- invlogit((landings_troph[,640:649]) + landings_troph[,5])
landings_planktivore <- invlogit((landings_troph[,765:774]) + landings_troph[,6])

model_data <- na.omit(landings_models[,c("grazer","invertivore","microphage",
                                 "omnivore","piscivore","planktivore",
                                 "geographic","Fisher.Name")])
ppc_dens_overlay(model_data$grazer,  posterior_predict(landings_troph_model)[1:100,,1])
ppc_dens_overlay(model_data$invertivore,  posterior_predict(landings_troph_model)[1:100,,2])
ppc_dens_overlay(model_data$microphage,  posterior_predict(landings_troph_model)[1:100,,3])
ppc_dens_overlay(model_data$omnivore,  posterior_predict(landings_troph_model)[1:100,,4])
ppc_dens_overlay(model_data$piscivore,  posterior_predict(landings_troph_model)[1:100,,5])
ppc_dens_overlay(model_data$planktivore,  posterior_predict(landings_troph_model)[1:100,,6])

landings_troph <- cbind(apply(landings_grazer,2,median), apply(landings_microphage,2,median),
                        apply(landings_planktivore,2,median), apply(landings_omnivore,2,median), 
                        apply(landings_invertivore,2,median), apply(landings_piscivore,2,median))
colnames(landings_troph) <- c("1_grazer", "2_microphage", "3_planktivore",
                               "4_omnivore", "5_invertivore", "6_piscivore")
rowSums(landings_troph)

troph_test <- as.data.frame(landings_troph - reef_troph)
troph_test <- stack(troph_test)
troph_test$geographic <- rep(landings_raw_avg$geographic)

#################################
## SCALE TO 1 FOR STACKED BARS ##
#################################

reef_troph_scaled <- make_relative(reef_troph)
rowSums(reef_troph_scaled)

landings_troph_scaled <- make_relative(landings_troph)
rowSums(landings_troph_scaled)

############################################
## DO STACKED BAR CHARTS PER ZONE BY EACH ##
############################################

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


#######################################
## TRAIT SPACE FIGURE WITH CENTROIDS ##
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

diet_cols <- kovesi.rainbow(length(unique(species_info$Diet)))

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

plot(PC1, PC2, cex=0, cex.lab=1.2, cex.axis=1.2)

points(PC1_centroid_reefs, PC2_centroid_reefs,
       pch=21, col=1, bg="blue", cex=2)

points(PC1_centroid_landings, PC2_centroid_landings,
       pch=21, col=1, bg="red", cex=2)

legend("topright", legend=c("Reef Observations", "Market Landings"),
       pch=19, col=c("blue","red"), pt.cex=1.5)

mtext("D", font=2, cex=1.5, adj=-0.1, line=0.5)




##############################
## MAP FUNCTIONAL DIVERSITY ##
##############################

########
# entropy #
########

graphics.off()
par(mfrow=c(1,2))

## REEF 
map_color <- variablecol(reef_entropy, col=brewer.oranges(length(reef_entropy)),
                         clim=range(c(reef_entropy, landings_entropy)))
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = reef_entropy,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.oranges(length(reef_entropy)),
          clim=range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Rao's Quadratic Entropy")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = reef_entropy,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(reef_entropy/max(c(reef_entropy,landings_entropy)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_entropy, col=brewer.oranges(length(landings_entropy)),
                         clim=range(c(reef_entropy, landings_entropy)))
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = landings_entropy,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.oranges(length(landings_entropy)),
          clim=range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Rao's Quadratic Entropy")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = landings_entropy,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_entropy/max(c(reef_entropy,landings_entropy)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Fisheries Landings")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)

## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(landings_entropy),
                      diff = reef_entropy - landings_entropy)
di_diff$sign <- ifelse(di_diff$diff < 0, "negative","positive")

ggplot(di_diff,aes(geographic,diff)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green"))


########
# FEVE #
########

graphics.off()
par(mfrow=c(1,2))


## REEF 
map_color <- variablecol(reef_feve, col=brewer.blues(length(reef_feve)),
                         clim=range(c(reef_feve, landings_feve)))
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = reef_feve,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.blues(length(reef_feve)),
          clim=range(c(reef_feve, landings_feve)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Functional Evenness")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = reef_feve,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(reef_feve/max(c(reef_feve,landings_feve)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_feve, col=brewer.blues(length(landings_feve)),
                         clim=range(c(reef_feve, landings_feve)))
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = landings_feve,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.blues(length(landings_feve)),
          clim=range(c(reef_feve, landings_feve)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Functional Evenness")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = landings_feve,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_feve/max(c(reef_feve,landings_feve)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Fisheries Landings")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(landings_feve),
                      diff = reef_feve - landings_feve)
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
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = reef_di,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.purples(length(reef_di)),
          clim=range(c(reef_di, landings_di)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1.5,"Proportion Distinct Species")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = reef_di,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(reef_di/max(c(reef_di,landings_di)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_di, col=brewer.purples(length(landings_di)),
                         clim=range(c(reef_di, landings_di)))
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = landings_di,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.purples(length(landings_di)),
          clim=range(c(reef_di, landings_di)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Proportion Distinct Species")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = landings_di,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_di/max(c(reef_di,landings_di)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Fisheries Landings")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


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
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = reef_vuln,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.reds(length(reef_vuln)),
          clim=range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
#mtext(side=4,line=1,"Proportion Vulnerable Species")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = reef_vuln,
          add=TRUE,colkey = FALSE, xlab.5=NA,ylab=NA,
          cex=(reef_vuln/max(c(reef_vuln,landings_vuln)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Reef Observations")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- variablecol(landings_vuln, col=brewer.reds(length(landings_vuln)),
                         clim=range(c(reef_vuln, landings_vuln)))
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=19, colvar = landings_vuln,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.255,5.375),
          col = brewer.reds(length(landings_vuln)),
          clim=range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
#mtext(side=4,line=1.5,"Proportion Vulnerable Species")
plot(FSM, xlim=range(reef_data$X),
     ylim=range(reef_data$Y), col="grey",add=TRUE)
scatter2D(geo_coords$X, geo_coords$Y, 
          pch=21, colvar = landings_vuln,
          add=TRUE,colkey = FALSE, xlab=NA,ylab=NA,
          cex=(landings_vuln/max(c(reef_vuln,landings_vuln)))*4,
          #cex=4,
          col = 1, bg=map_color)
title("Fisheries Landings")
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)

## DIFFERNECE BARPLOT
di_diff <- data.frame(geographic = names(landings_vuln),
                      diff = reef_vuln - landings_vuln)
di_diff$sign <- ifelse(di_diff$diff < 0, "negative","positive")

ggplot(di_diff,aes(geographic,diff)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green"))



###########################################################
## MODELS FOR DRIVERS OF REEF AND MARKET TRAIT DIVERSITY ##
###########################################################

benthic_site <- read.table("benthic_pcoa_site.txt")

benthic_geo <- read.table("benthic_pcoa_geo.txt")

######################
## IMPORT MSEC DATA ##
######################

#msec <- read.csv("msec_out.csv")
msec <- read.table("reef_drivers.txt")

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

reef_drivers <- merge(benthic_site, msec, by="site")
names(reef_drivers)[1] <- "Site"
reef_models <- merge(reef_models, reef_drivers, by="Site")
reef_models <- merge(reef_models, events_per_section, by="geographic")

################
# PC1 CENTROID #
################

reef_PC1_drivers_model <- stan_glmer(log(FIde_PC1+3) ~ scale(Axis.1) + scale(Axis.2) + scale(reef_area_5km) + 
                                       scale(fishing_events) + scale(log(total_gravity_pop)) + 
                                       (1 | geographic/Site), data=reef_models,
                                     family=gaussian, chains=4, iter=2000)
reef_PC1_drivers <- as.matrix(reef_PC1_drivers_model)
reef_PC1_drivers <- (reef_PC1_drivers[,2:6]) 
color_scheme_set("darkgray")
mcmc_intervals(reef_PC1_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)
model_data <- reef_models[,c("FIde_PC1","Axis.1","Axis.2","reef_area_5km","fishing_events",
                             "total_gravity_pop","geographic","Site")]

ppc_dens_overlay(log(reef_PC1_drivers_model$data$FIde_PC1+3), posterior_predict(reef_PC1_drivers_model,draws=100)[1:100,])
ppc_error_scatter_avg(log(reef_PC1_drivers_model$data$FIde_PC1+3), posterior_predict(reef_PC1_drivers_model,draws=100))

y <- log(reef_PC1_drivers_model$data$FIde_PC1+3)
yrep <- posterior_predict(reef_PC1_drivers_model,draws=100)
plot(colMeans(yrep), y)
summary(bayes_R2(reef_PC1_drivers_model))
r2_bayes(reef_PC1_drivers_model)
model_performance(reef_PC1_drivers_model)

##################
# entropy DIVERSITY #
##################

reef_entropy_drivers_model <- stan_glmer((FD_q1 +1)  ~ scale(Axis.1) + scale(Axis.2) + scale(reef_area_5km) + 
                                       scale(fishing_events) + scale(log(total_gravity_pop)) + 
                                       (1 | geographic/Site), data=reef_models,
                             family=Gamma(link="log"),chains=4, iter=2000)
reef_entropy_drivers <- as.matrix(reef_entropy_drivers_model)
reef_entropy_drivers <- (reef_entropy_drivers[,2:6]) 
mcmc_intervals(reef_entropy_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)
ppc_dens_overlay(reef_entropy_drivers_model$data$FD_q1+1, posterior_predict(reef_entropy_drivers_model,draws=100)[1:100,])

############
# EVENNESS #
############

reef_feve_drivers_model <- brm(FEve  ~ scale(Axis.1) + scale(Axis.2) + scale(reef_area_5km) + 
                                 scale(fishing_events) + scale(log(total_gravity_pop)) + 
                                 (1 | geographic/Site), data=reef_models,
                                      family=Beta(link = "logit", link_phi = "log"),chains=4, iter=3000,
                               control = list(adaptdelta=0.9))
reef_feve_drivers <- as.matrix(reef_feve_drivers_model)
reef_feve_drivers <- (reef_feve_drivers[,2:6]) 
mcmc_intervals(reef_feve_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)
ppc_dens_overlay(reef_feve_drivers_model$data$FEve, posterior_predict(reef_feve_drivers_model,draws=100)[1:100,])
ppc_error_scatter_avg(reef_feve_drivers_model$data$FEve, posterior_predict(reef_feve_drivers_model,draws=100))

y <- reef_feve_drivers_model$data$FEve
yrep <- posterior_predict(reef_feve_drivers_model,draws=100)
plot(colMeans(yrep), y)
r2_bayes(reef_feve_drivers_model)

#####################
## LANDINGS MODELs ##
#####################

# GRAVITY AS A GEOGRAPHIC MIDPOINT**

test <- aggregate(reef_models[,c(34:ncol(reef_models))], by=list(reef_models$geographic), FUN=mean,na.rm=TRUE)
colnames(test)[1] <- "geographic"

reef_intercepts <- data.frame(PC1_centroid_reefs, reef_entropy, reef_feve)
reef_intercepts$geographic <- rownames(reef_intercepts)

reef_intercepts <- merge(reef_intercepts, test, by="geographic")

landings_drivers <- merge(landings_models, reef_intercepts, by="geographic")

landings_drivers$Axis.1 <- with(benthic_geo,
                                        Axis.1[match(landings_drivers$geographic, geographic)])
landings_drivers$Axis.2 <- with(benthic_geo,
                                Axis.2[match(landings_drivers$geographic, geographic)])

# landings_drivers$fishing_events <- with(events_per_section,
#                                         fishing_events[match(landings_drivers$geographic, geographic)])
# landings_drivers$gravity <- with(test,
#                                         gravity[match(landings_drivers$geographic, geographic)])
# landings_drivers$wave_energy <- with(test,
#                                  wave_energy[match(landings_drivers$geographic, geographic)])

################
# PC1 CENTROID #
################

landings_PC1_drivers_model <- stan_glmer(FIde_PC1 ~ scale(PC1_centroid_reefs) + scale(fishing_events) + 
                                           scale(wave_energy) + scale(dist_market) +
                                           as.factor(moon_phase) + scale(wind) +
                                           as.factor(gear) + scale(num_fishers_lines) + 
                                           (1 | geographic) + (1 | Fisher.Name),
                                       family=gaussian, data=landings_drivers, chains=4, iter=2000)
landings_PC1_drivers <- as.matrix(landings_PC1_drivers_model)
landings_PC1_drivers <- (landings_PC1_drivers[,2:10]) 
mcmc_intervals(landings_PC1_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)
model_data <- landings_drivers[,c("FIde_PC1","PC1_centroid_reefs","PC1_centroid_reefs",
                                  "wave_energy","dist_market","moon_phase","wind",
                                  "gear","num_fishers_lines","geographic","Fisher.Name")]
model_data <- na.omit(model_data)
y <- model_data$FIde_PC1
yrep <- posterior_predict(landings_PC1_drivers_model,draws=100)
ppc_dens_overlay(y, yrep[1:100,])
ppc_stat(y, yrep, stat="mean")
plot(colMeans(yrep), y)
r2_bayes(landings_PC1_drivers_model)

##################
# entropy DIVERSITY #
##################

landings_entropy_drivers_model <- stan_glmer((FD_q1+1) ~ scale(reef_entropy+1) + scale(fishing_events) + 
                                            scale(wave_energy) + scale(dist_market) +
                                            as.factor(moon_phase) + scale(wind) +
                                            as.factor(gear) + scale(num_fishers_lines) + 
                                            (1 | geographic) + (1 | Fisher.Name),
                                       family=Gamma(link="log"), data=landings_drivers, chains=4, iter=2000)
landings_entropy_drivers <- as.matrix(landings_entropy_drivers_model)
landings_entropy_drivers <- (landings_entropy_drivers[,2:10])
mcmc_intervals(landings_entropy_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)
model_data <- landings_drivers[,c("FD_q1","reef_entropy","fishing_events",
                                  "wave_energy","dist_market","moon_phase","wind",
                                  "gear","num_fishers_lines","geographic","Fisher.Name")]
model_data <- na.omit(model_data)
y <- model_data$FIde_PC1
yrep <- posterior_predict(landings_PC1_drivers_model,draws=100)
ppc_dens_overlay(y, yrep[1:100,])
ppc_stat(y, yrep, stat="mean")
plot(colMeans(yrep), y)
r2_bayes(landings_PC1_drivers_model)

############
# EVENNESS #
############

landings_feve_drivers_model <- brm(FEve ~ scale(reef_feve) + scale(fishing_events) + 
                                     scale(wave_energy) + scale(dist_market) +
                                     as.factor(moon_phase) + scale(wind) +
                                     as.factor(gear) + scale(num_fishers_lines) + 
                                     (1 | geographic) + (1 | Fisher.Name),
                                        family=Beta(link = "logit", link_phi = "log"), data=landings_drivers, chains=4, iter=2000)
landings_feve_drivers <- as.matrix(landings_feve_drivers_model)
landings_feve_drivers <- (landings_feve_drivers[,2:10]) 
mcmc_intervals(landings_feve_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

r2_bayes(landings_feve_drivers_model)

##########################################
## INTERACTIONS AMONG DIVERSITY METRICS ##
##########################################

reef_intx_data <- data.frame(PC1=reef_FD$FIde_PC1, rao=reef_FD$reef_rao, entropy=reef_FD$FD_q1,
                             FEve = reef_FD$FEve, pref=reef_FD$reef_deviation,
                             geographic=reef_meta$geographic,
                             source=rep("reef"))
market_intx_data <- data.frame(PC1 = market_FD$FIde_PC1, rao=market_FD$market_rao,entropy=market_FD$FD_q1,
                               FEve =market_FD$FEve,pref=market_FD$market_deviation,
                               geographic=market_meta$geographic,
                               source=rep("market"))

intx_data <- rbind(reef_intx_data, market_intx_data)

# entropy & PC1
entropy_PC1_model <- brm((entropy+1) ~ PC1:source + (1 | geographic),
                      family=Gamma(link="log"), data=intx_data)
marginal_effects(entropy_PC1_model)
mean(PC1_centroid_reefs)
mean(PC1_centroid_landings)
mean(reef_entropy)
mean(landings_entropy)

# FEVE & PC1
FEve_PC1_model <- brm((FEve) ~ PC1:source + (1 | geographic),
                      family=Beta(link = "logit", link_phi = "log"), data=intx_data)
marginal_effects(FEve_PC1_model)

# entropy $ FEVE
entropy_feve_model <- brm((entropy+1) ~ FEve:source + (1 | geographic),
                       family=Gamma(link="log"), data=intx_data)
marginal_effects(entropy_feve_model)


########################################################
## RELATION OF MEAN DESIRABILITY TO DIVERSITY METRICS ##
########################################################

pref_model <- stan_glmer(market_deviation ~ scale(FIde_PC1) + scale(FEve) + scale(FD_q1) +
                                           (1 | geographic) + (1 | Fisher.Name),
                                         family=gaussian, data=landings_drivers, chains=4, iter=2000)
pref_drivers <- as.matrix(pref_model)
pref_drivers <- (pref_drivers[,2:4]) 
mcmc_intervals(pref_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)


############################################
# CPUE HIGHER OR LOWER FOR DESIRABLE FISH? #
############################################

CPUE <- brm(log(CPUE) ~ scale(market_deviation) +
              (1 | geographic) + (1 | Fisher.Name),
            family=gaussian, data=landings_drivers, chains=4, iter=2000)

marginal_effects(CPUE)
CPUE_drivers <- as.matrix(CPUE)
CPUE_drivers <- (CPUE_drivers[,1:2]) 
mcmc_intervals(CPUE_drivers, point_est = "median", prob = 0.5, prob_outer = 0.95,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)







