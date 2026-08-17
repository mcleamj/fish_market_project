


########################################################################
########################################################################
## SENSITIVITY TEST FOR REMOVING ALL NOCTURNAL AND OR non_reef SPECIES ##
########################################################################
########################################################################

##################################
## INSTALL AND LIBRARY PACKAGES ##
##################################

if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(ggh4x)){install.packages("ggh4x"); library(ggh4x)}
if(!require(FD)){install.packages("FD"); library(FD)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(funrar)){install.packages("funrar"); library(funrar)}
if(!require(rfishbase)){install.packages("rfishbase"); library(rfishbase)}
if(!require(mFD)){install.packages("mFD"); library(mFD)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}


# DEFINE non_reefS

non_reef <- c(
  "Elagatis bipinnulata", "Caesio caerulaurea", "Selar crumenophthalmus",
  "Ferdauia ferdau", "Caranx lugubris", "Ferdauia orthogrammus", "Caranx melampygus",
  "Aphareus furca",
  "Myripristis sp", "Myripristis adusta", "Sargocentron tiere", "Pempheris oualensis")

###################################################
## IMPORT SPECIES TRAITS AND TROPHIC INFORMATION ##
###################################################

traits <- read.table("data/clean traits 2025.txt") %>%
  tibble::rownames_to_column("species") %>%
  filter(! species %in% non_reef ) %>%
  arrange(species) %>%
  column_to_rownames("species")

species_info <- read.table("data/clean species info 2025.txt") %>%
  filter(! species %in% non_reef ) %>%
  arrange(species)

diet_cols <- kovesi.rainbow(length(unique(species_info$Diet)))

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

colnames(reef_fish) <- gsub("\\.", " ", colnames(reef_fish))

reef_fish <- reef_fish %>%
  select(-all_of(non_reef))

reef_log <- log10(reef_fish+1)

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

colnames(market_fish) <- gsub("\\.", " ", colnames(market_fish))

market_fish <- market_fish %>%
  select(-all_of(non_reef))

empty_catch <- which(rowSums(market_fish)==0)

market_fish <- market_fish[-empty_catch,]

market_meta <- market_meta[-empty_catch,]

market_log <- log10(market_fish+1)

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

vuln <- rfishbase::species(species_info$species, fields = c("Species", "Vulnerability"), version = "21.06")
vuln <- vuln %>% dplyr::rename(species = Species)
vuln <- merge(species_info[,c("species","family")], vuln, by="species",
              all=TRUE)
vuln$Vulnerability[vuln$species=="Ferdauia orthogrammus"] <- species("Carangoides orthogrammus", fields="Vulnerability", version = "21.06")$Vulnerability
vuln$Vulnerability[vuln$species=="Ferdauia ferdau"] <- species("Carangoides ferdau", fields="Vulnerability", version = "21.06")$Vulnerability
rownames(vuln) <- vuln$species
vuln$family <- NULL

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



# VULNERABLE
quars <- quantile(vuln$Vulnerability,probs=c(0.25,0.5,0.75))
vuln$quartile <- ifelse(vuln$Vulnerability <= quars[1],"Q1",
                        ifelse(vuln$Vulnerability > quars[1] & 
                                 vuln$Vulnerability <= quars[2],"Q2",
                               ifelse(vuln$Vulnerability > quars[2] & 
                                        vuln$Vulnerability <= quars[3],"Q3",
                                      ifelse(vuln$Vulnerability > quars[3],"Q4",NA))))


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
## INTERCEPT-ONLY HIERARCHICAL MODELS TO CALCULATE GEOGRAPHIC MEANS ##---------------------------------------------------------------------------
######################################################################

###############
## FOR REEFS ##----------------------------------------------------------------------------------------------------------------------------------
###############

reef_raw_avg <- aggregate(cbind(reef_FD,reef_proportion_di,reef_proportion_vuln), 
                          by=list(reef_meta$geographic),FUN=mean,na.rm=TRUE)
colnames(reef_raw_avg)[1] <- "geographic"

# TRANSFORMATION TO CONSTRAIN 0 AND 1 VALUES FOR BETA DISTRIBUTION MODEL
reef_models <- data.frame(reef_meta, reef_FD, reef_proportion_di, reef_proportion_vuln)
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


################## ##---------------------------------------------------------------------------------------------------------------------
##################
## FOR LANDINGS ##-------------------------------------------------------------------------------------------------------------------------
##################
################## ##----------------------------------------------------------------------------------------------------------------------

landings_raw_avg <- aggregate(cbind(market_FD,market_proportion_di,market_proportion_vuln),
                              by=list(market_meta$geographic),FUN=mean,na.rm=TRUE)
colnames(landings_raw_avg)[1] <- "geographic"

landings_models <- data.frame(market_meta, market_FD, market_proportion_di, market_proportion_vuln)
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

######################
## SAVE ALL METRICS ##
######################

reef_metrics <- data.frame(PC1_centroid_reefs, PC2_centroid_reefs, reef_entropy, reef_di, reef_vuln)
saveRDS(reef_metrics, "reef_metrics_no_non_reef.rds")

landings_metrics <- data.frame(PC1_centroid_landings, PC2_centroid_landings, landings_entropy, landings_di, landings_vuln)
saveRDS(landings_metrics, "landings_metrics_no_non_reef.rds")

######################################
## CORRELATION WITH ORIGINAL VALUES ##
######################################

reef_metrics_no_non_reef <- readRDS("reef_metrics_no_non_reef.rds")
reef_metrics_original <- readRDS("reef_metrics_original.rds")

landings_metrics_no_non_reef <- readRDS("landings_metrics_no_non_reef.rds")
landings_metrics_original <- readRDS("landings_metrics_original.rds")

all(rownames(reef_metrics_original) == rownames(reef_metrics_no_non_reef))
all(rownames(landings_metrics_no_non_reef) == rownames(landings_metrics_original))

# REEFS
reef_vars <- c("PC1_centroid_reefs",
               "PC2_centroid_reefs",
               "reef_entropy",
               "reef_di",
               "reef_vuln")

reef_cors <- sapply(reef_vars, function(v) {
  cor(
    reef_metrics_original[[v]],
    reef_metrics_no_non_reef[[v]],
    use = "complete.obs"
  )
})

reef_cors

range(reef_metrics_original$reef_di)
range(reef_metrics_no_non_reef$reef_di)

range(reef_metrics_no_non_reef$reef_vuln)
range(reef_metrics_original$reef_vuln)

# LANDINGS
landings_vars <- c("PC1_centroid_landings",
                   "PC2_centroid_landings",
                   "landings_entropy",
                   "landings_di",
                   "landings_vuln")

landings_cors <- sapply(landings_vars, function(v) {
  cor(
    landings_metrics_original[[v]],
    landings_metrics_no_non_reef[[v]],
    use = "complete.obs"
  )
})

landings_cors

range(landings_metrics_original$landings_di)
range(landings_metrics_no_non_reef$landings_di)

range(landings_metrics_original$landings_vuln)
range(landings_metrics_no_non_reef$landings_vuln)

########################
# CORRELATION PLOTS  ##
########################
plot_data <- rbind(
  landings_metrics_original %>%
    rownames_to_column("geographic") %>%
    mutate(source = "original"),
  landings_metrics_no_non_reef %>%
    rownames_to_column("geographic") %>%
    mutate(source = "non_reef")
)

df_long <- plot_data %>%
  pivot_longer(
    -c(geographic, source),
    names_to = "metric",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = value
  ) %>%
  filter(metric != "PC2_centroid_landings")

df_long$metric <- factor(
  df_long$metric,
  levels = c(
    "PC1_centroid_landings",
    "landings_entropy",
    "landings_di",
    "landings_vuln"
  ),
  labels = c(
    "Functional Identity (PC1)",
    "Functional Entropy",
    "Proportion Distinct",
    "Proportion Vulnerable"
  )
)

lims <- df_long %>%
  group_by(metric) %>%
  summarise(
    lo = min(c(original, non_reef), na.rm = TRUE),
    hi = max(c(original, non_reef), na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_long, aes(original, non_reef)) +
  geom_point(size = 3) +
  
  geom_text_repel(
    aes(label = geographic),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~metric, scales = "free") +
  
  facetted_pos_scales(
    x = purrr::map2(lims$lo, lims$hi,
                    ~scale_x_continuous(limits = c(.x, .y))),
    y = purrr::map2(lims$lo, lims$hi,
                    ~scale_y_continuous(limits = c(.x, .y)))
  )  +
  labs(
    x = "Original Data",
    y = "All Questionable Species Removed"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )




##########################
##  DIFFERENCE BAR PLOTS
##########################

diff_metrics <- data.frame(reef_metrics_no_non_reef - landings_metrics_no_non_reef) %>%
  select(-PC2_centroid_reefs) %>%
  rename("PC1 Centroid" = PC1_centroid_reefs,
         "Functional Entropy" = reef_entropy,
         "Proportional Functionally Distinct" = reef_di,
         "Proportion Fishing-Vulnerable" = reef_vuln) %>%
  rownames_to_column("geographic") %>%
  pivot_longer(
    cols = -geographic,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(sign =ifelse(value < 0, "negative","positive")) 

diff_metrics$metric <- factor(diff_metrics$metric,
                              levels = c("PC1 Centroid", "Functional Entropy",
                                         "Proportional Functionally Distinct",
                                         "Proportion Fishing-Vulnerable"))


ggplot(diff_metrics,aes(geographic,value)) +
  geom_bar(stat="identity",position="identity",aes(fill = sign))+
  scale_fill_manual(values=c(negative="red",positive="green")) +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 14)
    
  ) +
  labs(x=NULL, y="Difference Between Reef Observations and Landings")









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
     xlim=range(1*c(PC1_centroid_reefs, PC1_centroid_landings)),
     ylim=range(1*c(PC2_centroid_reefs, PC2_centroid_landings)))

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
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(PC1_centroid_reefs, col=brewer.blues(length(PC1_centroid_reefs)),
                         clim=range(c(PC1_centroid_reefs, PC1_centroid_landings))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = PC1_centroid_reefs,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.blues(length(PC1_centroid_reefs)),
          clim=range(c(PC1_centroid_reefs, PC1_centroid_landings)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"PC1 Centroid")
plot(map_poly, add=TRUE, col=map_poly$color)
title("Reef Observations")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(PC1_centroid_landings, col=brewer.blues(length(PC1_centroid_landings)),
                                            clim=range(c(PC1_centroid_reefs, PC1_centroid_landings))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = PC1_centroid_landings,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.blues(length(PC1_centroid_landings)),
          clim=range(c(PC1_centroid_reefs, PC1_centroid_landings)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"PC1 Centroid")
plot(map_poly, add=TRUE, col=map_poly$color)
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)
title("Fisheries Landings")


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
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(reef_entropy, col=brewer.oranges(length(reef_entropy)),
                                            clim=range(c(reef_entropy, landings_entropy))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = reef_entropy,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.oranges(length(reef_entropy)),
          clim=range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Functional Entropy")
plot(map_poly, add=TRUE, col=map_poly$color)
title("Reef Observations")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(landings_entropy, col=brewer.oranges(length(landings_entropy)),
                                            clim=range(c(reef_entropy, landings_entropy))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = landings_entropy,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.oranges(length(landings_entropy)),
          clim=range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Functional Entropy")
plot(map_poly, add=TRUE, col=map_poly$color)
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
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
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(reef_di, col=brewer.purples(length(reef_di)),
                                            clim=range(c(reef_di, landings_di))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = reef_di,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.purples(length(reef_di)),
          clim=range(c(reef_di, landings_di)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1.5,"Proportion Distinct Species")
plot(map_poly, add=TRUE, col=map_poly$color)
title("Reef Observations")
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)


## LANDINGS
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(landings_di, col=brewer.purples(length(landings_di)),
                                            clim=range(c(reef_di, landings_di))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = landings_di,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.purples(length(landings_di)),
          clim=range(c(reef_di, landings_di)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Proportion Distinct Species")
plot(map_poly, add=TRUE, col=map_poly$color)
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)
title("Fisheries Landings")


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
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(reef_vuln, col=brewer.reds(length(reef_vuln)),
                                            clim=range(c(reef_vuln, landings_vuln))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = reef_vuln,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.reds(length(reef_vuln)),
          clim=range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1,"Proportion Vulnerable Species")
plot(map_poly, add=TRUE, col=map_poly$color)
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
mtext("(a)",line=1,font=2,adj=-0.1,cex=1.75)
title("Reef Observations")


## LANDINGS
map_color <- data.frame(geographic = geo_coords$geographic,
                        color = variablecol(landings_vuln, col=brewer.reds(length(landings_vuln)),
                                            clim=range(c(landings_vuln, landings_vuln))))
map_poly <- kos_polys %>% left_join(map_color, by="geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch=19, colvar = landings_vuln,
          colkey = TRUE, cex=0,xlab=NA,ylab=NA,
          xlim=c(162.89,163.045),ylim=c(5.25,5.383),
          col = brewer.reds(length(landings_vuln)),
          clim=range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("lightblue",alpha=0.5))
mtext(side=4,line=1.5,"Proportion Vulnerable Species")
plot(map_poly, add=TRUE, col=map_poly$color)
plot(kos_shoreline, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey",add=TRUE)
plot(kos_reefs, xlim=range(reef_data$Lon),
     ylim=range(reef_data$Lat), col="grey90",add=TRUE)
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

# model_data <- reef_drivers[,c("fide_PC1","PC1_all","PC2_all","reef_area_5km","fishing_events",
#                               "tot_grav_pop","geographic","site","site_year")]
# 
# ppc_dens_overlay(log(reef_PC1_drivers_model$data$fide_PC1+3), posterior_predict(reef_PC1_drivers_model,draws=100)[1:100,])
# ppc_error_scatter_avg(log(reef_PC1_drivers_model$data$fide_PC1+3), posterior_predict(reef_PC1_drivers_model,draws=100))
# 
# y <- log(reef_PC1_drivers_model$data$fide_PC1+3)
# yrep <- posterior_predict(reef_PC1_drivers_model,draws=100)
# plot(colMeans(yrep), y)
# ppc_stat(y, yrep, stat="mean")

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

# ppc_dens_overlay(reef_entropy_drivers_model$data$FD_q1+1, posterior_predict(reef_entropy_drivers_model,draws=100)[1:100,])
# 
# y <- (reef_entropy_drivers_model$data$FD_q1 + 1)
# yrep <- posterior_predict(reef_entropy_drivers_model,draws=100)
# plot(colMeans(yrep), y)
# ppc_stat(y, yrep, stat="mean")

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
                                    #moon_phase + 
                                    #z_score_2sd(wind) +
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
         #'b_z_score_2sdwind',
         #'b_moon_phasemediummoon', 'b_moon_phasebigmoon',
         'b_gearShallowBottomFishing','spearguns', 'b_z_score_2sdnum_fishers_lines')
mcmc_intervals(landings_PC1_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

# model_data <- landings_drivers[,c("fide_PC1","PC1_centroid_reefs","fishing_events",
#                                   "wave_energy","tot_grav_pop","moon_phase","wind",
#                                   "gear","num_fishers_lines","geographic","Fisher.Name")]
# model_data <- na.omit(model_data)
# y <- model_data$fide_PC1
# yrep <- posterior_predict(landings_PC1_drivers_model,draws=100)
# ppc_dens_overlay(y, yrep[1:100,])
# ppc_stat(y, yrep, stat="mean")
# plot(colMeans(yrep), y)

#########################
# FUNCTIONAL ENTROPY  #
#########################

landings_entropy_drivers_model <- brm((FD_q1+1) ~ 
                                        z_score_2sd(reef_entropy+1) + 
                                        z_score_2sd(fishing_events) + 
                                        z_score_2sd(wave_energy) + 
                                        z_score_2sd(tot_grav_pop) +
                                        #moon_phase + 
                                        #z_score_2sd(wind) +
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
         #'b_z_score_2sdwind',
         #'b_moon_phasemediummoon', 'b_moon_phasebigmoon',
         'b_gearShallowBottomFishing','spearguns', 'b_z_score_2sdnum_fishers_lines')
mcmc_intervals(landings_entropy_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

# model_data <- landings_drivers[,c("FD_q1","reef_entropy","fishing_events",
#                                   "wave_energy","tot_grav_pop","moon_phase","wind",
#                                   "gear","num_fishers_lines","geographic","Fisher.Name")]
# model_data <- na.omit(model_data)
# y <- model_data$FD_q1+1
# yrep <- posterior_predict(landings_entropy_drivers_model,draws=100)
# ppc_dens_overlay(y, yrep[1:100,])
# ppc_stat(y, yrep, stat="mean")
# plot(colMeans(yrep), y)



