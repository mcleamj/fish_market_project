

#####################################################
#####################################################
## SENSITIVITY TEST FOR REMOVING NOCTURANL SPECIES ##
#####################################################
#####################################################

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

# DEFINE NOCTURNAL SPECIES

noctural <- c("Myripristis sp", "Myripristis adusta", 
              "Sargocentron tiere", "Pempheris oualensis")

###################################################
## IMPORT SPECIES TRAITS AND TROPHIC INFORMATION ##
###################################################

traits <- read.table("data/clean traits 2025.txt") %>%
  tibble::rownames_to_column("species") %>%
  filter(! species %in% noctural ) %>%
  arrange(species) %>%
  column_to_rownames("species")

species_info <- read.table("data/clean species info 2025.txt") %>%
  filter(! species %in% noctural ) %>%
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
  select(-all_of(noctural))

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
  select(-all_of(noctural))

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
saveRDS(reef_metrics, "reef_metrics_no_nocturnal.rds")

landings_metrics <- data.frame(PC1_centroid_landings, PC2_centroid_landings, landings_entropy, landings_di, landings_vuln)
saveRDS(landings_metrics, "landings_metrics_no_nocturnal.rds")

######################################
## CORRELATION WITH ORIGINAL VALUES ##
######################################

reef_metrics_no_nocturnal <- readRDS("reef_metrics_no_nocturnal.rds")
reef_metrics_original <- readRDS("reef_metrics_original.rds")

landings_metrics_no_nocturnal <- readRDS("landings_metrics_no_nocturnal.rds")
landings_metrics_original <- readRDS("landings_metrics_original.rds")

all(rownames(reef_metrics_original) == rownames(reef_metrics_no_nocturnal))
all(rownames(landings_metrics_no_nocturnal) == rownames(landings_metrics_original))

# REEFS
reef_vars <- c("PC1_centroid_reefs",
          "PC2_centroid_reefs",
          "reef_entropy",
          "reef_di",
          "reef_vuln")

reef_cors <- sapply(reef_vars, function(v) {
  cor(
    reef_metrics_original[[v]],
    reef_metrics_no_nocturnal[[v]],
    use = "complete.obs"
  )
})

reef_cors

range(reef_metrics_original$reef_di)
range(reef_metrics_no_nocturnal$reef_di)

range(reef_metrics_no_nocturnal$reef_vuln)
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
    landings_metrics_no_nocturnal[[v]],
    use = "complete.obs"
  )
})

landings_cors

range(landings_metrics_original$landings_di)
range(landings_metrics_no_nocturnal$landings_di)

range(landings_metrics_no_nocturnal$landings_vuln)
range(landings_metrics_original$landings_vuln)

########################
# CORRELATION PLOTS  ##
########################
plot_data <- rbind(
  landings_metrics_original %>%
    rownames_to_column("geographic") %>%
    mutate(source = "original"),
  landings_metrics_no_nocturnal %>%
    rownames_to_column("geographic") %>%
    mutate(source = "nocturnal")
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
    lo = min(c(original, nocturnal), na.rm = TRUE),
    hi = max(c(original, nocturnal), na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_long, aes(original, nocturnal)) +
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
    y = "Nocturnal Species Removed"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )


##########################
##  DIFFERENCE BAR PLOTS
##########################

diff_metrics <- data.frame(reef_metrics_no_nocturnal - landings_metrics_no_nocturnal) %>%
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


