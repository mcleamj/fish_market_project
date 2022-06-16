

library(vegan)
library(ape)

asinTransform <- function(p) { asin(sqrt(p)) }

#########################
## BENTHIC ORDINATIONS ##
#########################

kos_sites <- read.csv("data/KOS_spatial_biomass.csv")

##################
## AT SPC SCALE ##
##################

#################
## CORALS ONLY ##
#################
# 
# benthic_coral <- read.csv("Benthic_ordination_corals.csv")
# benthic_coral$site[which(nchar(benthic_coral$site)<6)] <- gsub(pattern = "KOS-", replacement ="KOS-0", benthic_coral$site[which(nchar(benthic_coral$site)<6)])
# benthic_meta <- benthic_coral[,1:9]
# benthic_meta$geographic <- with(kos_sites,
#                                 geographic[match(benthic_meta$site,
#                                                  Site)])
# corals <- benthic_coral[,10:ncol(benthic_coral)]
# missing_row <- which(rowSums(corals)==0)
# corals <- corals[,-which(colSums(corals)==0)]
# corals <- corals[-missing_row,]
# benthic_coral <- benthic_coral[-missing_row,]
# benthic_meta <- benthic_meta[-missing_row,]
# corals <- log10(corals+1)
# bray_corals <- vegdist(corals, method = "bray")
# coral_ord <- pcoa(bray_corals)
# biplot(coral_ord)
# 
# plot(coral_ord$vectors[,1], coral_ord$vectors[,2], pch=19, col=as.factor(benthic_meta$geographic))
# 
# scatter2D(coral_ord$vectors[,1], coral_ord$vectors[,2], pch=19, cex=1, 
#           colvar=corals$Porites..POR.)
# 
# scatter2D(coral_ord$vectors[,1], coral_ord$vectors[,2], pch=19, cex=1, 
#           colvar=corals$Acropora..ACROP.)
# 
# scatter2D(coral_ord$vectors[,1], coral_ord$vectors[,2], pch=19, cex=1, 
#           colvar=corals$Porites.massive..PORMAS.)
# 
# scatter2D(coral_ord$vectors[,1], coral_ord$vectors[,2], pch=19, cex=1, 
#           colvar=corals$Acropora.table..ACROTBL.)
# 
# vectors <- envfit(coral_ord$vectors[,1:2], corals)
# 
# plot(vectors)
# 
# spc_coords <- data.frame(coral_ord$vectors[,1:2])
# geo_coords <- aggregate(spc_coords, by=list(benthic_meta$geographic), FUN = mean,na.rm=TRUE)     
# names(geo_coords)[1] <- "geographic"
# 
# 
# ############################
# ##  DIRECTLY AT GEO LEVEL ##
# ############################
# 
# coral_geo <- aggregate(corals, by=list(benthic_meta$geographic), FUN = mean,na.rm=TRUE)
# coral_geo$Group.1 <- NULL
# 
# coral_geo <- log10(coral_geo+1)
# 
# geo_bray <- vegdist(coral_geo, method = "bray")
# 
# geo_ord <- pcoa(geo_bray)
# 
# biplot(geo_ord)
# 
# geo_coords_direct <- data.frame(geographic = geo_coords$geographic,
#                                   geo_ord$vectors[,1:2])
# 
# scatter2D(geo_ord$vectors[,1], geo_ord$vectors[,2], pch=19, cex=1, 
#           colvar=coral_geo$Heliopora..HELIO.)
# 
# vectors <- envfit(geo_ord$vectors[,1:2], coral_geo)
# 
# plot(vectors)
# 
# ######################
# ## HOW CORRELATED ? ##
# ######################
# 
# cor.test(geo_coords$Axis.1, geo_coords_direct$Axis.1)
# cor.test(geo_coords$Axis.2, geo_coords_direct$Axis.2) #Rotation?
# 
# write.table(geo_coords, "coral_pcoa.txt")
# 
# write.table(geo_coords_direct, "coral_pcoa.txt")



#######################
## AT TRANSECT SCALE ##
#######################

####################
## ALL SUBSTRATES ##
####################

benthic_all <- read.csv("data/Benthic_ordination_all.csv")
benthic_all$site[which(nchar(benthic_all$site)<6)] <- gsub(pattern = "KOS-", replacement ="KOS-0", benthic_all$site[which(nchar(benthic_all$site)<6)])
benthic_meta <- benthic_all[,1:9]
benthic_meta$geographic <- with(kos_sites,
                                geographic[match(benthic_meta$site,
                                                 Site)])
benthos <- benthic_all[,10:ncol(benthic_all)]
missing_row <- which(rowSums(benthos)==0)
benthos <- benthos[,-which(colSums(benthos)==0)]
#benthos <- benthos[-missing_row,]
#benthic_all <- benthic_all[-missing_row,]
#benthic_meta <- benthic_meta[-missing_row,]

# benthos <- log10(benthos+1)

# bray_benthos <- vegdist(benthos, method = "bray")
# benthos_ord <- pcoa(bray_benthos)
# biplot(benthos_ord)
# 
# plot(benthos_ord$vectors[,1], benthos_ord$vectors[,2], pch=19, col=as.factor(benthic_meta$geographic))
# 
# spc_coords <- data.frame(benthos_ord$vectors[,1:2])
# geo_coords <- aggregate(spc_coords, by=list(benthic_meta$geographic), FUN = mean,na.rm=TRUE)   
# names(geo_coords)[1] <- "geographic"
# 
# 
# scatter2D(benthos_ord$vectors[,1], benthos_ord$vectors[,2], pch=19, cex=1, 
#           colvar=benthos$Heliopora..HELIO.)
# 
# vectors <- envfit(benthos_ord$vectors[,1:2], benthos)
# 
# plot(vectors)



#############################
##  DIRECTLY AT SITE LEVEL ##
#############################

kos_sites <- read.csv("data/KOS_spatial_biomass.csv")
benthic_all <- read.csv("data/Benthic_ordination_all.csv")
benthic_all$site[which(nchar(benthic_all$site)<6)] <- gsub(pattern = "KOS-", replacement ="KOS-0", benthic_all$site[which(nchar(benthic_all$site)<6)])
benthic_meta <- benthic_all[,1:9]
benthic_meta$geographic <- with(kos_sites,
                                geographic[match(benthic_meta$site,
                                                 Site)])
benthos <- benthic_all[,10:ncol(benthic_all)]
missing_row <- which(rowSums(benthos)==0)
benthos <- benthos[,-which(colSums(benthos)==0)]
benthos_site <- aggregate(benthos, by=list(benthic_meta$site), FUN = mean,na.rm=TRUE)
site <- benthos_site$Group.1
benthos_site$Group.1 <- NULL

#benthos_site <- log10(benthos_site+1)
benthos_site <- asinTransform(benthos_site/100)

site_bray <- vegdist(benthos_site, method = "bray")

site_ord <- pcoa(site_bray)

site_ord <- prcomp(benthos_site, scale. = FALSE)
fviz_pca_var(site_ord, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

biplot(site_ord)

site_coords_direct <- data.frame(site = site,
                                site_ord$vectors[,1:2])

scatter2D(site_ord$vectors[,1], site_ord$vectors[,2], pch=19, cex=2, 
          colvar=benthos_site$Porites.massive..PORMAS.)

plot(site_coords_direct$Axis.1, site_coords_direct$Axis.2, pch=19,)

vectors <- envfit(site_ord$vectors[,1:2], benthos_site)

vectors_table <- data.frame(vectors$vectors$arrows, vectors$vectors$pvals, vectors$vectors$r)

plot(vectors)


## DEFINE "LOADINGS" ##
axis1 <- benthos_site*site_ord$vectors[,1]
axis2 <- benthos_site*site_ord$vectors[,2]
load <- data.frame(colMeans(axis1), colMeans(axis2))
colnames(load) <- c("Axis_1","Axis_2")

load_reduced <- subset(load, abs(load$Axis_1) >= quantile(abs(load$Axis_1),probs=0.85) |
                         abs(load$Axis_2) >= quantile(abs(load$Axis_2), probs=0.85))

top_taxa <- benthos_site[,colnames(benthos_site) %in% rownames(load_reduced)]

vectors_reduced <- envfit(site_ord$vectors[,1:2], top_taxa)

plot(site_coords_direct$Axis.1, site_coords_direct$Axis.2,cex=0,xlab="Axis 1", ylab="Axis 2",cex.lab=1.25)
text(site_coords_direct$Axis.1, site_coords_direct$Axis.2, site_coords_direct$site,cex=1.25)
plot(vectors_reduced)

write.table(site_coords_direct, "benthic_pcoa_site.txt")


############################
##  DIRECTLY AT GEO LEVEL ##
############################

benthos_geo <- aggregate(benthos, by=list(benthic_meta$geographic), FUN = mean,na.rm=TRUE)
geographic <- benthos_geo$Group.1
benthos_geo$Group.1 <- NULL

#benthos_geo <- log10(benthos_geo+1)
benthos_geo <- asinTransform(benthos_geo/100)

geo_bray <- vegdist(benthos_geo, method = "bray")

geo_ord <- pcoa(geo_bray)

biplot(geo_ord)

geo_coords_direct <- data.frame(geographic = geographic,
                                geo_ord$vectors[,1:2])

scatter2D(geo_ord$vectors[,1], geo_ord$vectors[,2], pch=19, cex=2, 
          colvar=benthos_geo$Porites.massive..PORMAS.)

plot(geo_coords_direct$Axis.1, geo_coords_direct$Axis.2, pch=19,)

vectors <- envfit(geo_ord$vectors[,1:2], benthos_geo)

vectors_table <- data.frame(vectors$vectors$arrows, vectors$vectors$pvals, vectors$vectors$r)

plot(vectors)

# geo_cmdscale <- cmdscale(geo_bray, k=2)
# 
# test <- add.spec.scores(geo_cmdscale,benthos_geo,method="cor.scores",multi=1,Rscale=F,scaling="1")
# ordiplot(test)
# test <- add.spec.scores(geo_cmdscale,benthos_geo,method="was.scores",multi=1,Rscale=F,scaling="1")
# ordiplot(test)
# test <- add.spec.scores(geo_cmdscale,benthos_geo,method="pcoa.scores",multi=1,Rscale=F,scaling="1")
# ordiplot(test)

## DEFINE "LOADINGS" ##
axis1 <- benthos_geo*geo_ord$vectors[,1]
axis2 <- benthos_geo*geo_ord$vectors[,2]
load <- data.frame(colMeans(axis1), colMeans(axis2))
colnames(load) <- c("Axis_1","Axis_2")

load_reduced <- subset(load, abs(load$Axis_1) >= quantile(abs(load$Axis_1),probs=0.85) |
                         abs(load$Axis_2) >= quantile(abs(load$Axis_2), probs=0.85))

top_taxa <- benthos_geo[,colnames(benthos_geo) %in% rownames(load_reduced)]

vectors_reduced <- envfit(geo_ord$vectors[,1:2], top_taxa)

plot(geo_coords_direct$Axis.1, geo_coords_direct$Axis.2,cex=0,xlab="Axis 1", ylab="Axis 2",cex.lab=1.25)
text(geo_coords_direct$Axis.1, geo_coords_direct$Axis.2, geo_coords_direct$geographic,cex=1.25)
plot(vectors_reduced)

write.table(geo_coords_direct, "benthic_pcoa_geo.txt")

######################
## HOW CORRELATED ? ##
######################

cor.test(geo_coords$Axis.1, geo_coords_direct$Axis.1)
cor.test(geo_coords$Axis.2, geo_coords_direct$Axis.2) #Rotation?


write.table(geo_coords, "benthic_pcoa.txt")

write.table(geo_coords_direct, "benthic_pcoa.txt")


###############################
## TRY WITH PCA ON HELLINGER ##
###############################

benthos_geo <- aggregate(benthos, by=list(benthic_meta$geographic), FUN = mean,na.rm=TRUE)
benthos_geo$Group.1 <- NULL

benthos_hellinger <- decostand(benthos_geo, method = "hellinger")
geo_PCA <- prcomp(benthos_geo)
biplot(geo_PCA)
load <- data.frame(geo_PCA$rotation)





