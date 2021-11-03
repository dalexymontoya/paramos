setwd("")
library(betapart)
library(vegan)
library(corrplot)
library(usdm)
lapply(c("sp", "raster", "RgoogleMaps", "ggmap", "rgdal", "dplyr", "data.table", "vegan", "usdm"), library, character.only = TRUE)

#Read ocurrence of plants in paramo complexes
occ_all = read.csv("ocurrencias_consolidado.csv") # 148330 observations
occ_all = occ_all[which(occ_all$taxonRank == "SPECIES"),] # Removing 48879 determined to genus
sector = aggregate(occ_all$CODSECTOR~occ_all$CODCOMPLEJ, FUN = unique) #Sector refers to the mountain range in which the complex is located

length(unique(occ_all$scientificName)) # 4731 species of Paramo flora

clim = read.csv("predictors_paramos.csv") # A set of variables for each complex extracted usin

#Creation of abundance matrix
abd = table(occ_all$scientificName, occ_all$CODCOMPLEJ)
abd = t(abd)

#Transformation to presence - absence data
pa = decostand(abd, method = "pa")
rich = specnumber(pa) #richness per each complex
rch_sp = rich[order(rich, decreasing = T)] #Sorting the data
barplot(rch_sp, las= 2, main = "Riqueza de complejo de páramos en Colombia") #Visualizing richness on each complex
barplot(tail(rch_sp,4), las= 2, main = "Páramos con menos especies") #Richness of the paramos with the fewest species


#Removing paramo complexes with few species
rownames(pa)
id = as.numeric(which(apply(pa, 2, sum) == 0)) # Identifying species unrecorded
pa = pa[-c(9, 18, 25, 36),-id] # Removing row of paramo complexes with few species and columns of the species unrecorded
par_rich = specnumber(pa)


# Pcoa and envfit all species ---------------------------------------------------------
pcoa = capscale(pa~1, distance="bray")
summary(pcoa)

## get species and site scores
spp.sc <- scores(pcoa, display = "species", shrink = shrink)
site.sc <- scores(pcoa, display = "sites")

## work out the ranges
ylim <- range(spp.sc[,2], site.sc[,2])
xlim <- range(spp.sc[,1], site.sc[,1])

cols <- c("yellow","blue","green","white", "gray")
cols = rainbow(5)

#ordination plot
plot(site.sc, xlim = xlim, ylim = ylim, type = "n", asp = 1,
     xlab = "PCoA1", ylab = "PCoA2", main = "PCoA sorensen")

## add the row
colrs = cols[sector$`occ_all$CODSECTOR`]
fit_env <- envfit(ord=pcoa, env=clim[-c(9, 18, 25, 36),])
plot(fit_env, p.max = 0.05, col = "gray")
text(site.sc, col = colrs[-c(9, 18, 25, 36)], labels = rownames(pa), cex = 1)
text(-2.5,1, paste("stress value =", round(spe.nmds$stress,3), " "), cex = 0.7)
legend("bottomright",levels(sector$`occ_all$CODSECTOR`), col = cols, cex = 1, bty = "n", pch = 16, y.intersp = 0.7, x.intersp = 0.4)
legend("bottomright",levels(sector$`occ_all$CODSECTOR`), cex = 1, bty = "n", pch = 1, y.intersp = 0.7, x.intersp = 0.4)



# Beta diversity all species ----------------------------------------------
# Calculating simpson dissimilarity between paramo complexes using betapart

#Decomposition of beta diversity computed for multiple sites
simp_diss = beta.multi(pa, index.family="sorensen")
simp_diss #high species turnover

#Decomposition of beta diversity computed between pairs
simp_pairwise = beta.pair(pa, index.family = "sorensen")
simp_pairwise
corrplot(as.matrix(simp_pairwise$beta.sim), type = "upper", order="hclust")



# Linear model to predict richness of all species -------------------------
clim_ok = clim[-c(9, 18, 25, 36),]
rownames(clim_ok)
clim_ok$rich = as.numeric(par_rich)
vif_par = vifstep(clim_ok[,c(2:55)]) #compute variance inflation factor to identify which variables are less correlated
vif_par

attach(clim_ok)
model_paramo = lm(log10(par_rich)~Annual_range_sd+Precp_dry_m_sd+Precp_seasonality_sd+Precp_warm_Q_sd+Precp_cold_Q_sd+prec_mean+sr_mean+sr_sd+vapr_sd+ws_mean+elev_med+elev_min+slope_mean+slope_sd+age_dev)
summary(model_paramo)
mpar_step = step(model_paramo)
summary(mpar_step)

m_par = lm(log10(par_rich) ~ Annual_range_sd + Precp_dry_m_sd + Precp_seasonality_sd + 
             Precp_warm_Q_sd + Precp_cold_Q_sd + prec_mean + sr_mean + 
             sr_sd + vapr_sd + ws_mean + elev_med + elev_min + slope_mean + 
             slope_sd + age_dev)
summary(m_par)
par(mfrow=c(2,2))
plot(m_par)


# Pcoa and envfit unique species ---------------------------------------------------------
id2 = as.numeric(which(apply(pa, 2, sum) == 1))
pa_spu = pa[,id2]
pcoa = capscale(pa_spu~1, distance="bray")
summary(pcoa)

## get species and site scores
spp.sc <- scores(pcoa, display = "species", shrink = shrink)
site.sc <- scores(pcoa, display = "sites")

## work out the ranges
ylim <- range(spp.sc[,2], site.sc[,2])
xlim <- range(spp.sc[,1], site.sc[,1])

cols <- c("yellow","blue","green","white", "gray")
cols = rainbow(5)

#ordination plot
plot(site.sc, xlim = xlim, ylim = ylim, type = "n", asp = 1,
     xlab = "PCoA1", ylab = "PCoA2", main = "PCoA sorensen")

## add the row
colrs = cols[sector$`occ_all$CODSECTOR`]
fit_env <- envfit(ord=pcoa, env=clim[-c(9, 18, 25, 36),])
plot(fit_env, p.max = 0.05, col = "gray")
text(site.sc, col = colrs[-c(9, 18, 25, 36)], labels = rownames(pa), cex = 1)
text(-4,1, paste("stress value =", round(spe.nmds$stress,3), " "), cex = 0.7)
legend("bottomright",levels(sector$`occ_all$CODSECTOR`), col = cols, cex = 1, bty = "n", pch = 16, y.intersp = 0.7, x.intersp = 0.4)
legend("bottomright",levels(sector$`occ_all$CODSECTOR`), cex = 1, bty = "n", pch = 1, y.intersp = 0.7, x.intersp = 0.4)


# Linear model to predict richness of unique species -------------------------
par_rich = specnumber(pa_spu[-c(9, 18, 25, 36),])
clim_ok = clim[-c(9, 18, 25, 36),]
rownames(clim_ok)
clim_ok$rich = as.numeric(par_rich)
vif_par = vifstep(clim_ok[,c(2:55)]) #compute variance inflation factor to identify which variables are less correlated
vif_par

attach(clim_ok)
model_paramo = lm(log10(par_rich)~Annual_range_sd+Precp_dry_m_sd+Precp_seasonality_sd+Precp_warm_Q_sd+Precp_cold_Q_sd+prec_mean+sr_mean+sr_sd+vapr_sd+ws_mean+elev_med+elev_min+slope_mean+slope_sd+age_dev)
summary(model_paramo)
mpar_step = step(model_paramo) # Selection of a model based on AIC
summary(mpar_step)

m_par = lm(log10(par_rich) ~ Annual_range_sd + sr_mean + ws_mean + 
             elev_med + slope_mean)
summary(m_par)
par(mfrow=c(2,2))
plot(m_par)

