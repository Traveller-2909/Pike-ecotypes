library(tidyverse)
library(vegan)
library(indicspecies)


# Load Data ---------------------------------------------------------------
oto <- read.csv2("Full Concentration Dataset.csv", sep = ",", stringsAsFactors = FALSE, dec = ".")
oto$Site <- as.factor(oto$Site)
oto$Fish <- as.factor(oto$Fish)

# Data Cleanup ------------------------------------------------------------

oto %>% 
    select(-Pb,-Y, -Ca) %>% 
    filter(Distanz_um > 0) %>% 
    mutate(Mn = ifelse(Mn < 0, 0, Mn)) %>% 
    mutate(Cu = ifelse(Cu < 0, 0, Cu)) %>% 
    group_by(Fish) %>% 
    summarise(site = first(Site),
              Mg_m = mean(Mg),
              Mn_m = mean(Mn),
              Cu_m = mean(Cu),
              Zn_m = mean(Zn),
              Sr_m = mean(Sr),
              Ba_m = mean(Ba)) -> pike

# as matrix, only data columns
pike_m <- pike[,3:ncol(pike)]
pike_m <- as.matrix(pike_m)
mode(pike_m) = "numeric"


# NMDS with averages ------------------------------------------------------
#initial plot with Bray distance measure
nmds = metaMDS(pike, distance = "bray")
plot(nmds) 

#Extract nmds scores
data.scores = as.data.frame(scores(nmds))

#add columns to data frame 
data.scores$ID <- pike$Fish
data.scores$Site = pike$site

head(data.scores)

#Plotting with NMDS scores
pplot <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 3, aes(colour = Site))+ 
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) +  
    annotate("text", x = 0, y = 0.24, label = paste0("Anosim R value = 0.4539",
             "; Anosim p < 0.001"), color = "#F8766D", fontface = 2, hjust = 0, vjust = 1)+
    labs(x = "NMDS1", colour = "Site", y = "NMDS2")  
pplot

# Statistical Test
ano = anosim(pike_m, pike$site, distance = "bray", permutations = 9999)
summary(ano)

# Species indicator test
pike_df = pike[,3:ncol(pike)]
psite = pike$site
inv = multipatt(pike_df, psite, func = "r.g", control = how(nperm=9999))
summary(inv)

#PCA (Timo)---------------------------------------------------------------------
#PCA visualization
library(factoextra)

#Compute PCA for all values
pike_pca <- prcomp(pike_m, scale = TRUE)

#Eigenvalue plot
fviz_eig(pike_pca)

#PCA visual for indviduals
groups <- as.factor(pike$site)
fviz(pike_pca, "ind", color = groups)

#PCA visual for variables
fviz_pca_var(pike_pca, col.var = "contrib", #color by contribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) #avoids text overlapping
#Individuals and variables
fviz_pca_biplot(pike_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
#Eigenvalues
eig.val <- get_eigenvalue(pike_pca)
eig.val             

#Results for variables
res.var <- get_pca_var(pike_pca)
res.var$coord                 #coordinates
res.var$contrib               #contribution to PCs
res.var$cos2                  #Quality of representation

#Same for individuals
res.ind <- get_pca_ind(pike_pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

#PCA similar to Engstedt et al 2014
oto %>% 
    select(-Pb,-Y, -Ca) %>% 
    filter(Distanz_um > 50 & Distanz_um < 250) %>% 
    mutate(Mn = ifelse(Mn < 0, 0, Mn)) %>% 
    mutate(Cu = ifelse(Cu < 0, 0, Cu)) %>% 
    group_by(Fish) %>% 
    summarise(site = first(Site),
              Mg_m = mean(Mg),
              Mn_m = mean(Mn),
              Cu_m = mean(Cu),
              Zn_m = mean(Zn),
              Sr_m = mean(Sr),
              Ba_m = mean(Ba)) -> pike_young

# as matrix, only data columns
pike_j <- pike_young[,3:ncol(pike_young)]
pike_j <- as.matrix(pike_j)
mode(pike_j) = "numeric"


#Compute PCA for first 200 um
pikej_pca <- prcomp(pike_j, scale = TRUE)

#Eigenvalue plot
fviz_eig(pikej_pca)

#PCA visual for indviduals
groups <- as.factor(pike_young$site)
fviz(pikej_pca, "ind", color = groups)

#PCA visual for variables
fviz_pca_var(pikej_pca, col.var = "contrib", #color by contribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) #avoids text overlapping
#Individuals and variables
fviz_pca_biplot(pikej_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
#Anosim for younglings

anoj = anosim(pike_j, pike_young$site, distance = "bray", permutations = 9999)
summary(anoj)

#Same analysis without Peenestrom-----------------------------------------------
oto %>% 
    select(-Pb,-Y, -Ca) %>% 
    filter(Distanz_um > 0 & Site != "Peene") %>% 
    mutate(Mn = ifelse(Mn < 0, 0, Mn)) %>% 
    mutate(Cu = ifelse(Cu < 0, 0, Cu)) %>% 
    group_by(Fish) %>% 
    summarise(site = first(Site),
              Mg_m = mean(Mg),
              Mn_m = mean(Mn),
              Cu_m = mean(Cu),
              Zn_m = mean(Zn),
              Sr_m = mean(Sr),
              Ba_m = mean(Ba)) -> pike_bodden

# as matrix, only data columns
pike_bm <- pike_bodden[,3:ncol(pike_bodden)]
pike_bm <- as.matrix(pike_bm)
mode(pike_bm) = "numeric"


# NMDS with averages ------------------------------------------------------
#initial plot with Bray distance measure
bnmds = metaMDS(pike_bm, distance = "bray")
plot(bnmds) 

#Extract nmds scores
data.scores = as.data.frame(scores(bnmds))

#add columns to data frame 
data.scores$ID <- pike_bodden$Fish
data.scores$Site = pike_bodden$site

head(data.scores)

#Plotting with NMDS scores
pplot_b <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 3, aes(colour = Site))+ 
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) +  
    annotate("text", x = -0.1, y = 0.24, label = paste0("Anosim R value = 0.071",
                                                     "; Anosim p < 0.1"), color = "#F8766D", fontface = 2, hjust = 0, vjust = 1)+
    labs(x = "NMDS1", colour = "Site", y = "NMDS2") 
pplot_b

# Statistical Test
anob = anosim(pike_bm, pike_bodden$site, distance = "bray", permutations = 9999)
summary(anob)

#Compute PCA for Bodden
pikeb_pca <- prcomp(pike_bm, scale = TRUE)

#Eigenvalue plot
fviz_eig(pikeb_pca)

#PCA visual for indviduals
groups <- as.factor(pike_bodden$site)
fviz(pikeb_pca, "ind", color = groups)

#PCA visual for variables
fviz_pca_var(pikeb_pca, col.var = "contrib", #color by contribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) #avoids text overlapping

# Species indicator test
pike_bdf = pike_bodden[,3:ncol(pike_bodden)]
psite = pike_bodden$site
binv = multipatt(pike_bdf, psite, func = "r.g", control = how(nperm=9999))
summary(binv)
