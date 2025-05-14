#Distribution of Mortierella spp. in Pennsylvania soybean agroecosystems 
#Authors: Dr. Michelle Paukett, and Dr. Paul Esker
#Department of Plant Pathology & Environmental Microbiology, The Pennsylvania State University, University Park, PA 16802, U.S.A.

#####Citations----
#Bui, An. (2020). Vegan cheat sheet. RPubs, <https://rpubs.com/an-bui/vegan-cheat-sheet>.

#Müller K (2020). _here: A Simpler Way to Find Your Files_. R package version 1.0.1,
#  <https://CRAN.R-project.org/package=here>.

#Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos
#  P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Carvalho
#  G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux
#  B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T,
#  Stier A, Ter Braak C, Weedon J (2022). _vegan: Community Ecology Package_. R package
#  version 2.6-4, <https://CRAN.R-project.org/package=vegan>.

#R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R
#  Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.

#Simpson G (2019). _ggvegan: 'ggplot2' Plots for the 'vegan' Package_. R package version 0.1-0.

#Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A,
#  Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J,
#  Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H
#  (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686.
#  doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

#####Notes----
#Analysis completed on 125 Mortierella spp. isolates
#Sites are ordered 1-29 in all files, even if no actual site column listed
#Analysis was completed at the regional level (Group / County Location - West, North, or South)

#Sections include----

#####Load necessary packages-----

library(here)                                                                    #File referencing/path-making for projects
library(tidyverse)                                                               #Data manipulation, includes dyplr & ggplot2 
library(vegan)                                                                   #Community ecology, primary package for distribution analysis
library(ggvegan)                                                                 #ggplot-based plots for vegan

#####Load necessary files----

mspecies <- read.csv("Mortierellapresence.csv")                                  #Species abundance per site
paenv <- read.csv("PASitesMetadata.csv")                                         #Site environmental metadata
mortr <- read.csv("MortierellaMetadata.csv")                                     #Mortierella spp. metadata

#####Create Tables & Data Frames for Analysis ----

site_type <- paenv %>% 
  select(group, texture, county)                              #Create an environmental metadata table

ov <- as_tibble(mortr)                                                           #Convert Mortierella metadata to a tibble 

#####Basic Species Presence----

regtotals <- ov %>%
  group_by(County.Location) %>%
  count(Seq.Result)%>%
  mutate(Percent = n/sum(n)*100)                                                 #Creates a table showing spp count & percentage by region

ggplot(regtotals, aes(fill=Seq.Result, y=Percent, x=County.Location)) + 
  geom_col(position="fill") + 
  labs(x="PA Region", y="Percent of Species Identified", fill="Species") + 
  scale_y_continuous(labels = scales::percent) + theme_bw() + 
  theme(text=element_text(size=9)) + 
  theme(legend.text = element_text(face = "italic")) + 
  scale_x_discrete(labels=c("N" = "North", "S" = "Southeast", "W" = "West"))     #Make a stacked bar graph

#####Species Richness----

sppr <- specnumber(mspecies)                                                     #Calculate & save # of species at each site
sppr_aov <- aov(sppr ~ group, data = site_type)                                 #Run anova for significance
summary(sppr_aov)                                                                #Show results in console
capture.output(sppr_aov, file = "sppr_aov.txt")                                  #Save results in text file

tky_sppr <- TukeyHSD(sppr_aov)                                                   #Run Tukey test to compare across treatments
capture.output(tky_sppr, file = "sppr_tukey.txt")                                #Save results in text file

#####Plot Species Richness---- 

sppr_df <- sppr %>% 
  enframe()                                                                      #Make a data frame of species richness

sppr_df_env <- cbind(sppr_df, site_type[1:3])                                    #Add environmental variables from metadata to data frame

gsppr_df_env_df <- sppr_df_env %>% 
  group_by(group) %>% 
  summarize(mean = round(mean(value), 2),
            err = sd(value)/sqrt(length(value))) %>% 
  dplyr::mutate(label = "mean") %>% 
  unite("mean_label", label, mean, sep = " = ", remove = FALSE)                  #Create data frame with SE & Mean for plot labels

plot_spprg2 <- ggplot(gsppr_df_env_df, aes(x = group, y = mean, fill = group)) +
  geom_col(color = "black") + update_geom_defaults("text", list(size=3))+
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err), width = 0.5) +
  geom_text(aes(x = group, y = mean + err + 0.07, label = mean_label)) +
  labs(x = "Region",
       y = "Number of Species per Region")  + theme_bw(base_size=9) +   
  scale_x_discrete(labels=c("N" = "North", "SE" = "Southeast", "SW" = "West"))   #Make a boxplot of species richness
plot_spprg2                                                                      #Plot boxplot

#####Calculate Shannon Diversity----

shannondiv <- diversity(mspecies)                                                #Calculate Shannon diversity
head(shannondiv)                                                                 #Show results in console

sppdiv_aov <- aov(shannondiv ~ group, data = site_type)                         #Run anova on Shannon Diversity
summary(sppdiv_aov)                                                              #Show results in console
capture.output(sppdiv_aov, file = "sppdiv_aov.txt")                              #Save anova results in text file

tky_sppdiv <- TukeyHSD(sppdiv_aov)                                               #Tukey's test to compare across treatments
capture.output(tky_sppdiv, file = "sppdiv_tky.txt")                              #Save Tukey's test results in text file

#####Plot Shannon Diversity----

shandiv_df <- shannondiv %>% 
  enframe()                                                                      #Make a data frame of diversity

shandiv_df_env <- cbind(shandiv_df, site_type[1:3])                              #Add environmental variables from metadata to data frame

Gdiv_plot_df <- shandiv_df_env %>% 
  group_by(group) %>% 
  summarize(mean = round(mean(value), 2),
            err = sd(value)/sqrt(length(value))) %>% 
  dplyr::mutate(label = "mean") %>% 
  unite("mean_label", label, mean, sep = " = ", remove = FALSE)                  #Create data frame with SE & Mean for plot labels

Gplot_shandiv <- ggplot(Gdiv_plot_df, aes(x = group, y = mean, fill = group)) +
  geom_col(color = "black") + update_geom_defaults("text", list(size=3))+
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err), width = 0.5) +
  geom_text(aes(x = group, y = mean + err + 0.07, label = mean_label)) +
  labs(x = "Region",
       y = "Shannon Diversity")   + theme_bw(base_size=9) +                     
  scale_x_discrete(labels=c("N" = "North", "SE" = "Southeast", "SW" = "West"))   #Plot boxplot of Shannon Diversity
Gplot_shandiv                                                                    #Show boxplot

#####PERMANOVA----

div_perm <- adonis2(mspecies ~ group, data = paenv)                              #Run PERMANOVA 
div_perm                                                                         #Show results in console
capture.output(div_perm, file = "div_perm.txt")                                  #Save results in text file

#####PCA (Principal Components Analysis)----

mPCA <- rda(mspecies)                                                            #Run PCA
mPCA                                                                             #Show results in console
capture.output(mPCA, file = "mPCA.txt")                                          #Save results in text file 
                                                                                 #Note: Manually add PC values from this for labels in plot below

PCA_biplot <- autoplot(mPCA)                                                     #Easy PCA autoplot using ggvegan
PCA_biplot                                                                       #Show plot 

                                                                                 #For customized PCA plot using ggplot2 
PCA_fortify <- fortify(mPCA)                                                     #Makes a data frame from which elements can be extracted

PCA_fort_sites <- PCA_fortify %>% 
  filter(Score == "sites")                                                       #Extract site coordinates (points)

PCA_fort_sites_env <- cbind(PCA_fort_sites, site_type[1:3])                      #Add environmental variables from metadata to data frame

PCA_fort_species <- PCA_fortify %>% 
  filter(Score == "species")                                                     #Extract species coordinates (arrows)

GPCA_fortify_plot <- ggplot() +
  geom_point(data = PCA_fort_sites_env, aes(x = PC1, y = PC2, col = group)) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = PCA_fort_species, aes(x = 0, xend = PC1, y = 0,
               yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = PCA_fort_species, aes(x = PC1, y = PC2, label = Label))+
  labs(x = "PC1 (3.45%)", y = "PC2 (1.36%)") + theme_bw(base_size=9)             #Plot ordination
                                                                                 #Note: PC %'s for axes were manually added!
GPCA_fortify_plot                                                                #Show ordination plot 

#####Non-metric Multidimensional Scaling (NMDS)----

mort_NMDS <- metaMDS(mspecies)                                                   #Run NMDS
mort_NMDS                                                                        #Show results in console
capture.output(mort_NMDS, file = "mort_NMDS.txt")                                #Save results text file

stressplot(mort_NMDS)                                                            #Plot NMDS stress 
plot(mort_NMDS)                                                                  #Show stress plot

plot_df <- scores(mort_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site")                                                     #Create data frame for NMDS plot
plot_df                                                                          #Show data frame in console

k <- cbind(plot_df, site_type[1:3])                                              #Add environmental variables from metadata to data frame

Gplot_nmds <- ggplot(k, aes(x = NMDS1, y = NMDS2, color = group, 
  shape = group)) + geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(linewidth = 2, size = 1) +
  theme_bw() + theme(text=element_text(size=9))                                  #Plot NMDS
Gplot_nmds                                                                       #Show NMDS plot 

#####Envfit-NMDS with Species----

fit <- envfit(mort_NMDS, mspecies, perm = 999)                                   #Add species into NMDS ordination & identify significance
fit                                                                              #Show results in console

fit_pvals <- fit$vectors$pvals %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::rename("pvals" = ".")                                                   #Extract p-values for each species

fit_spp <- fit %>% 
  scores(., display = "vectors") %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  full_join(., fit_pvals, by = "species") %>% 
  filter(pvals < .05)                                                            #Extract coordinates for species, only keep species with p-val < 0.05***

Gnmds_plot_new <- ggplot(k, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() +
  geom_point(aes(color = group, shape = group), size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = group)) +
  geom_segment(data = fit_spp, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), col = "black") +
  geom_text(data = fit_spp, aes(label = species)) + theme_bw() + 
  theme(text=element_text(size=9))                                               #Plot NMDS with significant species
Gnmds_plot_new                                                                   #Show species NMDS plot 

#####Redundancy analysis (RDA)----

rdatop <- rda(mspecies~elevation+temp+lat+long+snow+watervol+rain, data=paenv)   #Run RDA
rdatop                                                                           #Show RDA results in console
autoplot(rdatop, arrows = TRUE) + theme_bw() + theme(text=element_text(size=9))  #Autoplot RDA

rda_perm <- adonis2(mspecies ~ elevation+temp+lat+long+snow+watervol+rain, 
                   data = paenv)                                                 #Run PERMANOVA
rda_perm                                                                          #Show PERMANOVA results in console
capture.output(rda_perm , file = "rda_perm.txt")                                   #Save results in text file

#Canonical Correspondence Analysis (CCA)----

mortCCA <- cca(mspecies ~ elevation+temp+snow+lat+long+watervol+rain,
               data = paenv)                                                     #Run CCA
mortCCA                                                                          #Show CCA results in console
capture.output(mortCCA, file = "mortCCA.txt")                                    #Save results in text file

plot(mortCCA)                                                                    #Autoplot CCA

ccavectors <- as.matrix(scores(mortCCA, display = "bp", 
  scaling = "species")*7.627807) %>% 
  as.data.frame()                                                                #Get environmental vectors (arrows) for plot

site_data <- scores(mortCCA, display = "sites") %>% 
  as.data.frame() %>%  
  rownames_to_column("site")                                                     #Get site coordinates for plot

L <- cbind(site_data, site_type[1:3])                                            #Add environmental variables from metadata to data frame

species_data <- scores(mortCCA, display = "species") %>% 
  as.data.frame()                                                                #Get species coordinates for plot

Gplot_cca <- ggplot(L) +
  geom_point(aes(x = CCA1, y = CCA2, color = group),
  shape = 19, size = 2, alpha = 0.7) + coord_fixed() +
  geom_segment(data = ccavectors, aes(x = 0, y = 0, xend = CCA1,
  yend = CCA2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = species_data, aes(x = CCA1, y = CCA2), 
  shape = 17, size = 2, color = "slateblue") +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-4, 5.2)) +
  geom_text(data = ccavectors, aes(x = CCA1, y = CCA2, 
  label = rownames(ccavectors)), nudge_x = 0.3, nudge_y = 0.3) +
  theme_bw() + theme(text=element_text(size=9))                                  #Plot CCA 
Gplot_cca                                                                        #Show CCA plot

Gplot_ccaspp <- Gplot_cca +
  geom_text(data = species_data, aes(x = CCA1, y = CCA2, 
  label = rownames(species_data)), nudge_x = 0.3, nudge_y = 0.3)                 #Add spp names to plot
Gplot_ccaspp                                                                     #Show CCA plot with spp labels