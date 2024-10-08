
library(tidyverse)
library(jagsUI)
library(data.table)
library(grid)
library(ggpubr)
library(bayesplot)
library(HDInterval)
library(ggridges)

# set working directory
setwd("C:/Users/s444224/OneDrive - University of Canberra - STAFF/Honours 2023/aaa_Analysis/Occupancy modeling_")

# read in qPCR data
pcr <- read.csv("C:/Users/s444224/OneDrive - University of Canberra - STAFF/Honours 2023/aaa_Data/Final analysis csv data/PCR_results_cleaned.csv")
glimpse(pcr)

# exclude control sites (sites 31 - 35)
pcr <- pcr %>% filter(site < 31)

# label pcr reps and get unique sample values
pcr <- pcr |>
  group_by(site, sample) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  mutate(sample = paste(site, sample))

glimpse(pcr)

table(pcr$sample)

# data for the model
indat <- pcr |>
  dplyr::select(site, sample, rep, pa_edna) |>
  pivot_wider(values_from = pa_edna, names_from = rep, names_prefix = "rep") |>
  arrange(site, sample)

indat

y <- as.matrix(dplyr::select(indat, rep1, rep2, rep3))

# sites
site <- indat$site
N_site <- max(site)

# number of samples
N_samples <- nrow(indat)

# number of pcr reps
N_pcr_reps <- 3

# latent parameter for whether DNA is present
# at a site = z
# if DNA is detected at a site z = 1, if not z = NA
zz <- pcr |>
  group_by(site) |>
  summarise(n.pos = sum(pa_edna))

z <- ifelse(zz$n.pos > 0, 1, NA)

# in a sample = a
# if DNA is detected in a sample a = 1, if not a = NA
a <- ifelse(rowSums(y) > 0, 1, NA)

## Load in covariates
cov <- read.csv("C:/Users/s444224/OneDrive - University of Canberra - STAFF/Honours 2023/aaa_Data/Final analysis csv data/CT_occupancyModelling_covariates.csv")

cov <- cov %>% filter(site < 31) # remove the control sites

# function for standardising based on standard deviations
# ***note*** this is blocked out to facilitate making easier predictions
# it needs to be unblocked to generate nice figure for
sc <- function(x){
  (x-mean(x, na.rm = T))/(2*sd(x, na.rm = T))
}

table(is.na(cov$toad_density)) # Bolgers I lost data sheet so no time surveyed - no toads seen though so density set at 0 instead of NA
dens <- cov$toad_density %>% sc()

tads <- cov$tads

samp_vol <- cov$samp_vol %>% sc()

perim <- cov$perim %>% sc()

# ==============================================================================
# Occupancy model
# ==============================================================================

# JAGS code
mod <- "model {
# y is a binary variable taking values 1 (the jkth PCR replicate amplified)
#                                   or 0 (the jkth PCR replicate failed to amplify)

# loop through the samples & pcr reps
    for(j in 1:N_samples) {
      a[j] ~ dbern(z[site[j]] * prob.samp[j])
      
      logit(prob.samp[j]) <- beta0 + beta1 * tads[j] + beta2 * dens[j] + beta3 * samp_vol[j] + beta4 * perim[j]
      
      for(k in 1:N_pcr_reps) {
        y[j, k] ~ dbern(a[j] * prob.pcr)
      }
    }

# independent non-informative priors for each site
    for(i in 1:N_site) {
      z[i] ~ dbern(psi)
    }
    
# derived parameters
    N_occ <- sum(z[])
    
# priors
    psi ~ dunif(0, 1)
    prob.pcr ~ dunif(0, 1)
    beta0 ~ dnorm(0, 0.001)
    beta1 ~ dnorm(0, 0.001)
    beta2 ~ dnorm(0, 0.001)
    beta3 ~ dnorm(0, 0.001)
    beta4 ~ dnorm(0, 0.001)
}"

  # write model
  write(mod, "model.txt")

# run the model in JAGS with three chains
  mod.jags <- jags(model = "model.txt",
              data = list(y = y, site = site, N_samples = N_samples, N_pcr_reps = N_pcr_reps,
                          N_site = N_site, z = z, a = a, dens = dens, tads = tads, 
                          samp_vol = samp_vol, perim = perim),
#              inits = function() list(z = rep(1, N_site), a = rep(1, N_samples)),
              param = c("psi", "prob.samp", "prob.pcr", "z", "N_occ", 
                        "sensitivity", "beta0", "beta1", "beta2", "beta3", "beta4"),
              n.chains = 3,
              n.iter =50000,
              n.burnin = 20000,
              parallel = TRUE)


  #===============================================================================
  # ============================================ Get param estimates out
  #===============================================================================
  msdf <- as.data.frame(mod.jags$summary) %>% rownames_to_column()
  
  # Get posterior distribution values out of model
  mdf <- as.data.frame(as.matrix(mod.jags$samples)) %>%
    select(beta1, beta2, beta3, beta4)
  
  psi <- msdf %>% filter(rowname == "psi")
  prob.pcr <- msdf %>% filter(rowname == "prob.pcr")
  b0 <- msdf %>% filter(rowname == "beta0")
  b1 <- msdf %>% filter(rowname == "beta1") # note that tadpoles dist is skewed so using median
  b1.med <- median(mdf$beta1)
  b2 <- msdf %>% filter(rowname == "beta2")
  b3 <- msdf %>% filter(rowname == "beta3")
  b4 <- msdf %>% filter(rowname == "beta4")
  
  prob.samp_out <- mod.jags$summary[2:151, ] # mean and CRI  
  

  
#===============================================================================  
#===============================================================================
#===================== Figure 3
#===============================================================================
#===============================================================================
  


names(mdf) <- c("Tadpoles", "Toad encounter rate", "Volume sampled", "Waterbody perimeter")

mdf1 <- mdf %>% pivot_longer(names_to = "factor", cols = 1:4)
mdf1$factor <- as.character(mdf1$factor)
mdf1$factor <- factor(mdf1$factor,levels = c("Tadpoles", "Toad encounter rate", "Volume sampled", "Waterbody perimeter"))

# Posterior distribution plot for all factors excluding tadpoles
post1 <- ggplot() + 
  stat_density_ridges(mdf1 %>% filter(factor != "Tadpoles"), 
                      mapping = aes(x = value, y = 0, fill = after_stat(quantile)),
                      geom = "density_ridges_gradient",
                      calc_ecdf = TRUE,
                      quantiles = c(0.025, 0.975)) +
  
  stat_summary(mdf1 %>% filter(factor != "Tadpoles"), 
               mapping = aes(x = value,xintercept = after_stat(x), y = 0), 
               colour = "grey20", 
               fun = mean, 
               geom = "vline", 
               orientation = "y", 
               linewidth = .8) +
  
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.6) + 
  
  scale_fill_manual(values = c("transparent", "grey80", "transparent"), guide = "none") +
  scale_x_continuous(limits = c(-5, 5), breaks = c(-5,-2.5,0,2.5,5)) + 
  xlab("Standardised parameter value") +
  facet_wrap(~factor, scales = "free_y", ncol = 1) +

  theme_classic(5) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0),
        plot.background = element_blank(),
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 


# Posterior distribution plot for tadpoles
post2 <- ggplot() + 
  stat_density_ridges(mdf1 %>% filter(factor == "Tadpoles"), mapping = aes(x = value, y = 0, fill = after_stat(quantile)),
                      geom = "density_ridges_gradient",
                      calc_ecdf = TRUE,
                      quantiles = c(0.025, 0.975)) +
  stat_summary(mdf1 %>% filter(factor == "Tadpoles"), mapping = aes(x = value, 
                                                                    xintercept = after_stat(x), y = 0), colour = "grey20", 
               fun = median, geom = "vline", orientation = "y", linewidth = .8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", size = 0.6) + 
  scale_fill_manual(values = c("transparent", "grey80", "transparent"), guide = "none") +
  scale_x_continuous(limits = c(-150,150), breaks = c(-150,-75,0,75,150)) + 
  xlab("Value") +
  facet_wrap(~factor, scales = "free_y", ncol = 1) +
  
  
  theme_classic(5) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0),
        plot.background = element_blank(),
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) # +

p.dists <- ggarrange(post2 + rremove("x.title"), post1, ncol = 1, heights = c(1.1,3)) +
  theme(plot.margin = margin(5.5,0, 5.5, 15))

p.dists <- 
  annotate_figure(p.dists,
                  left = textGrob("Density",
                                  rot = 90, vjust = 1,
                                  gp = gpar(fontsize = 10))) 

p.dists


#===============================================================================
# Calculate survey sensitivity for each site (for histograms in fig 3)
#===============================================================================

# read in site level data
dat <- read.csv("C:/Users/s444224/OneDrive - University of Canberra - STAFF/Honours 2023/aaa_Data/Final analysis csv data/CT_allSurvey_siteLevel.csv")

glimpse(dat)

# select only relevant
dat <- dat %>% 
  select(site, mean_samp_vol, perim, tads, toad_density, survey_length, history) %>% 
  filter(site < 31)

# calculate survey sensitivity 
prob.samp.df <- as.data.frame(prob.samp_out)

cov$prob.samp <- prob.samp.df %>% pull(mean)
cov$zzz <- ((1 - cov$prob.samp) + cov$prob.samp * ((1 - prob.pcr %>% pull(mean)) ^ 3))

surveySens.df <- cov %>% 
  group_by(site) %>% 
  summarise(yyy = prod(zzz))  %>% 
  full_join(dat, by = "site")

surveySens.df$sens <- 1 - surveySens.df$yyy

# column for false negs with either method (1 = eDNA, 2 = visual)
surveySens.df$f.negs <- "both"
surveySens.df$f.negs[c(11,22,24,30)] <- "only visual"
surveySens.df$f.negs[c(1,6,8,25)] <- "only eDNA"

surveySens.df$history <- as.factor(surveySens.df$history)
levels(surveySens.df$history) <- c("Long invaded", "Recently invaded")

# ============== plot histogram 
sens.hist <- ggplot(data = surveySens.df, aes(x = sens))+
  geom_histogram(aes(fill = fct_rev(as.factor(f.negs))), colour = "black", bins = 10, alpha = 0.8, size = 0.3)+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.01,NA), oob = scales::oob_keep)+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  scale_fill_manual(values = c("white", "black", "#999999"))+
  labs(x = "Absolute site survey sensitivity", y = "Frequency")+ 
  guides(fill = guide_legend(title = "Detected with:")) +
  facet_wrap(history ~ .,dir = "v")+
  theme_classic(5)+
  theme(text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.35,0.85),
        plot.margin = margin(30,5.5,80,5.5),
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0))

sens.hist

# Generate figure 3
fig3 <- ggarrange(p.dists, sens.hist, labels = c("A", "B"), font.label = list(size = 10, face = "plain"),  ncol = 2, widths = c(1.3,1))
fig3
ggsave(plot = fig3, "Figure3.jpeg", width = 15, height = 15, units = "cm", dpi = 600)


# # calculate the number of samples the required for 0.95 sensitivity (for Table 4)
# surveySens.df$prob.samp <- b0 %>% pull(mean) + 
#   (b1 %>% pull(mean) * 0) + 
#   (b2 %>% pull(mean) * surveySens.df$toad_density) + 
#   (b3 %>% pull(mean) * surveySens.df$mean_samp_vol) + 
#   (b4 %>% pull(mean) * surveySens.df$perim)
# 
# surveySens.df$theta <-
#   exp(surveySens.df$prob.samp)/(1+exp(surveySens.df$prob.samp)) 
# 
# surveySens.df$n.samp95 <- log(0.05)/
#   (log( (1 - surveySens.df$theta) + 
#           surveySens.df$theta * ((1 - prob.pcr %>% pull(mean)) ^ n.pcrReps))) 
# 
# # filter for only false negative sites
# surveySens.df %>% filter(site %in% c(11,22,24,30))






# ==============================================================================
# ===== Number of samples vs PCR reps for increasing sensitivity (Appendix 1)
# ============================================================================== 
prob.samp <- c()
out <- c()
out0 <- c()
sim.sens <- 1 - seq(0, 1, 0.01)

# simulate a range of qPCR rep numbers 
n.pcr <- c(1, 2, 3, 6, 9)

for (i in n.pcr) {
  
  # calculate the prob.samp
    prob.samp <- b0 %>% pull(mean) + 
    (b1.med * 0) + 
    (b2 %>% pull(mean) * mean(dens)) + 
    (b3 %>% pull(mean) * mean(samp_vol)) + 
    (b4 %>% pull(mean) * mean(perim))
  # inverse logit transform
  theta <-exp(prob.samp)/(1+exp(prob.samp)) 
  # number of samples needed 
  n.samp95 <- log(sim.sens)/
    (log( (1 - theta) + 
            theta * ((1 - prob.pcr %>% pull(mean)) ^ i))) # i is pcr rep number
  # create df
  out0 <- data.frame(n.samp = n.samp95,
                     n.pcr = i,
                     sens = 1- sim.sens)
  # final df
  out <- rbind(out, out0)
}

out$n.pcr <- factor(out$n.pcr)
levels(out$n.pcr) <- c( "1", "2", "3", "6", "9")

rep.p0 <- ggplot(out, aes(y = sens, x = n.samp))+
  geom_line(aes(colour = n.pcr))+
  labs(x = "No. of samples", y = "Sensitivity", colour = "No. PCR reps")+
  theme_classic(18)+
  theme(legend.position = "bottom", text = element_text(size = 10))

# ============================================================== 
prob.samp <- c()
out <- c()
out0 <- c()
sim.sens <- 1 - seq(0, 1, 0.01)

n.pcr <- c(1, 2, 3, 6, 9)


for (i in n.pcr) {
  
  # calculate the mean prob.samp
  prob.samp <- b0 %>% pull(mean) + 
    (b1.med * 1) + 
    (b2 %>% pull(mean) * mean(dens)) + 
    (b3 %>% pull(mean) * mean(samp_vol)) + 
    (b4 %>% pull(mean) * mean(perim))

  # inverse logit transform
  theta <-exp(prob.samp)/(1+exp(prob.samp)) 
  
  # number of samples needed 
  # rearranged from Furlan et al. 2019 equation
  n.samp95 <- log(sim.sens)/
    (log( (1 - theta) + 
            theta * ((1 - prob.pcr %>% pull(mean)) ^ i))) # 3 is pcr rep number
  
  #
  out0 <- data.frame(n.samp = n.samp95,
                     n.pcr = i,
                     sens = 1- sim.sens)
  
  out <- rbind(out, out0)
}


out$n.pcr <- factor(out$n.pcr)
levels(out$n.pcr) <- c( "1", "2", "3", "6", "9")

rep.p1 <- ggplot(out, aes(y = sens, x = n.samp))+
  geom_line(aes(colour = n.pcr))+
  labs(x = "No. of samples", y = "Sensitivity", colour = "No. PCR reps")+
  theme_classic(18)+
  theme(legend.position = "bottom", text = element_text(size = 10))

#  arrange
ggarrange(rep.p1, rep.p0 + rremove("y.text") + rremove("y.ticks") + rremove("y.title"), 
          common.legend = T, legend = "bottom", labels = c("A", "B"),
          font.label = list(size = 10)) 

# save
ggsave("Samps vs reps plot.jpeg", height = 9, width = 14, units = "cm")





#===================================================================
#===================================================================
#======== Sample scheme plot (Figure 4 & Appendix 2)
#===================================================================
#===================================================================

# FUNCTION for estimating the number of samples required for 95% sensitivity 

samp.scheme_function <- function(beta0, beta1, beta2, beta3, beta4, prob.pcr.out) {
  
  # simulate a range of sampling conditions to estimate n_samp for 95
  # densities between 0 and 10 (toads seen per minute)
  
  d <- seq(0, 10, 0.01) #1:1001
  
  df <- c() # empty df for densities waterbody sizes (perimeter)
  df <- data.frame(dens = rep(d, 5), 
                   size = rep(c(10, 50, 100, 250, 500), each = 1001))
  
  # volumes per sample
  vols <- c(10, 50, 100, 250, 500)
  
  # empty objects for loop
  prob.samp0 <- c()
  out <- c()
  out0 <- c()
  n.pcrReps = 3 
  
  # for loop to calculate number of samples required given survey conditions
  for (i in vols) {
    
    # calculate the mean prob.samp
    prob.samp0 <- 
      beta0 +
      (beta1 * 0) +
      (beta2 * df$dens) + 
      (beta3 * i) + 
      (beta4 * df$size)
    
    
    # inverse logit transform
    theta0 <-exp(prob.samp0)/(1+exp(prob.samp0)) 
    
    # number of samples needed 
    # rearranged from Furlan et al. 2019 equation
    n.samp95 <- log(0.05)/
      (log( (1 - theta0) + 
              theta0 * ((1 - prob.pcr.out) ^ n.pcrReps))) # 3 is pcr rep number
    
    #
    out0 <- data.frame(size = df$size, 
                       dens = df$dens,
                       n = n.samp95,
                       vol = i,
                       tads = 0)
    
    out <- rbind(out, out0)
  }
  return(out)
  
} 

# use function to get predictions for mean
out <- samp.scheme_function(beta0 = b0 %>% pull(mean),
                            beta1 = b1.med,
                            beta2 = b2 %>% pull(mean),
                            beta3 = b3 %>% pull(mean),
                            beta4 = b4 %>% pull(mean),
                            prob.pcr.out = prob.pcr %>% pull(mean))

#=============================== Plot 
# vol labels
out$vol <- as.factor(out$vol)
levels(out$vol) <- c("10",  "50",  "100", "250", "500")

# perim labels
out$size <- as.factor(out$size)
levels(out$size) <- c("10",  "50",  "100", "250", "500")

# plot for no tad situations
samp.scheme <- ggplot(out, aes(x = dens, y = n))+
  geom_line(colour = "darkred")+
  labs(x = "Toad density (count per survey minute)",
       y = "Sample no. required for 95% sensitivity")+
  scale_y_continuous(n.breaks = 6, expand = c(0,0), limits = c(0,NA))+
  scale_x_continuous(n.breaks = 3, expand = c(0,0), limits = c(0,10))+
  geom_hline(yintercept = 5, linetype = "dotted")+
  facet_grid(as.factor(size) ~ as.factor(vol), scales = "free") +
  theme_light() +
  theme(panel.spacing = unit(0.6, "lines"), text = element_text(size = 10))
  

# arrange for export and add common axis labels
samp.schemeNoCrI <- ggarrange(samp.scheme + rremove("xlab") + rremove("ylab"))

samp.schemeNoCrI <- 
  annotate_figure(samp.schemeNoCrI,
                  bottom = textGrob("Toad encounter rate (count / survey minute)",
                                    vjust = 0, 
                                    gp = gpar(fontsize = 10)),
                  left = textGrob("No. of samples required for 95% confidence in detection",
                                  rot = 90, vjust = 1,
                                  gp = gpar(fontsize = 10)),
                  top = textGrob("Sample volume (ml)", 
                                 vjust = 1, 
                                 gp = gpar(fontsize = 10)),
                  right = textGrob("Waterbody perimeter (m)", 
                                   rot = 270, 
                                   vjust = 1, 
                                   gp = gpar(fontsize = 10))) 

samp.schemeNoCrI

# save
ggsave("CT No tads 95perc sample scheme from perimeter no CrI.jpeg", width = 16, height = 18, units = "cm")


#===============================================================================
#============== Create separate sampling scheme plot including CrI (Appendix 2)
#===============================================================================

# get predictions 
# lower 95% CrI 
l95.out <- samp.scheme_function(beta0 = b0 %>% pull(`2.5%`),
                                beta1 = b1 %>% pull(`2.5%`),
                                beta2 = b2 %>% pull(`2.5%`),
                                beta3 = b3 %>% pull(`2.5%`),
                                beta4 = b4 %>% pull(`2.5%`),
                                prob.pcr.out = prob.pcr %>% pull(`2.5%`))
# upper 95% CrI
u95.out <- samp.scheme_function(beta0 = b0 %>% pull(`97.5%`),
                                beta1 = b1 %>% pull(`97.5%`),
                                beta2 = b2 %>% pull(`97.5%`),
                                beta3 = b3 %>% pull(`97.5%`),
                                beta4 = b4 %>% pull(`97.5%`),
                                prob.pcr.out = prob.pcr %>% pull(`97.5%`))
# lower 50% CrI
l50.out <- samp.scheme_function(beta0 = b0 %>% pull(`25%`),
                                beta1 = b1 %>% pull(`25%`),
                                beta2 = b2 %>% pull(`25%`),
                                beta3 = b3 %>% pull(`25%`),
                                beta4 = b4 %>% pull(`25%`),
                                prob.pcr.out = prob.pcr %>% pull(`25%`))
# upper 59% CrI
u50.out <-samp.scheme_function(beta0 = b0 %>% pull(`75%`),
                               beta1 = b1 %>% pull(`75%`),
                               beta2 = b2 %>% pull(`75%`),
                               beta3 = b3 %>% pull(`75%`),
                               beta4 = b4 %>% pull(`75%`),
                               prob.pcr.out = prob.pcr %>% pull(`75%`))

samp.schemeCrI <- samp.scheme +
  geom_ribbon(aes(ymin = l95.out$n, ymax = u95.out$n), alpha = 0.1)+
  geom_ribbon(aes(ymin = l50.out$n, ymax = u50.out$n), alpha = 0.3)

samp.schemeCrI <- ggarrange(samp.schemeCrI + rremove("xlab") + rremove("ylab"))

samp.schemeCrI <- 
  annotate_figure(samp.schemeCrI, 
                  bottom = textGrob("Toad encounter rate (count / survey minute)",
                                    vjust = 0, 
                                    gp = gpar(fontsize = 10)),
                  left = textGrob("No. of samples required for 95% confidence in detection",
                                  rot = 90, 
                                  vjust = 1, 
                                  gp = gpar(fontsize = 10)),
                  top = textGrob("Sample volume (ml)", 
                                 vjust = 1, 
                                 gp = gpar(fontsize = 10)),
                  right = textGrob("Waterbody perimeter (m)", 
                                   rot = 270, 
                                   vjust = 1, 
                                   gp = gpar(fontsize = 10))) 

samp.schemeCrI

ggsave("CT No tads 95perc sample scheme from perimeter WITH CrI.jpeg", width = 16, height = 18, units = "cm")


p <- NCmisc::list.functions.in.file("C:/Users/s444224/OneDrive - University of Canberra - STAFF/Honours 2023/aaa_Analysis/R scripts_/2024_CaneToad_eDNA_Manuscript_OccupancyModelling_Fig3to6.R")
summary(p)
