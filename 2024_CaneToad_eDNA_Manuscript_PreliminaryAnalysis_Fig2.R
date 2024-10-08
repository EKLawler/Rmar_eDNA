
#===============================================================================
#=============== Script to generate figure 2 
#===============================================================================

library(tidyverse)
library(ggpubr)
library(grid)
library(ggpattern)

setwd("")

# read in site level data
df <- read.csv("CT_allSurvey_siteLevel.csv")

# read in sample level data
cov <- read.csv("CT_occupancyModelling_covariates.csv")

# exclude control sites for analysis and remove excess column
df <- df %>% filter(site < 31) %>% select(-c(X))

#-------------------------------------------------------------------------------
#-------------------- Histograms for factors
#-------------------------------------------------------------------------------

#-------------------- Toad density

df$history <- factor(df$history)
levels(df$history) <- c("Long invaded", "Recently invaded")

dplot1 <- ggplot(df, aes(x = toad_density, pattern = as.factor(pa_vis)))+
  geom_histogram_pattern(breaks = c(-1, 0,0.000000001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
                         fill = "grey80", 
                         colour = "black", 
                         linewidth = 0.3, 
                         pattern_density = 0.05)+
  scale_x_continuous(breaks = c(0,5,10), limits = c(-1,10))+
  scale_y_continuous(limits = c(0,9), breaks = c(0,4,8))+
  scale_pattern_manual(values = c("stripe", "none"))+
  scale_fill_brewer(guide = guide_legend(override.aes = list(pattern = c("stripe"))))+
  labs(x = "Toad encounter rate (count / survey minute)", y = "Frequency")+
  ggthemes::theme_base()+
  theme(text = element_text(size = 10), 
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0),
        plot.background = element_blank())+
  facet_wrap(~history)

#-------------------- Perimeter histogram

p.hist <- ggplot(df, aes(x = perim))+
  geom_histogram(bins = 20, colour = "black", fill = "grey80", linewidth = 0.3)+
  labs(x = "Waterbody perimeter (m)", y = "Frequency") + 
  scale_y_continuous(limits = c(0, 9), breaks = c(0,4,8))+
  scale_x_continuous(limits = c(0, 630), oob = scales::oob_keep,  breaks = c(0,300,600))+
  ggthemes::theme_base()+
  theme(text = element_text(size = 10),
        plot.background = element_blank())

#-------------------- Sample volume histogram

v.hist <- ggplot(df, aes(x = mean_samp_vol))+
  geom_histogram(bins = 20, colour = "black", fill = "grey80", size = 0.3)+
  scale_y_continuous(limits = c(0, 8), breaks = c(0,4,8))+
  scale_x_continuous(limits = c(0, 525), oob = scales::oob_keep, breaks = c(0,250,500))+
  labs(x = "Mean volume sampled (ml)", y = "Frequency")+
  ggthemes::theme_base()+
  theme(text = element_text(size = 10),
        plot.background = element_blank())


#--------------------- Arrange into one plot
all.hists <- ggarrange(dplot1 + rremove("ylab"),
          v.hist + rremove("ylab"), 
          p.hist + rremove("ylab"),
          ncol = 1, nrow = 3, 
          labels = c("B", "D", "F"), 
          hjust = 0.8, 
          font.label = list(size = 10, face = "plain"))

all.hists <- annotate_figure(all.hists,
                             left = textGrob("Frequency", 
                                             rot = 90, gp = gpar(fontsize = 10)))

all.hists
# save
ggsave("Perim and samp vol hists.jpeg", height = 15, width = 10, units = "cm")


#-------------------------------------------------------------------------------
#-------------------- GLMERs for each factor against the proportion of reps amplified
#-------------------------------------------------------------------------------

od <- 1:nrow(df) # to account for dispersion

# toad density + tadpoles
m1 <- lme4::glmer(cbind(npos,nneg) ~ toad_density + tads + (1|od), 
                  data = df,
                  family = "binomial" )

# sample volume + tadpoles
m2 <- lme4::glmer(cbind(npos,nneg) ~ mean_samp_vol + tads + (1|od), 
                  data = df,
                  family = "binomial" )

# perimeter + tadpoles
m3 <- lme4::glmer(cbind(npos,nneg) ~ perim + tads + (1|od), 
                  data = df,
                  family = "binomial" )

summary(m1) 
summary(m2)
summary(m3)

#-------------------------------------------------------------------------------
#-------------------- Plot GLMER
#-------------------------------------------------------------------------------

#--------------------- Toad density plot

# plot points 
dens_p <- ggplot(df, aes(x = toad_density, 
                              y = rep_prob_detect, 
                              colour = as.factor(tads)))+
  geom_point()+
  scale_color_manual(values = c("red", "blue"))+
  xlab( "Toad encounter rate (count / survey minute)")+
  ylab("Proportion of PCR replicates amplified")+
  labs(colour = "Tadpoles")+
  scale_x_continuous(expand = c(0,0), limits = c(-0.5,10.5), breaks = c(0,5,10))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.5,1)) +
  theme_classic(5)+
  theme(legend.position = "bottom",
        text = element_text(size = 10))

# create model to fit lines for plot 
dens_mod <- (glm(cbind(npos,nneg) ~ toad_density + tads, 
                 data = df,
                 family = "binomial")) 

# create data for tadpoles present and absent 
tads1 <- df %>% filter(tads == 1)
tads0 <- df %>% filter(tads == 0)

xx1 <- seq(0, max(tads1$toad_density, na.rm = T), by = 0.01)
pred_df1 <- data.frame(toad_density = xx1, tads = 1)
pred_df1$pred <- predict(dens_mod, newdata = pred_df1, type = "response")

xx0 <- seq(0, max(tads0$toad_density, na.rm = T), by = 0.01)
pred_df0 <- data.frame(toad_density = xx0, tads = 0)
pred_df0$pred <- predict(dens_mod, newdata = pred_df0, type = "response")

# add lines to point plot
dens_p <- dens_p + 
  geom_line(data = pred_df1, mapping = aes(x = toad_density, y = pred),
            colour = "blue")+
  geom_line(pred_df0, mapping = aes(x = toad_density, y = pred),
            colour = "red") 


#--------------------- Sample volume plot

vol_p <- ggplot(df, aes(x = mean_samp_vol, 
                          y = rep_prob_detect, 
                          colour = as.factor(tads)))+
  geom_point()+
  scale_color_manual(values = c("red", "blue"))+
  xlab("Mean volume sampled (ml)")+
  ylab("Proportion of PCR replicates amplified")+
  labs(colour = "Tadpoles seen")+
  scale_x_continuous(expand = c(0,0), limits = c(-25,525), breaks = c(0,250,500))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.5,1)) +
  theme_classic(5)+
  theme(legend.position = "bottom",
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

# create model to fit lines for plot
sampVol_mod <- glm(cbind(npos,nneg) ~ mean_samp_vol + tads, 
            data = df,
            family = "binomial") 

# create data for tadpoles present and absent
xx1 <- seq(0, max(tads1$mean_samp_vol, na.rm = T), by = 1)
pred_df1 <- data.frame(mean_samp_vol = xx1, tads = 1)
pred_df1$pred <- predict(sampVol_mod, newdata = pred_df1, type = "response")

xx0 <- seq(0, max(tads0$mean_samp_vol, na.rm = T), by = 0.01)
pred_df0 <- data.frame(mean_samp_vol = xx0, tads = 0)
pred_df0$pred <- predict(sampVol_mod, newdata = pred_df0, type = "response")

# add lines to point plot
vol_p <- vol_p + 
  geom_line(data = pred_df1, mapping = aes(x = mean_samp_vol, y = pred),
            colour = "blue")+
  geom_line(pred_df0, mapping = aes(x = mean_samp_vol, y = pred),
            colour = "red")

#--------------------- Perimeter plot

perim_p <-  ggplot(df, aes(x = perim, 
                       y = rep_prob_detect, 
                       colour = as.factor(tads)))+
  geom_point()+
  theme_classic(5)+
  scale_color_manual(values = c("red", "blue"))+
  xlab("Waterbody perimeter (m)")+
  ylab("Proportion of PCR replicates amplified")+
  labs(colour = "Tadpoles seen")+
  scale_x_continuous(expand = c(0,0), limits = c(-30,630), breaks = c(0,300,600))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.5,1)) +
  theme(legend.position = "bottom",
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

# create model to fit lines for plot
perim_mod <- glm(cbind(npos,nneg) ~ perim + tads, 
                   data = df,
                   family = "binomial") 

# create data for tadpoles present and absent 
xx1 <- seq(min(tads1$perim, na.rm = T), max(tads1$perim, na.rm = T), by = 1)
pred_df1 <- data.frame(perim = xx1, tads = 1)
pred_df1$pred <- predict(perim_mod, newdata = pred_df1, type = "response")
xx0 <- seq(min(tads0$perim, na.rm = T), max(tads0$perim, na.rm = T), by = 0.01)
pred_df0 <- data.frame(perim = xx0, tads = 0)
pred_df0$pred <- predict(perim_mod, newdata = pred_df0, type = "response")

# add lines to point plot
perim_p <- perim_p + 
  geom_line(pred_df1, mapping = aes(x = perim, y = pred),
            colour = "blue")+
  geom_line(pred_df0, mapping = aes(x = perim, y = pred),
            colour = "red")

# ------------ Arrange all into one figure

effect_fig <- ggarrange(dens_p + rremove("ylab") + font("xy.text", size = 10),
                        vol_p + rremove("ylab") + font("xy.text", size = 10),
                        perim_p + rremove("ylab")+ font("xy.text", size = 10), 
          nrow = 3, ncol = 1,
          align = "hv", labels = c("A", "C", "E"), hjust = 0.8, 
          font.label = list(size = 10, face = "plain"), common.legend = T, legend = "bottom") + 
  theme(plot.margin = margin(10, 5.5, 10, 5.5))

# common x-label
effect_fig <- annotate_figure(effect_fig, 
                left = textGrob("Proportion of qPCR replicates amplified", 
                                rot = 90, gp = gpar(fontsize = 10)))

effect_fig

#-------------------------------------------------------------------------------
#-------------------- Combine GLMER plots with histograms for figure 2
#-------------------------------------------------------------------------------


survCond <- ggarrange(effect_fig, all.hists, ncol = 2, nrow = 1, align = "hv")
survCond

ggsave(plot = survCond, "Survey conditions fig.jpeg", height = 17, width = 16, units = "cm", dpi = 600) 

p <- NCmisc::list.functions.in.file("C:/Users/s444224/OneDrive - University of Canberra - STAFF/Honours 2023/aaa_Analysis/R scripts_/2024_CaneToad_eDNA_Manuscript_PreliminaryAnalysis_Fig2.R")
summary(p)
