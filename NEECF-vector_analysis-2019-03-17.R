#==========================================================
# NEECF vetor analysis
# 
# First created by Xin Jing, 1/4/2019
# Last modified: 3/17/2019
#
# R version 3.5.0 (2018-04-23) -- "Joy in Playing"
# Contacts: Xin Jing <jingxin0123@gmail.com>
#==========================================================


rm(list = ls())

# load library
library(plyr)
library(tidyverse)
# library(devtools)
# install_github("pascal-niklaus/pascal/pascal")
library(pascal)
library(lme4)
library(lmerTest)
library(cowplot)


#==========================================================
# read data
paths <- "D:/My project/enzyme/Vector_analysis/data/NEECF_enzyme_data.csv"
enz.neecf <- read.csv(paths)  

# clean the data
enz.neecf.clean <- enz.neecf %>%
  select(site, site2, latitude, MAT, MAP, forest.type, 
         treat, N, P, new.block, new.plot, soil.layer, 
         STC, STN, STP, pH, SM, BG, CB, NAG, LAP, PHOS) %>%
  filter(soil.layer != "O-horizon") %>%
  mutate(site = factor(site, 
                       levels = c("JFL", "WYS", "GNJ", "DLS", "WY", "GH")),
         site2 = factor(site2, 
                        levels = c("JFL1", "JFL2", "WYS", "GNJ", "DLS1", "DLS2", "WY", "GH")),
         lat = latitude,
         fb = factor(new.block),
         fp = factor(new.plot),
         fd = factor(soil.layer, levels = c("0-10cm", "10-20cm", "20-40cm", "40-60cm")),
         AP = PHOS) %>%
  select(site, site2, lat, MAT, MAP, treat, N, P, fb, fp, fd, STC, STN, STP, pH, SM,
         BG, CB, NAG, LAP, AP) %>% 
  droplevels()

# convert the negative values of soil enzyme activity into zero
enz.vars <- c("BG", "CB", "NAG", "LAP", "AP")
enz.neecf.clean[names(enz.neecf.clean) 
                %in% enz.vars][enz.neecf.clean[names(enz.neecf.clean) 
                                               %in% enz.vars] < 0] <- 0

# convert "NaN" into "NA"
enz.neecf.clean[mapply(is.nan, enz.neecf.clean)] <- NA


#==========================================================
# calculate vector length and angle
enz.neecf.clean <- enz.neecf.clean %>% 
  mutate(y = (BG + CB)  / (BG + CB + NAG + LAP),
         x = (BG + CB) / (BG + CB + AP)) %>%
  mutate(vecL = sqrt(x^2 + y^2),
         vecA = atan2(y, x)*180/pi) %>% 
  na.omit() %>% 
  filter(vecL != 0 & vecA != 0)


#==========================================================
# Effects of treatment, site and depth on vector length and angle

# select variables
df <- enz.neecf.clean %>% 
  select(site2, lat, MAT, MAP, pH, SM, STC, STN, STP, treat, 
         fb, fp, fd, vecL, vecA) %>% 
  filter(treat %in% c("CK", "N50", "N100")) %>% 
  droplevels()

# linear-mixed effect model ---
# by vector lengths
fitL.lmer <- lmer(vecL ~ site2 * fd * treat - site2 : fd :treat + 
                    (1|fb/fp),
                  data = df)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")
# summarise site means of vecL
df %>% 
  group_by(site2) %>% 
  dplyr::summarise(., n = n(),
            mu = mean(vecL),
            se = sd(vecL) / sqrt(n()))

# by vector angles
fitA.lmer <- lmer(vecA ~ site2 * fd * treat - site2 : treat : fd + 
                    (1|fb/fp),
                  data = df)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")
df %>% 
  group_by(site2, fd) %>% 
  dplyr::summarise(., n = n(),
            mu = mean(vecA),
            se = sd(vecL) / sqrt(n()))

# ggplot
p.vecL <- df %>% 
  ggplot(aes(x = site2, y = vecL, color = fd)) +
  geom_boxplot() +
  # geom_jitter(shape = 1, width = 0.15, color = "grey") +
  labs(x = "Site", y = "Microbial C limitation") +
  theme_classic(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 0.5))
p.vecA <- df %>% 
  ggplot(aes(x = site2, y = vecA, color = fd)) +
  geom_boxplot() +
  geom_hline(yintercept = 45, lty = 2, color = "gray70") +
  labs(x = "Site", y = "Microbial N/P limitation") +
  theme_classic(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.20, 0.35),
        legend.text = element_text(size = 10.5),
        legend.title = element_text(size = 12.5),
        axis.line = element_line(size = 0.5)) +
  guides(color = guide_legend("Soil depth"))
plot_grid(p.vecL, p.vecA, ncol = 1, labels = c("a)", "b)"), align = "v")
ggsave("./outputs/VectorL&A_overview.pdf", width = 5.5, height = 7.0)


# Vector length VS degree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df %>% 
  group_by(fd) %>%
  dplyr::summarise(n())

df %>% 
  split(.$fd) %>% 
  map(~ lm(vecA ~ vecL, data = .x)) %>% 
  map(summary)

ggplot(df, aes(x = vecL, y = vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3, 
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  labs(x = "Microbial C limitation",
       y = "Microbial N/P limitation") +
  facet_grid(~ fd) +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')
ggsave("./outputs/microbial_CNP_limitation.pdf", width = 7.0, height = 3.5)


# Vector length VS degree across latitude ~~~~~~~~~~~~~~~~~
# aov.ftest to corect F-test
ylim.vecL <- ddply(df, .(fd), function(x) {
    max(x$vecL, na.rm = T) * 1.20
  }) %>% 
  mutate(ypos = max(V1)) %>% 
  ungroup() %>% 
  data.frame()
ylim.vecA <- ddply(df, .(fd), function(x) {
  max(x$vecA, na.rm = T) * 1.20
}) %>% 
  mutate(ypos = max(V1)) %>% 
  ungroup() %>% 
  data.frame()

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ lat + site2, data = x)
  aov.ftest(mod, lat ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                               ddply(df, .(fd), function(x) {
                                 fit <- lm(vecL ~ lat, data = x)
                                 summary(fit)$adj.r.squared
                               })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(35, 4)
p.vecL.lat <- df %>% 
  ggplot(aes(lat, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  facet_grid(~ fd) +
  geom_text(data = enz.df.aovftest, 
            aes(x = xpos, y = ylim.vecL$ypos, label = labels, 
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  labs(x = expression("Latitude ("*degree*N*")"), 
       y = "Soil microbial C limitation") +
  ylim(c(0, 1.65)) +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ lat + site2, data = x)
  aov.ftest(mod, lat ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ lat, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(35, 4)
p.vecA.lat <- df %>% 
  ggplot(aes(lat, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  facet_grid(~ fd) +
  geom_text(data = enz.df.aovftest, 
            aes(x = xpos, y = ylim.vecA$ypos, label = labels, 
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  labs(x = expression("Latitude ("*degree*N*")"), 
       y = "Soil microbial N/P limitation") +
  ylim(c(20, 110)) +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.lat, p.vecA.lat, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_latitude.pdf", width = 8.0, height = 7.0)


# Vector length VS degree across climate gradient ~~~~~~~~~

# by MAT
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ MAT + site2, data = x)
  aov.ftest(mod, MAT ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ MAT, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(8, 4)
p.vecL.mat <- df %>% 
  ggplot(aes(MAT, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  facet_grid(~ fd) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  labs(x = expression("Mean annual temperature ("*degree*C*")"), 
       y = "Soil microbial C limitation") +
  ylim(c(0, 1.65)) +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ MAT + site2, data = x)
  aov.ftest(mod, MAT ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ MAT, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(8, 4)
p.vecA.mat <- df %>% 
  ggplot(aes(MAT, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  facet_grid(~ fd) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  labs(x = expression("Mean annual temperature ("*degree*C*")"), 
       y = "Soil microbial N/P limitation") +
  ylim(c(20, 110)) +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.mat, p.vecA.mat, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_MAT.pdf", width = 8.0, height = 7.0)


# by MAP
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ MAP + site2, data = x)
  aov.ftest(mod, MAP ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ MAP, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(1400, 4)
p.vecL.map <- df %>% 
  ggplot(aes(MAP, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(0, 1.65)) +
  facet_grid(~ fd) +
  labs(x = "Mean annual precipitation (mm)", 
       y = "Soil microbial C limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ MAP + site2, data = x)
  aov.ftest(mod, MAP ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ MAP, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(1400, 4)
p.vecA.map <- df %>% 
  ggplot(aes(MAP, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 1),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(20, 110)) +
  facet_grid(~ fd) +
  labs(x = "Mean annual precipitation (mm)", 
       y = "Soil microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.map, p.vecA.map, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_MAP.pdf", width = 8.0, height = 7.0)


# Vector length VS degree across soil gradient ~~~~~~~~~~~~

# by STC
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ STC + site2, data = x)
  aov.ftest(mod, STC ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ STC, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(160, 4)
p.vecL.stc <- df %>% 
  ggplot(aes(STC, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(2, 2, 1, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(0, 1.65)) +
  facet_grid(~ fd) +
  labs(x = expression("Soil total C (mg g"^-1*")"), 
       y = "Soil microbial C limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ STC + site2, data = x)
  aov.ftest(mod, STC ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ STC, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(160, 4)
p.vecA.stc <- df %>% 
  ggplot(aes(STC, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(2, 2, 2, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(20, 110)) +
  facet_grid(~ fd) +
  labs(x = expression("Soil total C (mg g"^-1*")"), 
       y = "Soil microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.stc, p.vecA.stc, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_STC.pdf", width = 8.0, height = 7.0)


# by STN
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ STN + site2, data = x)
  aov.ftest(mod, STN ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ STN, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(7, 4)
p.vecL.stn <- df %>% 
  ggplot(aes(STN, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(2, 1, 1, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(0, 1.65)) +
  facet_grid(~ fd) +
  labs(x = expression("Soil total N (mg g"^-1*")"), 
       y = "Soil microbial C limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ STN + site2, data = x)
  aov.ftest(mod, STN ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ STN, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(7, 4)
p.vecA.stn <- df %>% 
  ggplot(aes(STN, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(2, 2, 2, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(20, 110)) +
  facet_grid(~ fd) +
  labs(x = expression("Soil total N (mg g"^-1*")"), 
       y = "Soil microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.stn, p.vecA.stn, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_STN.pdf", width = 8.0, height = 7.0)

# by STP
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ STP + site2, data = x)
  aov.ftest(mod, STP ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ STP, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(0.7, 4)
p.vecL.stp <- df %>% 
  ggplot(aes(STP, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(0, 1.65)) +
  facet_grid(~ fd) +
  labs(x = expression("Soil total P (mg g"^-1*")"), 
       y = "Soil microbial C limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ STP + site2, data = x)
  aov.ftest(mod, STP ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ STP, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(0.7, 4)
p.vecA.stp <- df %>% 
  ggplot(aes(STP, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(2, 2, 2, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(20, 110)) +
  facet_grid(~ fd) +
  labs(x = expression("Soil total P (mg g"^-1*")"), 
       y = "Soil microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.stp, p.vecA.stp, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_STP.pdf", width = 8.0, height = 7.0)

# by pH
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ pH + site2, data = x)
  aov.ftest(mod, pH ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ pH, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(5.4, 4)
p.vecL.ph <- df %>% 
  ggplot(aes(pH, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 2, 1, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(0, 1.65)) +
  facet_grid(~ fd) +
  labs(x = "Soil pH (unitless)", 
       y = "Soil microbial C limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ pH + site2, data = x)
  aov.ftest(mod, pH ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ pH, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(5.4, 4)
p.vecA.ph <- df %>% 
  ggplot(aes(pH, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 1),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(20, 110)) +
  facet_grid(~ fd) +
  labs(x = "Soil pH (unitless)", 
       y = "Soil microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.ph, p.vecA.ph, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_pH.pdf", width = 8.0, height = 7.0)

# by soil moisture
enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecL ~ SM + site2, data = x)
  aov.ftest(mod, SM ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecL ~ SM, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(40, 4)
p.vecL.sm <- df %>% 
  ggplot(aes(SM, vecL, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 2, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecL$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(0, 1.65)) +
  facet_grid(~ fd) +
  labs(x = "Soil moisture (%)", 
       y = "Soil microbial C limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

enz.df.aovftest <- ddply(df, .(fd), function(x) {
  mod <- aov(vecA ~ SM + site2, data = x)
  aov.ftest(mod, SM ~ site2, table = TRUE)[c(7, 8)]
})
enz.df.aovftest <- cbind(enz.df.aovftest,
                         ddply(df, .(fd), function(x) {
                           fit <- lm(vecA ~ SM, data = x)
                           summary(fit)$adj.r.squared
                         })$V1)
names(enz.df.aovftest)[4] <- "R"
enz.df.aovftest$labels <- with(enz.df.aovftest,
                               paste("R2 = ", round(R, 2), "; P =", P))
enz.df.aovftest$xpos <- rep(40, 4)
p.vecA.sm <- df %>% 
  ggplot(aes(SM, vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = TRUE, size = 0.3,
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 2, 2, 2),
                        guide = FALSE) +
  geom_text(data = enz.df.aovftest,
            aes(x = xpos, y = ylim.vecA$ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  ylim(c(20, 110)) +
  facet_grid(~ fd) +
  labs(x = "Soil moisture (%)", 
       y = "Soil microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none')

plot_grid(p.vecL.sm, p.vecA.sm, ncol = 1, align = "v")
ggsave("./outputs/VectorL&A_SM.pdf", width = 8.0, height = 7.0)

###########################################################
# site by site ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Jianfengling
df.JFL1 <- enz.neecf.clean %>% 
  select(site2, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site2 == "JFL1") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ N * P * fd - N:P:fd + (1|fb/fp),
                  data = df.JFL1)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ N * P * fd - N:P:fd + (1|fb/fp),
                  data = df.JFL1)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")
df.JFL1 %>% 
  group_by(fd, P) %>% 
  dplyr::summarise(., n = n(),
            mu = mean(vecA),
            se = sd(vecL) / sqrt(n()))

p.JFL1 <- ggplot(df.JFL1, aes(x = fd, y = vecA, color = P)) +
  geom_boxplot() +
  labs(title = "JFL1", x = "Depth", y = "Soil microbial N/P limitation") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

df.JFL2 <- enz.neecf.clean %>% 
  select(site2, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site2 == "JFL2") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ N * P * fd - N:P:fd + (1|fb/fp),
                  data = df.JFL2)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ N * P * fd - N:P:fd + (1|fb/fp),
                  data = df.JFL2)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")

p.JFL2 <- ggplot(df.JFL2, aes(x = fd, y = vecA, color = N)) +
  geom_boxplot() +
  labs(title = "JFL2", x = "Depth", y = "Soil microbial N/P limitation") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

# Wuyishan
df.WYS <- enz.neecf.clean %>% 
  select(site, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site == "WYS") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ treat * fd + (1|fb/fp),
                  data = df.WYS)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ treat * fd + (1|fb/fp),
                  data = df.WYS)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")

# Guniujiang
df.GNJ <- enz.neecf.clean %>% 
  select(site, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site == "GNJ") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ N + P + fd +
                    N:fd + P:fd + (1|fb/fp),
                  data = df.GNJ)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")
df.GNJ %>% 
  group_by(fd, N) %>% 
  dplyr::summarise(n = n(),
            mu = mean(vecL),
            se = sd(vecL) / sqrt(n()))
p.GNJ <- ggplot(df.GNJ, aes(x = fd, y = vecL, color = N)) +
  geom_boxplot() +
  labs(title = "GNJ", x = "Depth", y = "Soil microbial C limitation") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

fitA.lmer <- lmer(vecA ~ N + P + fd +
                    N:fd + P:fd + (1|fb/fp),
                  data = df.GNJ)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")
ggplot(df.GNJ, aes(x = fd, y = vecA)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector angle \n(Degrees)") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

# Donglingshan
df.DLS <- enz.neecf.clean %>% 
  select(site2, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site2 == "DLS1") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ treat * fd + (1|fb/fp),
                  data = df.DLS)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ treat * fd + (1|fb/fp),
                  data = df.DLS)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")

df.DLS <- enz.neecf.clean %>% 
  select(site2, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site2 == "DLS2") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ treat * fd + (1|fb/fp),
                  data = df.DLS)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ treat * fd + (1|fb/fp),
                  data = df.DLS)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")

# Wuying
df.WY <- enz.neecf.clean %>% 
  select(site, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site == "WY") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ treat * fd + (1|fb/fp),
                  data = df.WY)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ treat * fd + (1|fb/fp),
                  data = df.WY)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")
ggplot(df.WY, aes(x = fd, y = vecA)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector angle \n(Degrees)") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

# Genhe
df.GH <- enz.neecf.clean %>% 
  select(site, treat, N, P, fb, fp, fd, vecL, vecA) %>% 
  filter(site == "GH") %>% 
  mutate(N = factor(N),
         P = factor(P)) %>% 
  droplevels()
fitL.lmer <- lmer(vecL ~ treat * fd + (1|fb/fp),
                  data = df.GH)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")

fitA.lmer <- lmer(vecA ~ treat * fd + (1|fb/fp),
                  data = df.GH)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")

plot_grid(p.JFL1, p.JFL2, p.GNJ, ncol = 2, align = "v")
ggsave("./outputs/VectorL&A_sig.pdf", width = 10.5, height = 6.5)

sessionInfo()
###########################################################
#                  End of Script                          #
###########################################################