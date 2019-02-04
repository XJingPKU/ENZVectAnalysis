#==========================================================
# NEECF vetor analysis
# 
# First created by Xin Jing, 1/4/2019
# Last modified:
#
# R version 3.5.0 (2018-04-23) -- "Joy in Playing"
# Contacts: Xin Jing <jingxin0123@gmail.com>
#==========================================================


rm(list = ls())

# load library
library(tidyverse)
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
  select(site2, treat, fb, fp, fd, vecL, vecA) %>% 
  filter(treat %in% c("CK", "N50", "N100")) %>% 
  droplevels()

# linear-mixed effect model
fitL.lmer <- lmer(vecL ~ site2 * fd * treat - site2 : fd :treat + 
                  (1|fb/fp),
            data = df)
anova(fitL.lmer, type = "I", ddf = "Kenward-Roger")
fitA.lmer <- lmer(vecA ~ site2 * fd * treat - site2 : treat : fd + 
                    (1|fb/fp),
             data = df)
anova(fitA.lmer, type = "I", ddf = "Kenward-Roger")

# ggplot
p.vecL <- df %>% 
  ggplot(aes(x = site2, y = vecL)) +
  geom_boxplot(outlier.color = "white") +
  geom_jitter(shape = 1, width = 0.15, color = "grey") +
  labs(x = "Site", y = "Microbial C limitation") +
  theme_classic(base_size = 14.5) +
  theme(panel.grid = element_blank(),
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
plot_grid(p.vecL, p.vecA, ncol = 1, labels = c("a)", "b)"))
ggsave("./outputs/VectorL&A_overview.pdf", width = 5.5, height = 7.0)


# Vector length VS degree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df %>% 
  group_by(fd) %>%
  summarise(n())

df %>% 
  split(.$fd) %>% 
  map(~ lm(vecA ~ vecL, data = .x)) %>% 
  map(summary)

ggplot(df, aes(x = vecL, y = vecA, color = fd)) +
  geom_point(shape = 1, alpha = 0.8) +
  geom_smooth(method = 'lm', se = FALSE, size = 0.3, 
              aes(lty = fd)) +
  scale_linetype_manual(values = c(1, 1, 1, 2),
                        guide = FALSE) +
  labs(x = "Microbial C limitation",
       y = "Microbial N/P limitation") +
  theme_classic(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.text = element_text(size = 9.5),
        legend.title = element_text(size = 11.5)) +
  guides(color = guide_legend(title = "Soil depth"))
ggsave("./outputs/microbial_CNP_limitation.pdf", width = 5.0, height = 3.5)


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

ggplot(df.JFL1, aes(x = fd, y = vecA, color = P)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector angle \n(Degrees)") +
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

ggplot(df.JFL2, aes(x = fd, y = vecA, color = N)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector angle \n(Degrees)") +
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
ggplot(df.WYS, aes(x = fd, y = vecA)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector angle \n(Degrees)") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

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
ggplot(df.GNJ, aes(x = fd, y = vecL, color = N)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector length \n(Unitless)") +
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
ggplot(df.GH, aes(x = fd, y = vecA)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Vector angle \n(Degrees)") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

sessionInfo()
###########################################################
#                  End of Script                          #
###########################################################