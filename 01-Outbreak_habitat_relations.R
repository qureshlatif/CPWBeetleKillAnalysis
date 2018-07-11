require(dplyr)
require(stringr)
require(ggplot2)
require(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

dat.LP <- Cov.LP %>% tbl_df %>%
  mutate(Point = row.names(Cov.LP)) %>%
  select(Point, gridIndex:Rd_dens1km)

dat.SF <- Cov.SF %>% tbl_df %>%
  mutate(Point = row.names(Cov.SF)) %>%
  select(Point, gridIndex:Rd_dens1km)

#### Scatterplots - DeadConif VS YSO ####
p.LP <- ggplot(data = dat.LP, aes(x = YSO, y = DeadConif)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 1.05, label = "Lodgpole pine stratum", size = 4)

p.SF <- ggplot(data = dat.SF, aes(x = YSO, y = DeadConif)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 1.05, label = "Spruce-fir stratum", size = 4)

p <- ggdraw() + 
  draw_plot(p.LP, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p.SF, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Proportion conifer dead", "Years since outbreak"),
                  x = c(0.03, 0.4),
                  y = c(0.2, 0.03),
                  size = c(15, 15),
                  angle = c(90, 0),
                  hjust = c(0, 0),
                  vjust = c(0, 0))

save_plot("figure_DeadConif_VS_YSO.tiff", p, ncol = 2, nrow = 1, dpi = 200)

## Compare DeadConif between outbreak vs non-outbreak grids ##
#dat.LP %>%
#  group_by(gridIndex) %>%
#  summarise(DConif = mean(DeadConif, na.rm = T), DConif_max = max(DeadConif, na.rm = T), YSI = mean(YSO, na.rm = T)) %>%
#  ungroup %>%
#  mutate(Outbreak = !is.na(YSI)) %>%
#  group_by(Outbreak) %>%
#  summarise(DConif_mn = mean(DConif, na.rm = T), DConif_SD = sd(DConif, na.rm = T), DConif_max = max(DConif, na.rm = T)) %>%
#  View

#dat.SF %>%
#  group_by(gridIndex) %>%
#  summarise(DConif = mean(DeadConif, na.rm = T), DConif_max = max(DeadConif, na.rm = T), YSI = mean(YSO, na.rm = T)) %>%
#  ungroup %>%
#  mutate(Outbreak = !is.na(YSI)) %>%
#  group_by(Outbreak) %>%
#  summarise(DConif_mn = mean(DConif, na.rm = T), DConif_SD = sd(DConif, na.rm = T), DConif_max = max(DConif, na.rm = T)) %>%
#  View

#### Vegetation (mehanistic) VS outbreak metrics ####
r <- cor(dat.LP$YSI, dat.LP$CanCov, use = "complete") %>% round(digits = 3)
p. <- ggplot(dat.LP, aes(x = YSI, y = CanCov)) +
  geom_point(alpha = 0.3, size = 5) +
  geom_smooth() +
  xlab("Years since outbreak") + ylab("Canopy Cover") +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  annotate("text", x = 15, y = 60, label = str_c("r = ", r), size = 10)

p <- ggplot(dat.LP, aes(x = DeadConif_RCov, y = CanCov)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
