library(Rmisc)
library(tidyverse)
library(stringr)
library(scales)
library(grid)
library(ggpubr)
library(MASS)
library(lmerTest)
library(lme4)
library(stats)
library(modelr)
library(plotrix)
library(mgcv)
library(hexbin)
library(formattable)

rm(list=ls())

d.sim1a <- read.csv("results_data/simulation1a.csv") %>% 
  mutate(Prime = 'DO Prime',
         Cond = ifelse(Cond=='same', 'Same', 'Different'))
d.sim1b <- read.csv("results_data/simulation1b.csv") %>% 
  mutate(Prime = 'PO Prime',
         Cond = ifelse(Cond=='same', 'Same', 'Different'))
d.sim1c.DO <- read.csv("results_data/simulation1c.csv") %>%
  mutate(Cond = 'Control',
         Prime = 'DO Prime')
d.sim1c.PO <- read.csv("results_data/simulation1c.csv") %>%
  mutate(Cond = 'Control',
         Prime = 'PO Prime')

d.sim2a.f1 <- read.csv("results_data/simulation2a_f1.csv") %>% 
  mutate(Prime = 'DO Prime', Fillers = '1 Batch',
         Cond = ifelse(Cond=='same', 'Same', 'Different'))
d.sim2b.f1 <- read.csv("results_data/simulation2b_f1.csv") %>% 
  mutate(Prime = 'PO Prime', Fillers = '1 Batch',
         Cond = ifelse(Cond=='same', 'Same', 'Different'))
d.sim2c.f1.DO <- read.csv("results_data/simulation2c_f1.csv") %>%
  mutate(Cond = 'Control', Fillers = '1 Batch', Prime = 'DO Prime')
d.sim2c.f1.PO <- read.csv("results_data/simulation2c_f1.csv") %>%
  mutate(Cond = 'Control', Fillers = '1 Batch', Prime = 'PO Prime')

d.sim2a.f2 <- read.csv("results_data/simulation2a_f2.csv") %>% 
  mutate(Prime = 'DO Prime', Fillers = '2 Batch',
         Cond = ifelse(Cond=='same', 'Same', 'Different'))
d.sim2b.f2 <- read.csv("results_data/simulation2b_f2.csv") %>% 
  mutate(Prime = 'PO Prime', Fillers = '2 Batch',
         Cond = ifelse(Cond=='same', 'Same', 'Different'))
d.sim2c.f2.DO <- read.csv("results_data/simulation2c_f2.csv") %>%
  mutate(Cond = 'Control', Fillers = '2 Batch', Prime = 'DO Prime')
d.sim2c.f2.PO <- read.csv("results_data/simulation2c_f2.csv") %>%
  mutate(Cond = 'Control', Fillers = '2 Batch', Prime = 'PO Prime')

d.sim1 <- d.sim1a %>% rbind(d.sim1b) %>% rbind(d.sim1c.DO) %>% rbind(d.sim1c.PO)
d.sim2 <- d.sim2a.f1 %>% rbind(d.sim2b.f1) %>% rbind(d.sim2c.f1.DO) %>% rbind(d.sim2c.f1.PO) %>%
  rbind(d.sim2a.f2) %>% rbind(d.sim2b.f2) %>% rbind(d.sim2c.f2.DO) %>% rbind(d.sim2c.f2.PO)


################# Plots #################

p.bar.sim1 <- d.sim1 %>% 
  mutate(Cond = ifelse(Cond=='Control', 'Prior', Cond),
         Cond = factor(Cond, levels = c('Prior', 'Different', 'Same'))) %>%
  summarySE(measurevar="P_DO", groupvars=c("Prime","Cond")) %>%
  ggplot(aes(x=Cond, y=P_DO, fill=Cond)) +
  facet_grid(~Prime) +
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(P_DO, 3)),
            colour = "white", size = 3.5,
            position=position_stack(vjust=0.5)) +
  scale_fill_manual(values = c("Gray90","Gray65", "Gray40")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=10), 
        legend.position="none", 
        axis.line=element_line(color="black"), axis.title.x = element_blank()) +
  ylab("probability of DO")
p.bar.sim1

p.EffectSize.sim1 <- d.sim1 %>% 
  mutate(Cond = ifelse(Cond=='Control', 'Prior', Cond),
         Cond = factor(Cond, levels = c('Prior', 'Different', 'Same'))) %>%
  summarySE(measurevar="P_DO", groupvars=c("Prime","Item","Cond")) %>%
  dplyr::select(Item, Prime, Cond, P_DO) %>%
  reshape(idvar=c("Item", "Prime"), timevar="Cond", direction="wide") %>%
  mutate(Different = abs(log(P_DO.Different/(1-P_DO.Different)) - log(P_DO.Prior/(1-P_DO.Prior))),
         Same = abs(log(P_DO.Same/(1-P_DO.Same)) - log(P_DO.Prior/(1-P_DO.Prior)))) %>%
  dplyr::select(Item, Prime, Different, Same) %>%
  gather(Cond, EffectSize, Different:Same, factor_key=TRUE) %>%
  summarySE(measurevar="EffectSize", groupvars=c("Prime","Cond")) %>%
  ggplot(aes(x=Cond, y=EffectSize, fill=Cond)) +
  facet_grid(~Prime) +
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(EffectSize, 3)),
            colour = "white", size = 3.5,
            position=position_stack(vjust=0.5)) +
  scale_fill_manual(values = c("Gray65", "Gray40")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=10), 
        legend.position="none", 
        axis.title.x = element_blank(),
        axis.line=element_line(color="black")) +
  ylab("effect size (log odds)")
p.EffectSize.sim1

p.EffectSize.sim2 <- d.sim1 %>%
  mutate(Fillers = '0 Batch') %>%
  rbind(d.sim2) %>% 
  summarySE(measurevar="P_DO", groupvars=c("Prime","Item","Cond","Fillers")) %>%
  dplyr::select(Item, Prime, Cond, Fillers, P_DO) %>%
  reshape(idvar=c("Item", "Prime", "Fillers"), timevar="Cond", direction="wide") %>%
  mutate(Different = abs(log(P_DO.Different/(1-P_DO.Different)) - log(P_DO.Control/(1-P_DO.Control))),
         Same = abs(log(P_DO.Same/(1-P_DO.Same)) - log(P_DO.Control/(1-P_DO.Control)))) %>%
  dplyr::select(Item, Prime, Fillers, Different, Same) %>%
  gather(Cond, EffectSize, Different:Same, factor_key=TRUE) %>%
  summarySE(measurevar="EffectSize", groupvars=c("Prime","Cond","Fillers")) %>%
  ggplot(aes(x=Fillers, y=EffectSize, group=Cond)) +
  facet_grid(~Prime) +
  geom_text(aes(label=round(EffectSize, 3)), nudge_y=0.015) +
  geom_point(aes(shape=Cond), size=I(2.5)) +
  geom_line(aes(linetype=Cond),size=I(0.6)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=10), 
        legend.position="bottom", 
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.line=element_line(color="black")) +
  ylab("priming effect size (log odds)")
p.EffectSize.sim2

p.CogSci.sim1 <- ggarrange(p.bar.sim1, p.EffectSize.sim1,
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1, widths=c(2.5, 2), heights = c(1, 1))
pdf("plots/plot.CogSci.sim1.pdf", width = 8, height = 2)
p.CogSci.sim1
dev.off()

pdf("plots/plot.CogSci.sim2.pdf", width = 4.4, height = 3.2)
p.EffectSize.sim2
dev.off()






