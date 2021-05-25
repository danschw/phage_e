
library(here)
library(tidyverse)

d <- read_csv(here("data/gc-content.csv"), trim_ws = T)

d %>% 
  mutate(`Host GC(%)` = as.numeric(`Host GC(%)`)) %>% 
  ggplot(aes(`Phage GC(%)`, `Host GC(%)`))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(color = Environment), size = 2)+
  geom_text(aes(label = label), nudge_x = -1, nudge_y = -1)+
  xlim(30,70)+
  ylim(30,70)+
  theme_classic()+
  ggsave(here("plots/gc.png"), width = 5,height = 4)


d %>% 
  ggplot(aes(`Host (genus)`, `Phage GC(%)`))+
    geom_jitter(aes(fill = Environment, size = `# tRNA genes`),shape=21, width = .2, height = 0)+
  geom_text(aes(label = label), nudge_x = -0.1)+
  scale_fill_manual(values = c("marine"="white", "freshwater"="grey"))+  
  theme_bw()+
  coord_flip()+
  ggsave(here("plots/gc-phage.png"), width = 6,height = 4)
