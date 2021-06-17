
library(here)
library(tidyverse)
library(cowplot)

d <- read_csv(here("data/gc-content-s2.csv"), trim_ws = T)
d.host <- read_csv(here("data/host-gc-range.csv"), trim_ws = T)
d <-d %>% 
  left_join(., d.host, by = c("Host" = "Host.genus") )

# p1 <- d %>%
#   mutate(`Host GC(%)` = as.numeric(`Host GC(%)`)) %>% 
#   filter(!is.na(`Host GC(%)`)) %>% 
#   mutate(label  = if_else(phage.name == "Cr-LKS3", "Cr-LKS3", "")) %>% 
#   ggplot(aes(`Host GC(%)`, `Phage GC(%)`))+
#   geom_abline(slope = 1, intercept = 5, linetype = 2, color = "grey")+
#   geom_abline(slope = 1, intercept = -5, linetype = 2, color = "grey")+
#   geom_abline(slope = 1, intercept = 0)+
#   geom_point(aes(fill = Environment, shape = phage.type), size = 2, stroke = 1)+
#   geom_text(aes(label = label), nudge_x = -1, nudge_y = -1)+
#   xlim(30,70)+
#   ylim(30,70)+
#   theme_cowplot()+panel_border(color = "black")+
#   scale_shape_manual(values = c(24,25, 21))+
#   scale_fill_grey()
#   
# 
# p1 + ggsave(here("plots/gc.png"), width = 6,height = 4)
# 
# 
# p2 <- d %>% 
#   filter(!is.na(Environment)) %>% 
#   mutate(label  = if_else(phage.name == "Cr-LKS3", "Cr-LKS3", "")) %>% 
#   ggplot(aes(`Host`, `Phage GC(%)`))+
#   geom_segment(aes(y = GC.low, yend = GC.high, xend = Host), color = "pink", size = 4)+
#     geom_jitter(aes(fill = Environment, size = `# tRNA genes`),shape=21, width = .2, height = 0)+
#   geom_text(aes(label = label), nudge_x = -0.3)+
#   # scale_fill_manual(values = c("marine"="white", "freshwater"="grey"))+  
#   theme_cowplot()+panel_border(color = "black")+
#   coord_flip()+
#   scale_fill_grey()
# 
# p2+
#   ggsave(here("plots/gc-phage.png"), width = 6,height = 4)
# 
# 
# plot_grid(p1,p2, labels = letters, ncol = 1)+
#   ggsave(here("plots/gc-phage-2panels.png"), width = 6,height = 6)

#-------------------------#


p3 <- d %>%
  mutate(`Host GC(%)` = as.numeric(`Host GC(%)`)) %>%
  filter(!is.na(`Host GC(%)`)) %>% 
  mutate(label  = if_else(phage.name == "Cr-LKS3", "Cr-LKS3", "")) %>% 
  ggplot(aes(`Host GC(%)`, `Phage GC(%)`))+
  geom_abline(slope = 1, intercept = 5, linetype = 2, color = "grey")+
  geom_abline(slope = 1, intercept = -5, linetype = 2, color = "grey")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(fill = Environment), size = 2, stroke = 1, shape = 21)+
  geom_text(aes(label = label), nudge_x = -1, nudge_y = -3.5)+
  xlim(30,70)+
  ylim(30,70)+
  theme_cowplot()+panel_border(color = "black")+
  facet_wrap(~phage.type)+
  scale_fill_grey()


# p3 + ggsave(here("plots/gc-panels.png"), width = 8,height = 4)
# 
# 
# p4 <- d %>% 
#   filter(!str_detect(phage.type,"Unch")) %>% 
#   mutate(label  = if_else(phage.name == "Cr-LKS3", "Cr-LKS3", "")) %>% 
#   ggplot(aes(`Host`, `Phage GC(%)`))+
#   geom_jitter(aes(fill = Environment, size = `# tRNA genes`),shape=21, width = .2, height = 0)+
#   geom_text(aes(label = label), nudge_x = -0.3, nudge_y = -1)+
#   theme_cowplot()+panel_border(color = "black")+
#   coord_flip()+
#   facet_wrap(~phage.type)+
#     scale_fill_grey()
# 
# p4+
#   ggsave(here("plots/gc-phage-panels.png"), width = 8,height = 4)
# 
# 
# 
# plot_grid(p3,p4, labels = letters, ncol = 1)+
#   ggsave(here("plots/gc-phage-6panels.png"), width = 8,height = 6)


p5 <- d %>% 
  filter(!str_detect(phage.type,"Unch")) %>% 
  mutate(label  = if_else(phage.name == "Cr-LKS3", "Cr-LKS3", "")) %>% 
  ggplot(aes(`Host`, `Phage GC(%)`))+
  geom_segment(aes(y = GC.low, yend = GC.high, xend = Host), color = "pink", size = 4)+
  geom_jitter(aes(fill = Environment, size = `# tRNA genes`),shape=21, width = .2, height = 0)+
  geom_text(aes(label = label), nudge_x = -0.3, nudge_y = -2)+
  theme_cowplot()+panel_border(color = "black")+
  coord_flip()+
  facet_wrap(~fct_rev(phage.type), ncol = 2, scales = "free_y", as.table = F)+
  scale_fill_grey()+
  theme(legend.position = "right",
        legend.box="horizontal",
        legend.margin=margin(2,5,2,5),
        # legend.background = element_blank(),
        # legend.box.background = element_rect(color = "grey90"))
  )
# p5

legend = get_legend(p5)

p <- plot_grid(p3+theme(legend.position = "none"), 
          # p5,
          p5+theme(legend.position = "none"),
          labels = letters, ncol = 1,
          rel_heights = c(3,5))+
  draw_grob(legend,.6, 0, 0, 1, vjust = 0)

  p + ggsave( here("plots/gc-phage-6panels-b.png"), width = 8,height = 8)




#export to pptx using officer and rvg
library (officer)
library(rvg)

read_pptx() %>%
  add_slide(layout = , master = "Office Theme") %>%
  ph_with(dml(ggobj = p), location = ph_location(type = "body",
                                                      left = 0, top = 0, width = 7, height = 7)) %>%
  print(target = here("plots/gc-plot.pptx"))

