library(here)
library(tidyverse)
library(cowplot)

readBlast <- function(hits)
{
  fields <- c("queryId", "subjectId", "percIdentity", "alnLength", "mismatchCount", "gapOpenCount", "queryStart", "queryEnd", "subjectStart", "subjectEnd", "eVal", "bitScore")
  
  #read in Blast results
  blast <- read_tsv (hits, col_names = F)
    names (blast) <-fields
  
  # #filter out hits shorter than 50
  # blast <- subset (blast, blast$alnLength>49)
  return (blast)
}

#### list sra accesion by study#####
studies <- list.files(here("data/recruitment/sra_acc/"), 
                pattern = ".txt", full.names = TRUE)
d.studies <- tibble()

for (f in studies){
  d.studies <- read_tsv(f, col_names = F) %>% 
    rename(acc=X1) %>% 
    mutate(study = str_remove(f, ".txt")) %>% 
    mutate(study = str_remove(study, ".*sra_acc/")) %>% 
    bind_rows(d.studies,.)
}

#### Collect blast results#####
blast_hits <- list.files(here("data/recruitment/tblastn/"), 
                         pattern = ".tsv", full.names = TRUE)

d <- tibble()

for (f in blast_hits){
  d <- readBlast(f) %>% 
    bind_rows(d,.)
}

# add study
d <- d %>% 
  separate(subjectId, into = c("acc", "read", "pair"), sep = "\\.", remove = F) %>% 
  left_join(., d.studies, by = "acc" )

# # make numerical gene column
# d <- d %>% 
  
# d %>% 
#   ggplot(aes(alnLength, percIdentity)) +
#   geom_point()

p <- d %>% 
  # slice_head(n = 10000) %>% 
  ggplot(aes(parse_number(queryId), percIdentity)) +
  geom_jitter(shape=21, alpha = 0.5)+
  facet_grid(study ~ .)+
  theme_classic()+
  panel_border(color = "black", size = 1)+
  xlab("phage E gene number")+
  ylab("%ID (tblastn)")+
  scale_x_continuous(breaks = seq(0,72,2))
  # theme(axis.text.x = element_text(angle = 90))

ggsave(here("plots/recruit-prelim.png"), plot = p, width = 11, height = 8)  

