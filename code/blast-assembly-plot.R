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

# #### list sra accesion by study#####
# studies <- list.files(here("data/assembly/tblastn"),
#                 pattern = ".tsv", full.names = TRUE)
# d.studies <- tibble()
# 
# for (f in studies){
#   d.studies <- read_tsv(f, col_names = F) %>%
#     rename(acc=X1) %>%
#     mutate(study = str_remove(f, ".txt")) %>%
#     mutate(study = str_remove(study, ".*sra_acc/")) %>%
#     bind_rows(d.studies,.)
# }
# 
# #focus on Kinneret
# d.studies <- d.studies %>% 
#   filter(str_detect(study, "Kinneret"))
#### Collect blast results#####
blast_hits <- list.files(path = here("data/assembly/tblastn/"), 
                         pattern = ".tsv", full.names = TRUE)

d <- tibble()

for (f in blast_hits){
  d <- readBlast(f) %>% 
    mutate(sra = str_extract(f,"SRR[0-9]*")) %>% 
    bind_rows(d,.)
}


#multi-hit scaffolds
scaf <- d %>% 
  group_by(subjectId) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>%
  filter(n>=4)



# p <-
d2 <-   d %>% 
  filter(subjectId %in% scaf$subjectId) %>%
  mutate(sID = str_remove(subjectId, "NODE_") %>% 
           str_remove( "length_") %>% 
           str_remove( "cov_")) %>% 
  separate(sID, into = c("node", "length", "coverage"), sep = "_")

l_node <- d2 %>% 
  group_by(node, sra, length) %>% 
  summarise() %>% 
  mutate(x = 1,y = 1,xend = as.numeric(length))


d2%>%
  ggplot() +
  geom_segment(data = l_node, aes(x = x, xend = xend, y=y, yend = y),
               color = "pink", size = 5)+
    geom_segment(aes(x = subjectStart, y = parse_number(queryId),
                     xend= subjectEnd, yend = parse_number(queryId),
                     color = percIdentity), size = 3)+

  facet_wrap(~sra + paste("node",node), strip.position = "right", ncol = 2)+
  theme_classic()+
  panel_border(color = "black", size = 1)+
  ylim(1,NA)+
  ylab("phage Cr-LKS3 gene number")+
  xlab("scafold position (bp)")+
  # scale_x_continuous(breaks = seq(0,72,2))+
  scale_color_viridis_b()+
  # theme(axis.text.x = element_text(angle = 90))

ggsave(here("plots/Dalton_scaffolds.png"),  width = 8, height = 6)  

#########################
# Export results
dir.create(here("data/assembly/Dalton_scaffolds"))

write_csv(d2, here("data/assembly/Dalton_scaffolds", "tblastn_hits_scaffolds.csv"))

# extract scaffold sequences
setwd(here("data/assembly/Dalton_scaffolds"))

#scaffolds to extract
scaf <- d2 %>% 
  select(subjectId, node, sra) %>% 
  distinct()
  left_join(scaf, .)
  


for (i in 1:nrow(scaf)){
  slate.path <- paste0("/N/slate/danschw/phage_e/lm/assembly/",scaf$sra[i],"/scaffolds.fasta")
  prl.cmd <- 
    paste0("perl -ne \'if(/^>(\\S+)/){$c=grep{/^$1$/}qw(", scaf$subjectId[i],
           ")}print if $c\' ",
          slate.path," > ",
          scaf$sra[i],"_",scaf$node[i], ".fasta")
  system(prl.cmd)
}

