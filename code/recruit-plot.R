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


blast_hits <- list.files(here("data/recruitment/tblastn/"), 
                         pattern = ".tsv", full.names = TRUE)

d <- tibble()

for (f in blast_hits){
  d <- readBlast(f) %>% 
    bind_rows(d,.)
}

d %>% 
  ggplot(aes(alnLength, percIdentity)) +
  geom_point()


x <- list.files(here("data/recruitment/sra_acc/"), 
                pattern = ".txt", full.names = TRUE)
d.acc <- tibble()

for (f in x){
   d.acc <- read_tsv(f, col_names = F) %>% 
    mutate(study = str_remove(f, ".txt")) %>% 
     mutate(study = str_remove(study, ".*sra_acc/")) %>% 
     bind_rows(d.acc,.)
}
