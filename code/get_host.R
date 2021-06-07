# setwd("~/GitHub/phage_e/")
library(here)
library(tidyverse)

d.phage <- read_csv(here("data/phage-gc-content.csv"), trim_ws = T)

#looking at all   non-siphos
d.phage <- d.phage %>% 
  filter ( Family != "Siphoviridae" )

#import virus-host data
vh.db <- read_tsv(here("data","virushostdb.tsv") ) 

# extracting from vh.db a mapping-table between refseq and taxid
# there are many viruses that have multiple refseq ids (=multiple sequences)
# and there are viruses that have multiple hosts (= multiple rows)
map_taxid_refseq <- vh.db %>% 
  select(virus.tax.id, refseq.id) %>% 
  # get all sequence IDs
  separate(refseq.id, into = as.character(1:200), sep = ",",fill = "right") %>% 
  pivot_longer(-virus.tax.id, names_to = "num", values_to = "refseq.id") %>% 
  filter(! is.na(refseq.id)) %>% 
  select(-num) %>% 
  # discard replicated ID pairs
  distinct() %>% 
  mutate(refseq.id = trimws(refseq.id))


# join data
d.phage.join <- d.phage %>% 
  # add taxid
  left_join(., map_taxid_refseq, by = c("Accession"="refseq.id")) %>% 
  # use taxid to add vh-db data
  left_join(., vh.db, by ="virus.tax.id") %>% 
  arrange(host.name,env.source)

write_csv(d.phage.join,
          here("data/phage-host-marine-nonsipho.csv"))

# tests
# # extra rows due to viruses with multiple hosts
# look <- d.phage.join %>%
#   filter(!is.na(virus.tax.id)) %>%
#   filter(duplicated(virus.tax.id))
# # reciprocal presence
# d.phage$Accession %in% d.phage.join$Accession
# d.phage.join$Accession %in% d.phage$Accession
# 
# look <- d.phage.join %>%
#   group_by(virus.tax.id) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n)) 
