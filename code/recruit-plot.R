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

#focus on Kinneret
d.studies <- d.studies %>% 
  filter(str_detect(study, "Kinneret"))
#### Collect blast results#####
blast_hits <- list.files(path = here("data/recruitment/tblastn-Kinneret/"), 
                         pattern = ".tsv", full.names = TRUE)

d <- tibble()

for (f in blast_hits){
  d <- readBlast(f) %>% 
    mutate(sra = str_extract(f,"SRR[0-9]*")) %>% 
    bind_rows(d,.)
}

# add study
d <- d %>% 
  separate(subjectId, into = c("acc", "read", "pair"), sep = "\\.", remove = F) %>% 
  left_join(., d.studies, by = "acc" )

#sra meta data
d.meta <- read_csv(here("data/recruitment/Kinneret_SraRunTable.csv")) 

d <- d.meta %>% 
  select(Run, collection_date, geo_loc_name, isolation_source) %>% 
  left_join(d, ., by = c("sra"="Run")) 

# # make numerical gene column
# d <- d %>% 
  
# d %>%
#   ggplot(aes(alnLength, percIdentity)) +
#   geom_point()

p <- d %>% 
  # filter(sra == "SRR8088701") %>% 
  mutate(sra = str_remove(sra, "SRR")) %>% 
  ggplot(aes(parse_number(queryId), percIdentity)) +
  # geom_jitter(shape=21, alpha = 0.5)+
  geom_bin2d(bins = 72)+
  facet_wrap(~sra, strip.position = "right", ncol = 2 )+
  theme_classic()+
  panel_border(color = "black", size = 1)+
  xlab("phage E gene number")+
  ylab("%ID (tblastn)")+
  scale_x_continuous(breaks = seq(0,72,2))+
  scale_fill_viridis_b()
  # theme(axis.text.x = element_text(angle = 90))

ggsave(here("plots/recruit-Kinneret-prelim.png"), plot = p, width = 11, height = 8.5)  

p <- d %>% 
  filter(str_detect(geo_loc_name, "Dalton")) %>% 
  # filter(sra == "SRR8088692"|sra == "SRR8088693") %>%
  mutate(isolation_source = str_replace(isolation_source, "micro meter filter size", "μm")) %>%
  ggplot(aes(parse_number(queryId), percIdentity)) +
  # geom_jitter(shape=21, alpha = 0.5)+
  geom_bin2d(bins = 72)+
  facet_wrap(~sra + collection_date + isolation_source , strip.position = "right", ncol = 1 )+
  theme_classic()+
  panel_border(color = "black", size = 1)+
  xlab("phage E gene number")+
  ylab("%ID (tblastn)")+
  scale_x_continuous(breaks = seq(0,72,2))+
  scale_fill_viridis_b()
# theme(axis.text.x = element_text(angle = 90))

ggsave(here("plots/recruit-Dalton-prelim.png"), plot = p, width = 11, height = 8.5)  


################################
# better plots
library(seqinr)
faa <- read.fasta(here("data/recruitment/phageE_proteins.faa"))
d.faa <- tibble(gene = getName(faa),
                aa.length = getLength(faa),
                gene.start = NA,
                gene.end = NA,
                gene.middle = NA)


# make a concatenated protein "genome"

# assign first row
d.faa$gene.start[1] <- 1
d.faa$gene.end[1] <- d.faa$aa.length[1]
d.faa$gene.middle[1] <- 0.5*d.faa$aa.length[1]

for (i in 2:nrow(d.faa)){
  d.faa$gene.start[i] <- d.faa$gene.end[i-1]+1
  d.faa$gene.end[i] <- d.faa$gene.start[i]+d.faa$aa.length[i]-1
  d.faa$gene.middle[i] <- d.faa$gene.start[i]+0.5*d.faa$aa.length[i]
}

d <- d %>% 
  # fix error on name of gene_64 for compatbility
  mutate(queryId = str_replace(queryId, "gene_64.*", "gene_64")) %>% 
  left_join(., d.faa, by = c("queryId" = "gene"))

d <- d %>% 
  mutate(aln.start = gene.start + queryStart -1) %>% 
  mutate(aln.end = gene.start + queryEnd -1)

#background
d.bg <- 
  d.faa %>% 
  slice(which(row_number() %% 2 == 1))

p <- d %>% 
  mutate(sra = str_sub(sra, -3,-1)) %>% 
  mutate(geo_loc_name = str_remove(geo_loc_name, "Israel: ")) %>%
  ggplot() +
  geom_rect(data = d.bg, aes(xmin = gene.start, xmax = gene.end,
            ymin = 29, ymax = 100), fill = "grey90")+
  geom_segment(aes(x=aln.start,y=percIdentity, xend = aln.end, 
                   yend = percIdentity, color = geo_loc_name))+
  facet_wrap(~sra , strip.position = "right", ncol = 2 )+
  theme_classic(base_size = 10)+
  panel_border(color = "black", size = 1)+
  xlab("phage E concatenated proteins (gene #)")+
  ylab("%ID (tblastn)")+
  scale_x_continuous(breaks = d.faa$gene.middle, 
                     labels = str_remove(d.faa$gene, "gene_"), expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), breaks = c(50,100))+
  scale_color_brewer(type = "qual", palette = "Dark2")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(here("plots/recruit-IsraPonds.png"), plot = p, width = 11, height = 8.5) 

p <- d %>% 
  # filter(str_detect(geo_loc_name, "Dalton")) %>% 
  filter(sra == "SRR8088692"|sra == "SRR8088693") %>%
  mutate(isolation_source = str_replace(isolation_source, "micro meter filter size", "μm")) %>%
  ggplot() +
  geom_rect(data = d.bg, aes(xmin = gene.start, xmax = gene.end,
                             ymin = 29, ymax = 100), fill = "grey90")+
  geom_segment(aes(x=aln.start,y=percIdentity, xend = aln.end, yend = percIdentity), color = "blue")+
  facet_wrap(~sra + collection_date + isolation_source , strip.position = "right", ncol = 1 )+
  theme_classic()+
  panel_border(color = "black", size = 1)+
  xlab("phage E concatenated proteins (gene #)")+
  ylab("%ID (tblastn)")+
  scale_x_continuous(breaks = d.faa$gene.middle, 
                     labels = str_remove(d.faa$gene, "gene_"), expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), breaks = c(50,100))+
  scale_color_brewer(type = "qual", palette = "Dark2")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(here("plots/recruit-Dalton.png"), plot = p, width = 11, height = 6) 
