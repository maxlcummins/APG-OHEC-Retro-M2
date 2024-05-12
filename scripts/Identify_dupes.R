#Read in our packages
library(tidyverse)
library(vroom)
library(ggplot2)
library(xml2)
library(ggrepel)
library(caret)
library(here)

#Define our not in function
`%nin%` <- Negate(`%in%`)
genometa <- vroom("delims/genometa_n5631.txt", show_col_types = FALSE)

#Read in our tree
tree <- ggtree::read.tree("analysis/pad.nwk")

#Remove quotes from sample names
tree$tip.label <- gsub("'", "", tree$tip.label)

#Create a list of samples in our cohort
sample_names <- tree$tip.label %>% as.data.frame() %>% rename('name' = '.')


aussie_apec <- genometa %>% filter(Lab_Contact == "University of Technology, Sydney") %>% pull(name)

ANU_dupes <- genometa %>% filter(grepl("Australian National", Lab_Contact)) %>% pull(name)


potential_others <- cgMLST_dists2 %>% filter(`Sample 1` %nin% c(ANU_dupes, aussie_apec), `Sample 2` %nin% c(ANU_dupes, aussie_apec)) %>% filter(value == 0) %>% filter(`Sample 1` != `Sample 2`) %>% select(`Sample 1`, Revised_Source_Niche.x, value, Revised_Source_Niche.y, `Sample 2`)

potential_dupes <- potential_others %>% filter(grepl("SRR|ERR", `Sample 1`)) %>% filter(!grepl("SRR|ERR", `Sample 2`))
potential_dupes <- potential_others %>% filter(grepl("SRR|ERR", `Sample 2`)) %>% filter(!grepl("SRR|ERR", `Sample 1`)) %>% bind_rows(potential_dupes)


## None of the below turned out to be duplicates based on curation of metadata
#`Sample 1`   Revised_Source_Niche.x value Revised_Source_Niche.y `Sample 2`
#AusGEM_00746 Human                      0 Human                  ERR5355599
#AusGEM_00760 Human                      0 Human                  SRR7912422
#AusGEM_01980 Human                      0 Human                  SRR4067495
#TEA_GNB_2846 Human                      0 Human                  SRR7912397
#ERR5355599 Human                      0 Human                  AusGEM_00746
#SRR4067495 Human                      0 Human                  AusGEM_01980
#SRR7912397 Human                      0 Human                  TEA_GNB_2846
#SRR7912422 Human                      0 Human                  AusGEM_00760