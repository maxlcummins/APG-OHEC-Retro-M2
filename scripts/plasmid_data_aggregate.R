
#Set WD to project dir
setwd(here::here())

#Define paths to plasmid data
plasmid_abritamr_path <- "analysis/mobsuite_ref_plasmid_pipelord/mobsuite_refs_abritamr_resistance.txt"

plasmid_IncF_RST_path <- "analysis/mobsuite_ref_plasmid_pipelord/mobsuite_refs_IncF_RST.txt"

plasmid_genotype_path <- "analysis/mobsuite_ref_plasmid_pipelord/mobsuite_refs_genotype.txt"

plasmid_pMLST_path <- "analysis/mobsuite_ref_plasmid_pipelord/mobsuite_refs_pMLST.txt"

plasmid_IncF_RST_path <- "analysis/mobsuite_ref_plasmid_pipelord/mobsuite_refs_IncF_RST.txt"

#Read in plasmid clusters
clusters <- vroom("delims/mobsuite/clusters.txt", show_col_types = FALSE)

#Create clusters of plasmids
potential_plasmids <- clusters %>%
        #Group by primary cluster ID
        group_by(primary_cluster_id) %>%
        #Create a list of quartiles of plasmid sizes for each plasmid in a cluster
        summarise(q = list(quantile(size))) %>%
        unnest_wider(q, names_repair = ~paste0('size_Q_', sub('%', '', .))) %>%
        rename("primary_cluster_id" = size_Q_primary_cluster_id) %>%
        #Bind this back to our original table
        left_join(clusters, by = "primary_cluster_id") %>%
        #Remove plasmids which arent within Q2-Q3 (Central 50% of distribution)
        filter(size >= size_Q_25 & size <= size_Q_75) %>%
        #Add the count of rep types for a given primary cluster ID
        group_by(primary_cluster_id) %>%
        add_count(`rep_type(s)`) %>%
        #filter(primary_cluster_id %in% mobtyper_results$primary_cluster_id) %>%
        #Identify the most common replicon types within the remaning plasmids
        mutate(Majority_replicon_types = `rep_type(s)`[n == max(n)][1]) %>%
        dplyr::select(Majority_replicon_types, everything())


#Read in abricate data for tea dataset
abricateR::abricateR(abricate_in = plasmid_genotype_path, output = "plasmid_abricate", writecsv = FALSE)

#Rename our abricate data
plasmid_df <- plasmid_abricate_simple_summary_N90L90

#Add a prefix for the abricate columns
plasmid_df <- plasmid_df %>% rename_with(~gsub("^", "ABRICATE..", .x)) %>% rename("name" = "ABRICATE..name") %>% rename("ColV" = "ABRICATE..ColV")

#Read in pMLST data and IncF RST data
plas_pMLST <- vroom(plasmid_pMLST_path, show_col_types = FALSE)
plas_IncF_RST <- vroom(plasmid_IncF_RST_path, show_col_types = FALSE)

#Make pMLST table wide instead of long
plas_pMLST <- plas_pMLST %>% dplyr::select(name, pMLST_scheme, Sequence_type) %>% pivot_wider(names_from = "pMLST_scheme", values_from = Sequence_type, values_fill = "None")

#Combine both genotypic and abricate data
plasmid_df <- full_join(plasmid_df, plas_pMLST, by = "name") %>% full_join(plas_IncF_RST, by = "name") %>% rename('sample_id' = name)

#Join to original clusters
clusters <- clusters %>% full_join(plasmid_df, by = "sample_id")

#Read in plasmid clusters
plasmid_abritamr <- vroom("analysis/mobsuite_ref_plasmid_pipelord/mobsuite_refs_abritamr_resistance.txt", show_col_types = FALSE)

#Rename name column
plasmid_abritamr<- plasmid_abritamr %>% rename('sample_id' = name)

#Create a list of samples in our tree
sample_names <- clusters %>% select(sample_id)

#Filter out samples not in the tree, add in those missing from dataframe
plasmid_abritamr <- left_join(sample_names, plasmid_abritamr)

#Prepend plasmid_abritamr cols with ABRITAMR
colnames(plasmid_abritamr) <- gsub("^", "ABRITAMR.", colnames(plasmid_abritamr))

#Get rid of spaces from column names
colnames(plasmid_abritamr) <- gsub(" ", "_", colnames(plasmid_abritamr))

#Clean up 'sample_id' column sample_id for joining
plasmid_abritamr <- plasmid_abritamr %>% rename("sample_id" = "ABRITAMR.sample_id")

#Duplicate the table to determine CIA status
CIA_table <- plasmid_abritamr

#Select CIA gene columns
CIA_table <- CIA_table %>% select(sample_id, matches('Quinolone'), matches('[^not_]ESBL'), matches('[^or_]Carbapenemase'), matches('Colistin'), matches('Macrolide'), matches('Lincosamide'))



#Replace cell contents with 1 for cells which are not empty
CIA_table <- CIA_table %>%
        column_to_rownames('sample_id') %>%
        mutate(across(everything(),
                      ~str_replace(.x,
                                   ".*",
                                   "1"))) %>%
        mutate(across(everything(),
                      ~replace_na(.x,
                                  "0"))) %>%
        mutate(across(everything(),
                      ~as.numeric(.x)))

#Sum the number of CIA resistance genes and determine if greater than 1 (i.e. CIA resistant)
CIA_res <- CIA_table %>%
        mutate(CIA_res_status = rowSums(CIA_table)) %>%
        select(CIA_res_status)  %>%
        rownames_to_column('sample_id') %>%
        full_join(sample_names, by = "sample_id") %>%
        mutate(CIA_res_status = if_else(CIA_res_status >= 1, 1, 0)) %>%
        mutate(across(everything(), ~ replace_na(.x, 0)))

#Create our long format table, binding in CIA res and MDR/gene count data from above
long_format_plasmid_abritamr_processed <- plasmid_abritamr %>%
        pivot_longer(2:ncol(plasmid_abritamr),
                     names_to = "resistance",
                     values_drop_na = TRUE) %>%
        mutate(class_includes_CIA_res = str_detect(resistance,
                                                   'Quinolone|[^not_]ESBL|[^or_]Carbapenemase|Colistin|Macrolide|Lincosamide')) %>%
        rename('genes' = value) %>%
        mutate(resistance_genes_present = 1) %>%
        select(sample_id, resistance, class_includes_CIA_res, resistance_genes_present, genes) %>%
        mutate(resistance = str_replace(resistance, 'ABRITAMR\\.', '')) %>%
        #inner_join(MDR_status, by = 'sample_id') %>% 
        inner_join(CIA_res, by = 'sample_id') %>%
        mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
        group_by(sample_id) %>% 
        add_count() %>%
        rename("Plasmid_CIA_status" = `CIA_res_status`, "Plasmid_count_of_AMR_genes" = n)

clusters <- long_format_plasmid_abritamr_processed %>%
        select(sample_id, Plasmid_CIA_status, Plasmid_count_of_AMR_genes) %>%
        unique() %>%
        full_join(clusters) %>%
        mutate(across(.cols = c(Plasmid_CIA_status, Plasmid_count_of_AMR_genes), ~ replace_na(.x, 0))) %>%
        mutate(across(matches('ABRICATE'), ~ replace_na(.x, 0)))

cluster_VFs <- clusters %>% select(sample_id, matches("vfdb")) %>%
        column_to_rownames("sample_id") %>%
        mutate(Plasmid_VF_count = rowSums(.)) %>%
        rownames_to_column("sample_id") %>%
        full_join(sample_names, by = "sample_id") %>%
        select(sample_id, Plasmid_VF_count) %>%
        mutate(across(everything(), ~ replace_na(.x, 0)))

clusters <- left_join(clusters, cluster_VFs)

clusters <- clusters %>% select(sample_id, primary_cluster_id, secondary_cluster_id, predicted_mobility, Plasmid_CIA_status, Plasmid_count_of_AMR_genes, Plasmid_VF_count, `IncF RST`, any_of(colnames(plas_pMLST)), everything())

#long_format_plasmid_abritamr_processed <- full_join(long_format_plasmid_abritamr_processed, sample_names) %>% mutate(Strain_MDR_status = replace_na(Strain_MDR_status,0)) %>% mutate(Strain_count_of_AMR_gene_classes = replace_na(Strain_count_of_AMR_gene_classes,0)) %>% mutate(Strain_CIA_status = replace_na(Strain_CIA_status,0)) %>% mutate(genes = replace_na(genes, ''))

#Write to file
#long_format_plasmid_abritamr_processed%>% write_delim("delims/mobsuite_long_format_plasmid_cluster_data.txt", delim = "\t")
clusters %>% write_delim("delims/plasmid_cluster_data.txt", delim = "\t")
