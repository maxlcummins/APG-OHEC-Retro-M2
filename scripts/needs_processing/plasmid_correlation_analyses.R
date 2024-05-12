plasmid_correlation <- function(cluster_id){
        if(!exists("clusters")){
                if(file.exists("delims/plasmid_cluster_data.txt")){
                        clusters <- vroom("delims/plasmid_cluster_data.txt")     
                }else{source("scripts/plasmid_data_aggregate.R")}
        }        

        #Read in genometa
        genometa <- vroom("delims/genometa_n5631.txt", show_col_types = FALSE)
        
        #Read in our tree
        tree <- ggtree::read.tree("analysis/TEA/cgmlst_coreugate/pad.nwk")
        
        #Remove quotes from sample names
        tree$tip.label <- gsub("'", "", tree$tip.label)
        
        #Create a list of samples in our cohort
        sample_names <- tree$tip.label %>% as.data.frame() %>% rename('name' = '.')
        
        #Read in our mobtyper data
        mobtyper_df <- vroom("analysis/TEA/mobtyper_results_concatenated.txt", show_col_types = FALSE)
        
        #Filter in our mobtyper data and rename the sample ID column to name
        mobtyper_df <- mobtyper_df %>% rename('name' = sample_id) %>% mutate(name = str_replace(name, ":.*", "")) %>% filter(name %in% tree$tip.label)
        
        
        if(!exists("TEA_abricate_co_occurence_N90L90")){
                abricateR::abricateR(abricate_in = "/Users/131785/Dropbox/PostDoc/Austrakka/AusTrakka_R/analysis/TEA/genotype.txt", output = "TEA_abricate", writecsv = FALSE)}
        
        
        mobtyper_contig_wise <- read_delim("analysis/TEA/mobtyper_results.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
        
        
        AA176_APG_names <- mobtyper_df %>% filter(primary_cluster_id == cluster_id) %>% pull(name)
        
        #Read in our mobtyper data
        mobtyper_contig_wise <- mobtyper_contig_wise %>% rename('name' = sample_id) %>% filter(name %in% AA176_APG_names)
        
        small_mobtyper_hits <- mobtyper_contig_wise %>% select(name, contig_id, primary_cluster_id) %>% mutate(contig_id = str_replace(contig_id, " .*", ""))
        
        
        
        co_occurence_df <- TEA_abricate_co_occurence_N90L90 %>% rename('contig_id' = SEQUENCE) %>% inner_join(small_mobtyper_hits, by = c("name", "contig_id"))
        
        co_occurence_df_wide <- co_occurence_df %>%
                group_by(name, primary_cluster_id) %>%
                summarise(same_scaff = paste0(same_scaff, collapse = " ")) %>%
                #filter(primary_cluster_id %in% decat_known_PPCs_known_rep_wide$primary_cluster_id, primary_cluster_id != is.na(primary_cluster_id)) %>%
                mutate(same_scaff2 = strsplit(same_scaff, " ")) %>%
                unnest(same_scaff2) %>% 
                arrange(same_scaff2) %>%  
                pivot_wider(names_from = same_scaff2,
                            names_repair = "universal", values_from = same_scaff2, 
                            values_fill = 0, values_fn = length)
        
        AA176_scaffolds <- co_occurence_df_wide %>% filter(primary_cluster_id == cluster_id)
        
        AA176_plasdb <- clusters %>% filter(primary_cluster_id == cluster_id)
        
        AA176_plasdb <- AA176_plasdb %>% select(sample_id, primary_cluster_id, matches("vfdb"), matches("EC_custom"))
        
        AA176_plasdb <- AA176_plasdb %>% mutate(sample_id = paste(sample_id, primary_cluster_id)) %>% column_to_rownames('sample_id') %>% select(-primary_cluster_id)

        AA176_plasdb <- AA176_plasdb %>% mutate(across(everything(), ~ replace_na(., 0))) %>% select(where(~ is.numeric(.x) && sum(.x) !=0 )) 
        
        assign(x = "AA176_plasdb", value = AA176_plasdb, envir = globalenv())
                
        colnames(AA176_plasdb) <- gsub("ABRICATE..", "", colnames(AA176_plasdb))
        
        AA176_scaffolds <- AA176_scaffolds %>% column_to_rownames('name') 
        
        AA176_scaffolds <- AA176_scaffolds %>% select(any_of(colnames(AA176_plasdb)))
        
        AA176_plasdb <- AA176_plasdb %>% select(any_of(colnames(AA176_scaffolds)))
        
        plasmid_matrix <- AA176_plasdb
        genome_matrix <- AA176_scaffolds
        
        # Calculate Jaccard similarity between two binary vectors
        jaccard_similarity <- function(vector1, vector2) {
                intersect_count <- sum(vector1 == 1 & vector2 == 1)
                union_count <- sum(vector1 == 1 | vector2 == 1)
                return (intersect_count / union_count)
        }
        
        # Initialize an empty matrix to hold similarity scores
        num_plasmids <- nrow(plasmid_matrix)
        num_genomes <- nrow(genome_matrix)
        similarity_matrix <- matrix(0, nrow=num_plasmids, ncol=num_genomes)
        
        # Populate the similarity_matrix with Jaccard similarity scores
        for(i in 1:num_plasmids) {
                print(i)
                for(j in 1:num_genomes) {
                        similarity_matrix[i, j] <- jaccard_similarity(plasmid_matrix[i, ], genome_matrix[j, ])
                }
        }
        
        colnames(similarity_matrix) <- rownames(genome_matrix)
        rownames(similarity_matrix) <- rownames(plasmid_matrix)
        
        plas_heatmap<- pheatmap(similarity_matrix, fontsize_col = 2)
        
        plotme <- apply(X = similarity_matrix, FUN = max, MARGIN = 2) %>% as.data.frame()
        
        plotme <- plotme %>% rename('score' = '.') %>% mutate(test = 1)
        
        plotme %>% ggplot() + geom_boxplot(aes(x = score, y = test))
        
        plas_mean <- plotme$score[complete.cases(plotme$score)] %>% mean()
        plas_median <- plotme$score[complete.cases(plotme$score)] %>% median()
        
        assign(x = "cluster_id", value = cluster_id, envir = globalenv())
        assign(x = "plas_mean", value = plas_mean, envir = globalenv())
        assign(x = "plas_median", value = plas_median, envir = globalenv())
}

test <- data.frame("primary_cluster_id" = NA, 'mean_score' = NA, median_score = NA)

for(i in c("AA319", "AA179", "AA175", "AA176", "AA172", "AA315", "AA336", "AA337", "AA171", "AA297")){
#for(i in c("AA747", "AA176", "AA172", "AA315", "AA336", "AA337", "AA171", "AA297")){
        plasmid_correlation(i)
        df <- data.frame("primary_cluster_id" = cluster_id, 'mean_score' = plas_mean, median_score = plas_median)
        test <- bind_rows(test, df)
        
}
