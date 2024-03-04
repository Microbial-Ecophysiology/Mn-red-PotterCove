# function for modifying asv tables with only taxa above certain threshold

## correct script, previous script did not sum all ASVs per taxon, but checked if any ASV is above threshold

## info input
# table: input table, long format
# abundance threshold: taxa above which threshold in at least one sample should be kept separate, in %
# Abundance: name of column which has abundance numbers, not in %
# Sample: name of column which has unique name for every sequenced sample

# function to change taxon names if below set threshold
sort_abundant_taxa <- function(table, abundance_threshold, Abundance, Sample) {

  # check if rank names are correct
  if (identical(tail(colnames(phylo_melt), n=6), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))) {
    list_taxa_above_threshold <- function(input_table, rank, threshold, abundance_column, sample_column) {
      df1 <- input_table %>% 
        group_by({{rank}}, {{sample_column}}) %>%                                          ## group table by rank (Genus, Family, etc.) and by sample
        summarise(sum_tax_abundance = sum({{abundance_column}}), .groups = "keep") %>%     ## sum all ASVs of the same taxonomic assignment together, but separately per sample
        ungroup() %>%                                                                      ## ASVs now gone, next check across all samples if taxon is above threshold in at least one
        group_by({{rank}}) %>% 
        summarise(max_abundance = max(sum_tax_abundance), .groups = "keep") %>% 
        filter(max_abundance > threshold/100 & is.na({{rank}}) == FALSE) %>% 
        arrange(desc(max_abundance))
      out_vector <- pull(df1, {{rank}})
      return(out_vector)
    }
    
    ## create vectors for each rank which are above threshold
    genus_thr <- list_taxa_above_threshold(input_table = table,
                                           rank = Genus,
                                           threshold = abundance_threshold,
                                           abundance_column = Abundance,
                                           sample_column = Sample)
    
    family_thr <- list_taxa_above_threshold(input_table = table,
                                            rank = Family,
                                            threshold = abundance_threshold,
                                            abundance_column = Abundance,
                                            sample_column = Sample)
    
    order_thr <- list_taxa_above_threshold(input_table = table,
                                           rank = Order,
                                           threshold = abundance_threshold,
                                           abundance_column = Abundance,
                                           sample_column = Sample)
    
    class_thr <- list_taxa_above_threshold(input_table = table,
                                           rank = Class,
                                           threshold = abundance_threshold,
                                           abundance_column = Abundance,
                                           sample_column = Sample)
    
    phylum_thr <- list_taxa_above_threshold(input_table = table,
                                            rank = Phylum,
                                            threshold = abundance_threshold,
                                            abundance_column = Abundance,
                                            sample_column = Sample)
    
    ## rename taxa, if taxa is on that rank in all samples below threshold, it will be named after other of the next higher rank
    pf1 <- table %>% 
      mutate(Genus_mod = if_else(Genus %in% genus_thr, Genus,
                             if_else(Family %in% family_thr, paste0("other_f_", Family, "_<", abundance_threshold, "%"), 
                                     if_else(Order %in% order_thr, paste0("other_o_", Order, "_<", abundance_threshold, "%"), 
                                             if_else(Class %in% class_thr, paste0("other_c_", Class, "_<", abundance_threshold, "%"), 
                                                     if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                                             paste0("other_", Kingdom, "_<", abundance_threshold, "%")))))),
         Family_mod = if_else(Family %in% family_thr, Family, 
                              if_else(Order %in% order_thr, paste0("other_o_", Order, "_<", abundance_threshold, "%"), 
                                      if_else(Class %in% class_thr, paste0("other_c_", Class, "_<", abundance_threshold, "%"), 
                                              if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                                      paste0("other_", Kingdom, "_<", abundance_threshold, "%"))))),
         Order_mod = if_else(Order %in% order_thr, Order, 
                             if_else(Class %in% class_thr, paste0("other_c_", Class, "_<", abundance_threshold, "%"), 
                                     if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                             paste0("other_", Kingdom, "_<", abundance_threshold, "%")))),
         Class_mod = if_else(Class %in% class_thr, Class, 
                             if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                     paste0("other_", Kingdom, "_<", abundance_threshold, "%"))),
         Phylum_mod = if_else(Phylum %in% phylum_thr, Phylum, 
                              paste0("other_", Kingdom, "_<", abundance_threshold, "%"))
         ) %>% 
      arrange(Phylum_mod, Class_mod, Order_mod, Family_mod, Genus_mod)
    
    # create list which contains individual vectors for each taxon rank
    output <- vector("list", 6)
    output[[1]] <- genus_thr
    output[[2]] <- family_thr
    output[[3]] <- order_thr
    output[[4]] <- class_thr
    output[[5]] <- phylum_thr
    # actual long output table
    output[[6]] <- pf1
    
    # give the elements of the list names which contains information of the threshold used
    names(output) <- c(paste0("genera_abv_", abundance_threshold, "_perc"),
                   paste0("families_abv_", abundance_threshold, "_perc"),
                   paste0("orders_abv_", abundance_threshold, "_perc"),
                   paste0("classes_abv_", abundance_threshold, "_perc"),
                   paste0("phyla_abv_", abundance_threshold, "_perc"),
                   paste0("ASV_table_taxa_abv_", abundance_threshold, "_perc"))
    
    # output list
    return(output)
  }
  else {print("Error: rank names have to be changed to 'Kingdom', 'Phylum', 'Class', 'Order', 'Family',  'Genus'")}
  
}