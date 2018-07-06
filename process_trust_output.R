# hardcoding this so this routine can run on the cluster. it takes some time
setwd("/rstudio-common/dbortone/projects/bgv/Optimize_Diversity_Metrics_CRSV1371")
source("setup/init.R")
my_core_number = 32
library(Biostrings)
library(vegan)

#dataset_name = "trust_on_rnaseq_data"
#dataset_name = "trust_extended_on_rnaseq_data"
#dataset_name = "trust_cannonical_on_rnaseq_data"
dataset_name = "trust_cannonical_extended_on_rnaseq_data"


if(grepl("extended", dataset_name)){
  my_paths =list.files("/home/dbortone/scratch/Optimize_Diversity_Metrics_CRSV1371/trust_rnaseq/output", 
                       full.names = T, recursive = T, pattern = "-ALL-extended.txt$")
} else {
  my_paths =list.files("/home/dbortone/scratch/Optimize_Diversity_Metrics_CRSV1371/trust_rnaseq/output", 
                       full.names = T, recursive = T, pattern = "-ALL.txt$")
}

if(grepl("cannonical", dataset_name)){
  require_cannonical = TRUE
} else {
  require_cannonical = FALSE
}

output_path = get_formatted_matrix_path(dataset_name)

source("diversity_metric_functions.R")

my_chains = c("TRA", "TRB")#, "TRD", "TRG")
my_range = 1:length(my_paths)

sample_output = mclapply(my_paths[my_range], function(my_path){
  output_list = list()
  path_components = strsplit(my_path, "/")[[1]]
  my_sample = path_components[length(path_components)-1]
  output_list["Sample_ID"] = my_sample
  
  my_df = data.table::fread(my_path, header = T)
  
  # names(my_df) = c("Cell_Count", "V_Allele_1", "J_Allele_1", "ntCDR3_1", "RNA_1", "DNA_1","V_Allele_2", "J_Allele_2", "ntCDR3_2", "RNA_2", "DNA_2")
  my_df = data.frame( 
    Cell_Count = my_df$contig_reads_count,
    V_Allele = gsub(">", "", my_df$Vgene, fixed = T),
    J_Allele = gsub(">", "", my_df$Jgene, fixed = T),
    ntCDR3 = my_df$cdr3dna,
    aaCDR3 = my_df$cdr3aa
  )
  my_df$V_Allele = my_df$V_Allele %>% as.character()
  my_df$J_Allele = my_df$J_Allele %>% as.character()
  my_df$ntCDR3 = my_df$ntCDR3 %>% as.character()
  my_df$aaCDR3 = my_df$aaCDR3 %>% as.character()
  
  
  my_df$Chain = substr(paste0(my_df$V_Allele,my_df$J_Allele),1,3)
  
  if(require_cannonical){
    my_df$aaCDR3_length = nchar(my_df$aaCDR3)
    my_df$is_cannonical = TRUE
    my_df$is_cannonical[my_df$aaCDR3_length < 5] = FALSE
    my_df$is_cannonical[!grepl("^C", my_df$aaCDR3)] = FALSE
    my_df$is_cannonical[!grepl("F$", my_df$aaCDR3)] = FALSE
    my_df = my_df[my_df$is_cannonical, ]
  }
  
  # need to drop the last 9 nt since this is th fGXG motif
  # my_df$ntCDR3 = sapply(my_df$ntCDR3, function(x){ substring(x, 1, (nchar(x)-9))})
  
  # nt sequences can be duplicated by tcrer if you seperate the heavy from the light chain and if you drop the last GXG
  # lets add these together
  summed_counts = tapply(my_df$Cell_Count, my_df$ntCDR3, sum)
  
  my_df = my_df[!duplicated(my_df$ntCDR3), ]
  my_df$Cell_Count = summed_counts[my_df$ntCDR3]
  # mixcer data diversities are done by ntCDR3 so we'll do the same here
  # still lets translate to aaCDR3 to make sure there aren't any with stop codons
  # my_df$aaCDR3 = mclapply(my_df$ntCDR3, function(x){translate(DNAString(x)) %>% as.character()}, mc.cores = my_core_number) %>% unlist()
  # 
  # # drop sequences with stop codons. i filter out of frames and stop codons with my formatted mixcr exportClones --filter-out-of-frames \ --filter-stops call
  my_df =  my_df[!grepl("*", my_df$aaCDR3, fixed = T), ]
  
  my_df = my_df[order(as.numeric(as.character(my_df$Cell_Count)), decreasing = T), ]
  total_abundance = sum(my_df$Cell_Count, na.rm = T)
  # for each chain do our own diversity metrics and write out the vdjtool matrices
  for(my_chain in my_chains){
    
    # make individual chain subdat
    chain_file_name = paste0(my_sample, "_", my_chain, ".tsv")
    chain_df = my_df[grepl(paste0("^", my_chain), my_df$Chain), ]
    chain_counts = chain_df$Cell_Count %>% as.character() %>% as.numeric()
    chain_abundance = sum(chain_counts, na.rm = T)
    chain_richness = sum(chain_counts > 0)
    chain_chao1 = estimateR(chain_counts)["S.chao1"] %>% as.numeric()
    
    output_list[paste0(my_chain,"_Abundance")] = chain_abundance
    output_list[paste0(my_chain,"_Richness")] = chain_richness
    output_list[paste0(my_chain,"_Fraction")] = chain_abundance/total_abundance
    output_list[paste0(my_chain,"_d05")] = dXX_index(chain_counts, 0.05)
    output_list[paste0(my_chain,"_d10")] = dXX_index(chain_counts, 0.10)
    output_list[paste0(my_chain,"_d15")] = dXX_index(chain_counts, 0.15)
    output_list[paste0(my_chain,"_d20")] = dXX_index(chain_counts, 0.20)
    output_list[paste0(my_chain,"_d25")] = dXX_index(chain_counts, 0.25)
    output_list[paste0(my_chain,"_d30")] = dXX_index(chain_counts, 0.30)
    output_list[paste0(my_chain,"_d35")] = dXX_index(chain_counts, 0.35)
    output_list[paste0(my_chain,"_d40")] = dXX_index(chain_counts, 0.40)
    output_list[paste0(my_chain,"_d45")] = dXX_index(chain_counts, 0.45)
    output_list[paste0(my_chain,"_d50")] = dXX_index(chain_counts, 0.50)
    output_list[paste0(my_chain,"_d75")] = dXX_index(chain_counts, 0.75)
    output_list[paste0(my_chain,"_n10n20")] = n1_v_n2_index(chain_counts, 0.1, 0.2)
    output_list[paste0(my_chain,"_n25n50")] = n1_v_n2_index(chain_counts, 0.25, 0.5)
    output_list[paste0(my_chain,"_Shannon_Entropy")] = shannon_entropy(chain_counts)
    output_list[paste0(my_chain,"_Evenness")] = evenness(chain_counts)
    output_list[paste0(my_chain,"_Inv_Simpson")] = inv_simpson(chain_counts)
    output_list[paste0(my_chain,"_Chao1")] = chain_chao1
    output_list[paste0(my_chain,"_Richness_v_Chao1")] =  chain_richness / chain_chao1
    output_list[paste0(my_chain,"_NonSingletons")] =  sum(chain_counts > 1) / chain_richness 
  }
  return(output_list)
}, mc.cores = my_core_number)

sample_ouput_df = rbindlist(sample_output, use.names=TRUE, fill=TRUE) %>% as.data.frame()

log2_transform_cols = names(sample_ouput_df)[grepl("_Abundance$|_Richness$|_Inv_Simpson$|TR._Chao1$", names(sample_ouput_df))]

for (log2_transform_col in log2_transform_cols){
  sample_ouput_df[, paste0("Log2_", log2_transform_col)] = log2(sample_ouput_df[ , log2_transform_col] + 1)
}

fwrite(sample_ouput_df, output_path, sep = "\t")

