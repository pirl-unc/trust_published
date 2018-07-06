#USE_EXISTING_VDJCA = "true" # "true" or "false"
BATCH_NAME = "a"
MY_THREADS = 6
SBATCH_MEM_PER_CPU = "2g" # careful this might crash a step that can only use one thread, but needs all the memory
# SBATCH_MEM = "24g"
SH_HEADER = "#!/bin/bash"
MY_SPECIES = "hsa"
MY_TRUST_IMAGE = "dockerreg.bioinf.unc.edu:5000/trust_3.0.1:1"
MY_BOWTIE2_IMAGE = "dockerreg.bioinf.unc.edu:5000/bowtie2_2.3.3.1:1"
TIME_LIMIT = "UNLIMITED"
library(pryr)


scratch_dir = file.path("/home/dbortone/scratch", "Optimize_Diversity_Metrics_CRSV1371") # this is where the commands and fastqs will go

fastq_dir = file.path(scratch_dir, "stig_rnaseq", "output")
fastq_paths = list.dirs(fastq_dir, full.names = T) # pattern = paste0(BATCH_NAME,"_....$")
fastq_paths = fastq_paths[grepl(paste0(BATCH_NAME,"_....$"), fastq_paths)]


# check for repertoires
completed_reps = list.files(file.path(scratch_dir, "stig_rnaseq", "output"), pattern = "a_.....statistics.csv", recursive = T)
missing_reps = c()
for (this_num in 1:1000){
  my_search = str_pad(this_num, width = 4, side = c("left"), pad = "0")
  found_files = grep(paste0("^a_", my_search), completed_reps)
  if(length(found_files) == 0){missing_reps = c(missing_reps, paste0("a_", my_search))}
}

print(missing_reps)
# for(missing_rep in missing_reps){
#   remove_dir = file.path(scratch_dir, "stig_rnaseq", "output", missing_rep)
#   unlink(remove_dir, recursive = T) 
# }


base_names = basename(fastq_paths)
sample_df = data.frame(
  Sample_ID = base_names, 
  Sample_Folder = base_names, 
  Fastq_Directory = fastq_paths,
  Fastq_R1_File = paste0(base_names, "_R1.degraded.fastq.gz"),
  Fastq_R2_File = paste0(base_names, "_R2.degraded.fastq.gz")#,
  #Chain = "TCR", # from stig: The default value is a 27-mer that anchors on the reverse strand in EX1 of the beta chain C-region
  #RNA_Seq = "false"
  )

# sample_df = sample_df[grepl("^a_", sample_df$Sample_ID), ]
# sample_df = sample_df[sample_df$Sample_ID %ni% missing_reps, ]


batch_run_dir = file.path(scratch_dir, "trust_rnaseq")
batch_command_dir = file.path(batch_run_dir, "commands")
dir.create(batch_run_dir, showWarnings = F)
dir.create(batch_command_dir, showWarnings = F)


output_dir = file.path(batch_run_dir, "output")
dir.create(output_dir, showWarnings = F)


# writeLines(sample_df$Sample_ID %>% as.character(), file.path(batch_run_dir, paste0(BATCH_NAME, "_samples.txt")))


add_c = function(...){
  new_text = paste0(...)
  note_envir = pryr::where('my_commands')
  assign('my_commands', c(get('my_commands', envir = note_envir), new_text), envir = note_envir)
}


my_qcommands = c(SH_HEADER)
add_q = function(...){
  new_text = paste0(...)
  note_envir = pryr::where('my_qcommands')
  assign('my_qcommands', c(get('my_qcommands', envir = note_envir), new_text), envir = note_envir)
}


# rerun_samples = c("p106_Mid_PBMC_CD19_BCR", "p106_Pre_PBMC_CD19_BCR", "p107_Mid_PBMC_CD19_BCR", "p107_Pre_PBMC_CD19_BCR", "p102_Pre_PBMC_CD3_TCR")
# sample_df = sample_df[sample_df$Sample_ID %in% rerun_samples, ]
for(row_index in 1:(nrow(sample_df))){
  my_commands = c()

  my_sample = sample_df$Sample_Folder[row_index]
  
  input_path_1 = file.path(sample_df$Fastq_Directory[row_index], sample_df$Fastq_R1_File[row_index])
  input_path_2 = file.path(sample_df$Fastq_Directory[row_index], sample_df$Fastq_R2_File[row_index])
  my_output_dir = file.path(output_dir,  my_sample)
  my_sam_file_path = file.path(my_output_dir,  paste0(my_sample, ".sam"))
  my_unsorted_bam_file_path = file.path(my_output_dir,  paste0("unsorted_", my_sample, ".bam"))
  my_bam_file_path = file.path(my_output_dir,  paste0(my_sample, ".bam"))
  
  add_c(SH_HEADER)
  add_c("#SBATCH --job-name ", sample_df$Sample_ID[row_index] )
  add_c("#SBATCH --partition docker")
  add_c("#SBATCH --time ", TIME_LIMIT)
  add_c("#SBATCH --mincpus ", MY_THREADS)
  add_c("#SBATCH --mem-per-cpu ", SBATCH_MEM_PER_CPU)
  #add_c("#SBATCH --exclude r820-docker-2-1.local")
  add_c("#SBATCH --output ", file.path(my_output_dir, "_output.txt"))
  add_c("#SBATCH --error ", file.path(my_output_dir, "_error.txt"))
  add_c("")
  add_c("")
  add_c("# Run Bowtie2")
  add_c("docker run --rm=true \\")
  add_c("  -v /datastore:/datastore:shared \\")
  add_c("  -v /home/dbortone:/home/dbortone \\")
  add_c("  ", MY_BOWTIE2_IMAGE, " bash -c \\")
  add_c("  'bowtie2 \\")
  add_c("    --threads ", MY_THREADS, " \\")
  add_c("    -x /datastore/nextgenout2/share/labs/imgf/ref/bowtie2-build/bowtie2-hg38 \\")
  add_c("    -1 ", input_path_1, " \\")
  add_c("    -2 ", input_path_2, " \\")
  add_c("    -S ", my_sam_file_path, "'")
  add_c("")
  add_c("")
  add_c("# Convert samfile to sorted bam and index it")
  add_c("samtools view -@ ", MY_THREADS, " -S -b ", my_sam_file_path, " > ", my_unsorted_bam_file_path)
  add_c("samtools sort -@ ", MY_THREADS, " -o ", my_bam_file_path, " ", my_unsorted_bam_file_path)
  add_c("samtools index -@ ", MY_THREADS, " ", my_bam_file_path)
  add_c("")
  add_c("")
  add_c("# Run TRUST")
  add_c("docker run --rm=true \\")
  add_c("  -v /datastore:/datastore:shared \\")
  add_c("  -v /home/dbortone:/home/dbortone \\")
  add_c("  ", MY_TRUST_IMAGE, " bash -c \\")
  add_c("  'trust \\")
  add_c("    -H \\") # do heavy and light chains 
  add_c("    -E \\") # do extension step, creates a second file 
  add_c("    --CoreN=", MY_THREADS, " \\")
  add_c("    --genome=hg38 \\")
  add_c("    -f ", my_bam_file_path, " \\")
  add_c("    -o ", my_output_dir, "/'")
  add_c("")
  add_c("")
  add_c("# Clean up files")
  add_c("find ", my_output_dir, " -type f \\")
  add_c("  -not -name '*-TCR-ALL-extended.txt' \\")
  add_c("  -not -name '*-TCR-ALL.txt' \\")
  add_c("  -not -name '*_error.txt' \\")
  add_c("  -not -name '*_output.txt' \\")
  add_c("  -delete")
  add_c("")
  
  this_file = paste0(my_sample, ".job")
  
  # make directory and modify permissions
  add_q("mkdir -p ", my_output_dir)
  add_q("sbatch ", file.path(batch_command_dir, this_file))
  
  writeLines(text = my_commands, con = file.path(batch_command_dir, this_file))
}
writeLines(text = my_qcommands, con = file.path(batch_command_dir, paste0("_", BATCH_NAME, "_master.sh")))
