library(pould)

perform_LD_analysis <- function(input_file, output_file) {
  #cat(input_file)
  df <- read.table(input_file, sep=",", header=TRUE)
  result <- cALD(df, inPhase=FALSE, verbose=FALSE)
  write.table(result, file=output_file, sep=",", row.names=FALSE, col.names=TRUE)
}

process_all_files <- function(input_folder, output_folder) {
  files <- list.files(input_folder, pattern="*.csv", full.names=TRUE)
  
  for (file in files) {
    output_file <- file.path(output_folder, basename(file))
    perform_LD_analysis(file, output_file)
  }
}

# Example usage:
#process_all_files("D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/kir_LD", 
 #                 "D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/kir_LD_result")

#process_all_files("D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_LD", 
  #                "D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_LD_result")

#process_all_files("D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD", 
 #                 "D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD_result")

#process_all_files("D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD", 
   #               "D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD_result")

process_all_files("D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD", 
                  "D:/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD_result")