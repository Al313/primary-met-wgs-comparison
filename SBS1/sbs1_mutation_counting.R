






# for (i in args[1]:args[1]){
# 
# 
#   sample_id <- metadata_included$sample_id[i]
# 
# 
#   cohort <- tolower(as.vector(metadata_included$cohort[i]))
#   cancer_type <- as.character(metadata_included$cancer_type[i])
# 
#   print(i)
#   print(sample_id)
# 
#   if (dir.exists("/hpc/cuppen/")){
#     if (metadata_included$cohort[i] == "HMF") {
#       path_to_vcf <- paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/", sample_id, "/purple/",sample_id,".purple.somatic.vcf.gz")
#     } else if (metadata_included$cohort[i] == "PCAWG") {
#       path_to_vcf <- paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/", sample_id, "-from-jar/purplesoft3.3/", sample_id, "T.purple.somatic.postprocessed.vcf.gz")
#     }
#   } else {
#     if (metadata_included$cohort[i] == "HMF") {
#       path_to_vcf <- paste0(local, "/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/", sample_id, "/purple/",sample_id,".purple.somatic.vcf.gz")
#     } else if (metadata_included$cohort[i] == "PCAWG") {
#       path_to_vcf <- paste0(local, "/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/", sample_id, "-from-jar/purplesoft3.3/", sample_id, "T.purple.somatic.postprocessed.vcf.gz")
#     }
#   }
# 
# 
# 
# 
# 
#   if (cohort == "pcawg") {
#     sample1 <- sample_id
#     sample2 <- paste0(sample_id, "T")
#   } else if (cohort == "hmf"){
#     if (substr(sample_id, nchar(sample_id),nchar(sample_id)) == "I" | substr(sample_id, nchar(sample_id),nchar(sample_id)) == "V"){
#       sample_id2 <- str_split(sample_id, pattern = "TI")[[1]][1]
#     } else {
#       sample_id2 <- substr(sample_id, 1,(nchar(sample_id)-1))
#     }
#     sample1 <- sample_id2
#     sample2 <- sample_id
#   }
# 
#   tryCatch(vcf <- variantsFromVcf(vcf.file = path_to_vcf,
#                                   merge.consecutive = T,
#                                   vcf.fields = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER","INFO", "FORMAT", paste0(sample1, "R"), sample2)), error=function(e) NULL)
# 
#   if(!is.null(vcf)){
#     if (nrow(vcf) > 0){
# 
#       vcf <- vcf[vcf$filter == "PASS",]
#       vaf <- extractSigsSnv(df=vcf[,c(1:2,4:11,3)], output='df')
#       vaf$chrom <- sapply(str_split(vaf$chrom, pattern = "r"), "[[", 2)
# 
#       vaf <- vaf[,c(1:2, 11, 3:10, 12:15)]
# 
#       if (cancer_type == "Breast cancer"){
#         vaf <- vaf[grepl('\\w\\[C>T\\]G', vaf$context),]
#         vaf <- vaf[substr(vaf$context, 1,1) !="T",]
#       } else if (cancer_type == "Skin melanoma"){
#         vaf <- vaf[grepl('\\w\\[C>T\\]G', vaf$context),]
#         vaf <- vaf[substr(vaf$context, 1,1) %notin% c("C", "T"),]
#       } else {
#         vaf <- vaf[grepl('\\w\\[C>T\\]G', vaf$context),]
#         vaf <- vaf[substr(vaf$context, 1,1) !="T",]
#       }
# 
#       vaf <- vaf[,1:11]
# 
#       colnames(vaf) <- toupper(colnames(vaf))
#       rownames(vaf) <- 1:nrow(vaf)
# 
#       selelcted_info_fields <- getInfoValues(vaf$INFO, keys = c("TNC", "SUBCL"))
#       selelcted_info_fields$SUBCL[is.na(selelcted_info_fields$SUBCL)] <- 0
#       selelcted_info_fields$SUBCL <- as.numeric(selelcted_info_fields$SUBCL)
#       selelcted_info_fields$clonality <- NA
#       selelcted_info_fields$clonality[which(selelcted_info_fields$SUBCL >= 0.8)] <- "subclonal"
#       selelcted_info_fields$clonality[which(selelcted_info_fields$SUBCL < 0.8)] <- "clonal"
# 
#       vaf <- cbind(vaf, selelcted_info_fields)
# 
#     }
#   }
# 
#   if (metadata_included$cohort[i] == "HMF") {
#     if (dir.exists("/hpc/cuppen/")){
#       file_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
#     } else {
#       file_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
#     }
#   } else if (metadata_included$cohort[i] == "PCAWG") {
#     if (dir.exists("/hpc/cuppen/")){
#       file_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
#     } else {
#       file_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
#     }
#   }
# 
#   combined_df <- tryCatch(read.csv(file = file_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
# 
# 
#   if (!is.null(combined_df)){
#     vaf$coorPos <- paste0(vaf$CHROM, "_", vaf$REF, "_", vaf$POS, "_", vaf$ALT)
#     combined_df$coorPos <- paste0(combined_df$CHROM, "_", combined_df$REF, "_", combined_df$POS, "_", combined_df$ALT)
#     vaf <- merge(vaf, combined_df[,c("coorPos", "timing_class")], by = "coorPos")
#     vaf <- vaf[,-1]
#   } else {
#     vaf$timing_class <- NA
#   }
#   
#   write.table(vaf, file = gzfile(paste0(wd, "r-objects/clock-like-vafs-age-investigation/", sample_id, ".purple.somatic.clocklike-subset.vcf.gz")), sep = "\t", row.names = F, quote = F)
# 
# }



###############################################################
# sbs_count <- data.frame()
# 
# for (i in 1:nrow(metadata_included)){
# 
#   sample_id <- metadata_included$sample_id[i]
#   print(i)
# 
#   clock_like_muts <- tryCatch(read.csv(file = paste0(wd, "r-objects/clock-like-vafs-age-investigation/", sample_id, ".purple.somatic.clocklike-subset.vcf.gz"), sep = "\t", header = T, stringsAsFactors = F), error=function(e) NULL)
#   if(is.null(clock_like_muts)){
#     clock_like_number <- NA
#     subclonal_count <- NA
#     timing_late_clonal_count <- NA
#     timing_early_clonal_count <- NA
#     timing_na_clonal_count <- NA
#     timing_subclonal_count <- NA
#   } else {
#     clock_like_number <- nrow(clock_like_muts)
#     subclonal_count <- nrow(clock_like_muts[clock_like_muts$clonality == "subclonal",])
#     if (all(is.na(clock_like_muts$timing_class))) {
#       timing_late_clonal_count <- NA
#       timing_early_clonal_count <- NA
#       timing_na_clonal_count <- NA
#       timing_subclonal_count <- NA
#     } else {
#       timing_late_clonal_count <- nrow(clock_like_muts[clock_like_muts$timing_class == "clonal [late]",])
#       timing_early_clonal_count <- nrow(clock_like_muts[clock_like_muts$timing_class == "clonal [early]",])
#       timing_na_clonal_count <- nrow(clock_like_muts[clock_like_muts$timing_class == "clonal [NA]",])
#       timing_subclonal_count <- nrow(clock_like_muts[clock_like_muts$timing_class == "subclonal",])
#     }
#   }
# 
#   sbs_count <- rbind(sbs_count, c(sample_id, clock_like_number, subclonal_count, timing_late_clonal_count, timing_early_clonal_count,
#                                   timing_na_clonal_count, timing_subclonal_count))
# 
# }
# 
# 
# colnames(sbs_count) <- c("sample_id", "clock_like_number", "subclonal_count", "timing_late_clonal_count", "timing_early_clonal_count",
#                          "timing_na_clonal_count", "timing_subclonal_count")
# 
# 
# write.table(sbs_count, file = paste0(wd, "r-objects/clock-like-count-age-investigation2.txt"), quote = F, row.names = F, sep = "\t")


