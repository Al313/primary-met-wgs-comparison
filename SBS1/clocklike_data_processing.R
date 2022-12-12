


## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

## Dependencies + defined values --------------------------------

source(paste0(base_dir, 'drivers/analysis/dna-rep-ann/final-update/code_for_github/source_script.R'))


## Sample metadata ================================
sample_metadata <- read.delim(
  paste0(base_dir,'/passengers/processed/metadata/main/metadata_20221111_1251.txt.gz'),
  na.strings=c('NA','#N/A')
)


# Prepare the metadata

sample_metadata <- sample_metadata[!sample_metadata$is_blacklisted,]

row.names(sample_metadata) <- 1:nrow(sample_metadata)


sample_metadata[,"cohort"] <- factor(sample_metadata[,"cohort"], levels = cohort_order)
sample_metadata[,"tissue_group"] <- factor(sample_metadata[,"tissue_group"], levels = tissue_group_order)
sample_metadata[,"cancer_type"] <- factor(sample_metadata[,"cancer_type"], levels = cancer_type_order)
sample_metadata[,"cancer_type_code"] <- factor(sample_metadata[,"cancer_type_code"],
                                               levels = cancer_type_code_order)

sample_metadata[,"gender"] <- factor(sample_metadata[,"gender"], levels = c("MALE", "FEMALE"))

sample_metadata <- sample_metadata

cont_mat <- read.csv(file = paste0(base_dir, "passengers/processed/sigs_denovo/extractions/12_fixed_seeds/sig_contrib/fit_lsq.post_processed/denovo_contribs.lsq.post_processed.txt.gz"), header = T, stringsAsFactors = F, sep = "\t")
cont_mat$sample_id <- row.names(cont_mat)


cont_mat$`SBS5/40` <- cont_mat$SBS5 + cont_mat$SBS40


## final remake of the updated_sig5_&_40_contribs_ploidy.tsv and complete1-5-40.tsv


sample_metadata <- merge(sample_metadata, cont_mat[,c("sample_id", "SBS5", "SBS40", "SBS5/40")], by = "sample_id")


ploidy_prim <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/SuppTable2_karyotype _karyotype_271022.xlsx"), sheet = 1))
ploidy_met <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/SuppTable2_karyotype _karyotype_271022.xlsx"), sheet = 2))

ploidy <- rbind(ploidy_prim, ploidy_met)
colnames(ploidy)[1] <- "sample_id_2"


sample_metadata <- merge(sample_metadata, ploidy[,c("sample_id_2", "genome")])

colnames(sample_metadata)
sample_metadata <- sample_metadata[,c(1:50, 54, 51:53)]

sample_metadata$SBS5_ploidy_corrected <- sample_metadata$SBS5/sample_metadata$genome
sample_metadata$SBS40_ploidy_corrected <- sample_metadata$SBS40/sample_metadata$genome
sample_metadata$`SBS5/40_ploidy_corrected` <- sample_metadata$`SBS5/40`/sample_metadata$genome


# write.table(metadata_included, file = paste0(wd, "r-objects/rebuttal/updated_sig5_&_40_contribs_ploidy.tsv"), sep = "\t", row.names = F, quote = F)


#########################################################################################

sample_metadata$subclonal_tmb <- rowSums(sample_metadata[,c("sbs_load.subclonal", "dbs_load.subclonal", "indel_load.subclonal")])

head(sample_metadata)

# metadata_included <- readRDS(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/processed_metadata_24082022.rds"))
# colnames(metadata_included)
# nrow(metadata_included)
# adding age info
# Read in and process raw Hartwig age information


hartwig_meta <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/hmf-clinical-metadata-20211112.tsv"), stringsAsFactors = F, header = T, sep = "\t")

hartwig_meta$biopsyDateProcessed <- NA
hartwig_meta$biopsyDateProcessed <- unlist(lapply(str_split(hartwig_meta$biopsyDate, pattern = "-"), "[[", 1))

hartwig_meta$age_at_biopsy <- NA

hartwig_meta[,c("birthYear","biopsyDateProcessed")] <- apply(hartwig_meta[,c("birthYear","biopsyDateProcessed")], as.numeric, MARGIN = 2)

hartwig_meta$age <- hartwig_meta$biopsyDateProcessed - hartwig_meta$birthYear


colnames(hartwig_meta)[2] <- "sample_id"

hartwig_meta <- hartwig_meta[,c(2,ncol(hartwig_meta))]




# Read in and process raw PCWG age information


pcawg_meta <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/Metadata/donor.info.followup.all_projects.csv"), header = T, stringsAsFactors = F, sep = "\t")


colnames(pcawg_meta)[c(1,9)] <- c("sample_id", "age")
pcawg_meta <- pcawg_meta[,c(1,9)]




# Merging the age info of both cohorts to the sample metadata object
age_meta_all <- rbind(hartwig_meta, pcawg_meta)
sample_metadata <- merge(sample_metadata, age_meta_all, by = "sample_id")

colnames(sample_metadata)
nrow(sbs_count)
#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@








sbs_count <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/clock-like-count-age-investigation2.txt"), header = T, stringsAsFactors = F, sep = "\t")
sbs_count <- sbs_count[sbs_count$sample_id %in% sample_metadata$sample_id,]

sbs_count <- merge(sbs_count, sample_metadata[,c("sample_id", "age", "subclonal_tmb")], all.x = T)
colnames(sbs_count) <- c("sample_id", "sbs1_count", "sbs1_purple_subclonal", "sbs1_clonal_late", "sbs1_clonal_early", "sbs1_clonal_na", "sbs1_mutationTimeR_subclonal", "age", "tmb_purple_subclonal")

#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@

cont_mat <- read.csv(file = paste0(base_dir, "passengers/processed/sigs_denovo/extractions/12_fixed_seeds/sig_contrib/fit_lsq.post_processed/denovo_contribs.lsq.post_processed.txt.gz"), header = T, stringsAsFactors = F, sep = "\t")
cont_mat$sample_id <- row.names(cont_mat)


cont_mat$`SBS5/40` <- cont_mat$SBS5 + cont_mat$SBS40

cont_mat <- cont_mat[,c("sample_id", "SBS5", "SBS40", "SBS5/40")]

colnames(cont_mat)[2:4] <- c("sbs5_exposure", "sbs40_exposure", "sbs5/40_exposure") 

sbs_count <- merge(sbs_count, cont_mat, by = "sample_id")

#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@

ploidy_prim <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/SuppTable2_karyotype _karyotype_271022.xlsx"), sheet = 1))
ploidy_met <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/SuppTable2_karyotype _karyotype_271022.xlsx"), sheet = 2))

ploidy <- rbind(ploidy_prim, ploidy_met)

colnames(ploidy)[1] <- "sample_id_2"


ploidy <- merge(ploidy, metadata_included[,c("sample_id_2", "sample_id")])

colnames(ploidy)[2] <- "ploidy"

sbs_count <- merge(sbs_count, ploidy[,c("sample_id", "ploidy")], by = "sample_id")

#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@


all_timing_sbs5 <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/all_timing_processed_sbs5.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)
all_timing_sbs40 <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/all_timing_processed_sbs40.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)

all_timing_sbs5[all_timing_sbs5$sbs5_count == 0,3:7] <- 0
all_timing_sbs40[all_timing_sbs40$sbs40_count == 0,3:7] <- 0


colnames(all_timing_sbs5) <- c("sample_id", "sbs5_count", "sbs5_purple_subclonal", "sbs5_clonal_late", "sbs5_clonal_early", "sbs5_clonal_na", "sbs5_mutationTimeR_subclonal")
colnames(all_timing_sbs40) <- c("sample_id", "sbs40_count", "sbs40_purple_subclonal", "sbs40_clonal_late", "sbs40_clonal_early", "sbs40_clonal_na", "sbs40_mutationTimeR_subclonal")
all_timing_sbs5_40 <-all_timing_sbs5[,2:7] + all_timing_sbs40[,2:7]
all_timing_sbs5_40 <- cbind(all_timing_sbs5[,1], all_timing_sbs5_40)
colnames(all_timing_sbs5_40) <- c("sample_id", "sbs5/40_count", "sbs5/40_purple_subclonal", "sbs5/40_clonal_late", "sbs5/40_clonal_early", "sbs5/40_clonal_na", "sbs5/40_mutationTimeR_subclonal")



sbs_count <- merge(merge(merge(sbs_count, all_timing_sbs5, by = "sample_id"), all_timing_sbs40, by = "sample_id"), all_timing_sbs5_40, by = "sample_id")



#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@

# metadata_included <- metadata_included[!is.na(metadata_included$subclonal_tmb),]
# metadata_included <- metadata_included[!is.na(metadata_included$clonal_tmb),]
# 
# 
# metadata_included <- metadata_included[metadata_included$subclonal_tmb != 0,]
# metadata_included <- metadata_included[metadata_included$clonal_tmb != 0,]



metadata_included$total_clonality_ratio <- metadata_included$clonal_tmb/metadata_included$total_tmb

metadata_included$total_subclonality_ratio <- metadata_included$subclonal_tmb/metadata_included$total_tmb

colnames(metadata_included)[49] <- "tmb_total"

sbs_count <- merge(sbs_count, metadata_included[,c("sample_id", "total_clonality_ratio", "total_subclonality_ratio", "tmb_total")], by = "sample_id", all.x = T)





sbs_count$sbs1_clonality_ratio <- (sbs_count$sbs1_count - sbs_count$sbs1_purple_subclonal)/sbs_count$sbs1_count
sbs_count$sbs1_normalized_clonality <- sbs_count$sbs1_clonality_ratio/sbs_count$total_clonality_ratio


sbs_count$sbs1_subclonality_ratio <- sbs_count$sbs1_purple_subclonal/sbs_count$sbs1_count
sbs_count$sbs1_normalized_subclonality <- sbs_count$sbs1_subclonality_ratio/sbs_count$total_subclonality_ratio




#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@


all_timing <- read.csv(file = paste0(wd, "r-objects/all_timing_processed_sbs1.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)


all_timing <- all_timing[!is.na(all_timing$all_late_clonal),]


all_timing$sbs1_timing_normalized1 <- (all_timing$cln_late_clonal/(all_timing$cln_late_clonal + all_timing$cln_early_clonal))/(all_timing$all_late_clonal/(all_timing$all_late_clonal + all_timing$all_early_clonal))



all_timing$sbs1_timing_normalized2 <- (all_timing$cln_late_clonal/(all_timing$cln_late_clonal + all_timing$cln_early_clonal + all_timing$cln_na_clonal))/(all_timing$all_late_clonal/(all_timing$all_late_clonal + all_timing$all_early_clonal + all_timing$all_na_clonal))



colnames(all_timing)[c(1:5, 12:13)] <- c("sample_id", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na", "tmb_mutationTimeR_subclonal", "sbs1_timing_normalized1", "sbs1_timing_normalized2")

sbs_count <- merge(sbs_count, all_timing[,c("sample_id", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na", "tmb_mutationTimeR_subclonal", "sbs1_timing_normalized1", "sbs1_timing_normalized2")], by = "sample_id", all.x = T)

#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@


sbs_count$sbs5_clonality_ratio <- (sbs_count$sbs5_count - sbs_count$sbs5_purple_subclonal)/sbs_count$sbs5_count
sbs_count$sbs5_normalized_clonality <- sbs_count$sbs5_clonality_ratio/sbs_count$total_clonality_ratio


sbs_count$sbs40_clonality_ratio <- (sbs_count$sbs40_count - sbs_count$sbs40_purple_subclonal)/sbs_count$sbs40_count
sbs_count$sbs40_normalized_clonality <- sbs_count$sbs40_clonality_ratio/sbs_count$total_clonality_ratio



sbs_count$`sbs5/40_clonality_ratio` <- (sbs_count$`sbs5/40_count` - sbs_count$`sbs5/40_purple_subclonal`)/sbs_count$`sbs5/40_count`
sbs_count$`sbs5/40_normalized_clonality` <- sbs_count$`sbs5/40_clonality_ratio`/sbs_count$total_clonality_ratio


#@#@#@#@#@#@#@@#@#@#@#@#@#@#@#@#@#@

all_timing_sbs5 <- read.csv(file = paste0(wd, "r-objects/rebuttal/all_timing_processed_sbs5.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)


all_timing <- all_timing[,1:5]



all_timing_sbs5 <- merge(all_timing_sbs5, all_timing, by = "sample_id")



## Normalization method 1


all_timing_sbs5$sbs5_timing_normalized1 <- (all_timing_sbs5$sbs5_late_clonal/(all_timing_sbs5$sbs5_late_clonal + all_timing_sbs5$sbs5_early_clonal))/(all_timing_sbs5$tmb_clonal_late/(all_timing_sbs5$tmb_clonal_late + all_timing_sbs5$tmb_clonal_early))

## Normalization method 2

all_timing_sbs5$sbs5_timing_normalized2 <- (all_timing_sbs5$sbs5_late_clonal/(all_timing_sbs5$sbs5_late_clonal + all_timing_sbs5$sbs5_early_clonal + all_timing_sbs5$sbs5_na_clonal))/(all_timing_sbs5$tmb_clonal_late/(all_timing_sbs5$tmb_clonal_late + all_timing_sbs5$tmb_clonal_early + all_timing_sbs5$tmb_clonal_na))

sbs_count <- merge(sbs_count, all_timing_sbs5[,c("sample_id", "sbs5_timing_normalized1", "sbs5_timing_normalized2")], all.x = T)



all_timing_sbs40 <- read.csv(file = paste0(wd, "r-objects/rebuttal/all_timing_processed_sbs40.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)



all_timing_sbs40 <- merge(all_timing_sbs40, all_timing, by = "sample_id")



## Normalization method 1


all_timing_sbs40$sbs40_timing_normalized1 <- (all_timing_sbs40$sbs40_late_clonal/(all_timing_sbs40$sbs40_late_clonal + all_timing_sbs40$sbs40_early_clonal))/(all_timing_sbs40$tmb_clonal_late/(all_timing_sbs40$tmb_clonal_late + all_timing_sbs40$tmb_clonal_early))

## Normalization method 2


all_timing_sbs40$sbs40_timing_normalized2 <- (all_timing_sbs40$sbs40_late_clonal/(all_timing_sbs40$sbs40_late_clonal + all_timing_sbs40$sbs40_early_clonal + all_timing_sbs40$sbs40_na_clonal))/(all_timing_sbs40$tmb_clonal_late/(all_timing_sbs40$tmb_clonal_late + all_timing_sbs40$tmb_clonal_early + all_timing_sbs40$tmb_clonal_na))

sbs_count <- merge(sbs_count, all_timing_sbs40[,c("sample_id", "sbs40_timing_normalized1", "sbs40_timing_normalized2")], all.x = T)





all_timing_sbs5 <- read.csv(file = gzfile(paste0(wd, "r-objects/rebuttal/all_timing_processed_sbs5.tsv.gz")), sep = "\t", stringsAsFactors = F, header = T)

all_timing_sbs40 <- read.csv(file = gzfile(paste0(wd, "r-objects/rebuttal/all_timing_processed_sbs40.tsv.gz")), sep = "\t", stringsAsFactors = F, header = T)

all_timing_sbs5_dummy <- all_timing_sbs5
colnames(all_timing_sbs5_dummy)[2:7] <- c("V1", "V2", "V3", "V4", "V5", "v6")


all_timing_sbs40_dummy <- all_timing_sbs40
colnames(all_timing_sbs40_dummy)[2:7] <- c("V1", "V2", "V3", "V4", "V5", "v6")


all_timing_sbs5_40 <- all_timing_sbs5_dummy[,2:7] + all_timing_sbs40_dummy[,2:7]
all_timing_sbs5_40$sample_id <- all_timing_sbs5$sample_id

all_timing_sbs5_40 <- all_timing_sbs5_40[,c(7,1:6)]

colnames(all_timing_sbs5_40) <- c("sample_id", "sbs540_count", "sbs540_subcl_subclonal", "sbs540_late_clonal",
                                  "sbs540_early_clonal", "sbs540_na_clonal", "sbs540_timing_subclonal")




all_timing_sbs5_40 <- merge(all_timing_sbs5_40, all_timing, by = "sample_id")



## Normalization method 1

all_timing_sbs5_40$sbs540_timing_normalized1 <- (all_timing_sbs5_40$sbs540_late_clonal/(all_timing_sbs5_40$sbs540_late_clonal + all_timing_sbs5_40$sbs540_early_clonal))/(all_timing_sbs5_40$tmb_clonal_late/(all_timing_sbs5_40$tmb_clonal_late + all_timing_sbs5_40$tmb_clonal_early))



## Normalization method 2



all_timing_sbs5_40$sbs540_timing_normalized2 <- (all_timing_sbs5_40$sbs540_late_clonal/(all_timing_sbs5_40$sbs540_late_clonal + all_timing_sbs5_40$sbs540_early_clonal + all_timing_sbs5_40$sbs540_na_clonal))/(all_timing_sbs5_40$tmb_clonal_late/(all_timing_sbs5_40$tmb_clonal_late + all_timing_sbs5_40$tmb_clonal_early + all_timing_sbs5_40$tmb_clonal_na))


colnames(all_timing_sbs5_40)[12:13] <- c("sbs5/40_timing_normalized1", "sbs5/40_timing_normalized2") 

sbs_count <- merge(sbs_count, all_timing_sbs5_40[,c("sample_id", "sbs5/40_timing_normalized1", "sbs5/40_timing_normalized2")], all.x = T)


#%%%%%%%%%%%
sbs_count <- sbs_count %>% mutate_all(~ifelse(is.nan(.), NA, .))


# sbs_count1 <- read.csv(file = paste0(wd, "r-objects/complete1-5-40.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# 
# 
# sbs_count1 <- sbs_count1[sbs_count1$sample_id %in% sbs_count$sample_id,]
# 
# 
# sbs_count$tmb_mutationTimeR_subclonal <- sbs_count1$tmb_mutationTimeR_subclonal
#%%%%%%%%%%%


sbs_count <- sbs_count[,c("sample_id", "sbs1_count", "sbs1_purple_subclonal", "sbs1_clonal_late", "sbs1_clonal_early", "sbs1_clonal_na", "sbs1_mutationTimeR_subclonal",
                          "sbs5_exposure", "sbs5_count", "sbs5_purple_subclonal", "sbs5_clonal_late", "sbs5_clonal_early", "sbs5_clonal_na", "sbs5_mutationTimeR_subclonal",
                          "sbs40_exposure", "sbs40_count", "sbs40_purple_subclonal", "sbs40_clonal_late", "sbs40_clonal_early", "sbs40_clonal_na", "sbs40_mutationTimeR_subclonal",
                          "sbs5/40_exposure", "sbs5/40_count", "sbs5/40_purple_subclonal", "sbs5/40_clonal_late", "sbs5/40_clonal_early", "sbs5/40_clonal_na",
                          "sbs5/40_mutationTimeR_subclonal", "tmb_total", "tmb_purple_subclonal", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na", "tmb_mutationTimeR_subclonal",
                          "age", "ploidy", "total_clonality_ratio", "total_subclonality_ratio", "sbs1_clonality_ratio", 
                          "sbs1_normalized_clonality", "sbs1_subclonality_ratio", "sbs1_normalized_subclonality", "sbs1_timing_normalized1",
                          "sbs1_timing_normalized2", "sbs5_clonality_ratio", "sbs5_normalized_clonality", "sbs40_clonality_ratio", "sbs40_normalized_clonality",
                          "sbs5/40_clonality_ratio", "sbs5/40_normalized_clonality", "sbs5_timing_normalized1", "sbs5_timing_normalized2", "sbs40_timing_normalized1",
                          "sbs40_timing_normalized2", "sbs5/40_timing_normalized1", "sbs5/40_timing_normalized2")]



c("sample_id", "sbs1_count", "sbs1_purple_subclonal", "sbs1_clonal_late", "sbs1_clonal_early", "sbs1_clonal_na", "sbs1_mutationTimeR_subclonal",
  "sbs5_exposure", "sbs5_count", "sbs5_purple_subclonal", "sbs5_clonal_late", "sbs5_clonal_early", "sbs5_clonal_na", "sbs5_mutationTimeR_subclonal",
  "sbs40_exposure", "sbs40_count", "sbs40_purple_subclonal", "sbs40_clonal_late", "sbs40_clonal_early", "sbs40_clonal_na", "sbs40_mutationTimeR_subclonal",
  "sbs5/40_exposure", "sbs5/40_count", "sbs5/40_purple_subclonal", "sbs5/40_clonal_late", "sbs5/40_clonal_early", "sbs5/40_clonal_na",
  "sbs5/40_mutationTimeR_subclonal", "tmb_total", "tmb_purple_subclonal", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na",
  "tmb_mutationTimeR_subclonal", "age", "ploidy", "total_clonality_ratio", "total_subclonality_ratio", "sbs1_clonality_ratio",
  "sbs1_normalized_clonality", "sbs1_subclonality_ratio", "sbs1_normalized_subclonality", "sbs1_timing_normalized1",
  "sbs1_timing_normalized2", "sbs5_clonality_ratio", "sbs5_normalized_clonality", "sbs40_clonality_ratio", "sbs40_normalized_clonality",
  "sbs5/40_clonality_ratio", "sbs5/40_normalized_clonality", "sbs5_timing_normalized1", "sbs5_timing_normalized2", "sbs40_timing_normalized1",
  "sbs40_timing_normalized2", "sbs5/40_timing_normalized1", "sbs5/40_timing_normalized2")[c("sample_id", "sbs1_count", "sbs1_purple_subclonal", "sbs1_clonal_late", "sbs1_clonal_early", "sbs1_clonal_na", "sbs1_mutationTimeR_subclonal",
                                                                                            "sbs5_exposure", "sbs5_count", "sbs5_purple_subclonal", "sbs5_clonal_late", "sbs5_clonal_early", "sbs5_clonal_na", "sbs5_mutationTimeR_subclonal",
                                                                                            "sbs40_exposure", "sbs40_count", "sbs40_purple_subclonal", "sbs40_clonal_late", "sbs40_clonal_early","sbs40_clonal_na", "sbs40_mutationTimeR_subclonal",
                                                                                            "sbs5/40_exposure", "sbs5/40_count", "sbs5/40_purple_subclonal", "sbs5/40_clonal_late", "sbs5/40_clonal_early", "sbs5/40_clonal_na",
                                                                                            "sbs5/40_mutationTimeR_subclonal", "tmb_total", "tmb_purple_subclonal", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na",
                                                                                            "tmb_mutationTimeR_subclonal", "age", "ploidy", "total_clonality_ratio", "total_subclonality_ratio", "sbs1_clonality_ratio",
                                                                                            "sbs1_normalized_clonality", "sbs1_subclonality_ratio", "sbs1_normalized_subclonality", "sbs1_timing_normalized1",
                                                                                            "sbs1_timing_normalized2", "sbs5_clonality_ratio", "sbs5_normalized_clonality", "sbs40_clonality_ratio", "sbs40_normalized_clonality",
                                                                                            "sbs5/40_clonality_ratio", "sbs5/40_normalized_clonality", "sbs5_timing_normalized1", "sbs5_timing_normalized2", "sbs40_timing_normalized1",
                                                                                            "sbs40_timing_normalized2", "sbs5/40_timing_normalized1", "sbs5/40_timing_normalized2") %notin% colnames(sbs_count)]





meta <- read.csv(file = paste0(wd, "external-files/SuppTable1_sample_metadata - metadata-v2-13112022.tsv"), sep = "\t", stringsAsFactors = F, header = T)


metadata_included <- meta[!meta$is_blacklisted,]

sbs_count <- sbs_count[sbs_count$sample_id %in% metadata_included$sample_id,]

write.table(sbs_count, file = paste0(wd, "r-objects/rebuttal/complete1-5-40.tsv"), sep = "\t", row.names = F, quote = F)





