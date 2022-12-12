

if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/final-update/"
  local <- ""
} else {
  local <- "/home/ali313/Documents/studies/master/umc-project"
  wd <- paste0(local,"/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/final-update/")
}




included_cancer_type <- c("Breast carcinoma", "Glioblastoma multiforme", "Colorectal carcinoma", "Esophageal carcinoma",
                          "Stomach carcinoma", "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
                          "Pancreas neuroendocrine", "Cervical carcinoma", "Ovarian serous adenocarcinoma", "Uterus carcinoma",
                          "Upper respiratory tract carcinoma", "Kidney renal clear cell carcinoma", "Lung adenocarcinoma", 
                          "Lung squamous cell carcinoma", "Diffuse large B-cell lymphoma", "Prostate carcinoma", 
                          "Skin melanoma", "Leiomyosarcoma", "Liposarcoma", "Thyroid carcinoma", "Bladder urothelial carcinoma")




included_cancer_type_code <- c("BRCA", "GBM", "COREAD", "ESCA", "STAD", "CHOL", "LIHC", "PAAD", "PANET", "CESC", "OV",
                               "UCEC", "HNSC", "KIRC", "LUAD", "LUSC", "DLBCL", "PRAD", "SKCM", "LMS", "LPS", "THCA", "BLCA")



cancer_type_palette <- c("#702963", "#808080", "#B4C424", "#98FB98", "#228B22", "#008080", "#000080", "#40E0D0", 
                         "#87CEEB", "#FFB6C1", "#FF00FF", "#F33A6A", "#36454F", "#6F8FAF", "#A95C68", "#C38E96", 
                         "#8A9A5B", "#8B0000", "#FFC000", "#80461B", "#E49B0F", "#C2B280", "#EE4B2B")



included_cancer_type_fig1 <- c("Glioblastoma multiforme", "Upper respiratory tract carcinoma", "Thyroid carcinoma",
                               "Lung adenocarcinoma", "Lung squamous cell carcinoma","Diffuse large B-cell lymphoma", "Breast carcinoma",
                               "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
                               "Pancreas neuroendocrine", "Colorectal carcinoma", "Esophageal carcinoma", "Stomach carcinoma", 
                               "Kidney renal clear cell carcinoma", "Cervical carcinoma", "Ovarian serous adenocarcinoma", "Uterus carcinoma", 
                               "Bladder urothelial carcinoma", "Prostate carcinoma", "Leiomyosarcoma", "Liposarcoma", "Skin melanoma")



new_order <- match(included_cancer_type_fig1, included_cancer_type)

included_cancer_type_fig1 <- append("Pan-cancer", included_cancer_type_fig1)

included_cancer_type_code_fig1 <- append("PAN", included_cancer_type_code[new_order])



included_tissue_group_fig1 <- c("All tissues", "CNS", "Head_and_neck", "Thyroid", "Lung", "Breast", "Lymphoid", "GI_dev", "GI_core", "Kidney", "Gyn",
                                "Urothelial", "Prostate", "Soft_tissue", "Skin")





color_palette_cohort <- c("#f58134", "#9966CC")
# for fig 1 numbers #E66711 and #7217C5
cancer_type_palette_fig1 <- cancer_type_palette[match(included_cancer_type_fig1, included_cancer_type)]




# metadata_included <- readRDS(file = paste0(wd, "r-objects/processed_metadata_24082022.rds"))
# row.names(metadata_included) <- 1:nrow(metadata_included)
# 
# 
# metadata_included[,"cohort"] <- factor(metadata_included[,"cohort"], levels = c("PCAWG", "Hartwig"))
# metadata_included[,"tissue_group"] <- factor(metadata_included[,"tissue_group"], levels = included_tissue_group_fig1)
# metadata_included[,"cancer_type"] <- factor(metadata_included[,"cancer_type"], levels = included_cancer_type_fig1)
# metadata_included[,"cancer_type_code"] <- factor(metadata_included[,"cancer_type_code"],
#                                                  levels = included_cancer_type_code)
# 
# metadata_included[,"gender"] <- factor(metadata_included[,"gender"], levels = c("MALE", "FEMALE"))
# 
# metadata_included <- metadata_included[metadata_included$sample_id %in% meta$sample_id,]


meta <- read.csv(file = paste0(wd, "external-files/SuppTable1_sample_metadata - metadata-v2-13112022.tsv"), sep = "\t", stringsAsFactors = F, header = T)


metadata_included <- meta[!meta$is_blacklisted,]

row.names(metadata_included) <- 1:nrow(metadata_included)


metadata_included[,"cohort"] <- factor(metadata_included[,"cohort"], levels = c("PCAWG", "Hartwig"))
metadata_included[,"tissue_group"] <- factor(metadata_included[,"tissue_group"], levels = included_tissue_group_fig1)
metadata_included[,"cancer_type"] <- factor(metadata_included[,"cancer_type"], levels = included_cancer_type_fig1)
metadata_included[,"cancer_type_code"] <- factor(metadata_included[,"cancer_type_code"],
                                                 levels = included_cancer_type_code_fig1)

metadata_included[,"gender"] <- factor(metadata_included[,"gender"], levels = c("MALE", "FEMALE"))

metadata_hmf <- metadata_included[metadata_included$cohort == "Hartwig",]
# Loading libraries
library(readxl)
library(tidyr)
library(ggplot2)
library(stringr)


sex_df <- data.frame(cancer_type = character(48), cohort = character(48), male = numeric(48), female = numeric(48))








for (i in 1:24){
  for (j in 1:2){
    
    cancer_type <- included_cancer_type_fig1[i]
    cohort <- c("PCAWG", "Hartwig")[j]
    
    
    sex_df[i+(i-1) + (j-1),1:2] <- c(cancer_type, cohort)
    if (i != 1){
      sex_df[i+(i-1) + (j-1),3:4] <- 100*c(nrow(metadata_included[metadata_included$cancer_type == cancer_type & metadata_included$cohort == cohort &  metadata_included$gender == "MALE",]),
                                           nrow(metadata_included[metadata_included$cancer_type == cancer_type & metadata_included$cohort == cohort &  metadata_included$gender == "FEMALE",]))/nrow(metadata_included[metadata_included$cancer_type == cancer_type & metadata_included$cohort == cohort,])
    } else {
      sex_df[i+(i-1) + (j-1),3:4] <- 100*c(nrow(metadata_included[metadata_included$cohort == cohort &  metadata_included$gender == "MALE",]),
                                           nrow(metadata_included[metadata_included$cohort == cohort &  metadata_included$gender == "FEMALE",]))/nrow(metadata_included[metadata_included$cohort == cohort,])
    }
  }
}


sex_df <- tidyr::gather(sex_df, key = "gender", value = "proportion", 3:4)
sex_df$cancer_type <- factor(sex_df$cancer_type, levels = included_cancer_type_fig1)
sex_df$gender <- factor(sex_df$gender, levels = c("male", "female"))
sex_df$cohort <- factor(sex_df$cohort, levels = c("Hartwig", "PCAWG"))

kk <- ggplot(sex_df, aes(x = cohort,y = proportion, fill = gender)) + facet_wrap(~cancer_type, nrow = 44) +
  geom_bar(stat="identity") +
  scale_color_manual(values = rev(color_palette_cohort), labels = c("Hartwig", "PCAWG")) +
  scale_fill_manual(values = c("#8498f8", "#fc7571")) +
  guides(fill=FALSE) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, by = 50), labels = c("0", "50", "100%")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25,face="bold")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))






pdf(file = "/home/ali313/Desktop/v2/sex.pdf", height = 16.75, width = 14)
print(kk)
dev.off()







################# Getting the age info

# metadata_included <- readRDS(file = paste0(wd, "r-objects/processed_metadata_13122021.rds"))
# row.names(metadata_included) <- 1:nrow(metadata_included)
# 
# metadata_included[,"cohort"] <- factor(metadata_included[,"cohort"], levels = c("PCAWG", "Hartwig"))
# metadata_included[,"tissue_group"] <- factor(metadata_included[,"tissue_group"], levels = included_tissue_group_fig1)
# metadata_included[,"cancer_type"] <- factor(metadata_included[,"cancer_type"], levels = included_cancer_type_fig1)
# metadata_included[,"cancer_type_code"] <- factor(metadata_included[,"cancer_type_code"],
#                                                  levels = included_cancer_type_code)
# 
# metadata_included[,"gender"] <- factor(metadata_included[,"gender"], levels = c("MALE", "FEMALE"))


pcawg_meta <- readxl::read_xlsx(path = paste0(wd, "external-files/pcawg-clinical-metadata-20211112.xlsx"))
pcawg_meta <- as.data.frame(pcawg_meta)



hmf_meta <- read.csv(file = paste0(wd, "external-files/hmf-clinical-metadata-20211112.tsv"), stringsAsFactors = F, header = T, sep = "\t")

hmf_meta$biopsyDateProcessed <- NA
hmf_meta$biopsyDateProcessed <- unlist(lapply(str_split(hmf_meta$biopsyDate, pattern = "-"), "[[", 1))

hmf_meta$age_of_donor <- NA

hmf_meta[,c("birthYear","biopsyDateProcessed")] <- apply(hmf_meta[,c("birthYear","biopsyDateProcessed")], as.numeric, MARGIN = 2)

hmf_meta$age_of_donor <- hmf_meta$biopsyDateProcessed - hmf_meta$birthYear

# bnn <- merge(metadata_included[metadata_included$cohort == "PCAWG",], pcawg_meta[,c("icgc_donor_id","donor_age_at_diagnosis")], by.x = "sample_id", by.y = "icgc_donor_id")
wgd_original <- read.csv(file = paste0(wd, "external-files/pcawg-original-wgd-timing.txt"), header = T, stringsAsFactors = F, sep = "\t")
# wgd_original <- wgd_original[!duplicated(wgd_original$icgc_donor_id),]

bnn <- merge(metadata_included[metadata_included$cohort == "PCAWG",], wgd_original[,c("icgc_donor_id","age")], by.x = "sample_id", by.y = "icgc_donor_id")
colnames(bnn)[51] <- "age"

ann <- merge(metadata_included[metadata_included$cohort == "Hartwig",], hmf_meta[,c("sampleId","age_of_donor")], by.x = "sample_id", by.y = "sampleId")
colnames(ann)[51] <- "age"


metadata_included <- rbind(bnn, ann)
table(metadata_included$age[metadata_included$cohort == "Hartwig"], useNA = "always")
table(metadata_included$age[metadata_included$cohort == "PCAWG"], useNA = "always")

metadata_included[,"cohort"] <- factor(metadata_included[,"cohort"], levels = c("PCAWG", "Hartwig"))
metadata_included[,"tissue_group"] <- factor(metadata_included[,"tissue_group"], levels = included_tissue_group_fig1)
metadata_included[,"cancer_type"] <- factor(metadata_included[,"cancer_type"], levels = included_cancer_type_fig1)
metadata_included[,"cancer_type_code"] <- factor(metadata_included[,"cancer_type_code"],
                                                 levels = included_cancer_type_code)

metadata_included <- metadata_included[!duplicated(metadata_included),]

metadata_included_pan <- metadata_included
metadata_included_pan$cancer_type <- "Pan-cancer"

metadata_included <- rbind(metadata_included, metadata_included_pan)

library(plyr)
cohort <- c("PCAWG", "Hartwig")
p_meds <- ddply(metadata_included, .(cancer_type, cohort), summarise, med = median(age, na.rm = T))

# mean(p_meds$med[p_meds$cohort == "Hartwig"], na.rm = T) - mean(p_meds$med[p_meds$cohort == "PCAWG"], na.rm = T)

jj <- ggplot(metadata_included, aes(x = age, color = cohort))+ facet_wrap(~ cancer_type, ncol = 1)+
  geom_density(aes(color = cohort), size = 3) +
  scale_color_manual(values = c("#E66711","#7217C5")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  geom_vline(data=p_meds, aes(xintercept=med, colour=cohort),
             linetype=1, size=5) +
  guides(color = F)


median(metadata_included$age[metadata_included$cancer_type == "Skin melanoma" & metadata_included$cohort == "PCAWG"], na.rm = T)
median(metadata_included$age[metadata_included$cancer_type == "Skin melanoma" & metadata_included$cohort == "HMF"], na.rm = T)

# max(metadata_included$age, na.rm = T)
# min(metadata_included$age, na.rm = T)

pdf(file = "/home/ali313/Desktop/v2/age.pdf", height = 16.75, width = 14)
print(jj)
dev.off()


# write.table(metadata_included, file = paste0(wd, "r-objects/meta_data_with_age.tsv"), sep = "\t", quote = F, row.names = F)
# 
# metadata_included <- read.csv(file = paste0(wd, "r-objects/meta_data_with_age.tsv"), stringsAsFactors = F, sep = "\t", header = T)



################# Getting the treatment info








library(dplyr)
library(magrittr)
library(tidyr)


tt <- read.csv(file = "/home/ali313/Desktop/info_treatment.tsv",
               sep = "\t", stringsAsFactors = F, header = T)

metadata_included$treatment_info <- F
metadata_included$treatment_info[metadata_included$sample_id %in% unique(tt$sample_id)] <- T

table(metadata_included$treatment_info)


treatment_info1 <- metadata_included %>%
  filter(cohort == "HMF") %>%
  group_by(cancer_type) %>%
  summarize(treated = sum(treatment_info), untreated = sum(!treatment_info)) %>%
  mutate(treatment_prop = treated / (treated + untreated))

treatment_info$cancer_type <- factor(treatment_info$cancer_type, levels = included_cancer_type_fig1)

treatment_info[,2:3] <- treatment_info[,2:3]/rowSums(treatment_info[,2:3])
treatment_info_tibb <- gather(treatment_info, key = "treatment_status", value = "Prop", 2:3)
treatment_info_tibb$dummy <- 1

treatment_info_tibb <- treatment_info_tibb[treatment_info_tibb$treatment_status != "untreated",]

str(treatment_info_tibb)
table(treatment_info_tibb$treatment_status)

levels(treatment_info_tibb$cancer_type)
treatment_info_tibb$treatment_status <- factor(treatment_info_tibb$treatment_status, levels = c("treated", "untreated"))
treatment_info_tibb <- rbind(c("Pan-cancer", mean(treatment_info_tibb$treatment_prop), "treated", mean(treatment_info_tibb$treatment_prop), 1), treatment_info_tibb)

treatment_info_tibb[,c(2,4)] <- apply(treatment_info_tibb[,c(2,4)], 2, as.numeric)

nn <- ggplot(treatment_info_tibb, aes(x = dummy,y = treatment_prop, fill = treatment_status)) + facet_wrap(~cancer_type, nrow = 23) +
  geom_bar(stat="identity", size = 0.1) +
  # scale_color_manual(values = rev(color_palette_cohort), labels = c("HMF", "PCAWG")) +
  scale_fill_manual(values = c("#212f3d")) +
  guides(fill=FALSE) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, by = 50), labels = c("0", "50", "100%")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25,face="bold")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.spacing.y = unit(2, "lines"))




pdf(file = "/home/ali313/Desktop/treatment3.pdf", height = 16.75, width = 14)
print(nn)
dev.off()




################# Getting the treatment info





# treatment_summary <- data.frame(treatment_mech = c("Aromatase inhibitor", "Selective ER modulator", "Anti-ER", "Pyrimidine (ant)agonist",
#                                                    "Platinum", "Folate antagonist", "Experimental", "Alkaloid",
#                                                    "Taxane", "Multikinase inhibitor", "GnRH (ant)agonist", "Anti-AR",
#                                                    "Glucocorticoid", "Alkylating", "Antifolate", "Anthracycline",
#                                                    "Anthracycline, alkylating", "Anti-HER2", "mTOR inhibitor", "Somatostatin analogue",
#                                                    "Unknown", "Anti-EGFR", "Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor", "Microtubule inhibitor",
#                                                    "Anti-CTLA-4", "Anti-VEGF", "BRAF inhibitor", "Anti-PD-1",
#                                                    "MEK inhibitor", "Vinca Alkaloid", "Folinic acid", "Radionuclide",
#                                                    "HIPEC", "Topoisomerase inhibitor", "CDK4/6 inhibitor or placebo", "Antitumor antibiotic",
#                                                    "Anti-PDGFR-?", "Progestogen", "Immune therapy or placebo", "PI3K inhibitor",
#                                                    "CDK4/6 inhibitor", "Experimental or placebo", "Immunomodulator", "PARP inhibitor",
#                                                    "FGFR inhibitor", "Monoclonal antibody", "Oncolytic virus", "ALK/ROS1 inhibitor",
#                                                    "ALK inhibitor", "Anaplastic lymphoma kinase inhibitor, tyrosine kinase inhibitor", "Pyrimidine (ant)agonist, anthracycline, alkylating", "Anti-PD-L1",
#                                                    "Pyrimidine (ant)agonist, other", "Bisphosphonate", "FAK inhibitor", "Anti-CD20,MABthera-therapy",
#                                                    "Androgen inhibitor", "Other", "Radiotherapy", "Immunetherapy",
#                                                    "Pyrimidine (ant)agonist, platinum", "PARP inhibitor or placebo", "", "Anti-RANK-L",
#                                                    "Allosteric inhibitor", "Folinic acid, pyrimidine (ant)agonist, platinum", "Pyrimidine (ant)agonist, folinic acid", "Alkylating, vinca alkaloid, alkylating, glucocorticoid",
#                                                    "Anthracycline, antitumor antibiotic, vinca alkaloid, alkylating", "Vinca Alkaloid, alkaloid, glucocorticoid, anthracycline", "Alkylating, alkaloid, pyrimidine (ant)agonist, alkylating", "Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid",
#                                                    "PIPAC", "Nucleoside-analogon", "Alkylating, antifolate, pyrimidine (ant)agonist", "Chemoradiotherapy",
#                                                    "Taxane, anthracycline, alkylating", "Anti-PD-L1 or placebo", "Pyrimidine (ant)agonist, anthracycline, alkylating, taxane", "Anti-SMO",
#                                                    "Platinum, taxane", "Anti-GITR", "Anti-mitotic vinca alkaloid", "Hyperthermia",
#                                                    "Anti-PDGFR-? or placebo", "Nuclear therapy or placebo", "Anti-PSMA", "Immunotherapy",
#                                                    "Folate antagonist, Platinum", "Vinca Alkaloid, antitumor antibiotic, alkylating", "Anti-PD-1 or placebo", "Alkylating, alkylating, vinca alkaloid",
#                                                    "Platinum, Pyrimidine (ant)agonist", "Radiofrequency ablation", "Orteronel or placebo", "Taxane, Bisphosphonate",
#                                                    "Anti-CD20, alkylating, anthracycline, vinca alkaloid, glucocorticoid", "Proteasome inhibitor", "Immune modulator", "olinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor, platinum",
#                                                    "Aromatase inhibitor or placebo", "PI3K inhibitor or placebo"),
#                                 treatment_group = c("Hormone_therapy", "Hormone_therapy", "Hormone_therapy", "Chemotherapy",
#                                                     "Chemotherapy", "Chemotherapy", "Other", "Chemotherapy",
#                                                     "Chemotherapy", "Targeted_therapy", "Hormone_therapy", "Hormone_therapy",
#                                                     "Other", "Chemotherapy", "Chemotherapy", "Chemotherapy",
#                                                     "Chemotherapy", "Targeted_therapy", "Targeted_therapy", "Hormone_therapy",
#                                                     "Unknown", "Targeted_therapy", "Chemotherapy", "Chemotherapy",
#                                                     "Immunotherapy", "Targeted_therapy", "Targeted_therapy", "Immunotherapy",
#                                                     "Targeted_therapy", "Chemotherapy", "Chemotherapy", "Radiotherapy",
#                                                     "Chemotherapy", "Chemotherapy", "Targeted_therapy", "Chemotherapy",
#                                                     "Targeted_therapy", "Hormone_therapy", "Immunotherapy", "Targeted_therapy",
#                                                     "Targeted_therapy", "Other", "Immunotherapy", "Targeted_therapy",
#                                                     "Targeted_therapy", "Targeted_therapy", "Immunotherapy", "Targeted_therapy",
#                                                     "Targeted_therapy", "Targeted_therapy", "Chemotherapy", "Immunotherapy",
#                                                     "Chemotherapy", "Other", "Targeted_therapy", "Targeted_therapy",
#                                                     "Hormone_therapy", "Other", "Radiotherapy", "Immunotherapy",
#                                                     "Chemotherapy", "Targeted_therapy", "", "Targeted_therapy",
#                                                     "Targeted_therapy", "Chemotherapy", "Chemotherapy", "Chemotherapy",
#                                                     "Chemotherapy", "Chemotherapy", "Chemotherapy", "Chemotherapy",
#                                                     "Chemotherapy", "Chemotherapy", "Chemotherapy", "Chemotherapy|Radiotherapy",
#                                                     "Chemotherapy", "Immunotherapy", "Chemotherapy", "Targeted_therapy",
#                                                     "Chemotherapy", "Targeted_therapy", "Chemotherapy", "Other",
#                                                     "Targeted_therapy", "Radiotherapy", "Targeted_therapy", "Immunotherapy",
#                                                     "Chemotherapy", "Chemotherapy", "Immunotherapy", "Chemotherapy",
#                                                     "Chemotherapy", "Radiotherapy", "Hormone_therapy", "Chemotherapy",
#                                                     "Chemotherapy|Targeted_therapy", "Targeted_therapy", "Immunotherapy", "Chemotherapy",
#                                                     "Hormone_therapy", "Targeted_therapy"))
# 
# 
# write.table(treatment_summary, file = paste0(wd, "r-objects/treatment-summary.tsv"), sep = "\t", quote = F, row.names = F)



treatment_summary <- read.csv(file = paste0(wd, "r-objects/treatment-summary.tsv"), sep = "\t", stringsAsFactors = F, header = T)




tt <- read.csv(file = "/home/ali313/Desktop/info_treatment.tsv",
               sep = "\t", stringsAsFactors = F, header = T)


tt <- merge(tt, treatment_summary, by.x = "mechanism", by.y = "treatment_mech", sort = F)

tt <- tt[,c("sample_id", "treatment_group")]




# bb <- aggregate(data=tt,treatment_group~.,FUN=paste,collapse="|")
# 
# 
# metadata_included$treatment <- F
# metadata_included$treatment[metadata_included$sample_id %in% bb$sample_id] <- T
# 
# metadata_included_hmf <- metadata_included[metadata_included$cohort == "Hartwig",]
# metadata_included_hmf$had_radiotherapy[is.na(metadata_included_hmf$had_radiotherapy)] <- "FALSE"
# metadata_included_hmf$had_radiotherapy <- as.logical(metadata_included_hmf$had_radiotherapy)
# 
# row.names(metadata_included_hmf) <- 1:nrow(metadata_included_hmf)
# 
# 
# metadata_included_hmf$Chemotherapy <- F
# metadata_included_hmf$Hormone_therapy <- F
# metadata_included_hmf$Targeted_therapy <- F
# metadata_included_hmf$Immunotherapy <- F
# metadata_included_hmf$Radiotherapy <- F
# 
# for (i in 1:nrow(metadata_included_hmf)){
#   sample <- metadata_included_hmf$sample_id[i]
#   
#   if (metadata_included_hmf$treatment[metadata_included_hmf$sample_id == sample]) {
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Chemotherapy")) {
#       metadata_included_hmf$Chemotherapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Hormone_therapy")) {
#       metadata_included_hmf$Hormone_therapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Targeted_therapy")) {
#       metadata_included_hmf$Targeted_therapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Immunotherapy")) {
#       metadata_included_hmf$Immunotherapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Radiotherapy")) {
#       metadata_included_hmf$Radiotherapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#   }
#   
#   
# }
# 
# 
# metadata_included_hmf$treatment[!metadata_included_hmf$treatment & metadata_included_hmf$had_radiotherapy] <- T

metadata_included_hmf <- metadata_included[metadata_included$cohort == "Hartwig",]

metadata_included_hmf1 <- metadata_included_hmf

metadata_included_hmf1$cancer_type <- "Pan-cancer"

metadata_included_hmf <- rbind(metadata_included_hmf, metadata_included_hmf1)



df <- data.frame()

for (cancer_type in unique(metadata_included_hmf$cancer_type)){
  tmp <- metadata_included_hmf[metadata_included_hmf$cancer_type == cancer_type,]
  df <- rbind(df, c(cancer_type, nrow(tmp), sum(tmp$treatment_info_available), sum(tmp$had_chemotherapy, na.rm = T), sum(tmp$had_hormone_therapy, na.rm = T), sum(tmp$had_targeted_therapy, na.rm = T), sum(tmp$had_immunotherapy, na.rm = T), sum(tmp$had_radiotherapy, na.rm = T)))
}


colnames(df) <- c("cancer_type", "tot_sample", "tot_treatment", "Chemotherapy", "Hormone_therapy", "Targeted_therapy", "Immunotherapy", "Radiotherapy")





df[,2:8] <- apply(df[,2:8], 2, as.numeric)

df[,3:8] <- 100*df[,3:8]/df[,2]




df_tibb <- gather(df, key = "treatment_group", value = "percent", 3:8)

df_tibb$cancer_type <- factor(df_tibb$cancer_type, levels = included_cancer_type_fig1)
df_tibb$treatment_group <- factor(df_tibb$treatment_group, levels = rev(c("tot_treatment", "Chemotherapy", "Radiotherapy", "Targeted_therapy", "Immunotherapy","Hormone_therapy")))



# sum(metadata_included$treatment_info_available[metadata_included$cohort == "Hartwig"])/nrow(metadata_included[metadata_included$cohort == "Hartwig",])


nn <- ggplot(df_tibb, aes(x = treatment_group,y = percent, fill = treatment_group)) + facet_wrap(~cancer_type, nrow = 24) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values = c("#cc0000", "#999999", "#67a9cf", "#7FC97F", "#ef8a62", "#000000")) +
  # guides(fill=FALSE) +
  ylim(c(0,100)) +
  # scale_y_continuous(breaks = c(0,50,100), labels = c("0", "50", "100")) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(axis.text=element_text(size=10, color = "black"),
        axis.title=element_text(size=25,face="bold")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.spacing.y = unit(1, "lines"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))




pdf(file = "/home/ali313/Desktop/v2/treatment.pdf", height = 17.75, width = 15)
print(nn)
dev.off()



png(file = "/home/ali313/Desktop/treatment.png", height = 1500, width = 1000)
print(nn)
dev.off()

################# treatment info for supp table


treatment_summary <- read.csv(file = paste0(wd, "r-objects/treatment-summary.tsv"), sep = "\t", stringsAsFactors = F, header = T)




tt <- read.csv(file = "/home/ali313/Desktop/info_treatment.tsv",
               sep = "\t", stringsAsFactors = F, header = T)


tt <- merge(tt, treatment_summary, by.x = "mechanism", by.y = "treatment_mech", sort = F)

tt <- tt[,c("sample_id", "cancer_type", "treatment_group")]

bb <- aggregate(data=tt,treatment_group~.,FUN=paste,collapse="|")


bb <- bb[bb$sample_id %in% metadata_included$sample_id,]

metadata_hmf <- metadata_included[metadata_included$cohort == "Hartwig",]


metadata_hmf[,"had_radiotherapy"] <- as.logical(metadata_hmf[,"had_radiotherapy"])


metadata_supp_data <- merge(metadata_hmf, bb[,c(1, 3)], by = "sample_id", all.x = T)

metadata_supp_data$Other <- F
metadata_supp_data$Chemotherapy <- F
metadata_supp_data$Hormone_therapy <- F
metadata_supp_data$Targeted_therapy <- F
metadata_supp_data$Immunotherapy <- F
metadata_supp_data$Radiotherapy <- F

for (i in 1:nrow(metadata_supp_data)){
  print(i)
  
  sample <- metadata_supp_data$sample_id[i]
  
  if (!is.na(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample])){
    if (str_detect(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample], pattern = "Other")) {
      metadata_supp_data$Other[metadata_supp_data$sample_id == sample] <- T
    }
    if (str_detect(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample], pattern = "Chemotherapy")) {
      metadata_supp_data$Chemotherapy[metadata_supp_data$sample_id == sample] <- T
    }
    
    if (str_detect(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample], pattern = "Hormone_therapy")) {
      metadata_supp_data$Hormone_therapy[metadata_supp_data$sample_id == sample] <- T
    }
    
    if (str_detect(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample], pattern = "Targeted_therapy")) {
      metadata_supp_data$Targeted_therapy[metadata_supp_data$sample_id == sample] <- T
    }
    
    if (str_detect(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample], pattern = "Immunotherapy")) {
      metadata_supp_data$Immunotherapy[metadata_supp_data$sample_id == sample] <- T
    }
    
    # if (str_detect(metadata_supp_data$treatment_group[metadata_supp_data$sample_id == sample], pattern = "Radiotherapy")) {
    #   metadata_supp_data$Radiotherapy[metadata_supp_data$sample_id == sample] <- T
    # }
  } else {metadata_supp_data[metadata_supp_data$sample_id == sample,c("Other", "Chemotherapy", "Hormone_therapy", "Targeted_therapy", "Immunotherapy")] <- c(NA, NA, NA, NA, NA)
  }
  
  if (!is.na(metadata_supp_data$had_radiotherapy[metadata_supp_data$sample_id == sample])){
    
    if (metadata_supp_data$had_radiotherapy[metadata_supp_data$sample_id == sample]) {
      metadata_supp_data$Radiotherapy[metadata_supp_data$sample_id == sample] <- T
    }
  } else {
    metadata_supp_data$Radiotherapy[metadata_supp_data$sample_id == sample] <- NA
  }
  
}




metadata_supp_data$treatment_info_available1 <- rowSums(metadata_supp_data[,52:57], na.rm = T) >= 1
metadata_supp_data$treatment_info_available1[is.na(metadata_supp_data$Chemotherapy) & metadata_supp_data$Radiotherapy == F] <- T
metadata_supp_data$treatment_info_available1[rowSums(metadata_supp_data[,52:57], na.rm = T) ==0 & !is.na(metadata_supp_data$Chemotherapy)] <- T


metadata_hmf_supp_table <- metadata_supp_data[,c(1,51:58)]

colnames(metadata_hmf_supp_table)[9] <- "treatment_info_available"

write.table(metadata_hmf_supp_table, file = paste0(wd, "/r-objects/treatment_supp_table_291122.tsv"), sep = "\t", row.names = F, quote = F)





################# Getting the biopsy site info


# this is fran's file. It has some errors!

gg <- read.csv(file = paste0("/home/ali313/Videos/met_location.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")

gg <- gg[gg$sample_id %in% metadata_included$sample_id,]


metadata_included <- merge(metadata_included, gg[,c("sample_id", "major_site_metastatic_simplified")], by = "sample_id", all = T)



cc_type <- included_cancer_type_fig1[12]
cc_type
table(metadata_included$major_site_metastatic_simplified[metadata_included$cancer_type == cc_type & metadata_included$cohort == "HMF"], useNA = "always")
sum(as.numeric(table(metadata_included$major_site_metastatic_simplified[metadata_included$cancer_type == cc_type & metadata_included$cohort == "HMF"])))


# getting the biopsy location data and have it ready

vv <- read.csv(file = paste0(local, "/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/metadata.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")

vv <- vv[vv$sampleId %in% metadata_included$sample_id,]

metadata_included <- merge(metadata_included, vv[,c("sampleId", "biopsySite")], by.x = "sample_id", by.y = "sampleId", all = T)



cc_type <- included_cancer_type_fig1[12]
cc_type
table(metadata_included$biopsySite[metadata_included$cancer_type == cc_type & metadata_included$cohort == "HMF"], useNA = "always")
sum(as.numeric(table(metadata_included$biopsySite[metadata_included$cancer_type == cc_type & metadata_included$cohort == "HMF"])))

metadata_included <- metadata_included[!is.na(metadata_included$biopsySite),]
metadata_included$simplified_biopsiy_site <- NA

for (i in 2:23){
  cc_type <- included_cancer_type_fig1[i]
  
  if ( i == 2){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite != "ovarium"] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "ovarium"] <- "Distant"
  }
  
  if ( i == 3){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Primary"] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  
  if ( i == 4){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 5){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary","Lung")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 6){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 7){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "Breast")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 8){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Primary"] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("null", "other")] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 9){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "Liver")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 10){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "pancreas")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 11){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 12){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$major_site_metastatic_simplified == "Colorectum"] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$major_site_metastatic_simplified %in% c("Unknown", "Other site")] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$major_site_metastatic_simplified == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 13){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "oesophagus")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 14){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "stomach", "Stomach")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 15){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "Kidney")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 16){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "cervix")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 17){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "Ovarium CA", "ovary")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 18){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 19){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("null", "unknown")] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 20){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite %in% c("Primary", "prostate")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 21){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 22){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 23){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$major_site_metastatic_simplified == "Skin/Subcutaneous"] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$major_site_metastatic_simplified == "Unknown"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & metadata_included$major_site_metastatic_simplified == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_type == cc_type & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
}

table(metadata_included$cancer_type, metadata_included$simplified_biopsiy_site)
metadata_included[metadata_included$cancer_type == "Upper respiratory tract cancer" & metadata_included$simplified_biopsiy_site == "Distant",]



# local groups:
## Glioblastoma multiforme: all except for = "ovarium"
## Upper respiratory tract cancer: "just above a. subclavia left side", "Lung", "left mandibula"
## Thyroid cancer: 
## Non small cell lung cancer: "Lung", "primary"
## Diffuse large B-cell lymphoma: 
## breast cancer: "primary", "breast"
## Cholangiocarcinoma: "primary"
## Hepatocellular carcinoma: "liver", primary"
## Pancreas carcinoma: "pancreas", "primary"
## Pancreas neuroendocrine: "liver
## Colorectum carcinoma: fran's file
## Esophagus cancer: "primary", "oesophagus"
## Stomach cancer: "primary", "stomach"
## Kidney clear cell carcinoma: "primary", "kidney"
## Cervix carcinoma: "cervix", "primary"
## Ovarian cancer: "ovary", "ovary ca", "primary
## Uterus carcinoma: "primary"
## Urothelial cancer: "primary"
## Prostate carcinoma: "primary", "prostate
## Leiomyosarcoma: 
## Liposarcoma:
## Skin melanoma: fran's file and lymph the main file

table(metadata_included$cancer_type[metadata_included$cohort == "Hartwig"], metadata_included$metastatic_location[metadata_included$cohort == "Hartwig"], useNA = "always")

biopsy_summary <- data.frame(cancer_type = row.names(table(metadata_included$cancer_type_code, metadata_included$metastatic_location))[-1],
                             distant = c(1, 19, 11, 60, 11, 11, 
                                         493, 56, 11, 89,
                                         30, 423, 98, 26,
                                         59, 18, 78, 20, 66, 
                                         200, 28, 10, 104),
                             local = c(63, 5, 0, 62, 14, 0, 
                                       61, 3, 41, 4,
                                       0, 21, 9, 5,
                                       24, 5, 3, 1, 5, 
                                       16, 0, 0, 81),
                             lymph = c(0, 7, 9, 33, 4, 7, 
                                       114, 2, 0, 2,
                                       0, 19, 20, 6,
                                       14, 6, 23, 6, 42, 
                                       152, 3, 0, 94),
                             unknown = c(0, 7, 2, 100, 5, 2, 
                                         108, 5, 0, 4,
                                         4, 158, 13, 6,
                                         13, 9, 0, 5, 15, 
                                         33, 16, 15, 24))


biopsy_summary <- cbind(biopsy_summary, total = rowSums(biopsy_summary[,2:5]))


biopsy_summary_com <- rbind(c("PAN", sum(biopsy_summary$local), sum(biopsy_summary$lymph), sum(biopsy_summary$distant), sum(biopsy_summary$unknown), sum(biopsy_summary$total)),biopsy_summary)

biopsy_summary_com[,2:6] <- apply(biopsy_summary_com[,2:6], 2, as.numeric)

biopsy_summary_com[,2:5] <- 100*biopsy_summary_com[,2:5]/biopsy_summary_com[,6]

biopsy_summary_tibb <- gather(biopsy_summary_com, key = "biopsy_site", value = "percent", 2:5)

biopsy_summary_tibb$cancer_type <- factor(biopsy_summary_tibb$cancer_type, levels = included_cancer_type_code_fig1)
biopsy_summary_tibb$biopsy_site <- factor(biopsy_summary_tibb$biopsy_site, levels = rev(c("local", "lymph", "distant", "unknown")))

biopsy_summary_tibb$dummy <- 1


# biopsy_summary$total == (biopsy_summary$local + biopsy_summary$lymph + biopsy_summary$distant + biopsy_summary$unknown)

metadata_included$metastatic_location

# sum(!is.na(metadata_included$metastatic_location[metadata_included$cohort == "Hartwig"]))/nrow(metadata_included[metadata_included$cohort == "Hartwig",])
# sum(metadata_included$metastatic_location[metadata_included$cohort == "Hartwig"] == "Distant", na.rm = T)/nrow(metadata_included[metadata_included$cohort == "Hartwig",])


biopsy_plot <- ggplot(biopsy_summary_tibb, aes(x = dummy,y = percent, fill = biopsy_site)) + facet_wrap(~cancer_type, nrow = 24) +
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  scale_fill_manual(values = rev(c("#8f6014", "#14368f", "#8f1460", "#99a3a4"))) +
  ylim(c(0,100)) +
  theme(axis.text.y = element_blank(), # y axis
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), # x axis
        strip.background = element_blank(), # facet boxes
        strip.text.x = element_blank(),
        axis.text=element_text(size=10, color = "black"),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing.y = unit(2.5, "lines"),
        axis.line.x = element_line(color="black", size = 0.5))


pdf(file = "/home/ali313/Desktop/v2/biopsy.pdf", height = 17.75, width = 15)
print(biopsy_plot)
dev.off()




png(file = "/home/ali313/Desktop/biopsy.png", height = 1500, width = 1000)
print(biopsy_plot)
dev.off()







################# biopsy info for supp table





metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_13122021.tsv"), sep = "\t", header = T, stringsAsFactors = F)

metadata_hmf <- metadata[metadata$cohort == "HMF",]


##############################################
## merge gg and vv




gg <- read.csv(file = paste0("/home/ali313/Videos/met_location.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")



vv <- read.csv(file = paste0(local, "/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/metadata.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")

gg_my <- gg[gg$sample_id %in% metadata_hmf$sample_id[metadata_hmf$cancer_type_code %in% c("SKCM", "COREAD")],]


head(gg_my)

for (ss in vv$sampleId){
  if (ss %in% gg_my$sample_id){
    vv$biopsySite[vv$sampleId == ss] <- gg_my$major_site_metastatic_simplified[gg_my$sample_id == ss]
  }
  
}

write.table(vv, file = paste0(wd, "external-files/biopsy_loc_info_raw.tsv"), sep = "\t", row.names = F, quote = F)


#########################################################333

metadata_hmf <- merge(metadata_hmf, gg[,c("sample_id", "major_site_metastatic_simplified")], by = "sample_id")


metadata_hmf <- merge(metadata_hmf, vv[,c("sampleId", "biopsySite")], by.x = "sample_id", by.y = "sampleId", all.x = T)



metadata_hmf$simplified_biopsiy_site <- NA

metadata_hmf$biopsySite[metadata_hmf$cancer_type == cc_type]


for (i in 1:22){
  cc_type <- included_cancer_type[i]
  
  if ( i == 1){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite != "ovarium"] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "ovarium"] <- "Distant"
  }
  
  if ( i == 2){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Primary"] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  
  if ( i == 3){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 4){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary","Lung")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 5){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary","Lung")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 7){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 8){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "Breast")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 9){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Primary"] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("null", "other")] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 10){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "Liver")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 11){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "pancreas")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 12){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 13){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$major_site_metastatic_simplified == "Colorectum"] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$major_site_metastatic_simplified %in% c("Unknown", "Other site")] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$major_site_metastatic_simplified == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 14){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "oesophagus")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 15){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "stomach", "Stomach")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 16){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "Kidney")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 17){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "cervix")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 18){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "Ovarium CA", "ovary")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 19){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 20){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("null", "unknown")] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 21){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite %in% c("Primary", "prostate")] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 22){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 23){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$biopsySite == "null"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 24){
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$major_site_metastatic_simplified == "Skin/Subcutaneous"] <- "Local"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$major_site_metastatic_simplified == "Unknown"] <- "Unknown"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & metadata_hmf$major_site_metastatic_simplified == "Lymph node"] <- "Lymph"
    metadata_hmf$simplified_biopsiy_site[metadata_hmf$cancer_type == cc_type & is.na(metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
}


table(metadata_hmf$cancer_type, metadata_hmf$simplified_biopsiy_site)

metadata_hmf$biopsySite[metadata_hmf$cancer_type %in% c("Colorectal carcinoma", "Skin melanoma")] <- metadata_hmf$major_site_metastatic_simplified[metadata_hmf$cancer_type %in% c("Colorectal carcinoma", "Skin melanoma")]

metadata_hmf_biopsy_supp_table <- metadata_hmf[,c(1,3, 9, 52:53)]


# colnames(metadata_hmf)[21:22] <- c("biopsySite", "simplified_biopsiy_site")
# 
# metadata_hmf_biopsy_supp_table <- metadata_hmf[,c(1,3,11,21:22)]


write.table(metadata_hmf_biopsy_supp_table, file = paste0(wd, "/r-objects/biposySite_supp_table_291122.tsv"), sep = "\t", row.names = F, quote = F)




################################## Automated script for making lollipops for mutations in coding regions


# treatment groups

# treatments <- c("Aromatase_inhibitor", "untreated",
#                 "Anti_EGFR", "untreated",
#                 "Immunotherapy", "untreated", 
#                 "Anti_AR__GnRH", "untreated",
#                 "Alkaloid", "untreated",
#                 "Microtubule_inhibitor", "untreated",
#                 "Selective_ER_modulator", "untreated")
# 
# cancer_type_codes <- c("BRCA", "BRCA", 
#                        "NSCLC", "NSCLC",
#                        "NSCLC", "NSCLC",
#                        "PRAD", "PRAD",
#                        "BRCA", "BRCA", 
#                        "BRCA", "BRCA", 
#                        "OV", "OV")
# genes_of_interest <- c("ESR1", "ESR1", 
#                        "EGFR", "EGFR",
#                        "SMAD4", "SMAD4",
#                        "AR", "AR",
#                        "TLR7", "TLR7", 
#                        "BNC1", "BNC1", 
#                        "KRAS", "KRAS")
# 
# protein_length <- c(595, 595,
#                1210, 1210,
#                552, 552,
#                920, 920)


treatments <- c("Anti_EGFR", "untreated",
                "Anti_AR__GnRH", "untreated",
                "Aromatase_inhibitor", "untreated",
                "Aromatase_inhibitor", "untreated")

cancer_type_codes <- c("LUAD", "LUAD",
                       "PRAD", "PRAD",
                       "BRCA", "BRCA",
                       "BRCA_ERpos", "BRCA_ERpos")

genes_of_interest <- c("EGFR", "EGFR",
                       "AR", "AR",
                       "ESR1", "ESR1",
                       "ESR1", "ESR1")

protein_length <- c(1210, 1210,
                    920, 920,
                    595, 595,
                    595, 595)



for (i in 1:8){
  
  print(i)
  
  treatment <- treatments[i]
  cancer_type_code <- cancer_type_codes[i]
  gene_of_interest <- genes_of_interest[i]
  
  
  mutation_input <- read.csv(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/dndscv/", cancer_type_code, "/", treatment, ".dndscv.annotmuts.tsv.gz"),
                             header = T, stringsAsFactors = F, sep = "\t")
  
  mutations_of_interet <- mutation_input[mutation_input$gene == gene_of_interest,]
  
  if (nrow(mutations_of_interet) > 0){
    
    # Fix the indel labels
    
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) > nchar(mutations_of_interet$mut) & str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Inframe Deletion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) > nchar(mutations_of_interet$mut) & !str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Frameshift Deletion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) < nchar(mutations_of_interet$mut) & str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Inframe Insertion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) < nchar(mutations_of_interet$mut) & !str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Frameshift Insertion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) == nchar(mutations_of_interet$mut) & nchar(mutations_of_interet$mut) > 1] <- "MNV"
    
    
    mutations_of_interet$value <- 1
    
    mutations_of_interet$aachange[mutations_of_interet$aachange == "."] <- mutations_of_interet$impact[mutations_of_interet$aachange == "."]
    
    
    mutations_of_interet <- mutations_of_interet[,c(3, 15, 17, 1:2, 4:14, 16)]
    
    
    colnames(mutations_of_interet)[1:2] <- c("coord", "category")
    
    mutations_of_interet$coord[mutations_of_interet$codonsub != "."] <- str_sub(mutations_of_interet$aachange[mutations_of_interet$codonsub != "."],2,-2)
    mutations_of_interet$coord[mutations_of_interet$codonsub == "."] <- ceiling(as.numeric(lapply(str_split(mutations_of_interet$ntchange[mutations_of_interet$codonsub == "."], pattern = "-"), "[[", 1))/3)
    
    # mutations_of_interet$category[mutations_of_interet$codonsub == "."] <- paste0(mutations_of_interet$category[mutations_of_interet$codonsub == "."], " (", mutations_of_interet$coord[mutations_of_interet$codonsub == "."], ")")
    
    dups_coord <- mutations_of_interet$coord[duplicated(paste(mutations_of_interet$coord, mutations_of_interet$category))]
    
    for (i in 1:length(unique(dups_coord))){
      dups_aachange <- mutations_of_interet$category[mutations_of_interet$coord == unique(dups_coord)[i]]
      for (j in 1:length(unique(dups_aachange)))
        mutations_of_interet$value[mutations_of_interet$coord == unique(dups_coord)[i] & mutations_of_interet$category == unique(dups_aachange)[j]] <- nrow(mutations_of_interet[mutations_of_interet$coord == unique(dups_coord)[i] & mutations_of_interet$category == unique(dups_aachange)[j],])
    }
    
    
    mutations_of_interet <- mutations_of_interet[!duplicated(paste(mutations_of_interet$coord, mutations_of_interet$category)),]
    
    
    # possible_categories <- union(possible_categories, unique(mutations_of_interet$category))
    
    # print(unique(mutations_of_interet$category))
    
  } else {
    mutations_of_interet$value <- numeric(0)
    
    mutations_of_interet <- mutations_of_interet[,c(3, 15, 17, 1:2, 4:14, 16)]
    
    colnames(mutations_of_interet)[1:2] <- c("coord", "category")
    
    mutations_of_interet$coord <- as.character(mutations_of_interet$coord)
    
  }
  
  
  write.table(mutations_of_interet, file = paste0(wd, "r-objects/lollipops/positions/20221123-request/", gene_of_interest, ".", treatment, ".", cancer_type_code, ".tsv"),
              row.names = F, quote = F, sep = "\t")
}




possible_categories <- c("Missense", "Frameshift Deletion", "Frameshift Insertion", "Synonymous", "MNV", "Stop_loss", "Essential_Splice","Nonsense")






### making the lollipops (lollli.R in Desktop directory)


library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


library(mutsneedle)
library(shiny)
library(stringr)


treatments <- c("Anti_EGFR", "untreated",
                "Anti_AR__GnRH", "untreated",
                "Aromatase_inhibitor", "untreated",
                "Aromatase_inhibitor", "untreated")

cancer_type_codes <- c("LUAD", "LUAD",
                       "PRAD", "PRAD",
                       "BRCA", "BRCA",
                       "BRCA_ERpos", "BRCA_ERpos")

genes_of_interest <- c("EGFR", "EGFR",
                       "AR", "AR",
                       "ESR1", "ESR1",
                       "ESR1", "ESR1")

protein_length <- c(1210, 1210,
                    920, 920,
                    595, 595,
                    595, 595)




i <- 3
treatment <- treatments[i]
cancer_type_code <- cancer_type_codes[i]
gene_of_interest <- genes_of_interest[i]
length_of_protein <- protein_length[i]

genes_id <- as.character(unlist(mget(x=genes_of_interest,envir=org.Hs.egALIAS2EG)))

gene_id <- genes_id[i]

mutations_of_interet <- read.csv(file = paste0(wd, "r-objects/lollipops/positions/20221123-request/", gene_of_interest, ".", treatment, ".", cancer_type_code, ".tsv"),
                                 header = T, stringsAsFactors = F, sep = "\t")


mutations_of_interet <- mutations_of_interet[mutations_of_interet$category != "Essential_Splice",]



mutations_of_interet$coord <- as.character(mutations_of_interet$coord)

regiondata <- exampleRegionData()


regiondata[1,] <- c(paste0("0-", length_of_protein), "Gene")

regiondata <- regiondata[1,]

if (i %in% c(5:8)){
  
  regiondata[2,] <- c("42-181", "Oest_recep")
  regiondata[3,] <- c("183-252", "zf_C4")
  regiondata[4,] <- c("332-531", "Hormone_recep")
  regiondata[5,] <- c("552-595", "ESR1_C")
} else if (i %in% c(1,2)){
  regiondata[2,] <- c("57-168", "Recep_L_domain")
  regiondata[3,] <- c("177-338", "Furin-like")
  regiondata[4,] <- c("361-481", "Recep_L_domain")
  regiondata[5,] <- c("505-637", "GF_recep_IV")
  regiondata[6,] <- c("712-968", "PK_Tyr_Ser-Thr")
} else if (i %in% c(3,4)){
  regiondata[2,] <- c("6-449", "Androgen_recep")
  regiondata[3,] <- c("558-627", "zf-C4")
  regiondata[4,] <- c("690-881", "Hormone_recep")
} # else if (i %in% c(...)){regiondata[2,] <- c("36-137", "MH1")
# regiondata[3,] <- c("321-530", "MH2")}




category_vec <- c("Nonsense", "Missense", "Synonymous", "Frameshift Deletion", "Inframe Deletion", "Frameshift Insertion", "Inframe Insertion",
                  "MNV", "Stop_loss")



for (k in 1:length(category_vec)){
  if (category_vec[k] %notin% mutations_of_interet$category){
    mutations_of_interet <- rbind(mutations_of_interet, c(5000, category_vec[k], 0))
  }
}



index <- order(match(mutations_of_interet$category, category_vec, nomatch = F))
index <- index[index != 0]

mutations_of_interet <- mutations_of_interet[index,]


j <- 0

if (i %in% c(2)){
  while (j < 19){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[2,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} else if (i %in% c(4)){
  while (j < 20){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[2,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} else if (i %in% c(6)){
  while (j < 46){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[2,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} else if (i %in% c(8)){
  while (j < 40){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[6,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} 


row.names(mutations_of_interet) <- 1:nrow(mutations_of_interet)

mutations_of_interet1 <- mutations_of_interet[1,]


mm <- shinyApp(
  ui = mutsneedleOutput("id",width=800,height=500),
  server = function(input, output) {
    output$id <- renderMutsneedle(
      mutsneedle(mutdata=mutations_of_interet1,domains=regiondata,
                 gene = gene_of_interest,
                 # colorMap = ccc,
                 xlab="Protein Location",
                 ylab="Mutation Frequency",
                 maxlength = length_of_protein)
    )
  }
)



paste0(wd, "r-objects/lollipops/figs/20221123-request/", gene_of_interest, ".", treatment, ".", cancer_type_code)

mm

table(mutations_of_interet$aachange)
nrow(mutations_of_interet)


##########################################################3
# COHORT OVERVIEW SUPPLEMENTARY


meta <- read.csv(file = paste0(wd, "external-files/SuppTable1_sample_metadata - metadata-v2-13112022.tsv"), sep = "\t", stringsAsFactors = F, header = T)


metadata_included <- meta[!meta$is_blacklisted,]
metadata_included <- metadata_included[!is.na(metadata_included$cancer_type),]
nrow(metadata_included)

row.names(metadata_included) <- 1:nrow(metadata_included)


metadata_included[,"cohort"] <- factor(metadata_included[,"cohort"], levels = c("PCAWG", "Hartwig"))
metadata_included[,"tissue_group"] <- factor(metadata_included[,"tissue_group"], levels = included_tissue_group_fig1)
metadata_included[,"cancer_type"] <- factor(metadata_included[,"cancer_type"], levels = included_cancer_type_fig1)
metadata_included[,"cancer_type_code"] <- factor(metadata_included[,"cancer_type_code"],
                                                 levels = included_cancer_type_code_fig1)

metadata_included[,"gender"] <- factor(metadata_included[,"gender"], levels = c("MALE", "FEMALE"))



table(metadata_included2$cancer_subtype)

# 
# metadata <- read.csv(file = paste0(wd, "external-files/SuppTable1_sample_metadata - metadata  - with subtypes.tsv"), sep = "\t", stringsAsFactors = F, header = T)
# 
# 
# metadata_included <- metadata[!metadata$is_blacklisted,]


metadata_included2 <- metadata_included[metadata_included$cancer_type_code %in% c("BRCA", "COREAD", "UCEC"),]


metadata_included2$cancer_subtype[is.na(metadata_included2$cancer_subtype) & metadata_included2$cancer_type_code == "BRCA"] <- "NA"
metadata_included2$cancer_subtype[is.na(metadata_included2$cancer_subtype) & metadata_included2$cancer_type_code == "COREAD"] <- "NA"
metadata_included2$cancer_subtype[is.na(metadata_included2$cancer_subtype) & metadata_included2$cancer_type_code == "UCEC"] <- "NA"


metadata_included2$cancer_subtype[metadata_included2$cancer_subtype == "MSS" & metadata_included2$cancer_type_code == "COREAD"] <- "COREAD_MSS"
metadata_included2$cancer_subtype[metadata_included2$cancer_subtype == "MSS" & metadata_included2$cancer_type_code == "UCEC"] <- "UCEC_MSS"



metadata_included2$cancer_subtype[metadata_included2$cancer_subtype == "MSI/POLE" & metadata_included2$cancer_type_code == "COREAD"] <- "COREAD_MSI/POLE"
metadata_included2$cancer_subtype[metadata_included2$cancer_subtype == "MSI/POLE" & metadata_included2$cancer_type_code == "UCEC"] <- "UCEC_MSI/POLE"




metadata_included3 <- metadata_included2


metadata_included3$cancer_subtype <- metadata_included3$cancer_type_code


metadata_included_final <- rbind(metadata_included2, metadata_included3)

metadata_included_final <- metadata_included_final[metadata_included_final$cancer_subtype != "NA",]


metadata_included_final$cancer_subtype <- factor(metadata_included_final$cancer_subtype, levels = c(
  "BRCA", "ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "TNBC",
  "COREAD", "COREAD_MSS", "COREAD_MSI/POLE",
  "UCEC", "UCEC_MSS", "UCEC_MSI/POLE"))


included_cancer_subtype_fig1 <- levels(metadata_included_final$cancer_subtype)



color_palette_cohort <- c("#f58134", "#9966CC")


row.names(metadata_included_final) <- 1:nrow(metadata_included_final)


metadata_included_final[,"cohort"] <- factor(metadata_included_final[,"cohort"], levels = c("PCAWG", "Hartwig"))

metadata_included_final[,"gender"] <- factor(metadata_included_final[,"gender"], levels = c("MALE", "FEMALE"))



sex_df <- data.frame(cancer_subtype = character(22), cohort = character(22), male = numeric(22), female = numeric(22))


for (i in 1:11){
  for (j in 1:2){
    
    cancer_subtype <- included_cancer_subtype_fig1[i]
    cohort <- c("PCAWG", "Hartwig")[j]
    
    
    sex_df[i+(i-1) + (j-1),1:2] <- c(cancer_subtype, cohort)
    if (i != 1){
      sex_df[i+(i-1) + (j-1),3:4] <- 100*c(nrow(metadata_included_final[metadata_included_final$cancer_subtype == cancer_subtype & metadata_included_final$cohort == cohort &  metadata_included_final$gender == "MALE",]),
                                           nrow(metadata_included_final[metadata_included_final$cancer_subtype == cancer_subtype & metadata_included_final$cohort == cohort &  metadata_included_final$gender == "FEMALE",]))/nrow(metadata_included_final[metadata_included_final$cancer_subtype == cancer_subtype & metadata_included_final$cohort == cohort,])
    } else {
      sex_df[i+(i-1) + (j-1),3:4] <- 100*c(nrow(metadata_included_final[metadata_included_final$cohort == cohort &  metadata_included_final$gender == "MALE",]),
                                           nrow(metadata_included_final[metadata_included_final$cohort == cohort &  metadata_included_final$gender == "FEMALE",]))/nrow(metadata_included_final[metadata_included_final$cohort == cohort,])
    }
  }
}




sex_df <- tidyr::gather(sex_df, key = "gender", value = "proportion", 3:4)
sex_df$cancer_subtype <- factor(sex_df$cancer_subtype, levels = included_cancer_subtype_fig1)
sex_df$gender <- factor(sex_df$gender, levels = c("male", "female"))
sex_df$cohort <- factor(sex_df$cohort, levels = c("Hartwig", "PCAWG"))
sex_df$proportion[is.nan(sex_df$proportion)] <- 0


kk <- ggplot(sex_df, aes(x = cohort,y = proportion, fill = gender)) + facet_wrap(~cancer_subtype, nrow = 32) +
  geom_bar(stat="identity") +
  scale_color_manual(values = rev(color_palette_cohort), labels = c("Hartwig", "PCAWG")) +
  scale_fill_manual(values = c("#8498f8", "#fc7571")) +
  guides(fill=FALSE) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, by = 50), labels = c("0", "50", "100%")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25,face="bold")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

kk





pdf(file = "/home/ali313/Desktop/v2/subtypes/sex.pdf", height = 16.75, width = 14)
print(kk)
dev.off()





################# Getting the age info

metadata_included_final


pcawg_meta <- readxl::read_xlsx(path = paste0(wd, "external-files/pcawg-clinical-metadata-20211112.xlsx"))
pcawg_meta <- as.data.frame(pcawg_meta)



hmf_meta <- read.csv(file = paste0(wd, "external-files/hmf-clinical-metadata-20211112.tsv"), stringsAsFactors = F, header = T, sep = "\t")

hmf_meta$biopsyDateProcessed <- NA
hmf_meta$biopsyDateProcessed <- unlist(lapply(str_split(hmf_meta$biopsyDate, pattern = "-"), "[[", 1))

hmf_meta$age_of_donor <- NA

hmf_meta[,c("birthYear","biopsyDateProcessed")] <- apply(hmf_meta[,c("birthYear","biopsyDateProcessed")], as.numeric, MARGIN = 2)

hmf_meta$age_of_donor <- hmf_meta$biopsyDateProcessed - hmf_meta$birthYear

# bnn <- merge(metadata_included[metadata_included$cohort == "PCAWG",], pcawg_meta[,c("icgc_donor_id","donor_age_at_diagnosis")], by.x = "sample_id", by.y = "icgc_donor_id")
wgd_original <- read.csv(file = paste0(wd, "external-files/pcawg-original-wgd-timing.txt"), header = T, stringsAsFactors = F, sep = "\t")
# wgd_original <- wgd_original[!duplicated(wgd_original$icgc_donor_id),]

bnn <- merge(metadata_included_final[metadata_included_final$cohort == "PCAWG",], wgd_original[,c("icgc_donor_id","age")], by.x = "sample_id", by.y = "icgc_donor_id")
colnames(bnn)[51] <- "age"

ann <- merge(metadata_included_final[metadata_included_final$cohort == "Hartwig",], hmf_meta[,c("sampleId","age_of_donor")], by.x = "sample_id", by.y = "sampleId")
colnames(ann)[51] <- "age"


metadata_included <- rbind(bnn, ann)
table(metadata_included$age[metadata_included$cohort == "Hartwig"], useNA = "always")
table(metadata_included$age[metadata_included$cohort == "PCAWG"], useNA = "always")


metadata_included[,"cohort"] <- factor(metadata_included[,"cohort"], levels = c("PCAWG", "Hartwig"))
metadata_included <- metadata_included[!duplicated(metadata_included),]



metadata_included$cancer_subtype <- factor(metadata_included$cancer_subtype, levels = c(
  "NSCLC", "LUAD", "LUSC", #"NOS", 
  "NSCLC_NA",
  "BRCA", "ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "TNBC","BRCA_NA",
  "COREAD", "COREAD_MSS", "COREAD_MSI/POLE",
  "UCEC", "UCEC_MSS", "UCEC_MSI/POLE"))


included_cancer_subtype_fig1 <- levels(metadata_included$cancer_subtype)



library(plyr)
p_meds <- ddply(metadata_included, .(cancer_subtype, cohort), summarise, med = median(age, na.rm = T))


levels(metadata_included$cancer_subtype)

jj <- ggplot(metadata_included, aes(x = age, color = cohort))+ facet_wrap(~ cancer_subtype, ncol = 1)+
  geom_density(aes(color = cohort), size = 3) +
  scale_color_manual(values = c("#E66711","#7217C5")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  geom_vline(data=p_meds, aes(xintercept=med, colour=cohort),
             linetype=1, size=5) +
  guides(color = F)



pdf(file = "/home/ali313/Desktop/v2/subtypes/age.pdf", height = 16.75, width = 14)
print(jj)
dev.off()




################# Getting the treatment info


metadata_included_final

treatment_summary <- read.csv(file = paste0(wd, "r-objects/treatment-summary.tsv"), sep = "\t", stringsAsFactors = F, header = T)




tt <- read.csv(file = "/home/ali313/Desktop/info_treatment.tsv",
               sep = "\t", stringsAsFactors = F, header = T)


tt <- merge(tt, treatment_summary, by.x = "mechanism", by.y = "treatment_mech", sort = F)

tt <- tt[,c("sample_id", "treatment_group")]


# 
# bb <- aggregate(data=tt,treatment_group~.,FUN=paste,collapse="|")
# 
# 
# metadata_included_final$treatment <- F
# metadata_included_final$treatment[metadata_included_final$sample_id %in% bb$sample_id] <- T
# 
# metadata_included_hmf <- metadata_included_final[metadata_included_final$cohort == "Hartwig",]
# metadata_included_hmf$had_radiotherapy[is.na(metadata_included_hmf$had_radiotherapy)] <- "FALSE"
# metadata_included_hmf$had_radiotherapy <- as.logical(metadata_included_hmf$had_radiotherapy)
# 
# row.names(metadata_included_hmf) <- 1:nrow(metadata_included_hmf)
# 
# 
# metadata_included_hmf$Chemotherapy <- F
# metadata_included_hmf$Hormone_therapy <- F
# metadata_included_hmf$Targeted_therapy <- F
# metadata_included_hmf$Immunotherapy <- F
# metadata_included_hmf$Radiotherapy <- F
# 
# for (i in 1:nrow(metadata_included_hmf)){
#   sample <- metadata_included_hmf$sample_id[i]
#   
#   if (metadata_included_hmf$treatment[metadata_included_hmf$sample_id == sample]) {
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Chemotherapy")) {
#       metadata_included_hmf$Chemotherapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Hormone_therapy")) {
#       metadata_included_hmf$Hormone_therapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Targeted_therapy")) {
#       metadata_included_hmf$Targeted_therapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Immunotherapy")) {
#       metadata_included_hmf$Immunotherapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#     
#     if (str_detect(bb$treatment_group[bb$sample_id == sample], pattern = "Radiotherapy")) {
#       metadata_included_hmf$Radiotherapy[metadata_included_hmf$sample_id == sample] <- T
#     }
#   }
#   
#   
# }
# 
# 
# metadata_included_hmf$treatment[!metadata_included_hmf$treatment & metadata_included_hmf$had_radiotherapy] <- T


metadata_included_hmf <- metadata_included[metadata_included$cohort == "Hartwig",]


metadata_included_hmf1 <- metadata_included_hmf

metadata_included_hmf1$cancer_type <- "Pan-cancer"

metadata_included_hmf <- rbind(metadata_included_hmf, metadata_included_hmf1)

df <- data.frame()


for (cancer_subtype in unique(metadata_included_hmf$cancer_subtype)){
  tmp <- metadata_included_hmf[metadata_included_hmf$cancer_subtype == cancer_subtype,]
  df <- rbind(df, c(cancer_subtype, nrow(tmp), sum(tmp$treatment_info_available), sum(tmp$had_chemotherapy, na.rm = T), sum(tmp$had_hormone_therapy, na.rm = T), sum(tmp$had_targeted_therapy, na.rm = T), sum(tmp$had_immunotherapy, na.rm = T), sum(tmp$had_radiotherapy, na.rm = T)))
}


colnames(df) <- c("cancer_subtype", "tot_sample", "tot_treatment", "Chemotherapy", "Hormone_therapy", "Targeted_therapy", "Immunotherapy", "Radiotherapy")



df[,2:8] <- apply(df[,2:8], 2, as.numeric)

df[,3:8] <- 100*df[,3:8]/df[,2]




df_tibb <- gather(df, key = "treatment_group", value = "percent", 3:8)

df_tibb$cancer_subtype <- factor(df_tibb$cancer_subtype, levels = included_cancer_subtype_fig1)
df_tibb$treatment_group <- factor(df_tibb$treatment_group, levels = rev(c("tot_treatment", "Chemotherapy", "Radiotherapy", "Targeted_therapy", "Immunotherapy","Hormone_therapy")))


nn <- ggplot(df_tibb, aes(x = treatment_group,y = percent, fill = treatment_group)) + facet_wrap(~cancer_subtype, nrow = 11) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values = c("#cc0000", "#999999", "#67a9cf", "#7FC97F", "#ef8a62", "#000000")) +
  # guides(fill=FALSE) +
  ylim(c(0,100)) +
  # scale_y_continuous(breaks = c(0,50,100), labels = c("0", "50", "100")) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(axis.text=element_text(size=10, color = "black"),
        axis.title=element_text(size=25,face="bold")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.spacing.y = unit(1, "lines"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))





pdf(file = "/home/ali313/Desktop/v2/subtypes/treatment.pdf", height = 17.75, width = 15)
print(nn)
dev.off()




################# Getting the biopsy site info
metadata_included_final


# getting the biopsy location data and have it ready

vv <- read.csv(file = paste0(local, "/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/metadata.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")

vv <- vv[vv$sampleId %in% metadata_included_final$sample_id,]

metadata_included <- merge(metadata_included_final, vv[,c("sampleId", "biopsySite")], by.x = "sample_id", by.y = "sampleId", all = T)


included_cancer_type_fig1
included_cancer_subtype_fig1

metadata_included <- metadata_included[!is.na(metadata_included$biopsySite),]
metadata_included$simplified_biopsiy_site <- NA


for (i in 1:16){
  cc_subtype <- included_cancer_subtype_fig1[i]
  
  
  if ( i %in% 1:4){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite %in% c("Primary","Lung")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i %in% 5:10){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite %in% c("Primary", "Breast")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  
  if ( i %in% 11:13){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$major_site_metastatic_simplified == "Colorectum"] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$major_site_metastatic_simplified %in% c("Unknown", "Other site")] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$major_site_metastatic_simplified == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
  
  if ( i %in% 14:16){
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite %in% c("Primary")] <- "Local"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite == "null"] <- "Unknown"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & metadata_included$biopsySite == "Lymph node"] <- "Lymph"
    metadata_included$simplified_biopsiy_site[metadata_included$cancer_subtype == cc_subtype & is.na(metadata_included$simplified_biopsiy_site)] <- "Distant"
  }
  
}

table(metadata_included$cancer_subtype, metadata_included$simplified_biopsiy_site)

metadata_included[metadata_included$cancer_subtype == cc_subtype,]



gg <- read.csv(file = paste0("/home/ali313/Videos/met_location.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")

gg <- gg[gg$sample_id %in% metadata_included$sample_id,]


metadata_included <- merge(metadata_included_final, gg[,c("sample_id", "major_site_metastatic_simplified")], by = "sample_id", all = T)


cancer_subtype <- "COREAD_MSI/POLE"
table(metadata_included$major_site_metastatic_simplified[metadata_included$cancer_subtype == cancer_subtype & metadata_included$cohort == "Hartwig"], useNA = "always")
sum(as.numeric(table(metadata_included$major_site_metastatic_simplified[metadata_included$cancer_subtype == cancer_subtype & metadata_included$cohort == "Hartwig"])))


table(metadata_included$cancer_subtype, metadata_included$major_site_metastatic_simplified)





# local groups:
## Glioblastoma multiforme: all except for = "ovarium"
## Upper respiratory tract cancer: "just above a. subclavia left side", "Lung", "left mandibula"
## Thyroid cancer: 
## Non small cell lung cancer: "Lung", "primary"
## Diffuse large B-cell lymphoma: 
## breast cancer: "primary", "breast"
## Cholangiocarcinoma: "primary"
## Hepatocellular carcinoma: "liver", primary"
## Pancreas carcinoma: "pancreas", "primary"
## Pancreas neuroendocrine: "liver
## Colorectum carcinoma: fran's file
## Esophagus cancer: "primary", "oesophagus"
## Stomach cancer: "primary", "stomach"
## Kidney clear cell carcinoma: "primary", "kidney"
## Cervix carcinoma: "cervix", "primary"
## Ovarian cancer: "ovary", "ovary ca", "primary
## Uterus carcinoma: "primary"
## Urothelial cancer: "primary"
## Prostate carcinoma: "primary", "prostate
## Leiomyosarcoma: 
## Liposarcoma:
## Skin melanoma: fran's file and lymph the main file

sum(as.matrix(table(metadata_included$cancer_subtype, metadata_included$major_site_metastatic_simplified))[13,])

585  - 19 - 15 - 137

# biopsy_summary <- data.frame(cancer_subtype = included_cancer_subtype_fig1,
#                              local = c(129, 52, 11, 7, 59, 61, 11, 35, 5,
#                                        9, 1, 21, 19, 2, 1, 0, 1),
#                              lymph = c(63, 26, 2, 1, 34, 114, 13, 61, 8,
#                                        22, 10, 19, 15, 4, 6, 4, 2),
#                              distant = c(147, 47, 10, 12, 78, 493, 52, 354, 18,
#                                          51, 18, 423, 407, 16, 21, 20, 1),
#                              unknown = c(166, 77, 4, 2, 83, 109, 6, 61, 3,
#                                          26, 13, 158, 144, 14, 5, 2, 3))

biopsy_summary <- data.frame(cancer_subtype = included_cancer_subtype_fig1,
                             local = c(76, 50, 11, 15, 61, 11, 35, 5,
                                       9, 1, 21, 19, 2, 1, 0, 1),
                             lymph = c(37, 26, 2, 9, 114, 13, 61, 8,
                                       22, 10, 19, 15, 4, 6, 4, 2),
                             distant = c(71, 46, 10, 15, 493, 52, 354, 18,
                                         51, 18, 431, 414, 17, 20, 19, 1),
                             unknown = c(105, 77, 4, 24, 108, 6, 61, 3,
                                         26, 12, 150, 137, 13, 5, 2, 3))


biopsy_summary$total <- rowSums(biopsy_summary[,2:5])




biopsy_summary[,2:6] <- apply(biopsy_summary[,2:6], 2, as.numeric)

biopsy_summary[,2:5] <- 100*biopsy_summary[,2:5]/biopsy_summary[,6]

biopsy_summary_tibb <- gather(biopsy_summary, key = "biopsy_site", value = "percent", 2:5)

biopsy_summary_tibb$cancer_subtype <- factor(biopsy_summary_tibb$cancer_subtype, levels = included_cancer_subtype_fig1)
biopsy_summary_tibb$biopsy_site <- factor(biopsy_summary_tibb$biopsy_site, levels = rev(c("local", "lymph", "distant", "unknown")))

biopsy_summary_tibb$dummy <- 1


# biopsy_summary$total == (biopsy_summary$local + biopsy_summary$lymph + biopsy_summary$distant + biopsy_summary$unknown)


biopsy_plot <- ggplot(biopsy_summary_tibb, aes(x = dummy,y = percent, fill = biopsy_site)) + facet_wrap(~cancer_subtype, nrow = 16) +
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  scale_fill_manual(values = rev(c("#8f6014", "#14368f", "#8f1460", "#99a3a4"))) +
  ylim(c(0,100)) +
  theme(axis.text.y = element_blank(), # y axis
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), # x axis
        strip.background = element_blank(), # facet boxes
        strip.text.x = element_blank(),
        axis.text=element_text(size=10, color = "black"),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing.y = unit(2.5, "lines"),
        axis.line.x = element_line(color="black", size = 0.5))


pdf(file = "/home/ali313/Desktop/v2/subtypes/biopsy.pdf", height = 17.75, width = 15)
print(biopsy_plot)
dev.off()


table(metadata_included_final$cancer_subtype, metadata_included_final$cohort)


