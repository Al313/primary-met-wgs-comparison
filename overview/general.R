







`%notin%` <- Negate(`%in%`)


## Loading basic packages --------------------------------
library(readxl)
library(tidyr)
library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(magrittr)


## Defining values --------------------------------

cancer_type_order <- c("Pan-cancer", "Glioblastoma multiforme", "Upper respiratory tract carcinoma", "Thyroid carcinoma",
                               "Lung adenocarcinoma", "Lung squamous cell carcinoma","Diffuse large B-cell lymphoma", "Breast carcinoma",
                               "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
                               "Pancreas neuroendocrine", "Colorectal carcinoma", "Esophageal carcinoma", "Stomach carcinoma", 
                               "Kidney renal clear cell carcinoma", "Cervical carcinoma", "Ovarian serous adenocarcinoma", "Uterus carcinoma", 
                               "Bladder urothelial carcinoma", "Prostate carcinoma", "Leiomyosarcoma", "Liposarcoma", "Skin melanoma")

cancer_type_code_order <- c("PAN", "GBM", "HNSC", "THCA", "LUAD", "LUSC", "DLBCL", "BRCA", "CHOL", "LIHC", "PAAD", "PANET", "COREAD", "ESCA",
                                    "STAD", "KIRC", "CESC", "OV", "UCEC", "BLCA", "PRAD", "LMS", "LPS", "SKCM")


tissue_group_order <- c("All tissues", "CNS", "Head_and_neck", "Thyroid", "Lung", "Breast", "Lymphoid", "GI_dev", "GI_core", "Kidney", "Gyn",
                                "Urothelial", "Prostate", "Soft_tissue", "Skin")


cohort_order <- c("PCAWG", "Hartwig")


cohort_color_palette <- c("#f58134", "#9966CC")


