## ---------------------------
##
## Script name:
##
## Purpose of script: Running pyLMM on UCLA and GigaMuga density imputed grid for CC mice
##
## Author: Hao He
##
## Date Created: 2023-06-14
##
## Copyright (c) Hao He, 2023
## Email: hao.he@jax.org
##
## ---------------------------
##
## Notes:
## 1, genotypes are in /projects/chesler-lab/phenome/metasoft/data/
## GM_marker_positions_with_MUSter_v2_genotypes.RDS
## UCLA1_positions_with_MUSter_v2_genotypes.RDS
## 2, the measure ids/names that need to be run for addiction and high-fat diet meta-analyses are here:
## /projects/chesler-lab/phenome/metasoft/manuscript/data/forHao
## high-fat_diet_studynames.RDS
## addiction_studynames.RDS
##
## ---------------------------

# 1, Set working directory -----------------------------------------------------------------------
setwd("/projects/csna/meta-analysis/")    # working directory

# 2, Load up the packages we will need -----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
library(qtl2) #(0.18 version)
library(foreach)
library(doParallel)
library(parallel)
library(abind)
library(tidyverse)
library(data.table)
library(stringr)
library(DOQTL)
source("code/meta_utils.R")
#rankz transform
rz.transform <- function(y){
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))
  return(rzT)
}

# 3, Create folder -----------------------------------------------------------------------
#create t2d_addiction in data folder for ucla grid
if (!dir.exists("data/t2d_addiction/ucla"))
{ dir.create("data/t2d_addiction/ucla", recursive = TRUE) } #create folder
#create t2d_addiction in data folder for gigamuga grid
if (!dir.exists("data/t2d_addiction/gigamuga"))
{ dir.create("data/t2d_addiction/gigamuga", recursive = TRUE) } #create folder

#create t2d_addiction in analysis folder
if (!dir.exists("analysis/t2d_addiction/ucla/scripts"))
{ dir.create("analysis/t2d_addiction/ucla/scripts", recursive = TRUE) } #create folder
if (!dir.exists("analysis/t2d_addiction/ucla/job"))
{ dir.create("analysis/t2d_addiction/ucla/job", recursive = TRUE) } #create folder

if (!dir.exists("analysis/t2d_addiction/gigamuga/scripts"))
{ dir.create("analysis/t2d_addiction/gigamuga/scripts", recursive = TRUE) } #create folder
if (!dir.exists("analysis/t2d_addiction/gigamuga/job"))
{ dir.create("analysis/t2d_addiction/gigamuga/job", recursive = TRUE) } #create folder

#create t2d_addiction output folder for ucla grid
if (!dir.exists("output/t2d_addiction/pylmm/ucla"))
{ dir.create("output/t2d_addiction/pylmm/ucla", recursive = TRUE) } #create folder
if (!dir.exists("output/t2d_addiction/gemma/ucla"))
{ dir.create("output/t2d_addiction/gemma/ucla", recursive = TRUE) } #create folder

#create t2d_addiction output folder for gigamuga grid
if (!dir.exists("output/t2d_addiction/pylmm/gigamuga"))
{ dir.create("output/t2d_addiction/pylmm/gigamuga", recursive = TRUE) } #create folder
if (!dir.exists("output/t2d_addiction/gemma/gigamuga"))
{ dir.create("output/t2d_addiction/gemma/gigamuga", recursive = TRUE) } #create folder

# 4, Read and format the study name -----------------------------------------------------------------------
hf_diet = readRDS("/projects/chesler-lab/phenome/metasoft/manuscript/data/forHao/high-fat_diet_studynames.RDS")
addiction = readRDS("/projects/chesler-lab/phenome/metasoft/manuscript/data/forHao/addiction_studynames.RDS")
study_meta = rbind(data.frame(study_name = hf_diet, disease = "t2d"),
                   data.frame(study_name = addiction, disease = "addiction"))

study_meta = study_meta %>%
  dplyr::mutate(cohort = case_when(
    str_detect(study_name, "csna03_cc") & !str_detect(study_name, "price") ~ "cc",
    str_detect(study_name, "csna03_do") ~ "do",
    str_detect(study_name, "price_rz")  ~ "price.bxd",
    str_detect(study_name, ".price") ~ "price.cc",
    TRUE ~ "mpd"
  ))

# 5, Work on CC study based on UCLA grid -----------------------------------------------------------------------
cc <- study_meta %>%
  dplyr::filter(cohort == "cc") %>%
  dplyr::mutate(trait = gsub("cocaine_csna03_cc_", "", study_name))

#cc genotypes ucla
cc.geno <- readRDS("/projects/chesler-lab/phenome/metasoft/data/UCLA1_positions_with_MUSter_v2_genotypes.RDS") %>%
  dplyr::select(1:4, all_of("C57BL/6J"), starts_with("CC")) %>%
  dplyr::rename(b_6 = 5) %>%
  dplyr::select(-all_of("CC057/Unc")) %>%
  dplyr::rename_with(., ~gsub("\\/.*", "", .x)) %>%
  dplyr::mutate(across(6:last_col(), ~ case_when(
    . == "N" ~ "0 0",
    . == b_6 ~ "1 1",
    . == "H" ~ "1 2",
    TRUE ~ "2 2"
  ))) %>%
  dplyr::select(-all_of("b_6"))

#generate map file, four columns, chr, rs, 0, pos in Base-pair coordinate
map <- cc.geno %>%
  dplyr::select(chr, rs, bp38) %>%
  mutate(cM = 0, .after = rs) %>%
  mutate(rs = case_when(
    (rs=="") ~ paste0(chr, "_", bp38),
    TRUE ~ as.character(rs)
  ))
write.table(map, file = paste0("data/t2d_addiction/ucla/", "cc.f.map"), #female
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(map, file = paste0("data/t2d_addiction/ucla/", "cc.m.map"), #male
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#cc.geno
cc.geno <- cc.geno %>%
  dplyr::select(5:last_col()) %>%
  t(.) %>%
  as_tibble(.name_repair = "universal") %>%
  mutate(FID = colnames(cc.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         IID = colnames(cc.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         pid = 0,
         mid = 0,
         phe = -9,
         .before = 1)

# cc phenotype
# cc phenotypes ivsa
cc.pheno.ivsa <- read_csv("data/pheno/Prj04_IVSA_preqc_12072020_CC.csv",
                          col_names = TRUE, na = c("", "NA", "N/A", "NaN")) %>%
  rename_with(~str_remove(.x, "\\_cc|\\_do")) %>%
  rename_with(~(.x = "Strain"),   contains("_Strain")) %>%
  rename_with(~(.x = "Mouse.ID"), contains("Subject")) %>%
  rename_with(~(.x = "Sex"),      contains("Sex")) %>%
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male",
    TRUE       ~ as.character(Sex))) %>% #make Sex coding consistent
  mutate(across(Mouse.ID, as.character)) %>%
  dplyr::select(Mouse.ID, Strain, Sex, 25:last_col()) %>%
  dplyr::mutate(Strain = gsub("\\/.*", "", Strain))

# cc phenotypes sens
cc.pheno.sens <- read_csv("data/pheno/Prj02_Sensitization_preqc_12282020_CC.csv",
                          col_names = TRUE, na = c("", "NA", "N/A", "NaN")) %>%
  rename_with(~str_remove(.x, "\\_cc|\\_do")) %>%
  rename_with(~(.x = "Strain"),   contains("_Strain")) %>%
  rename_with(~(.x = "Mouse.ID"), contains("Subject")) %>%
  rename_with(~(.x = "Sex"),      contains("Sex")) %>%
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male",
    TRUE       ~ as.character(Sex))) %>% #make Sex coding consistent
  mutate(across(Mouse.ID, as.character)) %>%
  dplyr::select(Mouse.ID, Strain, Sex, 8:last_col()) %>%
  dplyr::mutate(Strain = gsub("\\/.*", "", Strain))

#full join
cc.pheno <- full_join(cc.pheno.ivsa, cc.pheno.sens) %>%
  dplyr::select(1:3, unique(substr(cc$trait, 1, nchar(cc$trait) - 2))) %>%
  dplyr::mutate(across(AQ_SessionsToAcquisition:last_col(), rz.transform))

#left_join cc.pheno and cc.geno
df.pheno.geno <- left_join(cc.pheno, cc.geno, by = c("Strain" = "IID")) %>%
  mutate(IID = FID, .after = FID) %>%
  relocate(Sex, .after = mid) %>%
  arrange(Mouse.ID) %>%
  mutate(Sex = case_when(
    Sex == "Male" ~ 1,
    Sex == "Female" ~ 2
  )) %>%
  mutate(FID = Mouse.ID) %>% #pylmm cannot handle duplicated ids
  mutate(IID = Mouse.ID)

#ped file
ped <- df.pheno.geno %>%
  dplyr::select(FID:last_col())
write.table(ped[ped$Sex == 2, ], file = paste0("data/t2d_addiction/ucla/", "cc.f.ped"), #female
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(ped[ped$Sex == 1, ], file = paste0("data/t2d_addiction/ucla/", "cc.m.ped"), #male
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#process phenotype as plink format
pheno <- df.pheno.geno %>%
  dplyr::select(IID, FID, Sex, AQ_SessionsToAcquisition:DR_Inf_Total_1p8mgkg)
#create phenotype file for each sex, NA or -9 to denote missing values.
write.table(pheno[pheno$Sex == 2, -3], file = "data/t2d_addiction/ucla/cc.f.rz.pheno",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
write.table(pheno[pheno$Sex == 1, -3], file = "data/t2d_addiction/ucla/cc.m.rz.pheno",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

# 6.1, Pylmm on CC study based on UCLA grid -----------------------------------------------------------------------
#1 create bed/fam/bim files.
system("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/; /projects/compsci/vmp/USERS/heh/software/plink --noweb --file cc.f --geno 0.1 --mind 0.1 --make-bed --out cc.f") #include only SNPs with a 90% genotyping rate (10% missing) use
system("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/; /projects/compsci/vmp/USERS/heh/software/plink --noweb --file cc.m --geno 0.1 --mind 0.1 --make-bed --out cc.m") #include only SNPs with a 90% genotyping rate (10% missing) use

#2 calculated the kinship matrix file using plink2
system("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/; /projects/csna/csna_workflow/code/plink2 --bfile cc.f  --make-rel square --out cc.f.pylmm.kin")
system("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/; /projects/csna/csna_workflow/code/plink2 --bfile cc.m  --make-rel square --out cc.m.pylmm.kin")

# pylmmgwas job for each phenotype
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/scripts/"
job.para <- pheno %>%
  dplyr::select(-1:-3) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 0:(length(.)-1)) %>%
  dplyr::mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
                stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
# pylmm job
cc.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/cc.pylmm.job.sh"
#Write script
write("#Submit job", cc.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=64GB # memory pool for all cores",                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 1-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            job.para[i, "script.file"], append = T)
  #female
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile cc.f  --kfile cc.f.pylmm.kin.rel --phenofile cc.f.rz.pheno ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".f.pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  #male
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile cc.m  --kfile cc.m.pylmm.kin.rel --phenofile cc.m.rz.pheno ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".m.pylmm.out || true"),
        job.para[i, "script.file"], append = T)

  write(paste0("sbatch ", job.para[i, "script.file"]), cc.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            cc.pylmm.job, append = T)
}

# 6.2, GEMMA on CC study based on UCLA grid-----------------------------------------------------------------------
# job.para
job.para <- pheno %>%
  dplyr::select(-1:-3) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 1:(length(.)))

# gemma job
cc.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/cc.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         cc.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 cc.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      cc.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              cc.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      cc.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         cc.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/cc.gemma.job.stderr"), cc.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/cc.gemma.job.stdout"), cc.gemma.job, append = T)
write("module load singularity",                                             cc.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            cc.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #male
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile cc.m --pheno cc.m.rz.pheno --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out cc.m_p"),  cc.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.m_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/ucla/ -o cc.m_p"), cc.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.m_p",
               " -k ", "cc.m_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/ucla -o ",
               job.para[i,"pheno"], ".m"), cc.gemma.job, append = T)
  #female
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile cc.f --pheno cc.f.rz.pheno --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out cc.f_p"),  cc.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.f_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/ucla/ -o cc.f_p"), cc.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.f_p",
               " -k ", "cc.f_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/ucla -o ",
               job.para[i,"pheno"], ".f"), cc.gemma.job, append = T)
  write("\n", cc.gemma.job, append = T)
}

# 7, Work on CC study based on gigamuga grid -----------------------------------------------------------------------
#cc genotypes gigamuga
cc.geno <- readRDS("/projects/chesler-lab/phenome/metasoft/data/GM_marker_positions_with_MUSter_v2_genotypes.RDS") %>%
  dplyr::select(1:4, all_of("C57BL/6J"), starts_with("CC")) %>%
  dplyr::rename(b_6 = 5) %>%
  dplyr::select(-all_of("CC057/Unc")) %>%
  dplyr::rename_with(., ~gsub("\\/.*", "", .x)) %>%
  dplyr::mutate(across(6:last_col(), ~ case_when(
    . == "N" ~ "0 0",
    . == b_6 ~ "1 1",
    . == "H" ~ "1 2",
    TRUE ~ "2 2"
  ))) %>%
  dplyr::select(-all_of("b_6")) %>%
  dplyr::filter(chr != "M") %>%
  dplyr::distinct(., chr, rs, bp38, .keep_all = TRUE)

#generate map file, four columns, chr, rs, 0, pos in Base-pair coordinate
map <- cc.geno %>%
  dplyr::select(chr, rs, bp38) %>%
  mutate(cM = 0, .after = rs) %>%
  mutate(rs = case_when(
    (rs=="") ~ paste0(chr, "_", bp38),
    TRUE ~ as.character(rs)
  )) %>%
  dplyr::filter(chr != "M") %>%
  dplyr::distinct(., chr, rs, bp38, .keep_all = TRUE)
write.table(map, file = paste0("data/t2d_addiction/gigamuga/", "cc.f.map"), #female
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(map, file = paste0("data/t2d_addiction/gigamuga/", "cc.m.map"), #male
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#cc.geno
cc.geno <- cc.geno %>%
  dplyr::select(5:last_col()) %>%
  t(.) %>%
  as_tibble(.name_repair = "universal") %>%
  mutate(FID = colnames(cc.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         IID = colnames(cc.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         pid = 0,
         mid = 0,
         phe = -9,
         .before = 1)

# cc phenotype
# cc phenotypes ivsa
cc.pheno.ivsa <- readr::read_csv("data/pheno/Prj04_IVSA_preqc_12072020_CC.csv",
                          col_names = TRUE, na = c("", "NA", "N/A", "NaN")) %>%
  rename_with(~str_remove(.x, "\\_cc|\\_do")) %>%
  rename_with(~(.x = "Strain"),   contains("_Strain")) %>%
  rename_with(~(.x = "Mouse.ID"), contains("Subject")) %>%
  rename_with(~(.x = "Sex"),      contains("Sex")) %>%
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male",
    TRUE       ~ as.character(Sex))) %>% #make Sex coding consistent
  mutate(across(Mouse.ID, as.character)) %>%
  dplyr::select(Mouse.ID, Strain, Sex, 25:last_col()) %>%
  dplyr::mutate(Strain = gsub("\\/.*", "", Strain))

# cc phenotypes sens
cc.pheno.sens <- readr::read_csv("data/pheno/Prj02_Sensitization_preqc_12282020_CC.csv",
                          col_names = TRUE, na = c("", "NA", "N/A", "NaN")) %>%
  rename_with(~str_remove(.x, "\\_cc|\\_do")) %>%
  rename_with(~(.x = "Strain"),   contains("_Strain")) %>%
  rename_with(~(.x = "Mouse.ID"), contains("Subject")) %>%
  rename_with(~(.x = "Sex"),      contains("Sex")) %>%
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male",
    TRUE       ~ as.character(Sex))) %>% #make Sex coding consistent
  mutate(across(Mouse.ID, as.character)) %>%
  dplyr::select(Mouse.ID, Strain, Sex, 8:last_col()) %>%
  dplyr::mutate(Strain = gsub("\\/.*", "", Strain))

#full join
cc.pheno <- full_join(cc.pheno.ivsa, cc.pheno.sens) %>%
  dplyr::select(1:3, unique(substr(cc$trait, 1, nchar(cc$trait) - 2))) %>%
  dplyr::mutate(across(AQ_SessionsToAcquisition:last_col(), rz.transform))

#left_join cc.pheno and cc.geno
df.pheno.geno <- left_join(cc.pheno, cc.geno, by = c("Strain" = "IID")) %>%
  mutate(IID = FID, .after = FID) %>%
  relocate(Sex, .after = mid) %>%
  arrange(Mouse.ID) %>%
  mutate(Sex = case_when(
    Sex == "Male" ~ 1,
    Sex == "Female" ~ 2
  )) %>%
  mutate(FID = Mouse.ID) %>% #pylmm cannot handle duplicated ids
  mutate(IID = Mouse.ID)

#ped file
ped <- df.pheno.geno %>%
  dplyr::select(FID:last_col())
write.table(ped[ped$Sex == 2, ], file = paste0("data/t2d_addiction/gigamuga/", "cc.f.ped"), #female
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(ped[ped$Sex == 1, ], file = paste0("data/t2d_addiction/gigamuga/", "cc.m.ped"), #male
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#process phenotype as plink format
pheno <- df.pheno.geno %>%
  dplyr::select(IID, FID, Sex, AQ_SessionsToAcquisition:DR_Inf_Total_1p8mgkg)
#create phenotype file for each sex, NA or -9 to denote missing values.
write.table(pheno[pheno$Sex == 2, -3], file = "data/t2d_addiction/gigamuga/cc.f.rz.pheno",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
write.table(pheno[pheno$Sex == 1, -3], file = "data/t2d_addiction/gigamuga/cc.m.rz.pheno",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

# 8.1, Pylmm on CC study based on gigamuga grid -----------------------------------------------------------------------
#1 create bed/fam/bim files.
system("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/; /projects/compsci/vmp/USERS/heh/software/plink --noweb --file cc.f --geno 0.1 --mind 0.1 --make-bed --out cc.f") #include only SNPs with a 90% genotyping rate (10% missing) use
system("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/; /projects/compsci/vmp/USERS/heh/software/plink --noweb --file cc.m --geno 0.1 --mind 0.1 --make-bed --out cc.m") #include only SNPs with a 90% genotyping rate (10% missing) use

#2 calculated the kinship matrix file using plink2
system("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/; /projects/csna/csna_workflow/code/plink2 --bfile cc.f  --make-rel square --out cc.f.pylmm.kin")
system("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/; /projects/csna/csna_workflow/code/plink2 --bfile cc.m  --make-rel square --out cc.m.pylmm.kin")

# pylmmgwas job for each phenotype
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/scripts/"
job.para <- pheno %>%
  dplyr::select(-1:-3) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 0:(length(.)-1)) %>%
  dplyr::mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
                stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
# pylmm job
cc.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/cc.pylmm.job.sh"
#Write script
write("#Submit job", cc.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=64GB # memory pool for all cores",                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 1-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            job.para[i, "script.file"], append = T)
  #female
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile cc.f  --kfile cc.f.pylmm.kin.rel --phenofile cc.f.rz.pheno ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".f.pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  #male
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile cc.m  --kfile cc.m.pylmm.kin.rel --phenofile cc.m.rz.pheno ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".m.pylmm.out || true"),
        job.para[i, "script.file"], append = T)

  write(paste0("sbatch ", job.para[i, "script.file"]), cc.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            cc.pylmm.job, append = T)
}

# 8.2, GEMMA on CC study based on gigamuga grid-----------------------------------------------------------------------
# job.para
job.para <- pheno %>%
  dplyr::select(-1:-3) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 1:(length(.)))

# gemma job
cc.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/cc.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         cc.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 cc.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      cc.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              cc.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      cc.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         cc.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/cc.gemma.job.stderr"), cc.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/cc.gemma.job.stdout"), cc.gemma.job, append = T)
write("module load singularity",                                             cc.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            cc.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #male
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile cc.m --pheno cc.m.rz.pheno --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out cc.m_p"),  cc.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.m_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/ -o cc.m_p"), cc.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.m_p",
               " -k ", "cc.m_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/gigamuga -o ",
               job.para[i,"pheno"], ".m"), cc.gemma.job, append = T)
  #female
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile cc.f --pheno cc.f.rz.pheno --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out cc.f_p"),  cc.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.f_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/ -o cc.f_p"), cc.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile cc.f_p",
               " -k ", "cc.f_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/gigamuga -o ",
               job.para[i,"pheno"], ".f"), cc.gemma.job, append = T)
  write("\n", cc.gemma.job, append = T)
}
