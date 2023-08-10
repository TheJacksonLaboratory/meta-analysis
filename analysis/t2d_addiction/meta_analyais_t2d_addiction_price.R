## ---------------------------
##
## Script name:
##
## Purpose of script: Running pyLMM on UCLA and GigaMuga density imputed grid for price bxd mice
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

# 5, Work on price study based on UCLA grid -----------------------------------------------------------------------
price <- study_meta %>%
  dplyr::filter(cohort == "price.bxd") %>%
  dplyr::mutate(trait = gsub("cocaine_price_rz_", "", study_name))

#pheno
price.pheno <- read.csv("data/BXD_traits.csv", header = TRUE)
colnames(price.pheno) = price.pheno[6,]
price.pheno <- price.pheno[27:nrow(price.pheno), ] %>%
  dplyr::rename(., Strain = Symbol) %>%
  dplyr::filter(!str_detect(Strain, "_SE")) %>%
  dplyr::filter(!str_detect(Strain, "x"))
price.pheno[price.pheno == "x"] <- NA

#geno
price.geno <- readRDS("/projects/chesler-lab/phenome/metasoft/data/UCLA1_positions_with_MUSter_v2_genotypes.RDS") %>%
  dplyr::select(1:4, all_of("C57BL/6J"), starts_with("BXD")) %>%
  dplyr::rename(b_6 = 5)
name.idx = gsub("\\/.*", "", colnames(price.geno))

#subset
strains = intersect(unique(price.pheno$Strain), name.idx)
price.pheno <- price.pheno %>%
  dplyr::filter(Strain %in% strains)

index <- NULL
for(s in strains) {
  idx = which(name.idx == s)
  if(length(idx) == 1){
    index[s] = idx
  }else{
    index[s] <- ifelse(sum(price.geno[, idx[[1]]] == "") > sum(price.geno[, idx[[2]]] == ""),
                       idx[[2]], idx[[1]])
  }
}
price.geno <- price.geno[, c(1:5, index)]
colnames(price.geno) = gsub("\\/.*", "", colnames(price.geno))

#geno
price.geno <- price.geno %>%
  dplyr::mutate(across(6:last_col(), ~ case_when(
    . == "N" ~ "0 0",
    . == b_6 ~ "1 1",
    . == "H" ~ "1 2",
    TRUE ~ "2 2"
  ))) %>%
  dplyr::select(-all_of("b_6"))

#generate map file, four columns, chr, rs, 0, pos in Base-pair coordinate
map <- price.geno %>%
  dplyr::select(chr, rs, bp38) %>%
  mutate(cM = 0, .after = rs) %>%
  mutate(rs = case_when(
    (rs=="") ~ paste0(chr, "_", bp38),
    TRUE ~ as.character(rs)
  ))
write.table(map, file = paste0("data/t2d_addiction/ucla/", "price.map"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#price.geno
price.geno <- price.geno %>%
  dplyr::select(5:last_col()) %>%
  t(.) %>%
  as_tibble(.name_repair = "universal") %>%
  mutate(FID = colnames(price.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         IID = colnames(price.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         pid = 0,
         mid = 0,
         phe = -9,
         .before = 1)
#left_join price.pheno and price.geno
df.pheno.geno <- left_join(price.pheno, price.geno, by = c("Strain" = "IID")) %>%
  mutate(IID = FID, .after = FID) %>%
  mutate(Sex = -9, .after = mid) %>%
  mutate(FID = 1:nrow(price.pheno)) %>% #pylmm cannot handle duplicated ids
  mutate(IID = 1:nrow(price.pheno))

#ped file
ped <- df.pheno.geno %>%
  dplyr::select(FID:last_col())
write.table(ped, file = paste0("data/t2d_addiction/ucla/", "price.ped"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#process phenotype as plink format
pheno <- df.pheno.geno %>%
  dplyr::select(IID, FID, COC_IVSA_SES_TO_ACQ:COC_IVSA_INACTV_1.8) %>%
  dplyr::mutate(across(COC_IVSA_SES_TO_ACQ:COC_IVSA_INACTV_1.8, rz.transform))
#create phenotype file for NA or -9 to denote missing values.
write.table(pheno, file = "data/t2d_addiction/ucla/price.rz.pheno",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

# 6.1, Pylmm on price study based on UCLA grid -----------------------------------------------------------------------
# pylmmgwas job for price
#1 create bed/fam/bim files.
system("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/; /projects/compsci/vmp/USERS/heh/software/plink --noweb --file price --geno 0.1 --make-bed --out price")
#2 calculated the kinship matrix file using plink2
system("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/; /projects/csna/csna_workflow/code/plink2 --bfile price  --make-rel square --out price.pylmm.kin")

# pylmmgwas job for each phenotype
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/scripts/"
job.para <- pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 0:(length(.)-1)) %>%
  dplyr::mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
                stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
# pylmm job
price.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/price.pylmm.job.sh"
#Write script
write("#Submit job", price.pylmm.job, append = F)
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
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile price  --kfile price.pylmm.kin.rel --phenofile price.rz.pheno ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".price.pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  write(paste0("sbatch ", job.para[i, "script.file"]), price.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            price.pylmm.job, append = T)
}

# 6.2, GEMMA on price study based on UCLA grid-----------------------------------------------------------------------
# job.para
job.para <- pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 1:(length(.)))

# gemma job
price.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/price.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         price.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 price.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      price.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              price.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      price.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         price.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/price.gemma.job.stderr"), price.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/price.gemma.job.stdout"), price.gemma.job, append = T)
write("module load singularity",                                             price.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            price.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile price --pheno price.rz.pheno --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out price_p"),  price.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile price_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/ucla/ -o price_p"), price.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile price_p",
               " -k ", "price_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/ucla -o ",
               job.para[i,"pheno"]), price.gemma.job, append = T)
  write("\n", price.gemma.job, append = T)
}

# 7, Work on price study based on gigamuga grid -----------------------------------------------------------------------
#pheno
price.pheno <- read.csv("data/BXD_traits.csv", header = TRUE)
colnames(price.pheno) = price.pheno[6,]
price.pheno <- price.pheno[27:nrow(price.pheno), ] %>%
  dplyr::rename(., Strain = Symbol) %>%
  dplyr::filter(!str_detect(Strain, "_SE")) %>%
  dplyr::filter(!str_detect(Strain, "x"))
price.pheno[price.pheno == "x"] <- NA

#geno
price.geno <- readRDS("/projects/chesler-lab/phenome/metasoft/data/GM_marker_positions_with_MUSter_v2_genotypes.RDS") %>%
  dplyr::select(1:4, all_of("C57BL/6J"), starts_with("BXD")) %>%
  dplyr::rename(b_6 = 5) %>%
  dplyr::filter(chr != "M") %>%
  dplyr::distinct(., chr, rs, bp38, .keep_all = TRUE)
name.idx = gsub("\\/.*", "", colnames(price.geno))

#subset
strains = intersect(unique(price.pheno$Strain), name.idx)
price.pheno <- price.pheno %>%
  dplyr::filter(Strain %in% strains)

index <- NULL
for(s in strains) {
  idx = which(name.idx == s)
  if(length(idx) == 1){
    index[s] = idx
  }else{
    index[s] <- ifelse(sum(price.geno[, idx[[1]]] == "") > sum(price.geno[, idx[[2]]] == ""),
                       idx[[2]], idx[[1]])
  }
}
price.geno <- price.geno[, c(1:5, index)]
colnames(price.geno) = gsub("\\/.*", "", colnames(price.geno))

#geno
price.geno <- price.geno %>%
  dplyr::mutate(across(6:last_col(), ~ case_when(
    . == "N" ~ "0 0",
    . == b_6 ~ "1 1",
    . == "H" ~ "1 2",
    TRUE ~ "2 2"
  ))) %>%
  dplyr::select(-all_of("b_6"))

#generate map file, four columns, chr, rs, 0, pos in Base-pair coordinate
map <- price.geno %>%
  dplyr::select(chr, rs, bp38) %>%
  mutate(cM = 0, .after = rs) %>%
  mutate(rs = case_when(
    (rs=="") ~ paste0(chr, "_", bp38),
    TRUE ~ as.character(rs)
  )) %>%
  dplyr::filter(chr != "M") %>%
  dplyr::distinct(., chr, rs, bp38, .keep_all = TRUE)
write.table(map, file = paste0("data/t2d_addiction/gigamuga/", "price.map"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#price.geno
price.geno <- price.geno %>%
  dplyr::select(5:last_col()) %>%
  t(.) %>%
  as_tibble(.name_repair = "universal") %>%
  mutate(FID = colnames(price.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         IID = colnames(price.geno)[-1:-4], #str_replace_all(colnames(cc.geno)[-1:-6], "/GeniUnc|/Unc|/TauUnc", ""),
         pid = 0,
         mid = 0,
         phe = -9,
         .before = 1)
#left_join price.pheno and price.geno
df.pheno.geno <- left_join(price.pheno, price.geno, by = c("Strain" = "IID")) %>%
  mutate(IID = FID, .after = FID) %>%
  mutate(Sex = -9, .after = mid) %>%
  mutate(FID = 1:nrow(price.pheno)) %>% #pylmm cannot handle duplicated ids
  mutate(IID = 1:nrow(price.pheno))

#ped file
ped <- df.pheno.geno %>%
  dplyr::select(FID:last_col())
write.table(ped, file = paste0("data/t2d_addiction/gigamuga/", "price.ped"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#process phenotype as plink format
pheno <- df.pheno.geno %>%
  dplyr::select(IID, FID, COC_IVSA_SES_TO_ACQ:COC_IVSA_INACTV_1.8) %>%
  dplyr::mutate(across(COC_IVSA_SES_TO_ACQ:COC_IVSA_INACTV_1.8, rz.transform))
#create phenotype file for NA or -9 to denote missing values.
write.table(pheno, file = "data/t2d_addiction/gigamuga/price.rz.pheno",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

# 8.1, Pylmm on price study based on gigamuga grid -----------------------------------------------------------------------
# pylmmgwas job for price
#1 create bed/fam/bim files.
system("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/; /projects/compsci/vmp/USERS/heh/software/plink --noweb --file price --geno 0.1 --mind 0.1 --make-bed --out price")
#2 calculated the kinship matrix file using plink2
system("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/; /projects/csna/csna_workflow/code/plink2 --bfile price  --make-rel square --out price.pylmm.kin")

# pylmmgwas job for each phenotype
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/scripts/"
job.para <- pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 0:(length(.)-1)) %>%
  dplyr::mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
                stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
# pylmm job
price.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/price.pylmm.job.sh"
#Write script
write("#Submit job", price.pylmm.job, append = F)
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
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile price  --kfile price.pylmm.kin.rel --phenofile price.rz.pheno ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".price.pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  write(paste0("sbatch ", job.para[i, "script.file"]), price.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            price.pylmm.job, append = T)
}

# 8.2, GEMMA on price study based on gigamuga grid-----------------------------------------------------------------------
# job.para
job.para <- pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 1:(length(.)))

# gemma job
price.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/price.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         price.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 price.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      price.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              price.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      price.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         price.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/price.gemma.job.stderr"), price.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/price.gemma.job.stdout"), price.gemma.job, append = T)
write("module load singularity",                                             price.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            price.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile price --pheno price.rz.pheno --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out price_p"),  price.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile price_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/ -o price_p"), price.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile price_p",
               " -k ", "price_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/gigamuga -o ",
               job.para[i,"pheno"]), price.gemma.job, append = T)
  write("\n", price.gemma.job, append = T)
}
