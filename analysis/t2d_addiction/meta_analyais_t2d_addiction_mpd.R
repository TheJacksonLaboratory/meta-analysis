## ---------------------------
##
## Script name:
##
## Purpose of script: Running pyLMM on UCLA and GigaMuga density imputed grid for MPD study
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

# 4, Format genotype data into plink -----------------------------------------------------------------------
#ucla
generate_ped_map(genotype.file = "/projects/chesler-lab/phenome/metasoft/data/UCLA1_positions_with_MUSter_v2_genotypes.RDS",
                 ped.file = "data/t2d_addiction/ucla/ucla.ped",
                 map.file = "data/t2d_addiction/ucla/ucla.map")
#gigamuga
generate_ped_map(genotype.file = "/projects/chesler-lab/phenome/metasoft/data/GM_marker_positions_with_MUSter_v2_genotypes.RDS",
                 ped.file = "data/t2d_addiction/gigamuga/gigamuga.ped",
                 map.file = "data/t2d_addiction/gigamuga/gigamuga.map")

# 5, Read and format the study name -----------------------------------------------------------------------
HFD <- readxl::read_excel(path = "data/measurements_w_HFD.xlsx",
                          sheet = 2) %>%
  dplyr::filter(Analyze == 1) %>%
  dplyr::filter(paneldesc != "DO population") %>% #no DO population
  dplyr::select(1, 11)
HFD_both <- HFD %>%
  dplyr::filter(sextested == "both")
HFD_m_f <- HFD %>%
  dplyr::filter(sextested != "both") %>%
  dplyr::mutate(study_name = paste0(measnum,"_", sextested))

hf_diet = readRDS("/projects/chesler-lab/phenome/metasoft/manuscript/data/forHao/high-fat_diet_studynames.RDS")
addiction = readRDS("/projects/chesler-lab/phenome/metasoft/manuscript/data/forHao/addiction_studynames.RDS")

study_meta = rbind(data.frame(study_name = c(hf_diet, HFD_m_f$study_name,
                                             paste0(HFD_both$measnum, "_m"),
                                             paste0(HFD_both$measnum, "_f")
                                             ), disease = "t2d"),
                   data.frame(study_name = addiction, disease = "addiction"))

study_meta = study_meta %>%
  dplyr::mutate(cohort = case_when(
    str_detect(study_name, "csna03_cc") ~ "cc",
    str_detect(study_name, "csna03_do") ~ "do",
    str_detect(study_name, "price_rz")  ~ "price",
    TRUE ~ "mpd"
  ))

# 6, Work on MPD study -----------------------------------------------------------------------
mpd <- study_meta %>%
  dplyr::filter(cohort == "mpd") %>%
  dplyr::mutate(study_id = str_extract(study_name, "^[^_]*")) %>%
  # dplyr::filter(study_name %in% c(HFD_m_f$study_name,
  #                                paste0(HFD_both$measnum, "_m"),
  #                                paste0(HFD_both$measnum, "_f"))) %>%
  #dplyr::filter(!(study_name %in% hf_diet)) %>% #new added
  dplyr::select(-2) %>%
  dplyr::distinct()

#extract mpd measures
mpd.pheno <- extract_pheno_measure(url = "https://phenome.jax.org/", measure_id = unique(mpd$study_id))
mpd.pheno <- mpd.pheno %>%
  dplyr::mutate(sex = str_replace(sex, "FALSE", "f"))
#generate_pheno_plink for ucla
generate_pheno_plink(ped.file = "data/t2d_addiction/ucla/ucla.ped",
                     map.file = "data/t2d_addiction/ucla/ucla.map",
                     pheno   = mpd.pheno,
                     outdir = "data/t2d_addiction/ucla/",
                     ncore = 20)
#generate_pheno_plink for gigamuga
generate_pheno_plink(ped.file = "data/t2d_addiction/gigamuga/gigamuga.ped",
                     map.file = "data/t2d_addiction/gigamuga/gigamuga.map",
                     pheno   = mpd.pheno,
                     outdir = "data/t2d_addiction/gigamuga/",
                     ncore = 20)

# 7.1, Pylmm on MPD study based on UCLA grid-----------------------------------------------------------------------
# job.para
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/scripts/"
job.para <- mpd %>%
  dplyr::mutate(file = gsub("_", ".", study_name)) %>%
  dplyr::mutate(kin =  paste0(file, ".kin")) %>%
  dplyr::mutate(script.file = paste0(script.Dir, file,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, file,  ".pylmm.err"),
                stdout      = paste0(script.Dir, file,  ".pylmm.out")) %>%
  as.data.frame()

# calculated the kinship matrix file using plink2
mpd.pylmm.kinship.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.pylmm.kinship.job.sh"
#Write script
write("#!/bin/bash",                                                         mpd.pylmm.kinship.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 mpd.pylmm.kinship.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              mpd.pylmm.kinship.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.pylmm.kinship.job.stderr"), mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.pylmm.kinship.job.stdout"), mpd.pylmm.kinship.job, append = T)
write("module load singularity",                                             mpd.pylmm.kinship.job, append = T)
#Write into mpd.pylmm.kinship.job.sh
write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/", mpd.pylmm.kinship.job, append = T)
#bed/bim/fam
write(paste0("/projects/csna/meta-analysis/code/plink2 --pedmap ", job.para[, "file"],
             " --make-bed --geno 0.1 --mind 0.1 --out ", job.para[, "file"]), mpd.pylmm.kinship.job, append = T)
write(paste0("/projects/csna/meta-analysis/code/plink2 --pedmap ", job.para[, "file"],
             " --make-rel square --out ", job.para[, "kin"]), mpd.pylmm.kinship.job, append = T)

# pylmm job
mpd.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.pylmm.job.sh"
#Write script
write("#Submit job", mpd.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=64GB # memory pool for all cores",                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            job.para[i, "script.file"], append = T)
  #pylmm
  #zscore
  write(paste0("python /projects/csna/meta-analysis/code/pylmm/scripts/pylmmGWAS.py -v --bfile ",
               job.para[i,"file"], " --kfile ", job.para[i,"file"], ".kin.rel",
               " --phenofile ", job.para[i,"file"], ".pheno",
               " -p ", 0, " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"file"], ".pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  write(paste0("sbatch ", job.para[i, "script.file"]), mpd.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            mpd.pylmm.job, append = T)
}

# 7.2, GEMMA on MPD study based on UCLA grid-----------------------------------------------------------------------
# job.para
job.para <- mpd %>%
  dplyr::mutate(file = gsub("_", ".", study_name)) %>%
  as.data.frame()

# gemma job
mpd.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         mpd.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 mpd.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      mpd.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              mpd.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      mpd.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         mpd.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.gemma.job.stderr"), mpd.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/mpd.gemma.job.stdout"), mpd.gemma.job, append = T)
write("module load singularity",                                             mpd.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            mpd.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Generate bed file
  write(paste0("/projects/csna/meta-analysis/code/plink2 --pedmap ", job.para[i, "file"],
               " --make-bed --geno 0.1 --mind 0.1 --out ", job.para[i, "file"]),  mpd.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile ",
               job.para[i,"file"],
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/ucla/ -o ",
               job.para[i,"file"]), mpd.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile ",
               job.para[i,"file"], " -k ",
               job.para[i,"file"], ".cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/ucla -o ",
               job.para[i,"file"]), mpd.gemma.job, append = T)
  write("\n", mpd.gemma.job, append = T)
}

# 8.1, Pylmm on MPD study based on gigamuga grid-----------------------------------------------------------------------
# job.para
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/scripts/"
job.para <- mpd %>%
  dplyr::mutate(file = gsub("_", ".", study_name)) %>%
  dplyr::mutate(kin =  paste0(file, ".kin")) %>%
  dplyr::mutate(script.file = paste0(script.Dir, file,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, file,  ".pylmm.err"),
                stdout      = paste0(script.Dir, file,  ".pylmm.out")) %>%
  as.data.frame()

# calculated the kinship matrix file using plink2
mpd.pylmm.kinship.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.pylmm.kinship.job.sh"
#Write script
write("#!/bin/bash",                                                         mpd.pylmm.kinship.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 mpd.pylmm.kinship.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              mpd.pylmm.kinship.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.pylmm.kinship.job.stderr"), mpd.pylmm.kinship.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.pylmm.kinship.job.stdout"), mpd.pylmm.kinship.job, append = T)
write("module load singularity",                                             mpd.pylmm.kinship.job, append = T)
#Write into mpd.pylmm.kinship.job.sh
write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/", mpd.pylmm.kinship.job, append = T)
#bed/bim/fam
write(paste0("/projects/csna/meta-analysis/code/plink2 --pedmap ", job.para[, "file"],
             " --make-bed --geno 0.1 --mind 0.1 --out ", job.para[, "file"]), mpd.pylmm.kinship.job, append = T)
write(paste0("/projects/csna/meta-analysis/code/plink2 --pedmap ", job.para[, "file"],
             " --make-rel square --out ", job.para[, "kin"]), mpd.pylmm.kinship.job, append = T)

# pylmm job
mpd.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.pylmm.job.sh"
#Write script
write("#Submit job", mpd.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=64GB # memory pool for all cores",                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            job.para[i, "script.file"], append = T)
  #pylmm
  #zscore
  write(paste0("python /projects/csna/meta-analysis/code/pylmm/scripts/pylmmGWAS.py -v --bfile ",
               job.para[i,"file"], " --kfile ", job.para[i,"file"], ".kin.rel",
               " --phenofile ", job.para[i,"file"], ".pheno",
               " -p ", 0, " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"file"], ".pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  write(paste0("sbatch ", job.para[i, "script.file"]), mpd.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            mpd.pylmm.job, append = T)
}

# 8.2, GEMMA on MPD study based on gigamuga grid-----------------------------------------------------------------------
# job.para
job.para <- mpd %>%
  dplyr::mutate(file = gsub("_", ".", study_name)) %>%
  as.data.frame()

# gemma job
mpd.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         mpd.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 mpd.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      mpd.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              mpd.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      mpd.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         mpd.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.gemma.job.stderr"), mpd.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/mpd.gemma.job.stdout"), mpd.gemma.job, append = T)
write("module load singularity",                                             mpd.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            mpd.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Generate bed file
  write(paste0("/projects/csna/meta-analysis/code/plink2 --pedmap ", job.para[i, "file"],
               " --make-bed --geno 0.1 --mind 0.1 --out ", job.para[i, "file"]),  mpd.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile ",
               job.para[i,"file"],
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/ -o ",
               job.para[i,"file"]), mpd.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile ",
               job.para[i,"file"], " -k ",
               job.para[i,"file"], ".cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/gigamuga -o ",
               job.para[i,"file"]), mpd.gemma.job, append = T)
  write("\n", mpd.gemma.job, append = T)
}

