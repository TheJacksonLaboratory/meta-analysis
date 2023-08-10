## ---------------------------
##
## Script name:
##
## Purpose of script: Running pyLMM on UCLA and GigaMuga density imputed grid for DO mice
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

# 4, Interpolate DO GigaMUGA array genoprobs to UCLA grid (132277 grid)  -----------------------------------------------------------------------
if(FALSE){ #this step has been done by csna_workflow/analysis/scripts/11_pr_apr_132k.R and 11_pr_apr_132k.sh
  source("/projects/csna/csna_workflow/code/reconst_utils.R")
  load("/projects/csna/csna_workflow/data/Jackson_Lab_12_batches/gm_DO3173_qc_newid.RData")#gm_after_qc
  load("/projects/csna/csna_workflow/data/Jackson_Lab_12_batches/pr_DO3173.RData")#pr
  #x chr
  if(dim(pr$X)[2] == 44){
    #X chromosome should have 36 states
    pr$X <- pr$X[,(-37:-44),]
  }

  str(pr)
  print(paste0("Number of samples, ", dim(pr[[1]])[1]))
  #interpolate in chunks
  chunk_size = 167 #19*167 = 3173
  idx <- list()
  pr.132k <- list()
  for(chunk_number in 1:19){
    print(chunk_number)
    idx[[chunk_number]] <-  ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
    #Create a large 3D array with samples/states/SNPs in dimensions 1,2,3.
    #combine all chrs into one 3d array
    pr.3d <- do.call("abind",list(pr[idx[[chunk_number]],],along = 3))

    #gigamuga snps
    load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
    #define from
    from = GM_snps[intersect(dimnames(pr.3d)[[3]],GM_snps$marker),1:4]
    #define to
    grid <- read.csv("/projects/csna/csna_workflow/data/meta/snp_UCLA1_full.csv", header = TRUE, fill = TRUE, na.strings = "") %>%
      select(1:3) %>%
      mutate(rs = case_when(
        is.na(rs) ~ paste0(chr, "_", bp38),
        TRUE ~ as.character(rs)
      )) %>%
      mutate(bp38 = bp38/10^6)  %>%
      select(marker = rs,
             chr = chr,
             pos = bp38)
    rownames(grid) <- grid$marker
    to = grid[,1:3]
    #subset pr.3d to the marker
    sub.pr.3d <- pr.3d[,,from$marker]
    #interplote to 132k
    #interplote function
    interplote_fun <- function(i){
      z <- t(interpolate.markers(data = as.matrix(t(sub.pr.3d[i,,])), from = from, to = to))
      return(z)
    }
    #each chunk
    system.time({
      pr.132k[[chunk_number]] <- mclapply(X = 1:dim(sub.pr.3d)[[1]], FUN = interplote_fun, mc.preschedule = F, mc.cores = 10)
    })
    names(pr.132k[[chunk_number]]) <- dimnames(sub.pr.3d)[[1]]
    pr.132k[[chunk_number]] <- do.call("abind",list(pr.132k[[chunk_number]],along = 3))
    pr.132k[[chunk_number]] <- aperm(pr.132k[[chunk_number]], perm = c(3,1,2))
    print(str(pr.132k[[chunk_number]]))
  }
  pr.132k <- do.call("abind",list(pr.132k,along = 1))
  save(pr.132k, file = "/projects/csna/csna_workflow/data/Jackson_Lab_12_batches/pr_132k_DO3173.RData")
  print("pr_132k_DO3173 done")
}

# 5, Impute genotypes based on DO genoprobs of 132k and generate ped/map file from sanger SNP -----------------------------------------------------------------------
if(FALSE){#this step had been done by csna_workflow/analysis/scripts/13_impute_sanger_snp.R and 13_impute_sanger_snp.sh
  source("/projects/csna/csna_workflow/code/reconst_utils.R")
  source("/projects/csna/csna_workflow/code/predict_sanger_snp_geno.R")
  source("/projects/csna/csna_workflow/code/do2sanger.helper.R")

  # DO
  #load pr.132k
  load("/projects/csna/csna_workflow/data/Jackson_Lab_12_batches/pr_132k_DO3173.RData") #pr.132k

  #snps
  snp_UCLA1 <- read.csv("/projects/csna/csna_workflow/data/meta/snp_UCLA1_full.csv", header = TRUE, fill = TRUE, na.strings = "")
  snps <- snp_UCLA1 %>%
    dplyr::select(1:3) %>%
    mutate(rs = case_when(
      is.na(rs) ~ paste0(chr, "_", bp38),
      TRUE ~ as.character(rs)
    )) %>%
    mutate(bp38 = bp38/10^6)  %>%
    dplyr::select(marker = rs,
                  chr = chr,
                  pos = bp38)
  rownames(snps) <- snps$marker

  #predict_sanger_snp_geno
  sanger_snp_geno <- predict_sanger_snp_geno(pr = pr.132k,
                                             snps =  snps,
                                             out.dir = "/projects/csna/csna_workflow/output/meta/",
                                             output.file = "snp_geno_132k.txt",
                                             snp.file = "/projects/csna/csna_workflow/data/meta/cc.snps.NCBI38.txt.gz",
                                             cores = 20)
}

# 6, Read and format the study name -----------------------------------------------------------------------
hf_diet = readRDS("/projects/chesler-lab/phenome/metasoft/manuscript/data/forHao/high-fat_diet_studynames.RDS")
addiction = readRDS("/projects/chesler-lab/phenome/metasoft/manuscript/data/forHao/addiction_studynames.RDS")
study_meta = rbind(data.frame(study_name = hf_diet, disease = "t2d"),
                   data.frame(study_name = addiction, disease = "addiction"))

study_meta = study_meta %>%
  dplyr::mutate(cohort = case_when(
    str_detect(study_name, "csna03_cc") ~ "cc",
    str_detect(study_name, "csna03_do") ~ "do",
    str_detect(study_name, "price_rz")  ~ "price",
    TRUE ~ "mpd"
  ))

# 7, Work on DO study based on UCLA grid -----------------------------------------------------------------------
do <- study_meta %>%
  dplyr::filter(cohort == "do") %>%
  dplyr::mutate(trait = gsub("cocaine_csna03_do_|cocaine_csna03_do_ex_ri_|cocaine_csna03_do_sens_", "", study_name)) %>%
  dplyr::mutate(trait = gsub(".{2}$", "", trait, perl = TRUE))

# for ped and map file, we need to
#1, filter by sex, create bed/fam/bim files.
#2, calculated the kinship matrix file using pylmmKinship.py
#3, create phenotype file for each sex, NA or -9 to denote missing values.

#1 filter by sex, create bed/fam/bim files.
#cd /projects/csna/csna_workflow/output/meta
#/projects/compsci/vmp/USERS/heh/software/plink --noweb --file DO3173_snp_geno_ucla --filter-males --make-bed --out /projects/csna/meta-analysis/data/t2d_addiction/ucla/DO3173_snp_geno_ucla_male
#/projects/compsci/vmp/USERS/heh/software/plink --noweb --file DO3173_snp_geno_ucla --filter-females --make-bed --out /projects/csna/meta-analysis/data/t2d_addiction/ucla/DO3173_snp_geno_ucla_female

#2 calculated the kinship matrix file using pylmmKinship.py
# cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/
# python2 /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO3173_snp_geno_ucla_male DO3173_snp_geno_ucla_male.pylmm.kin
# python2 /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO3173_snp_geno_ucla_female DO3173_snp_geno_ucla_female.pylmm.kin

#3, create phenotype file for each sex, NA or -9 to denote missing values.
csna.allqtl.pheno <- read.csv("/projects/csna/csna_workflow/data/pheno/csna.allqtl.pheno.csv", header = T) %>%
  dplyr::mutate(MID = Mouse.ID,
                FID = Mouse.ID,
                .after = predict.sex) %>%
  dplyr::select(-1:-5) %>%
  dplyr::select(MID, FID, all_of(unique(do$trait)))

#male
df.male <- read.table("/projects/csna/meta-analysis/data/t2d_addiction/ucla/DO3173_snp_geno_ucla_male.fam") %>%
  dplyr::rename(MID = V1)
#male pheno
csna.allqtl.pheno.male <- semi_join(csna.allqtl.pheno, df.male) %>% #1512 male
  dplyr::mutate(across(ends_with(c(".raw", ".nooutlier")), rankZ))
write.table(csna.allqtl.pheno.male, file = "/projects/csna/meta-analysis/data/t2d_addiction/ucla/DO3173_snp_geno_ucla_male.rz.phenos", row.names = FALSE, col.names = FALSE, sep = " ")
#female pheno
csna.allqtl.pheno.female <- anti_join(csna.allqtl.pheno, df.male) %>% #1661 female
  dplyr::mutate(across(ends_with(c(".raw", ".nooutlier")), rankZ))
write.table(csna.allqtl.pheno.female, file = "/projects/csna/meta-analysis/data/t2d_addiction/ucla/DO3173_snp_geno_ucla_female.rz.phenos", row.names = FALSE, col.names = FALSE, sep = " ")

# 8.1, Pylmm on DO study based on UCLA grid -----------------------------------------------------------------------
# pylmmgwas job for each phenotype
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/scripts/"
job.para <- csna.allqtl.pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 0:(length(.)-1)) %>%
  dplyr::mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
                stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
# pylmm job
do.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/do.pylmm.job.sh"
#Write script
write("#Submit job", do.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=64GB # memory pool for all cores",                     job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 1-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            job.para[i, "script.file"], append = T)
  #male
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO3173_snp_geno_ucla_male  --kfile DO3173_snp_geno_ucla_male.pylmm.kin --phenofile DO3173_snp_geno_ucla_male.rz.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".m.pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  #female
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO3173_snp_geno_ucla_female  --kfile DO3173_snp_geno_ucla_female.pylmm.kin --phenofile DO3173_snp_geno_ucla_female.rz.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".f.pylmm.out || true"),
        job.para[i, "script.file"], append = T)

  write(paste0("sbatch ", job.para[i, "script.file"]), do.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            do.pylmm.job, append = T)
}

# 8.2, GEMMA on DO study based on UCLA grid-----------------------------------------------------------------------
# job.para
job.para <- csna.allqtl.pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 1:(length(.)))

# gemma job
do.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/do.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         do.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 do.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      do.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              do.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      do.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         do.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/do.gemma.job.stderr"), do.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/do.gemma.job.stdout"), do.gemma.job, append = T)
write("module load singularity",                                             do.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",            do.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #male
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile DO3173_snp_geno_ucla_male --pheno DO3173_snp_geno_ucla_male.rz.phenos --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out DO3173_snp_geno_ucla_male_p"),  do.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_ucla_male_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/ucla/ -o DO3173_snp_geno_ucla_male_p"), do.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_ucla_male_p",
               " -k ", "DO3173_snp_geno_ucla_male_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/ucla -o ",
               job.para[i,"pheno"], ".m"), do.gemma.job, append = T)
  #female
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile DO3173_snp_geno_ucla_female --pheno DO3173_snp_geno_ucla_female.rz.phenos --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out DO3173_snp_geno_ucla_female_p"),  do.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_ucla_female_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/ucla/ -o DO3173_snp_geno_ucla_female_p"), do.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_ucla_female_p",
               " -k ", "DO3173_snp_geno_ucla_female_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/ucla -o ",
               job.para[i,"pheno"], ".f"), do.gemma.job, append = T)
  write("\n", do.gemma.job, append = T)
}

# 9, Work on DO study based on gigamuga grid -----------------------------------------------------------------------

# for ped and map file, we need to
#1, filter by sex, create bed/fam/bim files.
#2, calculated the kinship matrix file using pylmmKinship.py
#3, create phenotype file for each sex, NA or -9 to denote missing values.

#1 filter by sex, create bed/fam/bim files.
#cd /projects/csna/csna_workflow/data/GCTA
#/projects/csna/meta-analysis/code/Plink2/plink2 --bfile 12_batches_QC --rm-dup force-first --filter-males --make-bed --out /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO3173_snp_geno_gigamuga_male
#/projects/csna/meta-analysis/code/Plink2/plink2 --bfile 12_batches_QC --rm-dup force-first --filter-females --make-bed --out /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO3173_snp_geno_gigamuga_female

#2 calculated the kinship matrix file using pylmmKinship.py
# cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/
# python2 /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO3173_snp_geno_gigamuga_male DO3173_snp_geno_gigamuga_male.pylmm.kin
# python2 /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO3173_snp_geno_gigamuga_female DO3173_snp_geno_gigamuga_female.pylmm.kin

#3, create phenotype file for each sex, NA or -9 to denote missing values.
csna.allqtl.pheno <- read.csv("/projects/csna/csna_workflow/data/pheno/csna.allqtl.pheno.csv", header = T) %>%
  dplyr::mutate(MID = Mouse.ID,
                FID = Mouse.ID,
                .after = predict.sex) %>%
  dplyr::select(-1:-5) %>%
  dplyr::select(MID, FID, all_of(unique(do$trait)))

#male
df.male <- read.table("/projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO3173_snp_geno_gigamuga_male.fam") %>%
  dplyr::rename(MID = V1)
#male pheno
csna.allqtl.pheno.male <- semi_join(csna.allqtl.pheno, df.male) %>% #1512 male
  dplyr::mutate(across(ends_with(c(".raw", ".nooutlier")), rankZ))
write.table(csna.allqtl.pheno.male, file = "/projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO3173_snp_geno_gigamuga_male.rz.phenos", row.names = FALSE, col.names = FALSE, sep = " ")
#female pheno
csna.allqtl.pheno.female <- anti_join(csna.allqtl.pheno, df.male) %>% #1661 female
  dplyr::mutate(across(ends_with(c(".raw", ".nooutlier")), rankZ))
write.table(csna.allqtl.pheno.female, file = "/projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO3173_snp_geno_gigamuga_female.rz.phenos", row.names = FALSE, col.names = FALSE, sep = " ")

# 10.1, Pylmm on DO study based on gigamuga grid -----------------------------------------------------------------------
# pylmmgwas job for each phenotype
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/scripts/"
job.para <- csna.allqtl.pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 0:(length(.)-1)) %>%
  dplyr::mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
                stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
                stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
# pylmm job
do.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/do.pylmm.job.sh"
#Write script
write("#Submit job", do.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=64GB # memory pool for all cores",                     job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 1-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            job.para[i, "script.file"], append = T)
  #male
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO3173_snp_geno_gigamuga_male  --kfile DO3173_snp_geno_gigamuga_male.pylmm.kin --phenofile DO3173_snp_geno_gigamuga_male.rz.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".m.pylmm.out || true"),
        job.para[i, "script.file"], append = T)
  #female
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO3173_snp_geno_gigamuga_female  --kfile DO3173_snp_geno_gigamuga_female.pylmm.kin --phenofile DO3173_snp_geno_gigamuga_female.rz.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".f.pylmm.out || true"),
        job.para[i, "script.file"], append = T)

  write(paste0("sbatch ", job.para[i, "script.file"]), do.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            do.pylmm.job, append = T)
}

# 10.2, GEMMA on DO study based on gigamuga grid-----------------------------------------------------------------------
# job.para
job.para <- csna.allqtl.pheno %>%
  dplyr::select(-1:-2) %>%
  colnames(.) %>%
  data.frame(pheno = .,
             p     = 1:(length(.)))

# gemma job
do.gemma.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/do.gemma.job.sh"
#Write script
write("#!/bin/bash",                                                         do.gemma.job, append = F)
write(paste0("#SBATCH -p compute -q batch"),                                 do.gemma.job, append = T)
write("#SBATCH -N 1 # number of nodes",                                      do.gemma.job, append = T)
write(paste0("#SBATCH -n 1 # number of cores"),                              do.gemma.job, append = T)
write("#SBATCH --mem=64GB # memory pool for all cores",                      do.gemma.job, append = T)
write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         do.gemma.job, append = T)
write(paste0("#SBATCH -e ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/do.gemma.job.stderr"), do.gemma.job, append = T)
write(paste0("#SBATCH -o ", "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/do.gemma.job.stdout"), do.gemma.job, append = T)
write("module load singularity",                                             do.gemma.job, append = T)
write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",            do.gemma.job, append = T)
for(i in 1:dim(job.para)[[1]]){#pheno
  #male
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile DO3173_snp_geno_gigamuga_male --pheno DO3173_snp_geno_gigamuga_male.rz.phenos --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out DO3173_snp_geno_gigamuga_male_p"),  do.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_gigamuga_male_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/ -o DO3173_snp_geno_gigamuga_male_p"), do.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_gigamuga_male_p",
               " -k ", "DO3173_snp_geno_gigamuga_male_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/gigamuga -o ",
               job.para[i,"pheno"], ".m"), do.gemma.job, append = T)
  #female
  #Generate bed file
  write(paste0("/projects/compsci/vmp/USERS/heh/software/plink --noweb --bfile DO3173_snp_geno_gigamuga_female --pheno DO3173_snp_geno_gigamuga_female.rz.phenos --mpheno ",
               job.para[i, "p"],
               " --make-bed --recode --out DO3173_snp_geno_gigamuga_female_p"),  do.gemma.job, append = T)
  #Generate relatedness matrix:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_gigamuga_female_p",
               " -gk -outdir /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/ -o DO3173_snp_geno_gigamuga_female_p"), do.gemma.job, append = T)
  #Perform the association test:"
  write(paste0("/projects/compsci/vmp/USERS/heh/software/gemma/gemma -bfile DO3173_snp_geno_gigamuga_female_p",
               " -k ", "DO3173_snp_geno_gigamuga_female_p.cXX.txt ",
               "-lmm 4 ",
               "-outdir /projects/csna/meta-analysis/output/t2d_addiction/gemma/gigamuga -o ",
               job.para[i,"pheno"], ".f"), do.gemma.job, append = T)
  write("\n", do.gemma.job, append = T)
}



#collect results
gm = readRDS("/projects/chesler-lab/phenome/metasoft/data/GM_marker_positions_with_MUSter_v2_genotypes.RDS") %>%
  dplyr::filter(chr != "M") %>%
  dplyr::mutate(rs = if_else(rs == "", paste0(chr, "_", bp38), rs))

map(colnames(csna.allqtl.pheno)[-1:-2], function(x){
  df.m <- read.table(file = paste0("output/t2d_addiction/pylmm/gigamuga0/", x, ".m.pylmm.out"),
                     header = TRUE) %>%
    dplyr::rename(marker = SNP_ID) %>%
    dplyr::left_join(gm[,3:4]) %>%
    dplyr::filter(!is.na(rs)) %>%
    dplyr::rename(SNP_ID = rs) %>%
    dplyr::relocate(SNP_ID, .before = 1) %>%
    dplyr::select(-2) %>%
    write.table(., file = paste0("output/t2d_addiction/pylmm/gigamuga/", x, ".m.pylmm.out"),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

  df.f <- read.table(file = paste0("output/t2d_addiction/pylmm/gigamuga0/", x, ".f.pylmm.out"),
                     header = TRUE) %>%
    dplyr::rename(marker = SNP_ID) %>%
    dplyr::left_join(gm[,3:4]) %>%
    dplyr::filter(!is.na(rs)) %>%
    dplyr::rename(SNP_ID = rs) %>%
    dplyr::relocate(SNP_ID, .before = 1) %>%
    dplyr::select(-2) %>%
    write.table(., file = paste0("output/t2d_addiction/pylmm/gigamuga/", x, ".f.pylmm.out"),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
})
