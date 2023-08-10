## ---------------------------
##
## Script name:
##
## Purpose of script: Running pyLMM on UCLA and GigaMuga density imputed grid for Attie DO
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

if(FALSE){
  # 1, get the cross2 object for 232_Attie_DO_Islets-----------------------------------------------------------------------
  #pull the 232_Attie_DO_Islets data from dodb and process to get the cross2 object
  attie <- read_cross2(file = "data/thecube/DO500_attie/Attie-232_Attie_DO_Islets-GigaMUGA.json")
  attie
  #Let’s omit markers without any genotype data and re-code the genotypes
  #so that the first allele is that which is most common among the founders.
  attie <- drop_nullmarkers(attie)
  #Dropping markers with no data
  for(chr in seq_along(attie$founder_geno)) {
    fg <- attie$founder_geno[[chr]]
    g <- attie$geno[[chr]]
    f1 <- colSums(fg==1)/colSums(fg != 0)

    fg[fg==0] <- NA
    g[g==0] <- NA

    fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
    g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]

    fg[is.na(fg)] <- 0
    g[is.na(g)] <- 0

    attie$founder_geno[[chr]] <- fg
    attie$geno[[chr]] <- g
  }
  attie
  #remove one duplicate AA-DO021f
  attie <- attie[ind_ids(attie)[ind_ids(attie) != "AA-DO021f"],]
  #get the new id
  new_ids <- str_remove(ind_ids(attie), "AA-")
  new_ids <- str_remove(new_ids, "f")
  new_ids <- str_remove(new_ids, "m")
  new_ids <- setNames(new_ids,
                      ind_ids(attie))
  attie_gm <- replace_ids(attie, new_ids)
  attie_gm <- attie_gm[, c(1:19, "X")]
  attie_gm
  save(attie_gm, file = "data/thecube/DO500_attie/gm_DOCUBE500.RData")

  # 2, calculate genoprobs by qtl2 -----------------------------------------------------------------------
  pr <- calc_genoprob(attie_gm, cores=20)
  save(pr, file = "data/thecube/DO500_attie/pr_DOCUBE500.RData")

  load("data/thecube/DO500_attie/pr_DOCUBE500.RData")
  # 3, interpolate to 132k grid by DOQTL -----------------------------------------------------------------------
  if(TRUE){
    #x chr
    if(dim(pr$X)[2] == 44){
      #X chromosome should have 36 states
      pr$X <- pr$X[,(-37:-44),]
    }

    str(pr)

    print(paste0("Number of samples, ", dim(pr[[1]])[1]))
    #interpolate in chunks
    chunk_size = 50
    idx <- list()
    pr.132k <- list()
    for(chunk_number in 1:10){
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
      grid <- read.csv("data/meta/snp_UCLA1_full.csv", header = TRUE, fill = TRUE, na.strings = "") %>%
        select(1:3) %>%
        mutate(rs = case_when(
          is.na(rs) ~ paste0(chr, "_", bp38),
          TRUE ~ as.character(rs)
        )) %>%
        mutate(bp38 = bp38/10^6)  %>%
        select(marker = rs,
               chr = chr,
               pos = bp38) %>%
        arrange(chr, pos)
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
        pr.132k[[chunk_number]] <- mclapply(X = 1:dim(sub.pr.3d)[[1]], FUN = interplote_fun, mc.preschedule = F, mc.cores = 20)
      })
      names(pr.132k[[chunk_number]]) <- dimnames(sub.pr.3d)[[1]]
      pr.132k[[chunk_number]] <- do.call("abind",list(pr.132k[[chunk_number]],along = 3))
      pr.132k[[chunk_number]] <- aperm(pr.132k[[chunk_number]], perm = c(3,1,2))
      print(str(pr.132k[[chunk_number]]))
    }
    pr.132k <- do.call("abind",list(pr.132k,along = 1))
    save(pr.132k, file = "data/thecube/DO500_attie/pr_132k_DOCUBE500.RData")
    print("pr_132k_DOCUBE500 done")
  }

  # 4, predict_sanger_snp_geno ---------------------------------------------------------------------
  #snps
  snps <- read.csv("data/meta/snp_UCLA1_full.csv", header = TRUE, fill = TRUE, na.strings = "") %>%
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

  # predict_sanger_snp_geno
  sanger_snp_geno <- predict_sanger_snp_geno(pr = pr.132k,
                                             snps =  snps,
                                             out.dir = "output/thecube/DO500_attie/",
                                             output.file = "snp_geno_132k.txt",
                                             snp.file = "data/meta/cc.snps.NCBI38.txt.gz",
                                             cores = 20)

  #sex df. Sex (1=male; 2=female; other=unknown)
  sex.df <- attie_gm$covar %>%
    rownames_to_column(var = "name") %>%
    dplyr::select(id    = name,
                  sex   = sex) %>%
    mutate(sex = case_when(
      sex == "M" ~ 1,
      sex == "F" ~ 2
    ))

  # 5, covert *.snp_geno_132k.txt into ped and map file ---------------------------------------------------------------------
  sanger_snp_geno_pedmap <- sanger_snp_geno_2pedmap(snp_geno.file = paste0("output/thecube/DO500_attie/chr", c(1:19, "X"), ".snp_geno_132k.txt"),
                                                    sex.df        = sex.df,
                                                    snps          = snps,
                                                    output.file   = "output/thecube/DO500_attie/DO500_attie_snp_geno_132k",
                                                    cores         = 10)
}

# 6, prepare files for pylmm ---------------------------------------------------------------------
# for ped and map file, we need to
#1, filter by sex, create bed/fam/bim files.
#2, calculated the kinship matrix file using pylmmKinship.py
#3, create phenotype file for each sex, NA or -9 to denote missing values.

#1 filter by sex, create bed/fam/bim files.
#cd /projects/csna/csna_workflow/output/thecube/DO500_attie
#/projects/compsci/USERS/heh/software/plink --noweb --file DO500_attie_snp_geno_132k --filter-males   --make-bed --out DO500_attie_snp_geno_132k_male
#/projects/compsci/USERS/heh/software/plink --noweb --file DO500_attie_snp_geno_132k --filter-females --make-bed --out DO500_attie_snp_geno_132k_female

#2 calculated the kinship matrix file using pylmmKinship.py
#python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO500_attie_snp_geno_132k_male   DO500_attie_snp_geno_132k_male.pylmm.kin
#python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO500_attie_snp_geno_132k_female DO500_attie_snp_geno_132k_female.pylmm.kin

system("cp /projects/csna/csna_workflow/output/thecube/DO500_attie/DO500_attie_snp_geno_132k_male* /projects/csna/meta-analysis/data/t2d_addiction/ucla/")
system("cp /projects/csna/csna_workflow/output/thecube/DO500_attie/DO500_attie_snp_geno_132k_female* /projects/csna/meta-analysis/data/t2d_addiction/ucla/")

#3, create phenotype file for each sex, NA or -9 to denote missing values.
attieDO <- readxl::read_excel(path = "data/measurements_w_HFD.xlsx",
                              sheet = 2) %>%
  dplyr::filter(Analyze == 1) %>%
  dplyr::filter(paneldesc == "DO population") #DO population

pheno <- extract_pheno_measure(url = "https://phenome.jax.org/", measure_id = unique(attieDO$measnum))
#change factors from long to wide format
#female
pheno_w_f <-  pheno %>%
  dplyr::select(c(1,3,5,7,11)) %>%
  dplyr::filter(sex == "f") %>%
  dplyr::select(c(-1, -3)) %>%
  pivot_wider(names_from = varname, values_from = zscore)
attieDO_pheno_f <- read.table("/projects/csna/csna_workflow/output/thecube/DO500_attie/DO500_attie_snp_geno_132k_female.fam") %>%
  dplyr::select(1:2) %>%
  left_join(., pheno_w_f,by = c("V1" = "animal_id"))
#female pheno
write.table(attieDO_pheno_f,
            file = "data/t2d_addiction/ucla/DO500_attie_snp_geno_132k_female.phenos",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

#male
pheno_w_m <-  pheno %>%
  dplyr::select(c(1,3,5,7,11)) %>%
  dplyr::filter(sex == "m") %>%
  dplyr::select(c(-1, -3)) %>%
  pivot_wider(names_from = varname, values_from = zscore)
attieDO_pheno_m <- read.table("/projects/csna/csna_workflow/output/thecube/DO500_attie/DO500_attie_snp_geno_132k_male.fam") %>%
  dplyr::select(1:2) %>%
  left_join(., pheno_w_m,by = c("V1" = "animal_id"))
#male pheno
write.table(attieDO_pheno_m,
            file = "data/t2d_addiction/ucla/DO500_attie_snp_geno_132k_male.phenos",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

# 7, pylmmgwas job for ucla--------------------------------------------------------------
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/scripts/"
job.para <- data.frame(pheno = colnames(attieDO_pheno_m)[-1:-2],
                       p     = 0:(length(colnames(attieDO_pheno_m)[-1:-2])-1)) %>%
  mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
         stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
         stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
AttieDO.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/ucla/job/AttieDO.pylmm.job.sh"
#Write script
write("#Submit job", AttieDO.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=8GB # memory pool for all cores",                       job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/ucla/",         job.para[i, "script.file"], append = T)
  #male
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO500_attie_snp_geno_132k_male --kfile DO500_attie_snp_geno_132k_male.pylmm.kin --phenofile DO500_attie_snp_geno_132k_male.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".m.pylmm.out || true"), job.para[i, "script.file"], append = T)
  #female
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO500_attie_snp_geno_132k_female --kfile DO500_attie_snp_geno_132k_female.pylmm.kin --phenofile DO500_attie_snp_geno_132k_female.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/ucla/", job.para[i,"pheno"], ".f.pylmm.out || true"), job.para[i, "script.file"], append = T)

  write(paste0("sbatch ", job.para[i, "script.file"]), AttieDO.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            AttieDO.pylmm.job, append = T)
}

# gigamuga ----------------------------------------------------------------
if(FALSE){
  source("/projects/csna/csna_workflow/code/predict_sanger_snp_geno.R")
  source("/projects/csna/csna_workflow/code/reconst_utils.R")
  source("/projects/csna/csna_workflow/code/impute_dogenome.R")

  # 1, get the cross2 object for 232_Attie_DO_Islets
  load("/projects/csna/csna_workflow/data/thecube/DO500_attie/gm_DOCUBE500.RData")
  load("/projects/csna/csna_workflow/data/thecube/DO500_attie/pr_DOCUBE500.RData")
  gigamuga <- readRDS("/projects/chesler-lab/phenome/metasoft/data/GM_marker_positions_with_MUSter_v2_genotypes.RDS")
  # 2, predict_sanger_snp_geno
  #snps
  snps <- gigamuga %>%
    dplyr::select(1:4) %>%
    mutate(rs = case_when(
      rs == "" ~ paste0(chr, "_", bp38),
      TRUE ~ as.character(rs)
    )) %>%
    mutate(bp38 = bp38/10^6)  %>%
    dplyr::select(marker,
                  chr = chr,
                  pos = bp38,
                  snp = rs) %>%
    dplyr::distinct(., marker, .keep_all = TRUE) %>%
    dplyr::filter(chr != "M")
  rownames(snps) <- snps$marker

  #x chr
  if(dim(pr$X)[2] == 44){
    #X chromosome should have 36 states
    pr$X <- pr$X[,(-37:-44),]
  }
  pr = do.call("abind",list(pr,along = 3))
  # predict_sanger_snp_geno
  sanger_snp_geno <- predict_sanger_snp_geno(pr =  pr,
                                             snps =  snps,
                                             out.dir = "data/t2d_addiction/gigamuga/",
                                             output.file = "DO500_attie_snp_geno_gigamuga.txt",
                                             snp.file = "/projects/csna/csna_workflow/data/meta/cc.snps.NCBI38.txt.gz",
                                             cores = 20)

  #sex df. Sex (1=male; 2=female; other=unknown)
  sex.df <- attie_gm$covar %>%
    rownames_to_column(var = "name") %>%
    dplyr::select(id    = name,
                  sex   = sex) %>%
    mutate(sex = case_when(
      sex == "M" ~ 1,
      sex == "F" ~ 2
    ))

  # 5, covert *.snp_geno_132k.txt into ped and map file ---------------------------------------------------------------------
  snps <- gigamuga %>%
    dplyr::select(1:4) %>%
    mutate(rs = case_when(
      rs == "" ~ paste0(chr, "_", bp38),
      TRUE ~ as.character(rs)
    )) %>%
    mutate(bp38 = bp38/10^6)  %>%
    dplyr::select(marker = rs,
                  chr = chr,
                  pos = bp38,
                  snp = marker) %>%
    dplyr::distinct(., marker, .keep_all = TRUE) %>%
    dplyr::filter(chr != "M")
  rownames(snps) <- snps$marker
  sanger_snp_geno_pedmap <- sanger_snp_geno_2pedmap(snp_geno.file = paste0("data/t2d_addiction/gigamuga/chr", c(1:19, "X"), ".DO500_attie_snp_geno_gigamuga.txt"),
                                                    sex.df        = sex.df,
                                                    snps          = snps,
                                                    output.file   = "data/t2d_addiction/gigamuga/DO500_attie_snp_geno_gigamuga",
                                                    cores         = 20)
}
# 6, prepare files for pylmm ---------------------------------------------------------------------
# for ped and map file, we need to
#1, filter by sex, create bed/fam/bim files.
#2, calculated the kinship matrix file using pylmmKinship.py
#3, create phenotype file for each sex, NA or -9 to denote missing values.

#1 filter by sex, create bed/fam/bim files.
#cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga
#/projects/csna/meta-analysis/code/Plink2/plink2 --pedmap DO500_attie_snp_geno_gigamuga --rm-dup force-first --filter-males   --make-bed --out DO500_attie_snp_geno_gigamuga_male
#/projects/csna/meta-analysis/code/Plink2/plink2 --pedmap DO500_attie_snp_geno_gigamuga --rm-dup force-first --filter-females --make-bed --out DO500_attie_snp_geno_gigamuga_female

#2 calculated the kinship matrix file using pylmmKinship.py
#python2 /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO500_attie_snp_geno_gigamuga_male   DO500_attie_snp_geno_gigamuga_male.pylmm.kin
#python2 /projects/csna/csna_workflow/code/pylmm/scripts/pylmmKinship.py --bfile DO500_attie_snp_geno_gigamuga_female DO500_attie_snp_geno_gigamuga_female.pylmm.kin

#3, create phenotype file for each sex, NA or -9 to denote missing values.
attieDO <- readxl::read_excel(path = "data/measurements_w_HFD.xlsx",
                              sheet = 2) %>%
  dplyr::filter(Analyze == 1) %>%
  dplyr::filter(paneldesc == "DO population") #DO population

pheno <- extract_pheno_measure(url = "https://phenome.jax.org/", measure_id = unique(attieDO$measnum))
#change factors from long to wide format
#female
pheno_w_f <-  pheno %>%
  dplyr::select(c(1,3,5,7,11)) %>%
  dplyr::filter(sex == "f") %>%
  dplyr::select(c(-1, -3)) %>%
  pivot_wider(names_from = varname, values_from = zscore)
attieDO_pheno_f <- read.table("/projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO500_attie_snp_geno_gigamuga_female.fam") %>%
  dplyr::select(1:2) %>%
  left_join(., pheno_w_f,by = c("V1" = "animal_id"))
#female pheno
write.table(attieDO_pheno_f,
            file = "data/t2d_addiction/gigamuga/DO500_attie_snp_geno_gigamuga_female.phenos",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

#male
pheno_w_m <-  pheno %>%
  dplyr::select(c(1,3,5,7,11)) %>%
  dplyr::filter(sex == "m") %>%
  dplyr::select(c(-1, -3)) %>%
  pivot_wider(names_from = varname, values_from = zscore)
attieDO_pheno_m <- read.table("/projects/csna/meta-analysis/data/t2d_addiction/gigamuga/DO500_attie_snp_geno_gigamuga_male.fam") %>%
  dplyr::select(1:2) %>%
  left_join(., pheno_w_m,by = c("V1" = "animal_id"))
#male pheno
write.table(attieDO_pheno_m,
            file = "data/t2d_addiction/gigamuga/DO500_attie_snp_geno_gigamuga_male.phenos",
            row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

# 7, pylmmgwas job for --------------------------------------------------------------
script.Dir <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/scripts/"
job.para <- data.frame(pheno = colnames(attieDO_pheno_m)[-1:-2],
                       p     = 0:(length(colnames(attieDO_pheno_m)[-1:-2])-1)) %>%
  mutate(script.file = paste0(script.Dir, pheno,  ".pylmm.sh"),
         stderr      = paste0(script.Dir, pheno,  ".pylmm.err"),
         stdout      = paste0(script.Dir, pheno,  ".pylmm.out")) #name and path of script
AttieDO.pylmm.job <- "/projects/csna/meta-analysis/analysis/t2d_addiction/gigamuga/job/AttieDO.pylmm.job.sh"
#Write script
write("#Submit job", AttieDO.pylmm.job, append = F)
for(i in 1:dim(job.para)[[1]]){#pheno
  #Write script
  write("#!/bin/bash",                                                         job.para[i, "script.file"], append = F)
  write(paste0("#SBATCH -p compute -q batch"),                                 job.para[i, "script.file"], append = T)
  write("#SBATCH -N 1 # number of nodes",                                      job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -n 1 # number of cores"),                              job.para[i, "script.file"], append = T)
  write("#SBATCH --mem=8GB # memory pool for all cores",                       job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -t 3-00:00 # time (D-HH:MM)"),                         job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -o ", job.para[i, "stderr"]),                          job.para[i, "script.file"], append = T)
  write(paste0("#SBATCH -e ", job.para[i, "stdout"]),                          job.para[i, "script.file"], append = T)

  write("module load singularity",                                             job.para[i, "script.file"], append = T)
  write("cd /projects/csna/meta-analysis/data/t2d_addiction/gigamuga/",         job.para[i, "script.file"], append = T)
  #male
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO500_attie_snp_geno_gigamuga_male --kfile DO500_attie_snp_geno_gigamuga_male.pylmm.kin --phenofile DO500_attie_snp_geno_gigamuga_male.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".m.pylmm.out || true"), job.para[i, "script.file"], append = T)
  #female
  write(paste0("python /projects/csna/csna_workflow/code/pylmm/scripts/pylmmGWAS.py -v --bfile DO500_attie_snp_geno_gigamuga_female --kfile DO500_attie_snp_geno_gigamuga_female.pylmm.kin --phenofile DO500_attie_snp_geno_gigamuga_female.phenos ",
               "-p ", job.para[i,"p"], " /projects/csna/meta-analysis/output/t2d_addiction/pylmm/gigamuga/", job.para[i,"pheno"], ".f.pylmm.out || true"), job.para[i, "script.file"], append = T)

  write(paste0("sbatch ", job.para[i, "script.file"]), AttieDO.pylmm.job, append = T)
  write(paste0("sleep 3 "),                            AttieDO.pylmm.job, append = T)
}
