library(foreach)
library(doParallel)
library(parallel)
library(abind)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(DOQTL)

##########
# Helper function to get Sanger SNPs and place them in the DO samples.
do2sanger.helper = function(snp.file, chr, start, end, geno, keep) {
  sanger = get.variants(file = snp.file, chr = chr, start = start, end = end,
                        strains = c("129S1/SvImJ", "A/J", "C57BL/6J", "CAST/EiJ",
                                    "NOD/ShiLtJ", "NZO/HlLtJ", "PWK/PhJ", "WSB/EiJ"), quality = 1)
  if(length(sanger) > 0 && nrow(sanger) > 0) {
    sanger = sanger[,keep]
    geno.cols = (ncol(sanger) - 7):ncol(sanger)
    map = do.colors[,1:2]
    map = map[match(colnames(sanger)[geno.cols], map[,2]),]
    colnames(sanger)[geno.cols] = map[,1]

    # Here, we assume that the founders are homozygous. This may not be the case
    # and this is an area for improvement.
    sanger[,geno.cols] = apply(sanger[,geno.cols], 2, substring, 1, 1)
    hdr = sanger[,1:5]
    sanger = as.matrix(sanger[,-1:-5])
    # Turn the sanger SNPs into 0 for reference allele and 1 for alternate allele.
    sanger = matrix(as.numeric(sanger != hdr[,4]), nrow(sanger), ncol(sanger),
                    dimnames = dimnames(sanger))
    # Create a 2 x n matrix with founder haplotypes for each DO founder.
    samples = names(geno)
    geno = matrix(unlist(strsplit(geno, split = "")), nrow = 2, dimnames =
                    list(NULL, samples))
    # Convert the founder haplotypes to numeric genotypes.
    ##### NOTE: We are assuming bi-allelic SNPs here! This may not be the case
    #####       for all Sanger SNPs.
    geno = matrix(sanger[,geno[1,]] + sanger[,geno[2,]],
                  nrow(sanger), ncol(geno), dimnames = list(NULL, samples))
    # Code the major allele as 0 .
    num.zero = rowMeans(geno == 0, na.rm = TRUE)
    num.two  = rowMeans(geno == 2, na.rm = TRUE)
    swap.rows = which(num.two > num.zero)
    geno[swap.rows,] = 2 - geno[swap.rows,]
    return(cbind(hdr, geno))
  } else {
    return(NULL)
  } # else
} # do2sanger.helper()

#' Given a set of DO genotype probability generated/interpolated by qtl2 and the location of the Tabix indexed Sanger file,
#' impute the Sanger SNPs on to DO genomes.
#' @param pr - genotype probability. 3-D array, sample*36*snp
#' @param snps - Data.frame containing the marker locations. SNP ID, chromosome and Mb in columns 1 through 3, respectively.
#' @param out.dir - out directory
#' @param output.file - Character string to write the results to.
#' @param snp.file - Character string with path to a Tabix indexed SNP file.
#' @param cores - numeric for cores used
#' @return  - pmap list required by qtl2
#'
predict_sanger_snp_geno <- function(pr, snps, out.dir, output.file, snp.file, cores) {
  if (missing(snps)) {
    stop(paste("Please provide SNP locations that match the SNPs in the genotype probabilities"))
  }
  prsmth = t(pr[1,,])
  states = colnames(prsmth)
  snps = snps[snps[, 1] %in% rownames(prsmth), ]
  prsmth = prsmth[rownames(prsmth)[[3]] %in% snps[, 1], ]
  gt = matrix("", dim(pr)[[1]], nrow(prsmth), dimnames = list(dimnames(pr)[[1]], rownames(prsmth)))
  print("Calling DO genotypes...")
  for (i in 1:dim(pr)[[1]]) {
    gt[i, ] = states[apply(t(pr[i,,]), 1, which.max)]
  }
  snps = snps[snps[, 1] %in% rownames(prsmth), ]
  snps = split(x = snps, f = snps[, 2])
  sanger = get.variants(file = snp.file, chr = 1, start = 0,
                        end = 4, strains = c("129S1/SvImJ", "A/J", "C57BL/6J",
                                             "CAST/EiJ", "NOD/ShiLtJ", "NZO/HlLtJ", "PWK/PhJ",
                                             "WSB/EiJ"))
  keep = setdiff(1:ncol(sanger), grep("quality", colnames(sanger)))

  #setup parallel backend to use many processors
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  #log files
  writeLines(c(""), paste0(out.dir, "predict_sanger_snp_geno.log"))

  foreach(c=1:length(snps), .packages='DOQTL', .errorhandling = 'pass') %dopar% {
    source(file = "code/do2sanger.helper.R")
    #log
    sink(paste0(out.dir, "predict_sanger_snp_geno.log"), append=TRUE)
    #chr
    curr.chr = names(snps)[c]
    cat(paste("CHR", curr.chr, "\n"))
    num.snps = 0
    cat("Mapping Sanger SNPs to DO samples...\n")
    con = file(description = paste(out.dir, "chr", curr.chr, ".",
                                   output.file, sep = ""), open = "w")
    writeLines(paste("ID", "CHR", "GRC38.bp", "Ref", "Alt",
                     paste(rownames(gt), collapse = "\t"), sep = "\t"),
               con)
    gc()
    result = do2sanger.helper(snp.file = snp.file,
                              chr      = curr.chr,
                              start    = 0,
                              end      = mean(snps[[curr.chr]][1:2, 3]),
                              geno     = gt[, snps[[curr.chr]][1, 1]],
                              keep     = keep)
    if (!is.null(result)) {
      result = apply(result, 1, paste, collapse = "\t")
      writeLines(result, con)
      num.snps = num.snps + nrow(result)
    }
    for (i in 2:(nrow(snps[[curr.chr]]) - 1)) {
      if (i%%100 == 0) {
        cat(paste0(i, "\n"))
      }
      if (all(!is.na(snps[[curr.chr]][(i - 1):(i + 1), 3]))) {
        result = do2sanger.helper(snp.file = snp.file,
                                  chr      = curr.chr,
                                  start    = mean(snps[[curr.chr]][(i - 1):i, 3]),
                                  end      = mean(snps[[curr.chr]][i:(i + 1), 3]),
                                  geno     = gt[, snps[[curr.chr]][i,1]],
                                  keep     = keep)
        if (!is.null(result)) {
          result = apply(result, 1, paste, collapse = "\t")
          writeLines(result, con)
          num.snps = num.snps + nrow(result)
        }
      }
      gc()
    }
    if (all(!is.na(snps[[curr.chr]][(i - 1):(i + 1), 3]))) {
      result = do2sanger.helper(snp.file = snp.file,
                                chr      = curr.chr,
                                start    = mean(snps[[curr.chr]][(nrow(snps[[curr.chr]]) -1):nrow(snps[[curr.chr]]), 3]),
                                end      = 200,
                                geno     = gt[, snps[[curr.chr]][nrow(snps[[curr.chr]]), 1]],
                                keep     = keep)
      if (!is.null(result)) {
        result = apply(result, 1, paste, collapse = "\t")
        writeLines(result, con)
        num.snps = num.snps + nrow(result)
      }
    }
    close(con)
    gc()
    cat(paste("CHR", curr.chr, ":", num.snps, "\n"))
  }
  #stop cluster
  stopCluster(cl)
}

#' Given a set of all snp_geno.files location and name generated by predict_sanger_snp_geno, to generate the ped/map plink format files
#' @param snp_geno.files - all snp_geno.files location and name generated by predict_sanger_snp_geno
#' @param sex.df - a data.frame containing two columns "id" and "sex". (1=male; 2=female; other=unknown)
#' @param snps - Data.frame containing the marker locations. SNP ID, chromosome and Mb in columns 1 through 3, respectively.
#' @param output.file - Character string to write the results to. eg. "data/meta/do_geno_132k"
#' @param cores - numeric for cores used
#' @return  - .ped/.map file
#'
sanger_snp_geno_2pedmap <- function(snp_geno.file, sex.df, snps, output.file, cores) {
  #setup parallel backend to use many processors
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  #read, rbind and filter
  df <- foreach(f = 1:length(snp_geno.file), .combine = rbind, .packages = c("data.table", "tidyverse")) %dopar% {
    snp_geno <- fread(snp_geno.file[[f]])
    snp_geno <- snp_geno %>%
      mutate(ID = case_when(
        (ID == ".") ~ paste0(CHR, "_", GRC38.bp),
        TRUE ~ as.character(ID)
      )) %>%
      filter(ID %in% snps[, "marker"])
    snp_geno
  }
  #stop cluster
  stopCluster(cl)
  #generate map file, four columns, chr, rs, 0, pos in Base-pair coordinate
  map <- df %>%
    select(CHR, ID, GRC38.bp) %>%
    mutate(cM = 0, .after = ID)
  write.table(map, file = paste0(output.file, ".map"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  #genreate ped file
  ##SNPs coded as 0, 1 or 2 into "1 1", "1 2", "2 2" for ped file
  df2 <- df[,-1:-5]
  df2[df2 == 0] <- "1 1"
  df2[df2 == 1] <- "1 2"
  df2[df2 == 2] <- "2 2"
  ped <- df2 %>%
    t(.) %>%
    as_tibble(.name_repair = "universal") %>%
    mutate(fid = colnames(df)[-1:-5],
           iid = colnames(df)[-1:-5],
           pid = 0,
           mid = 0,
           phe = -9,
           .before = 1) %>%
    left_join(sex.df, by = c("fid" = "id")) %>%
    relocate(sex, .after = mid)
  write.table(ped, file = paste0(output.file, ".ped"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
}
