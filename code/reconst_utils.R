library(qtl2)
library(DOQTL)
library(jsonlite)

################################################################################
# Interpolate the genoprobs of a matrix from one set of markers to another.
# Daniel Gatti
# Dan.Gatti@jax.org
# Dec. 20, 2014
################################################################################
# Arguments: data: numeric matrix containing the genotype or haplotype
#                  probabilities. Named markers in rows and genotypes in
#                  columns.
#            from: marker information (names, chr, pos) for the markers to
#                  impute from. Must match the marker names in data.
#            to:  marker information (names, chr, pos) for the markers to
#                 impute to. May be higher or lower density than from.
interpolate.markers = function(data, from, to) {
  # Verify that the number of row in data and from matches.
  if(nrow(data) != nrow(from)) {
    stop(paste("Number of rows in data (", nrow(data), ") must equal the",
               "number or rows in from (", nrow(from) ,")."))
  } # if(nrow(data) != nrow((from))
  # Verify that the marker names in data and from are identical.
  if(any(rownames(data) != from[,1])) {
    stop(paste("The rownames in \'data\' do not equal the marker names",
               "in column 1 of \'from\'."))
  } # if(any(rownames(data) != from[,1]))
  # If no 'to' argument was provided, jsut return the data.
  if(missing(to) || length(to) == 0) {
    return(data)
  } # if(missing(to) || length(to) == 0)
  data = as.matrix(data)
  if(!is.numeric(data)) {
    stop(paste("\'data\' must be a numeric matrix."))
  } # if(!is.numeric(data))
  # Figure out which marker set has higher density.
  higher = from
  lower = to
  if(nrow(to) > nrow(from)) {
    higher = to
    lower  = from
  } # if(nrow(to) > nrow(from))
  # Impute the values by chromosome. Use the chromosomes in the smaller set.
  chr = intersect(unique(from[,2]), unique(to[,2]))
  from = from[from[,2] %in% chr,]
  to   = to[to[,2] %in% chr,]
  newdata = matrix(0, nrow(to), ncol(data), dimnames =
                     list(to[,1], colnames(data)))
  for(c in chr) {
    from.c = from[from[,2] == c,]
    to.c   = to[to[,2] == c,]
    data.c    = data[from.c[,1],]
    newdata.c = newdata[to.c[,1],]

    # Interior markers.
    to.prox = outer(from.c[,3], to.c[,3], "<=")
    to.dist = outer(from.c[,3], to.c[,3], ">")
    to.prox = sapply(apply(to.prox, 2, which), max)
    to.dist = sapply(apply(to.dist, 2, which), min)
    # Trim off the ends where the array grids don't overlap.
    rng.to = which(is.finite(to.prox) & is.finite(to.dist))
    to.prox = to.prox[rng.to]
    to.dist = to.dist[rng.to]
    denom1 = from.c[to.dist,3] - from.c[to.prox,3]
    num1 = to.c[rng.to,3] - from.c[to.prox,3]
    num2 = from.c[to.dist,3] - to.c[rng.to,3]
    newdata.c[rng.to,] = (num1 * data.c[to.prox,] + num2 * data.c[to.dist,]) /
      denom1
    # Start of chromosome.
    # If the first marker in 'from' is > the first marker in 'to', then
    # copy the values in from[1,] to to[1,].
    to.lt.from = which(to.c[,3] < from.c[1,3])
    if(length(to.lt.from) > 0) {
      newdata.c[to.lt.from,] = matrix(newdata.c[min(rng.to),], length(to.lt.from),
                                      ncol(newdata), byrow = T)
    } # if(length(to.lt.from) > 0)
    # End of chromosome.
    to.gt.from = which(to.c[,3] > from.c[nrow(from.c),3])
    if(length(to.gt.from) > 0) {
      newdata.c[to.gt.from,] = matrix(newdata.c[max(rng.to),], length(to.gt.from),
                                      ncol(newdata), byrow = T)
    } # if(length(to.gt.from) > 0)
    newdata[to[,2] == c,] = newdata.c
  } # for(c)
  newdata
} # interpolate.markers()

# Read the requested remote path and look for a *FinalReport and a
# Sample_Map.* file. Download the 2 files and return their names.
# remove_user: user name and server where the raw Neogen data is stored
#              in the form "user@server.jax.org"
# path: path to the requested data directory.
# scratch_dir: local scratch directory where files will be copied to.
# Returns:list with *FinalReport and Sample_Map filenames.
get_neogen_files = function(remote_user, path, scratch_dir) {

  # Get a directory listing from the remote host (churchilldev).
  remote.ls = system2(command = "ssh",
                      args = c(remote_user, "ls", path),
                      stdout = TRUE)

  # Get the FinalReport and Sample_Map file names.
  finalreport = remote.ls[grep("FinalReport\\.", remote.ls)]

  if(nchar(finalreport) == 0) {
    stop(paste("Can't find FinalReport in", path))
  }

  samplemap   = remote.ls[grep("Sample_Map", remote.ls)]
  if(nchar(samplemap) == 0) {
    stop(paste("Can't find Sample_Map in", path))
  }

  # Copy the files from churchilldev to a temp directory.
  remote.fr = paste0(path, finalreport)
  local.fr  = paste0(scratch_dir, finalreport)
  ret.code  = system2(command = "scp", args =
                        c(paste0(remote_user, ":", dQuote(sQuote(remote.fr))), scratch_dir))

  ret.code  = system2(command = "scp", args =
                        c(paste0(remote_user, ":", remote.fr), scratch_dir))

  if(ret.code != 0) {
    stop(paste("Could not retrieve", remote.fr))
  } # if(ret.code != 0)

  remote.sm = paste0(path, samplemap)
  local.sm  = paste0(scratch_dir, samplemap)
  ret.code  = system2(command = "scp", args =
                        c(paste0(remote_user, ":", remote.sm), scratch_dir))

  if(ret.code != 0) {
    stop(paste("Could not retrieve", remote.sm))
  } # if(ret.code != 0)

  return(list(local.fr, local.sm))

} # get_neogen_files()


# Read in the *FinalReport and Sample_Map files and return the
# allele calls and X and Y intensities.
# path: directory where the files are.
# markers: marker annotation file from UNC.
read_neogen = function(path, markers) {

  report_file = dir(path = path, pattern = "_FinalReport.zip$",
                    full.names = TRUE)

  report_file_hdr = unlist(strsplit(report_file, split = "/"))
  report_file_hdr = report_file_hdr[length(report_file_hdr)]
  #header = fread(paste("unzip -cq", report_file), nrows = 2, sep = "\t", skip = 0)
  header = read.table(unz(report_file, sub(pattern = "zip$", replacement = "txt", report_file_hdr)),
                      sep = "\t", skip = 1, nrows = 2)
  report = fread(paste("unzip -cq", report_file),
                 skip = 9, showProgress = TRUE, sep = "\t")

  samples = fread(paste0("unzip -cq ",  path, "Sample_Map.zip"))

  colnames(report) = gsub(" |-", "_", colnames(report))

  # It's critical to make Sample_ID a factor with the levels in the
  # order in the file. If you don't, dplyr helpfully sorts them
  # alphabetically and messes up the order of the sample IDs.
  geno = report %>%
    select(Sample_ID, SNP_Name, Allele1___Forward, Allele2___Forward) %>%
    mutate(Sample_ID = factor(Sample_ID, levels = unique(report$Sample_ID))) %>%
    unite(genotype, c("Allele1___Forward", "Allele2___Forward"), sep = "") %>%
    spread(Sample_ID, genotype)

  geno = geno[match(markers$marker, geno$SNP_Name),]
  rn = geno$SNP_Name
  geno = as.matrix(geno[,-1])
  dimnames(geno) = list(rn, samples$ID)

  x = report %>%
    select(Sample_ID, SNP_Name, X) %>%
    mutate(Sample_ID = factor(Sample_ID, levels = unique(report$Sample_ID))) %>%
    spread(Sample_ID, X)
  x = x[match(markers$marker, x$SNP_Name),]
  rn = x$SNP_Name
  x = as.matrix(x[,-1])
  dimnames(x) = list(rn, samples$ID)

  y = report %>%
    select(Sample_ID, SNP_Name, Y) %>%
    mutate(Sample_ID = factor(Sample_ID, levels = unique(report$Sample_ID))) %>%
    spread(Sample_ID, Y)
  y = y[match(markers$marker, y$SNP_Name),]
  rn = y$SNP_Name
  y = as.matrix(y[,-1])
  dimnames(y) = list(rn, samples$ID)

  run_date = strsplit(unlist(header[2,2]), " ")[[1]][1]
  samples = cbind(samples, run_date = run_date)

  return(list(geno = geno, x = x, y = y, samples = samples))

} # read_neogen()


# Write out the X and Y intensities for Chr X, Y & M for sex and Chr M & Y
# haplotyping.
write_geno_intensities = function(geno, x, y, out_file) {

  x_int_file = paste0(dest_dir, "x_chrxym.csv")
  y_int_file = paste0(dest_dir, "y_chrxym.csv")

  geno_chrxym = geno[grepl("[XYM]", attr(geno, "map")$chr), ]

  # Get intensities.
  x_chrxym = attr(geno_chrxym, "intensity")$x
  y_chrxym = attr(geno_chrxym, "intensity")$y

  append = file.exists(x_int_file)

  write.table(t(x_chrxym), file = x_int_file, row.names = TRUE,
              sep = ",", quote = FALSE, append = append, col.names = !append)
  write.table(t(y_chrxym), file = y_int_file, row.names = TRUE,
              sep = ",", quote = FALSE, append = append, col.names = !append)

} # write_geno_intensities()

##########
# Sex the samples and get the Chr Y & M haplogroups.
# x: numeric matrix containing the X intensities. Samples x markers.
# y: numeric matrix containing the Y intensities. Samples x markers.
# markers: data.frame containing marker annotation.
determine_sex_chry_m = function(x, y, markers) {

  sex = determine_sex(x, y, markers)

  ### TND: We need the founder intensities to do this.
  #  males = which(sex == "M")
  #  y_haplogroup = determine_chrm(x[males,], y[males,])

  #  m_haplogroup = determine_chrm(x, y)

  #  return(list(sex = sex, y_haplogroup = y_haplogroup, m_haplogroup = m_haplogroup))

  return(list(sex = sex))

} # determine_sex_chry_m()


determine_sex = function(x, y, markers) {

  chrx = markers$marker[which(markers$chr == "chrX")]
  chry = markers$marker[which(markers$chr == "chrY")]

  chrx_int = colMeans(x[chrx,] + y[chrx,], na.rm = T)
  chry_int = colMeans(x[chry,] + y[chry,], na.rm = T)
  df = data.frame(x = chrx_int, y = chry_int)

  mc = Mclust(data = df, G = 2, model = "VVV")

  # Determine which cluster is females by selecting the cluster with the
  # higher Chr X intensity mean.
  cluster_means = mc$parameters$mean
  female_cl = which.max(cluster_means["x",])

  sex = rep("M", length(mc$classification))
  names(sex) = colnames(x)
  sex[mc$classification == female_cl] = "F"

  return(sex)

} # determine_sex()

# Function to convert AA,CC,GG,TT,-- format to ACGTHN format.
convert_2letter_to1letter = function(data) {

  new_data = data

  new_data[new_data == "AA"] = "A"
  new_data[new_data == "CC"] = "C"
  new_data[new_data == "GG"] = "G"
  new_data[new_data == "TT"] = "T"
  new_data[new_data == "--"] = "N"

  new_data[nchar(new_data) == 2] = "H"

  stopifnot(unique(as.vector(new_data)) %in% c("A", "C", "G", "T", "H", "N"))

  return(new_data)

} # convert_2letter_to1letter()


##########
# Function to convert genotype data to 123 format.
convert_genotypes = function(data, founders) {

  stopifnot(nrow(data) == nrow(founders))
  stopifnot(rownames(data) == rownames(founders))

  data = cbind(data, founders)

  data = t(data)
  dn = dimnames(data)
  data = data.frame(data, stringsAsFactors = F)
  data = lapply(data, factor, levels = c("A", "C", "G", "T", "H", "N"))
  data = lapply(data, as.numeric)
  data = matrix(unlist(data), nrow = length(data[[1]]), ncol = length(data))

  # Record the location of the het calls.
  het = which(data == 5)
  # Set them to NA temporarily.
  data[het] = NA

  # Set '6' to NA, indicating no-call.
  data[data == 6] = NA

  # For each marker, set the lowest allele to 1 and the highest to 3.
  colrng = apply(data, 2, range, na.rm = T)
  colrng[is.infinite(colrng)] = NA

  # Set the '1' allele.
  one = matrix(colrng[1,], nrow = nrow(data), ncol = ncol(data), byrow = T)
  data[data == one] = 1
  rm(one)

  # NOTE: At this point, the '2' values are 'C' alleles, not het calls.

  # If there is only one allele, set it to NA in the column range matrix.
  one.allele = apply(colrng, 2, unique)
  one.allele = which(sapply(one.allele, length) == 1)
  colrng[,one.allele] = NA

  # Set the '3' allele.
  three = matrix(colrng[2,], nrow = nrow(data), ncol = ncol(data), byrow = T)
  data[data == three] = 3
  rm(three)

  # Set het calls = 2.
  data[het] = 2

  dimnames(data) = dn

  data = t(data)

  # Separate out founders.
  wh_founders = which(colnames(data) %in% colnames(founders))
  founders = data[,wh_founders]
  data     = data[,-wh_founders]

  return(list(data = data, founders = founders))

} # convert_genotypes()


# reorganize SNP intensity data (for one marker)
reorg_intensities <-
  function(intensities, marker=NULL)
  {
    if(!is.null(marker)) {
      marker <- as.character(marker)
      wh <- which(intensities[,1] == marker)
      if(sum(wh==0)) stop("Marker ", marker, " not found")
      intensities <- intensities[wh,]
    }

    result <- data.frame(X=as.numeric(intensities[1,-(1:2)]),
                         Y=as.numeric(intensities[2,-(1:2)]))
    rownames(result) <- colnames(intensities)[-(1:2)]

    result
  }

# combine intensities and genotypes
grab_gni <-
  function(marker=NULL, cross=cross,
           intensities_data=NULL, drop_bad_samples=TRUE)
  {
    intensities <- intensities_data

    if(!is.null(marker)) {
      marker <- as.character(marker)
      wh <- which(intensities[,1] == marker)
      if(sum(wh==0)) stop("Marker ", marker, " not found")
      intensities <- intensities[wh,]
    }
    int <- data.frame(X=as.numeric(intensities[1,-(1:2)]),
                      Y=as.numeric(intensities[2,-(1:2)]))
    rownames(int) <- colnames(intensities)[-(1:2)]

    marker <- as.character(marker)

    # get genotypes
    chr <- qtl2::find_markerpos(cross, marker)$chr
    g <- cross$geno[[chr]][,marker]
    gg <- setNames(rep(0, nrow(int)), rownames(int))
    gg[names(g)] <- g

    if(drop_bad_samples) {
      int <- int[names(g),,drop=FALSE]
      gg <- gg[names(g)]
    }

    cbind(int, g=gg)
  }

# load intensities and plot them, colored by genotype calls
#
# drop_bad_samples: if TRUE, don't plot points for samples that are not in the cross object
#                   (e.g., were omitted previously as being bad DNAs)
#
plot_intensities <-
  function(marker, cross=gm,
           intensities_data=NULL, drop_bad_samples=TRUE, geno=NULL, ...)
  {
    if(is.character(marker) || is.factor(marker)) {
      marker <- as.character(marker)
      gni <- grab_gni(marker, cross=gm,
                      intensities_data=intensities_data, drop_bad_samples=drop_bad_samples)
    } else {
      gni <- marker
    }

    if(!is.null(geno)) {
      if(is.logical(geno)) geno <- geno + 1 # FALSE/TRUE -> 1/2
      gni[names(geno),"g"] <- geno
    }

    internal_plot <-
      function(pch=21, bg=broman::brocolors("f2"),
               xlab="allele 1", ylab="allele 2",
               xlim=c(0, max(gni$X, na.rm=TRUE)),
               ylim=c(0, max(gni$Y, na.rm=TRUE)),
               ...)
      {
        grayplot(gni$X, gni$Y, pch=pch, bg=bg[match(gni$g, c(1:3,0))],
                 xlab=xlab, ylab=ylab,
                 xlim=xlim, ylim=ylim, ...)
      }

    internal_plot(...)

    invisible(gni)
  }


#' For a given DO mice allele props returned by qtl2, the function returns a dataframe for founder proportion for samples in generation g
#' @param probs - a list containing 20 arrays for 20 chromosome, each array is 3D one [subjects，allele，markers]#
#' @param gm - a cross2 object. probs is calcuated by using this gm object.
#' @param g - a integer, representing generation information, eg, 22
#' @return - a dataframe
#'
# NOTE: I'm using a lower case L as the beginning of 129 becuase
#       the tidy functions won't allow a number or a "_" at the
#       beginning of a variable name.
names(CCcolors) = c("A_J", "C57BL_6J", "l29S1_SvImJ", "NOD_ShiLtJ",
                    "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")
founder_prop_g = function(probs, gm, g) {
  #subset all the subjects in generation g
  probs <- probs[names(gm$cross_info[gm$cross_info == g,]),]
  #combine all chrs into one 3d array
  probs <- do.call("abind",list(probs,along = 3))
  # Load in the markers.
  load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
  #remove na chromosome
  GM_snps <- GM_snps[!is.na(GM_snps$chr),]
  GM_snps$chr = factor(GM_snps$chr, levels = c(1:19, "X"))
  #subset down to overlapped markers in Gigamuga and probs
  markers = GM_snps[intersect(dimnames(probs)[[3]],GM_snps$marker),1:4]
  # subset to the markers
  probs <- probs[,,markers$marker]
  #starting matrix
  mat = matrix(0, nrow = nrow(probs) * dim(probs)[3], ncol = 8,
               dimnames = list(rep(rownames(probs), dim(probs)[3]),
                               names(CCcolors)))
  for(i in 1:dim(probs)[3]) {
    st = (i - 1) * nrow(probs) + 1
    en = i * nrow(probs)
    mat[st:en,] = probs[,,i]
  } # for(i)
  fp <- data.frame(chr = rep(markers$chr, each = nrow(probs)),
                   pos = rep(markers$pos, each = nrow(probs)), mat)
  # Summarize founder proportion by chromosome.
  fp = fp %>% group_by(chr, pos) %>%
    summarize_all(mean) %>%
    gather(founder, prop, 3:10)
  fp$founder = factor(fp$founder, levels = names(CCcolors))
  fp$gen <- as.factor(g)
  return(fp)
}

#' For a given DO mice allele props returned by qtl2, the function returns a dataframe for founder proportion across all generations. It's used for visualization.
#' @param probs - a list containing 20 arrays for 20 chromosome, each array is 3D one [subjects，allele，markers]#.
#' @param gm - a cross2 object. probs is calcuated by using this gm object.
#' @param seq_g - a vector of integers, representing all generations ub gm object.
#' @return - a dataframe for founder proportion across all generations.
#'
founder_prop = function(probs, gm, cores) {
  fp <- mclapply(X = unique(gm$cross_info[,1]), function(i) {
    founder_prop_g(probs = probs, gm = gm, g = i)
  }, mc.cores = cores)
  #list names
  names(fp) <- unique(gm$cross_info[,1])
  #combine for all generations
  fp <- do.call(rbind.data.frame,fp)
  return(fp)
}

#' Given the lod matrix and marker grid, plot the heatmap
#' @param x         - lod matrix: number of phenotype * number of markers (like 69k grid)
#' @param markers   - marker datafame (eg. 69k grid)
#' @return plot     - QTL Heatmap
qtl.heatmap <- function(x, markers, ...){

  # Get Chr lengths and midpoints.
  map = split(markers[,3], markers[,2])
  chrlen = sapply(map, length)
  chrlen = chrlen[order(as.numeric(names(chrlen)))]
  chrsum = cumsum(chrlen)
  chrmid = c(1, chrsum[-length(chrsum)]) + chrlen / 2
  names(chrmid) = names(chrlen)

  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }

  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,0.25), heights=c(1,0.1))

  # YlOrRd
  ColorRamp   <- brewer.pal(n=length(Lst$zlim[1]:Lst$zlim[2]), "YlOrRd")
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Data Map
  par(mar = c(3,20,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max), ann = FALSE)
  # NOTE: we print every other analyte to make the plot legible.
  at = seq(1, ncol(t(x)), 2)
  mtext(side = 1, line = 1.5, at = chrmid, text = names(chrmid), cex = 0.75)
  abline(v = chrsum, col = "gray50")
  box()
  if( !is.null(title) ){
    title(main=title)
  }
  axis(side = 2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.8)

  # Color Scale
  par(mar = c(4,1.5,4,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp, axes = FALSE,
        xlab="",ylab="")
  axis(side = 4, las=2, cex.axis=0.8, lwd = 0, lwd.ticks = 1)

  layout(1)
}

#' Given DO mice qtl peak, gwas catelog, human phenotype keywords, plot the circos for integration
#' @param mm.peak        - a data.frame, generated by qtl2::find_peaks, lodcolumn must contain phenotype name
#' @param gwas.cat.file  - gwas catalog file with absolute path. Default:gwas_catalog_v1.0.2-associations_e100_r2020-12-15.tsv. This is the latest version
#' @param key1           - a character. Key word in human gwas. eg. dependence
#' @param key2           - a character. Key word in human gwas. eg. "anxiety|Anxiety". multiple words seperated by |
#' @param csv.file       - a character. Path and name for a csv, which contains mouse and human gene info
#' @param pdf.file       - a character. Path and name for a pdf, plotting circos
#' @return               - a csv contains mouse and human gene info; a plot of circos
#'
#'example:
#'load("~/projects/csna_workflow/qtl/circos.RData")
# mm.peak
# gwas.cat.file <- "~/projects/csna_workflow/gwas_catalog_v1.0.2-associations_e100_r2020-12-15.tsv"
# key1 = "dependence"
# key2 = "anxiety|Anxiety"
# csv.file = "~/projects/csna_workflow/out.csv"
# pdf.file = "~/projects/csna_workflow/out.pdf"
# circos_mouse_human(mm.peak, gwas.cat.file, key1, key2, csv.file, pdf.file)
circos_mouse_human <- function(mm.peak, gwas.cat.file, key1, key2, csv.file, pdf.file){
  #library
  library(circlize)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  library(biomaRt)
  library(org.Mm.eg.db)
  library(ComplexHeatmap)
  library(tidyverse)

  #human genome
  cytoband_hg = read.cytoband(species = "hg19")
  cytoband_df_hg = cytoband_hg$df
  cytoband_df_hg <- subset(cytoband_df_hg, V1!="chrY")
  windows <- floor(cytoband_hg$chr.len[-24] / 1000000)
  chroms <- names(windows)

  #gwas catelog
  gwas <- read_tsv(gwas.cat.file)

  # subset gwas based on key1 and key2
  key.gwas <- c(key1,key2) %>% map(function(x){
    gwas %>%
      rename_with(~(gsub("/", "_", .x, fixed = TRUE))) %>%
      rename_with(~(gsub("-", "_", .x, fixed = TRUE))) %>%
      filter(str_detect(DISEASE_TRAIT, x)) %>%
      filter(!is.na(CHR_ID)) %>%
      filter(!str_detect(CHR_ID, " x ")) %>%
      mutate(CHR_ID = case_when(str_detect(CHR_ID, ";") ~ word(CHR_ID, 1, sep = "\\;"),
                                TRUE                    ~  CHR_ID)) %>%
      mutate(chr = "chr") %>%
      unite(chr, c("chr", "CHR_ID"), sep = "") %>%
      mutate(end   = CHR_POS + 1,
             value = -log10(P_VALUE)) %>%
      dplyr::select(chr    = chr,
                    start  = CHR_POS,
                    end    = end,
                    marker = SNPS,
                    trait  = DISEASE_TRAIT,
                    value  = value,
                    gene   = MAPPED_GENE)
  })
  #for each key.gwas
  plotType <- list()
  for(k in 1:length(key.gwas)){
    myVect <- c()
    for(i in 1:length(chroms)) {
      print(chroms[i])
      toSUM <- which(key.gwas[[k]]$chr == chroms[i])
      print(toSUM)
      for(j in 1:windows[i]) {
        Upper <- j * 1000000
        Lower <- (j-1) * 1000000
        wSum <- ifelse(length(toSUM)!=0, length(which(key.gwas[[k]][toSUM,"start"] <= Upper & key.gwas[[k]][toSUM,"start"] > Lower)), 0)
        myVect <- c(myVect, wSum)
      }
    }
    intervalS <- (do.call(c, sapply(windows, function(x) seq(1:x)))) * 1000000
    plotType[[k]] <- data.frame(chr = (rep(names(windows), windows)),
                                start = intervalS - 1000000,
                                end=intervalS,
                                value = myVect,
                                stringsAsFactors = FALSE)
    # Due to plot area, limit maximum of number of SNPs to 20 in the plot.
    plotType[[k]][,4][plotType[[k]][,4] > 20] <- 20
  }

  # circlize plot----------------------------------------------------------------
  # Set up base of circle plot using tracks in circlize package:
  cytoband_hg = read.cytoband(species = "hg19")
  cytoband_df_hg = cytoband_hg$df
  cytoband_df_hg <- subset(cytoband_df_hg, V1!="chrY")

  cytoband_mm = read.cytoband(species = "mm10")
  cytoband_df_mm = cytoband_mm$df
  cytoband_df_mm <- cytoband_df_mm[nrow(cytoband_df_mm):1,]
  cytoband_df_mm <- subset(cytoband_df_mm, V1!="chrY")
  cytoband_df_mm$V1 <- paste0("mm10_", cytoband_df_mm$V1)

  xrange = c(cytoband_hg$chr.len[1:23], rev(cytoband_mm$chr.len)[-1])
  human_chr_index = 1:23
  mouse_chr_index = 24:43
  names(xrange)[mouse_chr_index] <-  paste0( "mm10_", names(xrange)[mouse_chr_index])

  sector.width = c(xrange[human_chr_index] / sum(xrange[human_chr_index]),
                   xrange[mouse_chr_index] / sum(xrange[mouse_chr_index]))

  Both_Genome <- rbind(cytoband_df_hg, cytoband_df_mm)
  Both_Genome$V1 <- factor(Both_Genome$V1, levels = unique(Both_Genome$V1))
  NAMES <- c(1:22, "X","X", 19:1)

  #Now read in clinical QTL
  clin_data <- mm.peak %>%
    transmute(chr       = chr,
              pos       = pos*1000000,
              CIStart   = ci_lo*1000000,
              CIEnd     = ci_hi*1000000,
              lod       = lod,
              lodcolumn = lodcolumn)

  #ChIPseeker package has nice function to get nearest gene. To use it for the closest gene
  #need a small region
  get_NearestGenes <- data.frame(clin_data) %>%
    mutate(start = pos - 1,
           end   = pos + 1) %>%
    dplyr::select(1,7,8,2:6) %>%
    mutate(chr = paste0("chr", chr)) # Need the chr prefix
  gData <- makeGRangesFromDataFrame(get_NearestGenes,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqnames.field = "chr",
                                    start.field="start",
                                    end.field="end")
  peakAnno <- annotatePeak(gData,
                           tssRegion=c(-1, 1),
                           TxDb=txdb,
                           annoDb="org.Mm.eg.db",
                           ignoreOverlap = FALSE,
                           overlap="all")
  toUse <- data.frame(peakAnno@anno)

  # Had to manually annotate a few, I just used UCSC genome browser and
  # for each QTL position looked for nearest gene that was not a predicted gene or
  # gene model because those will not have orthologs.
  # These were all the closest to the QTL and within the orginal QTL confidence interval.
  toUse <- toUse[,c(1, 6:9, 10, 19, 21)]

  ### Download ortholog information from biomaRt:
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  MMattributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene", "hsapiens_homolog_perc_id")
  orth.mouse = getBM(MMattributes, filters="with_hsapiens_homolog",values =TRUE, mart = mouse, bmHeader=FALSE)

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  HGattributes <- c("chromosome_name", "start_position", "end_position", "hgnc_symbol","ensembl_gene_id")
  human_convert <- getBM(HGattributes, "ensembl_gene_id", orth.mouse$hsapiens_homolog_ensembl_gene, mart= human)
  mouse_convert <- getBM(c("mgi_symbol","ensembl_gene_id"), "ensembl_gene_id", orth.mouse$ensembl_gene_id, mart= mouse)

  orth.1 <- merge(orth.mouse, mouse_convert, by="ensembl_gene_id", all=TRUE)
  orth.2 <- merge(orth.1, human_convert, by.x="hsapiens_homolog_ensembl_gene", by.y="ensembl_gene_id", all.x=TRUE)
  # will use orth.2

  ## Now merge back the ortholog info in.
  Hgenes <- orth.2[which(orth.2$mgi_symbol %in% toUse$SYMBOL),c(4:8, 3)]

  # Check any with no ortholog?
  setdiff(toUse$SYMBOL, Hgenes$mgi_symbol)
  subset(toUse, !(toUse$SYMBOL %in%Hgenes$mgi_symbol))

  # This one does have one that I found by searching manually,
  # Not sure why it didn't show up here. I've added it directly:
  # Hgenes <- rbind(Hgenes, data.frame(mgi_symbol=setdiff(toUse$SYMBOL, Hgenes$mgi_symbol),
  #                                    chromosome_name=c(22),
  #                                    start_position=c(22890123),
  #                                    end_position=c(22901768),
  #                                    hgnc_symbol=c("PRAME"),
  #                                    hsapiens_homolog_perc_id=c(NA)))

  Hgenes[,2] <- paste0("chr", Hgenes[,2])
  Hgenes <- Hgenes[,c(2:4,1,5,6)]
  names(Hgenes)[1:4] <- c("chr", "start", "end", "SYMBOL")

  toUse[,1]<-paste0("mm10_", toUse[,1])
  # merge orthologs and clinical qtl together
  ORTHO <- merge(Hgenes, toUse, by="SYMBOL", all=F)

  # Some mouse genes map to multiple orthologs in human
  # Remove duplicates
  # Choose based on percent homology or matching symbol name
  subset(ORTHO, ORTHO$SYMBOL %in% unique(ORTHO[duplicated(ORTHO[,7:9]),1]))
  ORTHO <- ORTHO[-c(7,8,36),]

  # The following is formatting for the plot and creating the supplement file:
  # We want the names to look nice for the plot:
  clin_data[,1] <- paste0("mm10_chr", clin_data[,1])
  colnames(clin_data) <- c("chr.mm", "qtl", "qtl.lowCI", "qtl.highCI",  "lod", "lodcolumn")
  clin_data$name.for.fig <- clin_data$lodcolumn

  # Additional manual formatting for the plot
  allHMLG_data <- ORTHO[order(ORTHO$seqnames, ORTHO$CIStart),]
  h_data <- allHMLG_data[,2:4]
  h_data[,2] <- round(as.numeric(h_data[,2]) / 1000000) #nearest Mbp
  h_data[,3] <- h_data[,2] + 1
  h_data[,2:3] <- h_data[,2:3] * 1000000
  colnames(h_data) <- colnames(ORTHO[,2:4])

  anno_data <- clin_data[,c(1,3,4)]
  colnames(anno_data) <- colnames(ORTHO[,2:4])
  anno_data <- rbind(anno_data, h_data)
  anno_data$value <- 1
  anno_data$value[1:nrow(clin_data)] <- clin_data[,7]
  anno_data$value[c(1+nrow(clin_data)):nrow(anno_data)] <- anno_data[c(1+nrow(clin_data)):nrow(anno_data),2] / 1000000
  anno_data$cols1 <- (c(anno_data$value[1:nrow(clin_data)], rep("black",nrow(anno_data) - nrow(clin_data))))
  anno_data <- anno_data[order(factor(anno_data[,1], levels = levels(Both_Genome$V1)), anno_data[,2], anno_data[,3]),]

  B_hg <- h_data[,1:3]
  B_mm <- allHMLG_data[,c(7,9,10)]

  ## This is to format the data as a supplementary file for the paper:
  toSupp <- allHMLG_data
  colnames(toSupp) <- c("mgi_symbol", "chr.hg", "start.hg", "end.hg", "hgnc_symbol", "percent_homology", "chr.mm", "qtl", "qtl.lowCI", "qtl.highCI", "lod", "lodcolumn", "distanceToTSS")
  toSupp <- merge(toSupp, clin_data, by=c("chr.mm", "qtl.lowCI", "qtl.highCI", "qtl"))
  toSupp[,4] <- round(as.numeric(toSupp[,4]) / 1000000)
  toSupp[,7] <- round(as.numeric(toSupp[,7]) / 1000000)
  toSupp[,6] <- paste0("hg19_", toSupp[,6])
  toSupp[,11] <- round(as.numeric(toSupp[,11]), 2)
  #better organizing/names for output
  toSupp <- toSupp[,c(16,1,4,11,5,9,6,7)]
  colnames(toSupp) <- c( "trait", "mm10_chr", "qtl_Mbp", "lod.score", "nearest_gene_mgi","syntenic_gene_hg","syntenic_hg19_chr", "syntenic_hg19_Mbp")
  clinData <- toSupp
  write.table(toSupp, file = csv.file, quote=F, row.names=F, sep=",")

  #circos plot
  set.seed(2020)
  colfunc <- colorRampPalette(c("paleturquoise1", "blue"))
  COL1 <- "#a8df49" #colfunc(24)[4]
  colfunc2 <- colorRampPalette(c("orange", "red"))
  COL2 <- colfunc2(21)[2]

  anno_data <- anno_data %>% mutate(value = case_when(
    str_detect(value, 'raw') ~ str_remove(value, '.raw'),
    TRUE ~ value
  ))

  pdf(pdf.file, width = 28, height = 26)
  circos.clear()
  circos.par(start.degree = 85, "gap.degree" = c(rep(1, 22), 10, rep(1, 19), 10),
             canvas.xlim=c(-.1,.1), canvas.ylim=c(-.95, .95), track.margin=c(.01, .01),
             cell.padding=c(.01, .01, .01, .01), track.height = 0.3, points.overflow.warning = FALSE)
  circos.initializeWithIdeogram(plotType=NULL, Both_Genome, sector.width = sector.width,
                                sort.chr = FALSE, track.height=.8)
  #first plotType
  circos.par(track.margin=c(.01, .01),
             cell.padding=c(.001, .001, .001, .001), track.height = .2)
  circos.genomicTrack(plotType[[1]], track.index=1, bg.border="white", ylim=c(-10,20),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, col = "brown3", type='l', lwd=4)
                      })
  circos.yaxis(side = "left", sector.index = "chr1", track.index=1, labels.cex=1.75, col="black", labels.col="black")
  #second plotType
  circos.genomicTrack(plotType[[2]], track.index=1, bg.border="white", ylim=c(-10,20),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, -value, col = "dodgerblue", type='l', lwd=3)
                      })
  circos.par(track.margin=c(.01, .01),
             cell.padding=c(.001, .001, .001, .001), track.height = .8)
  circos.genomicLabels(anno_data, labels.column = 4, cex=1.75, track.margin=c(0,0),
                       side="outside", padding=.1, col="black",
                       connection_height=convert_height(10, "mm"), line_lwd=4,
                       labels_height = (convert_height(1, "cm")))
  circos.par(track.margin=c(.01, .01),
             cell.padding=c(.001, .001, .001, .001), track.height = 0.1)

  COLORS <- c(rep(COL1, 23), rep(COL2, 20))
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = COLORS[CELL_META$sector.numeric.index], border="white")
    circos.text(mean(xlim), mean(ylim), NAMES[CELL_META$sector.numeric.index], col = "black", cex=1.75,
                facing="downward", niceFacing = TRUE)
  }, track.height = 0.25, bg.border = NA)

  circos.par("track.height" = 0.1, cell.padding = c(.1, .1, .1, .1))
  circos.genomicLink(B_mm, B_hg, col=gray.colors(nrow(B_mm), start = .1, end=.8, alpha=1))
  dev.off()
}

#correlation plot
my.plotcorr <- function (corr, outline = FALSE, col = "grey", upper.panel = c("ellipse", "number", "none"), lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"), digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...)
{
  # this is a modified version of the plotcorr function from the ellipse package
  # this prints numbers and ellipses on the same plot but upper.panel and lower.panel changes what is displayed
  # diag now specifies what to put in the diagonal (numbers, ellipses, nothing)
  # digits specifies the number of digits after the . to round to
  # unlike the original, this function will always print x_i by x_i correlation rather than being able to drop it
  # modified by Esteban Buz
  if (!require('ellipse', quietly = TRUE, character = TRUE)) {
    stop("Need the ellipse library")
  }
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ #diag behavior
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ #lower half of plot
      if (lower.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { #upper half of plot
      if (upper.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}

#' covert marker data.frame into pmap list required by qtl2
#' @param marker - marker data.frame containing chr, pos, snp. pos is marker position in Mbp
#' @return pmap - pmap list required by qtl2
#'
get_pmap <- function(marker){
  #group
  marker <- marker %>%
    group_by(chr) #group
  #list
  marker <- marker %>%
    group_map(~(
      .x %>%
        pull(pos) %>%
        setNames(., nm = .x$snp)
    )) %>%
    set_names(., nm = group_keys(marker)$chr) %>%
    .[c(1:19, "X")]
}

# retrieve.dat.measnum.ontology -------------------------------------------
library(jsonlite)
library(RCurl)
#outliers are outside of mean +/- 3SD
remove.outliers <- function(x, thres = 3, na.rm = TRUE, ...) {
  y <- x
  y[y <= mean(y,na.rm = na.rm) - thres * sd(y, na.rm = na.rm) | y >= mean(y,na.rm = na.rm) + thres * sd(y, na.rm = na.rm)] <- NA
  y
}
#rankz transform
rz.transform <- function(y){
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))
  return(rzT)
}

#' Given a ontology.id, retrieve DO csna dataset of measnum mapped to ontology.id,
#' @param ontology.id - a character. eg, "VT0010718"
#' @param geno.fam    - a data.frame, three column names must be FID, IID and sex
#' @param rm.outlier  - logical. TRUE/FALSE. whether to remove outlier or not
#' @param out.dir     - out directory
#' @return            - write out for male and female, respectively.
#'
retrieve.dat.measnum.ontology <- function(ontology.id, geno.fam, rm.outlier, out.dir){
  #retrieve data from api of MPD for ontology id
  api.url <- "http://bhmpd01lp.jax.org:82/api/pheno/measures_by_ontology/"
  url <- paste0(api.url, ontology.id)
  raw.dat <- tryCatch(RCurl::getURLContent(url),
                      HTTPError = function(e) {
                        return(paste("HTTP error: ", trimws(e$message)))
                      })
  #get mueasure id for the ontology.id
  dat <- jsonlite::fromJSON(raw.dat)
  #filter to ontology_measure_mappings, DO population w/par, CSNA03 on dat$measures
  df <- dat$measures %>%
    filter(measnum %in% dat$ontology_measure_mappings$measnum) %>%
    filter(paneldesc == "DO population w/par") %>%
    filter(projsym   == "CSNA03")
  #get DO population data set for each measnum
  df.measnum <- unique(df$measnum) %>%
    map(function(x){
      url <- paste0("http://bhmpd01lp.jax.org:82/api/pheno/animalvals/", x)
      raw.dat <- tryCatch(RCurl::getURLContent(url),
                          HTTPError = function(e) {
                            return(paste("HTTP error: ", trimws(e$message)))
                          })
      dat <- jsonlite::fromJSON(raw.dat)$animaldata %>%
        filter(strain == "J:DO")
      colnames(dat)[colnames(dat) == "value"] = as.character(unique(dat$varname))
      dplyr::select(dat, -c(measnum, varname, zscore))
    }) %>%
    reduce(full_join, by = c("animal_id", "animal_projid", "projsym", "sex", "stocknum", "strain", "strainid")) %>%
    mutate(IID = as.numeric(str_remove(animal_id, "_0")), .after = animal_id) %>%
    group_by(sex)
  #remove.outliers
  if(rm.outlier){
    df.measnum <- df.measnum %>%
      mutate(across(!c(animal_id, IID, animal_projid,
                       projsym, stocknum, strain, strainid),
                    remove.outliers))
  }
  #rankz transformation
  df.measnum <- df.measnum %>%
    mutate(across(!c(animal_id, IID, animal_projid,
                     projsym, stocknum, strain, strainid),
                  rz.transform)) %>%
    select(!c(animal_id, animal_projid,
              projsym, stocknum, strain, strainid)) %>%
    mutate(FID = IID, .before = IID) %>%
    right_join(geno.fam) #to those having genotype data in geno.fam

    #write out for male and female, respectively.
    #header
    write.table(df.measnum[1, -3][-1, ], file = paste0(out.dir, gsub(":", "", ontology.id), ".header.pheno"),
              row.names = FALSE, col.names = TRUE, sep = " ")
    write.table(df.measnum[df.measnum$sex == "m", -3], file = paste0(out.dir, gsub(":", "", ontology.id), ".m.pheno"),
                row.names = FALSE, col.names = FALSE, sep = " ")
    write.table(df.measnum[df.measnum$sex == "f", -3], file = paste0(out.dir, gsub(":", "", ontology.id), ".f.pheno"),
                row.names = FALSE, col.names = FALSE, sep = " ")
}
