library(InteractionSet)
library(dplyr)
library(stringr)
library(readr)

read_rao_huntley_data <- function(dir, resolution, mapq, chr, inter_chr = FALSE, 
                          read_expected = FALSE, seqinfo = NULL){
  if (length(dir) != 1) {stop("Please provide a single data directory.")}
  if (length(resolution) > 1) {stop("Can only read in a single resolution at a time.")}
  if (length(mapq) > 1) {stop("Please choose a single mapq value.")}
  
  bin_size <- get_bin_size(resolution)
  resolution <- paste0(resolution, "_resolution_intrachromosomal")
  
  if (mapq == 0){
    mapq <- "MAPQG0"
  } else if (mapq == 30){
    mapq <- "MAPQGE30"
  } else {
    stop("MAPQ should be one of 0 or 30")
  }
  
  message("Will read data for ", length(chr), " chromosomes.")
  
 islist <- lapply(chr, function(chr){
    d <- paste(dir, resolution, chr, mapq, "/", sep = "/")
    df_list <- read_data_from_dir(dir = d, bin_size = bin_size, read_expected = read_expected)
    
    iset <- make_gi_from_data(raw_df = df_list$raw, norm_df = df_list$norm, 
                              exp_df = df_list$exp,
                              chr, bin_size, seqinfo = seqinfo)
    return(iset)
  })
 
 if (inter_chr){
   
 }
 
 iset_all <- do.call("c", islist)
 return(iset_all)
 
}

get_bin_size <- function(res){
  pre <- as.numeric(str_extract(res, "[0-9]+"))
  multiplier <- str_to_lower(str_extract(res, "[aA-zZ]+"))
  lookup <- list(mb = 10^6, kb = 1000, bp = 1)
  if(!multiplier %in% names(lookup)) {
    stop("Unknown multiplier for bin size, should be: ",
         paste(names(lookup), collapse = ","), ". Found: ", multiplier)
    }
  multiplier <- lookup[[multiplier]]
  bin_size <- multiplier * pre
  return(bin_size)
}

read_data_from_dir <- function(dir, bin_size, read_expected = FALSE){
  message("Reading data from: ", dir)
  lf <- list.files(dir)

  obs_idx <- grepl("RAWobserved", lf)
  norm_idx <- grepl("norm", lf)
  exp_idx <- grepl("expected", lf)
  
  if (sum(obs_idx) > 1){stop("More than one 'RAWobserved' file found in this directory!")}
  if (sum(obs_idx) == 0){stop("No 'RAWobserved' file found in this directory!")}
  
  rawfile <- paste0(dir, lf[obs_idx])
  norm_files <- paste0(dir, lf[norm_idx])
  exp_files <- paste0(dir, lf[exp_idx])
  
  #read raw data
  raw_df <- read_tsv(rawfile, col_names = c("start1", "start2", "obs"), col_types = "iid")
  
  #read normalisation vectors
  norm_list <- lapply(norm_files, function(f){
    col_name <- strsplit(basename(f), "\\.")[[1]][2]
    read_tsv(f, col_names = col_name, col_types = "d")
  })
  norm_df <- bind_cols(norm_list)
  
  # connect start positions with bin numbers to link raw / norm data
  bin_df <- data_frame(bin_start = seq(from = 0, by = bin_size, length.out = nrow(norm_df)),
                       bin_num = seq_along(bin_start))
  
  raw_df <- raw_df %>%
    left_join(bin_df, by = c("start1" = "bin_start")) %>%
    dplyr::rename(bin_num1 = bin_num) %>%
    left_join(bin_df, by = c("start2" = "bin_start")) %>%
    dplyr::rename(bin_num2 = bin_num)
  
  norm_df <- bind_cols(bin_df, norm_df)
  
  if (read_expected){
    exp_list <- lapply(exp_files, function(f){
      col_name <- strsplit(basename(f), "\\.")[[1]][2]
      read_tsv(f, col_names = col_name, col_types = "d")
    })
    exp_df <- bind_cols(exp_list)
    exp_df$distance <- (1:nrow(exp_df)-1)*bin_size  
  } else {
    exp_df <- NULL
  }
  
  return(list(raw = raw_df, norm = norm_df, exp = exp_df))
}

make_gi_from_data <- function(raw_df, norm_df, exp_df = NULL, chr, bin_size, seqinfo = NULL){
  #get metadata columns for regions
  norm_cols <- names(norm_df)[-match("bin_start", names(norm_df))]
  #make regions (bins) GRanges
  regions <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = norm_df$bin_start, width = bin_size),
                    mcols = as.data.frame(norm_df[,norm_cols]),
                    seqinfo = seqinfo)
  names(mcols(regions)) <- gsub("mcols.", "", names(mcols(regions)))
  
  #make interactions and InteractionSet
  ints <- GInteractions(anchor1 = raw_df$bin_num1, anchor2 = raw_df$bin_num2,
                regions = regions)
  iset <- InteractionSet(list(RAWobserved = as.matrix(raw_df$obs)), ints)
  
  if (!is.null(exp_df)){
    exp_vecs <- names(exp_df)
    exp_vecs <- exp_vecs[-grep("distance", exp_vecs)]
    if (length(exp_vecs) == 0){
      warning("No expected vectors found, will skip calculation of expected matrices.")
    } else{
      dists <- pairdist(ints)
      idx <- match(dists, exp_df$distance)
      
      exp_list <- lapply(exp_vecs, function(exp){
        message("Calculating expected values for ", exp)
        as.matrix(exp_df[[exp]][idx])
      })
      names(exp_list) <- exp_vecs
      assays(iset) <- c(assays(iset), exp_list)
    }
  }
  return(iset)
}

normalise_hic <- function(iset){
  obs <- assays(iset)[["RAWobserved"]]
  if (is.null(obs)){ stop("No observed data to normalise!")}
  
  vecs <- names(mcols(regions(iset)))
  vecs <- vecs[grepl("norm", vecs)]
  if (length(vecs)==0) { stop("No normalisation vectors available!")}
  message("Available normalisation vectors: ", paste(vecs, collapse = ", "))
  
  
  assay_list <- lapply(vecs, function(v){
    message("Normalising using ", v)
    vec <- mcols(regions(iset))[[v]]
    norm_assay <- obs * vec[anchors(iset, type = "first", id = TRUE)] * 
      vec[anchors(iset, type = "second", id = TRUE)] 
  })
  names(assay_list) <- vecs
  assays(iset) <- c(assays(iset), assay_list)
  return(iset)
}

dir <- "~/Projects/genome_folding/rao_huntley_data/"
resolution <- "100kb"
mapq <- 30
chr <- "chr19"
chrs <- chr
read_expected <- TRUE

chrs <- paste0("chr", c(1:19, "X"))

iset <- read_rao_data("~/Projects/genome_folding/rao_huntley_data/", 
                      resolution = "100kb", 
                      mapq = 30, chr = chrs, 
                      read_expected = TRUE)

iset <- normalise_hic(iset)


