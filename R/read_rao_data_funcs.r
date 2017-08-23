######################## FUNCTIONS TO HELP READ RAO DATA ###############################
#' Read data from a directory
#'
#' Read data from the correct subdirectory of a given directory, depending on desired
#' resolution, mapq, etc.
#'
#' @param dir Data directory containing subdirectories with different resolutions.
#' Default: current working directory.
#' @param resolution Resolution to read in, e.g. 1Mb, 10kb. Default: 1Mb.
#' @param mapq MAPQ threshold to use. Only 0 and 30 are supported. Default: 30
#' @param chr Chromosomes to read in. Default: chr1.
#' @param seqinfo If you don't provide seqinfo, you will get warnings when
#' combining data from different chromosome. Default: NULL.
#' @param show_progress Whether to show progress bar for file reading. Default: TRUE.
#' @return An InteractionSet object containing data for the specified resolution
#' and chromosomes.
#'
#' @import InteractionSet
#'
#' @export
read_rao_huntley_data <- function(dir = ".", resolution = "1Mb", mapq = 30,
                                  chr = "chr1",read_expected = FALSE,
                                  seqinfo = NULL, show_progress = TRUE){
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
    df_list <- read_data_from_dir(dir = d, bin_size = bin_size,
                                  read_expected = read_expected,
                                  show_progress = show_progress)

    iset <- make_is_from_data(raw_df = df_list$raw, norm_df = df_list$norm,
                              exp_df = df_list$exp,
                              chr, bin_size, seqinfo = seqinfo)
    return(iset)
  })

  iset_all <- do.call("c", islist)
  return(iset_all)

}


#' Interchromosomal version of the main read function
#' Doesn't support expected reads!
#' @param dir Data directory containing subdirectories with different resolutions.
#' Default: current working directory.
#' @param resolution Resolution to read in, e.g. 1Mb, 10kb. Default: 1Mb.
#' @param mapq MAPQ threshold to use. Only 0 and 30 are supported. Default: 30
#' @param chr Chromosomes to read in. Default: c("chr1","chr2"). Only takes pairs of chromosomes, for more combinations
#' call the function with an apply-like approach.
#' @param show_progress Whether to show progress bar for file reading. Default: TRUE.
#' @return An InteractionSet object containing data for the specified resolution
#' and chromosomes.
#' @importFrom gtools mixedsort
read_rao_huntley_data_interchr <- function(dir = ".", resolution = "1Mb", mapq = 30,
                                           chr = c("chr1","chr2"), show_progress = TRUE){
  if ((length(chr) != 2) | (chr[1] == chr[2])) {stop("Please provide two different chromosomes in a vector!")}
  chr = gtools::mixedsort(chr) ## makes sure the lower chr number is first
  bin_size <- get_bin_size(resolution)
  resolution <- paste0(resolution, "_resolution_interchromosomal")

  if (mapq == 0){
    mapq <- "MAPQG0"
  } else if (mapq == 30){
    mapq <- "MAPQGE30"
  } else {
    stop("MAPQ should be one of 0 or 30")
  }

  message("Will read data for ", length(chr), " chromosomes.")

  ### Read interchromosome data and create InteractionSet
  d <- paste(dir, resolution, paste0(chr,collapse="_"), mapq, "/", sep = "/")
  df_list <- read_data_from_dir_interchr(dir = d, bin_size = bin_size,
                                         show_progress = show_progress)
  iset_inter <- make_is_from_data(raw_df = df_list$raw, norm_df = df_list$norm,
                                  chr = chr, bin_size = bin_size, seqinfo = seqinfo, interchr = TRUE)

  return(iset_inter)

}

#' Get bin size from string resolution
#'
#' Get bin size in numeric base pairs from string resolution. Mb, kb, or bp
#' are supported.
#'
#' @param res Resolution as a string
#' @return numeric bin size
#'
#' @importFrom stringr str_extract str_to_lower
#' @export
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

#' Read data from given directory into a list of data frames
#'
#' Read data from a given directory, for a single chromosome and resolution,
#' into a list of data frames.
#'
#' @param dir Directory with data.
#' @param bin_size Numeric bin size
#' @param read_expected Whether to read in vectors of expected interactions at
#' different distances. Default: FALSE.
#' @param show_progress Whether to show progress bar for file reading. Default: TRUE.
#' @return A list of data frames, with raw data, normalisation vectors, and
#' expected data.
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr bind_cols left_join rename "%>%" data_frame
#'
#' @export
read_data_from_dir <- function(dir, bin_size, read_expected = FALSE, show_progress = TRUE){
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
  raw_df <- read_tsv(rawfile, col_names = c("start1", "start2", "obs"),
                     col_types = "iid", progress = show_progress)

  #read normalisation vectors
  norm_list <- lapply(norm_files, function(f){
    col_name <- strsplit(basename(f), "\\.")[[1]][2]
    read_tsv(f, col_names = col_name, col_types = "d", progress = show_progress)
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
      read_tsv(f, col_names = col_name, col_types = "d", progress = show_progress)
    })
    exp_df <- bind_cols(exp_list)
    exp_df$distance <- (1:nrow(exp_df)-1)*bin_size
  } else {
    exp_df <- NULL
  }

  return(list(raw = raw_df, norm = norm_df, exp = exp_df))
}

#' Read data from given directory into a list of data frames
#'
#' Read data from a given directory, for two different chromosomes and resolution,
#' into a list of data frames.
#'
#' @param dir Directory with data.
#' @param bin_size Numeric bin size
#' @param show_progress Whether to show progress bar for file reading. Default: TRUE.
#' @return A list of data frames, with raw data, normalisation vectors.
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr bind_cols left_join rename "%>%" data_frame
#' @importFrom gtools mixedsort
#' @export
read_data_from_dir_interchr <- function(dir, bin_size, show_progress = TRUE){

  message("Reading data from: ", dir)
  lf <- list.files(dir)

  obs_idx <- grepl("RAWobserved", lf)
  norm_idx <- grepl("norm", lf)

  if (sum(obs_idx) > 1){stop("More than one 'RAWobserved' file found in this directory!")}
  if (sum(obs_idx) == 0){stop("No 'RAWobserved' file found in this directory!")}

  rawfile <- paste0(dir, lf[obs_idx])
  norm_files <- paste0(dir, lf[norm_idx])
  norm_files <- gtools::mixedsort(norm_files) ## sort so chrA is always above chrB
  ## if interchromosomal normalization vectors available, use them
  norm_ixs <- grep("norm",norm_files)
  if(length(norm_ixs) > 0){
    norm_files <- norm_files[norm_ixs]
    message("Normalization files found!")
  } else(stop("No normalization files found, stopping."))

  #read raw data
  raw_df <- read_tsv(rawfile, col_names = c("start1", "start2", "obs"),
                     col_types = "iid", progress = show_progress)
  #read normalisation vectors
  norm_list <- lapply(norm_files, function(f){
    col_name <- strsplit(basename(f), "\\.")[[1]][2]
    read_tsv(f, col_names = col_name, col_types = "d", progress = show_progress)
  })
  ## stack norm vectors on top of each other and use later offset indexing to track which chr it is
  INTERKRix = grep("INTERKRnorm",norm_files,fixed=T)
  INTERVCix = grep("INTERVCnorm",norm_files)
  GWKRix = grep("GWKRnorm",norm_files)
  GWVCix = grep("GWVCnorm",norm_files)
  KRix = grep(".KRnorm",norm_files,fixed=T)
  SQRTVCix = grep("SQRTVCnorm",norm_files)
  VCix = grep(".VCnorm",norm_files,fixed=T)
  
  INTERKR = rbind(norm_list[[INTERKRix[1]]],norm_list[[INTERKRix[2]]])
  INTERVC = rbind(norm_list[[INTERVCix[1]]],norm_list[[INTERVCix[2]]])
  GWKR = rbind(norm_list[[GWKRix[1]]],norm_list[[GWKRix[2]]])
  GWVC = rbind(norm_list[[GWVCix[1]]],norm_list[[GWVCix[2]]])
  KR = rbind(norm_list[[KRix[1]]],norm_list[[KRix[2]]])
  SQRTVC = rbind(norm_list[[SQRTVCix[1]]],norm_list[[SQRTVCix[2]]])
  VC = rbind(norm_list[[VCix[1]]],norm_list[[VCix[2]]])
  
  #norm_df = bind_cols(INTERKR,INTERVC,GWKR,GWVC,KR,SQRTVC,VC)
  norm_df = bind_cols(INTERKR,INTERVC,GWKR,GWVC)
  
  # connect start positions with bin numbers to link raw / norm data for chrA and chrB
  bin_df_A <- data_frame(bin_start = seq(from = 0, by = bin_size, length.out = nrow(norm_list[[INTERKRix[1]]])),
                         bin_num = seq_along(bin_start))
  bin_df_B <- data_frame(bin_start = seq(from = 0, by = bin_size, length.out = nrow(norm_list[[INTERKRix[2]]])),
                         bin_num = seq_along(bin_start))
  bin_df <- rbind(bin_df_A,bin_df_B)
  raw_df <- raw_df %>%
    left_join(bin_df_A, by = c("start1" = "bin_start")) %>%
    dplyr::rename(bin_num1 = bin_num) %>%
    left_join(bin_df_B, by = c("start2" = "bin_start")) %>%
    dplyr::rename(bin_num2 = bin_num)

  chr_col <- data.frame("Source"=c(rep("chrA",nrow(bin_df_A)),rep("chrB",nrow(bin_df_B)))) #to keep track of which bins belong where
  norm_df <- bind_cols(bin_df, norm_df, chr_col)

  return(list(raw = raw_df, norm = norm_df))
}

#' Make InteractionSet object from a list of data frames
#'
#' Makes an InteractionSet object from a list of data frames containing raw
#' counts between bins, normalisation vectors, and (optional) expected data.
#'
#' @param raw_df data.frame of raw data with bin start positions, sizes, bin numbers, and
#' observed counts.
#' @param norm_df data.frame of normalisation vectors for each bin.
#' @param exp_df data.frame of expected vectors, can be NULL to not calculate
#' expected interactions. Default: NULL.
#' @param chr Chromosome of data.
#' @param bin_size Numeric bin size in bp.
#' @param seqinfo If you don't provide seqinfo, you will get warnings when
#' combining data from different chromosome. Default: NULL.
#' @return An InteractionSet object with assays containing observed counts and
#' optionally expected counts. Normalisation vectors are returned as metadata
#' columns for the `regions()` of the object.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import InteractionSet
#'
#' @export
make_is_from_data <- function(raw_df, norm_df, exp_df = NULL, chr, bin_size, seqinfo = NULL, interchr = FALSE){
  #get metadata columns for regions
  norm_cols <- names(norm_df)[-grep("bin_start|bin_num|Source", names(norm_df))]
  #make regions (bins) GRanges
  chr_A_end <- 0 ## needed for offsetting the indices in GInteractions call
  if(interchr){
    chr_A_end <- tail(which(norm_df$Source == "chrA"),n=1) ## this is where bins reset
    regions <- GRanges(seqnames = c(rep(chr[1],chr_A_end),rep(chr[2] , nrow(norm_df) - chr_A_end)),
                       ranges = IRanges(start = norm_df$bin_start, width = bin_size),
                       mcols = as.data.frame(norm_df[,norm_cols]))
    names(mcols(regions)) <- gsub("mcols.", "", names(mcols(regions)))
  }
  else{
    regions <- GRanges(seqnames = chr,
                       ranges = IRanges(start = norm_df$bin_start, width = bin_size),
                       mcols = as.data.frame(norm_df[,norm_cols]),
                       seqinfo = seqinfo)
    names(mcols(regions)) <- gsub("mcols.", "", names(mcols(regions)))
  }


  #make interactions and InteractionSet, chr_A_end = 0 for intrachromosomal ints

  ints <- GInteractions(anchor1 = raw_df$bin_num1, anchor2 = raw_df$bin_num2+chr_A_end,
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

#' Normalise observed data by normalisation vectors
#'
#' Normalise observed assays by multiplying them by the appropriate values for
#' the interacting bins. Result = observed /(anchor 1 value * anchor 2 value)
#'
#' @param iset InteractionSet object with observed interactions as an assay and
#' normalisation vectors as metadata columns of `regions()`.
#' @param obs Name of assay with observed data. Default: 'RAWobserved'.
#' @param norm names of columns in `mcols(regions(iset))` to use for normalisation,
#' or NULL to use all columns matching 'norm'. Default: NULL.
#' @return An InteractionSet object with additional assays holding normalised data.
#'
#' @export
normalise_hic <- function(iset, obs = "RAWobserved", norm = NULL){
  obs <- assays(iset)[[obs]]
  if (is.null(obs)){ stop("No observed data to normalise!")}

  if (is.null(norm)){
    vecs <- names(mcols(regions(iset)))
    vecs <- vecs[grepl("norm", vecs)]
  } else {
    vecs <- norm
  }

  if (length(vecs)==0) { stop("No normalisation vectors available!")}
  message("Using normalisation vectors: ", paste(vecs, collapse = ", "))

  assay_list <- lapply(vecs, function(v){
    #message("Normalising using ", v)
    vec <- mcols(regions(iset))[[v]]

    norm_assay <- obs / (vec[anchors(iset, type = "first", id = TRUE)] *
                           vec[anchors(iset, type = "second", id = TRUE)])
  })
  names(assay_list) <- vecs
  assays(iset) <- c(assays(iset), assay_list)
  return(iset)
}
