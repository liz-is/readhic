# Read HiC contact matrices and normalisation vectors

Functions to read and process HiC contact matrices and normalisation vectors 
from Rao, Huntley et al 2014, as provided [on GEO](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525). 

You can specify the location of the data, the resolution to read in (as a string), 
and the chromosomes to read in. You can also specify whether you wish to read in data corresponding 
to expected interactions (if available). Data is returned as an InteractionSet object, with 
observed and expected data stored as `assays`. See the
[InteractionSet](http://www.bioconductor.org/packages/release/bioc/html/InteractionSet.html) 
and [GenomicInteractions](http://www.bioconductor.org/packages/release/bioc/html/GenomicInteractions.html) 
packages for more information and useful functions for manipulating these objects.

## Example usage

```
iset <- read_rao_huntley_data(dir = ".", resolution = "100kb",
  mapq = 30, chr = c("chr18", "chr19"), read_expected = TRUE)
# Will read data for 2 chromosomes.
# Reading data from: ./100kb_resolution_intrachromosomal/chr18/MAPQGE30//
# Calculating expected values for KRexpected
# Calculating expected values for RAWexpected
# Calculating expected values for SQRTVCexpected
# Calculating expected values for VCexpected
# Reading data from: ./100kb_resolution_intrachromosomal/chr19/MAPQGE30//
# Calculating expected values for KRexpected
# Calculating expected values for RAWexpected
# Calculating expected values for SQRTVCexpected
# Calculating expected values for VCexpected
# Warning message:
# In .Seqinfo.mergexy(x, y) :
#   The 2 combined objects have no sequence levels in common. (Use
#   suppressWarnings() to suppress this warning.)

iset
# class: InteractionSet 
# dim: 499539 1 
# metadata(0):
# assays(5): RAWobserved KRexpected RAWexpected SQRTVCexpected VCexpected
# rownames: NULL
# rowData names(0):
# colnames: NULL
# colData names(0):
# type: GInteractions
# regions: 1524

iset <- normalise_hic(iset)
# Using normalisation vectors: KRnorm, SQRTVCnorm, VCnorm

iset
# class: InteractionSet 
# dim: 499539 1 
# metadata(0):
# assays(8): RAWobserved KRexpected ... SQRTVCnorm VCnorm
# rownames: NULL
# rowData names(0):
# colnames: NULL
# colData names(0):
# type: GInteractions
# regions: 1524
```
