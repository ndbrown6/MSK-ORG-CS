#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
suppressPackageStartupMessages(library("gdata"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("superheat"))
suppressPackageStartupMessages(library("ggsignif"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggforce"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("preseqR"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("fuzzyjoin"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("ggnomics"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("dichromat"))
suppressPackageStartupMessages(library("ggpubr"))

url_manifest <- "../data/manifest.txt"
url_ihc_her2 <- "../data/IHC_Her2.txt"
url_ihc_trop2 <- "../data/IHC_Trop2.txt"
exclude <- "cs42"

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

'scale_range' <- function(x, center = TRUE, scale = TRUE) {
	x = as.vector(scale(x, center = center, scale = scale))
	return(invisible(x))
}
