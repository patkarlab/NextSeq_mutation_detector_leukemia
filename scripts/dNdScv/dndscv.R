library("dndscv")
library(optparse)

# Define options
option_list <- list(
  make_option(c("-f", "--infile"), help = "file with variants"),
  make_option(c("-g", "--genelist"), help = "list of genes"),
  make_option(c("-o", "--geneout"), help = "genewise output file"),
  make_option(c("-v", "--variantout"), help = "variant wise output file")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)

# Access the arguments
# cat("Infile:", opt$infile, "\n")
# cat("Genelist:", opt$genelist, "\n")

infile <- read.table(opt$infile, header = TRUE)
GeneList <- read.table(opt$genelist)

dnds = dndscv(infile, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

write.table(dnds$sel_cv, opt$geneout, row.names = FALSE, sep="\t")
write.table(dnds$annotmuts, opt$variantout, row.names = FALSE, sep="\t")