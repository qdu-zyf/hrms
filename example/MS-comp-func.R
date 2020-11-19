library(optparse)
library(hrms)

args <- commandArgs(trailingOnly=TRUE)

option_list <- list(  
  make_option(c("-i", "--abd_table"), type="character", help="Input abundance matrix table [Required]"),
  make_option(c("-o", "--distfile"), type="character", default='default', help="Output distance matrix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$abd_table)) stop('Please input a abundance matrix table')


abds <- read.table(opts$abd_table, header = T, row.names = 1)
abds <- as.matrix(abds)
dis <- CompFunc(abds)
write.table(dis, opts$distfile, col.names = NA, sep = '\t', quote = FALSE)

