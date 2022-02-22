ts <- read.csv(snakemake@input[[1]], header=T, sep = ',')
ts$p.adjust <-p.adjust(ts$p.val,method="fdr", n=length(ts$p.val))
write.csv(ts, file=snakemake@output[[1]], row.names = FALSE, quote=FALSE)

