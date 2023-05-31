

source("util.R")
args <- commandArgs(trailingOnly = TRUE)
target = args[1]
info(sprintf("Running fix QC on %s", target))
qc = read.delim(target)
ncol.orig = ncol(qc)
qc = qc[ , !(colnames(qc) %in% c("Sample.1","Sample.2","Sample.3","Sample.4","Unnamed..56"))]

info(sprintf("Removed %s bad columns", ncol.orig-ncol(qc)))

info(sprintf("Adding dummy 'Mouse' and 'Plate' columns"))

qc$Mouse = "dummy"
qc$Plate = "dummy"

target = gsub(".txt", "", target, fixed=T)
output = paste(target, "fixed.txt", sep="_")
info(sprintf("Writing output to %s", output))
write.table(qc, file=output, sep="\t", quote=F, row.names=F)
