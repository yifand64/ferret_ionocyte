# Visualise the distribution of various QC metrics (gene count, rRNA content, etc)
# across a set of samples. this is mostly barplots.


source("util.r")
#script("multiplot")
source("qc.r")
source("samples.r")

# cat(sprintf("Sourcing scripts..\n"))
# scripts = c("samples.r",
# 			"multiplot.r",
# 			"qc.r")
# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# this.script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.dir <- dirname(this.script.name)

# for (i in 1:length(scripts)) {
#   	script.name = scripts[[i]]
#   	other.name <- paste(sep="/", script.dir, script.name)
#     print(paste("Sourcing", other.name, "from", this.script.name))
#   	source(other.name)
    
# }


info("Generating summary qc plots..")

args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
qc_data_file = args[1]
if(n_args > 1){
output = args[2]
}else{output="summary.pdf"}
cat(sprintf("\"Loading qc data from  --> %s\"\n", qc_data_file))


# find sample names:
data <- read.delim(qc_data_file)
found = colnames(data)
print("Found columns: ")
print(found)

visualise.qc.summary(data)




