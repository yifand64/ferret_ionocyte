
library(ggplot2)
extract.field=function(string,field=1,delim="_") return(strsplit(string,delim)[[1]][field])


args <- commandArgs(trailingOnly = TRUE)
input = args[1]

jobs <- read.csv(input, sep="")
colnames(jobs) = c("JOBID", "USER", "STAT", "QUEUE", "FROM_HOST", "EXEC_HOST", "JOB_NAME", "SUBMIT_TIME", "C2")
jobs["JOB"]=rownames(jobs)
colnames(jobs) = c("USER", "STAT", "QUEUE", "FROM_HOST", "EXEC_HOST", "JOB_NAME", "SUBMIT_MOTNTH", "SUBMIT_DAY","SUBMIT_TIME", "JOBID")
jobs["COMPUTE_CORES"] =  unlist(lapply(as.character(jobs$EXEC_HOST),extract.field,1, "*"))
jobs$COMPUTE_CORES = gsub("v", 1, jobs$COMPUTE_CORES)
jobs$COMPUTE_CORES = gsub("d", 1, jobs$COMPUTE_CORES)
jobs$COMPUTE_CORES = gsub("o", 1, jobs$COMPUTE_CORES)
jobs$COMPUTE_CORES = gsub("n", 1, jobs$COMPUTE_CORES)
jobs$COMPUTE_CORES = gsub("a", 1, jobs$COMPUTE_CORES)
running = jobs[jobs$STAT == "RUN", ]
running = running[, c("USER", "COMPUTE_CORES")]


running$COMPUTE_CORES = as.numeric(running$COMPUTE_CORES)
df <- lapply(running, function(x) type.convert(as.character(x)))
by.user = aggregate(. ~ USER, df, sum)
print(by.user)
pdf("queue.usage.pdf");ggplot(by.user, aes(x=USER, y=COMPUTE_CORES, fill=USER)) + geom_bar(stat="identity", position="stack") + theme_bw() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE) + ggtitle(paste("Regevlab queue usage", Sys.time(), sep="\n")) ; dev.off()
