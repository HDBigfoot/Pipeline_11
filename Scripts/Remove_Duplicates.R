args <- commandArgs(trailingOnly = TRUE)

csv_file <- args[1]

dupedSnpDistance <- read.csv(csv_file, header = FALSE)

length <- nrow(dupedSnpDistance)
colnames(dupedSnpDistance) <- c("Sample1", "Sample2", "SNP.Distance")

Sample1 <- c()
Sample2 <- c()
SNP.Distance <- c()
iteration <- 1
stop_index <- length(unique(dupedSnpDistance$Sample1))
sample_number <- length(unique(dupedSnpDistance$Sample1))
current_index <- 1
while(current_index < nrow(dupedSnpDistance)){
  current_index <- current_index + iteration
  while(current_index <= stop_index){
    Sample1 <- append(Sample1, dupedSnpDistance$Sample1[current_index])
    Sample2 <- append(Sample2, dupedSnpDistance$Sample2[current_index])
    SNP.Distance <- append(SNP.Distance, dupedSnpDistance$SNP.Distance[current_index])
    current_index <- current_index + 1
  }
  iteration <- iteration + 1
  stop_index <- stop_index + sample_number
}

SnpDistanceDedup <- data.frame(Sample1, Sample2, SNP.Distance)

write.csv(SnpDistanceDedup, row.names = FALSE)

