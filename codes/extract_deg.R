extract <- function(gene,ext,list) {
  setwd("/Users/rishav_raj/Desktop/bioinformatics/geo/paper2/")
  #dir = "/Users/rishav_raj/Desktop/bioinformatics/geo/paper2/"+gene+".csv"
  geodata <- read.csv(gene,header=TRUE)
  dim(geodata)
  deg <- geodata[which(geodata$P.Val<0.05),]
  deg <- deg[,c("ID","P.Value","t","logFC","Gene.symbol","Gene.title","Gene.ID")]
  dim(deg)
  
  
  
  #deg<- deg[2]
  #saving the dataframe in the local
  #setwd("/Users/rishav_raj/Desktop/bioinformatics/geo/paper2/DEG")
  
  # Ignore row names/numbers
  #output <- "de_genes"+genes+".csv"
  write.table(deg,ext,sep=",", row.names=FALSE)
  
}


extract("GSE58545.csv","58545deg.csv")

