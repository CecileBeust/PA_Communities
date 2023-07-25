# R code to read gmt files, export gene sets and associated terms to csv files
# code from: https://bioinformatics.stackexchange.com/questions/5400/how-to-convert-data-in-gmt-format-to-dataframe

library(GSA)

# REACTOME GMT file
data <- GSA.read.gmt("../../_00_data/gmt_files/hsapiens.REAC.name.gmt")
len_vec = c() # create vector that will contain the length of genes at each position
len_vec[1] = 3
for(i in 1:length(data$genesets)){len_vec[i] <- c(length(data$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(data$geneset.names, len_vec),use.names = FALSE) # Now create a vector for all the pathways in the data 
desired_df <- as.data.frame(cbind(pathway_vec,unlist(data$genesets,use.names = FALSE))) # This gives your desired dataframe
colnames(desired_df) <- c("Pathway term", "Gene")
write.csv(desired_df, "REAC_genes.csv", row.names = FALSE)

# GO BP GMT file
data <- GSA.read.gmt("../../_00_data/gmt_files/hsapiens.GO:BP.name.gmt")
len_vec = c() # create vector that will contain the length of genes at each position
len_vec[1] = 3
for(i in 1:length(data$genesets)){len_vec[i] <- c(length(data$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(data$geneset.names, len_vec),use.names = FALSE) # Now create a vector for all the pathways in the data 
desired_df <- as.data.frame(cbind(pathway_vec,unlist(data$genesets,use.names = FALSE))) # This gives your desired dataframe
colnames(desired_df) <- c("GO:BP term", "Gene")
write.csv(desired_df, "GOBP_genes.csv", row.names = FALSE)


# GO CC GMT file
data <- GSA.read.gmt("../../_00_data/gmt_files/hsapiens.GO:CC.name.gmt")
len_vec = c() # create vector that will contain the length of genes at each position
len_vec[1] = 3
for(i in 1:length(data$genesets)){len_vec[i] <- c(length(data$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(data$geneset.names, len_vec),use.names = FALSE) # Now create a vector for all the pathways in the data 
desired_df <- as.data.frame(cbind(pathway_vec,unlist(data$genesets,use.names = FALSE))) # This gives your desired dataframe
colnames(desired_df) <- c("GO:CC term", "Gene")
write.csv(desired_df, "GOCC_genes.csv", row.names = FALSE)