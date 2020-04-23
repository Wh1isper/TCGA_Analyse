options(stringsAsFactors = F)
dir.create("SampleFiles")
filepath = dir (path = ".\\gdc_download",full.names = TRUE)
for (wd in filepath){
  files = dir(path = wd,pattern = 'gz$')
  fromfilepath = paste(wd,"\\",files,sep="")
  tofilepath = paste(".\\SampleFiles\\",files,sep="")
  file.copy(fromfilepath,tofilepath)
}

setwd(".\\SampleFiles")
countsFiles = dir(path = ".\\", pattern="gz$")
library(R.utils)
sapply(countsFiles, gunzip)

library(rjson)
metadata_json_File = fromJSON(file = '..\\metadata.cart.2020-04-13.json')
json_file_info = data.frame(filesName = c(), TCGA_Barcode = c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode = metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name = metadata_json_File[[i]][["file_name"]]
  json_file_info = rbind(json_file_info,data.frame(filesName = file_name, TCGA_Barcode = TCGA_Barcode))
}
rownames(json_file_info) = json_file_info[,1]
write.csv(json_file_info,file = "..\\json_file_info.csv")



filesName_To_TCGA_BarcodeFile = json_file_info[-1]
countsFileNames = dir(pattern = "counts$")

allSampleRawCounts = data.frame()
for (txtFile in countsFileNames){
  sampleCounts = read.table(txtFile,header = FALSE)
  rownames(sampleCounts) = sampleCounts[,1]
  sampleCounts = sampleCounts[-1]
  colnames(sampleCounts) = filesName_To_TCGA_BarcodeFile[paste(txtFile,".gz",sep = ""),]
  if (dim(allSampleRawCounts)[1]==0){
    allSampleRawCounts = sampleCounts
  }
  else{
    allSampleRawCounts = cbind(allSampleRawCounts,sampleCounts)
  }
}
write.csv(allSampleRawCounts,file = "..\\allSampleRawCounts.csv")
ensembl_id = substr(row.names(allSampleRawCounts),1,15)
rownames(allSampleRawCounts) = ensembl_id
write.csv(allSampleRawCounts,file = "..\\RawCounts.csv")



# ID转换
RawCounts = allSampleRawCounts
Ensembl_ID = data.frame(Ensembl_ID = row.names(RawCounts))
rownames(Ensembl_ID) = Ensembl_ID[,1]
RawCounts = cbind(Ensembl_ID,RawCounts)

get_map = function(input){
  if (is.character(input)){
    if(!file.exists(input)) stop("bad input file.")
    message("Treat input as file")
    input = data.table::fread(input,header = FALSE)
  } else {
    data.table::setDT(input)
  }
  
  input = input[input[[3]]=="gene",]
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  
  Ensembl_ID_To_Genename = data.frame(gene_id = gene_id, 
                                      gene_name = gene_name,
                                      stringsAsFactors = FALSE)
  return (Ensembl_ID_To_Genename)
}

Ensembl_ID_To_Genename = get_map("..\\gencode.v33lift37.annotation.gtf")
gtf_Ensembl_ID = substr(Ensembl_ID_To_Genename[,1],1,15)
Ensembl_ID_To_Genename = data.frame(gtf_Ensembl_ID, Ensembl_ID_To_Genename[,2])
colnames(Ensembl_ID_To_Genename) = c("Ensembl_ID","gene_id")
write.csv(Ensembl_ID_To_Genename,file = "..\\Ensembl_ID_To_Genename.csv")


mergeRawCount = merge(Ensembl_ID_To_Genename,RawCounts,by="Ensembl_ID")
mergeRawCount = mergeRawCount[order(mergeRawCount[,"gene_id"]),]
index = duplicated(mergeRawCount$gene_id)
mergeRawCount = mergeRawCount[!index,]
rownames(mergeRawCount) = mergeRawCount[,"gene_id"]
LUAD_Counts_expMatrix = mergeRawCount[,-c(1:2)]
write.csv(LUAD_Counts_expMatrix, file = "..\\LUAD_Counts_expMatrix.csv")
