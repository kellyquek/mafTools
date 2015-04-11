#' Basic statistics about MAF file.
#' 
#' Takes MAF file as input and gives intial impression of the data, such as number of mutations per sample, variant types, etc. Also plots barplot and boxplot of variant distributions.
#' 
#' @param maf_file input file in MAF format.
#' @param removeSilent removeSilent logical. Whether to discard silent mutations ("Silent","Intron","RNA","3'UTR"). Default is TRUE.
#' @import ggplot2
#' @import plyr
#' @import reshape
#' @return returns a list of data frames with number of mutations per sample, mutations classified according to variant type and variant classification.
#' @export

#' 
maf_stats = function(maf_file,removeSilent = TRUE){
  
  require(package = "RColorBrewer",quietly = T)
  require(package = "gridExtra",quietly = T)
  
  #Read file
  tot.muts = read.delim(file = maf_file,header = T,sep = "\t",stringsAsFactors = F,comment.char = "#")
  
  #Remove silent mutations
  if(removeSilent){
    tot.muts = tot.muts[!tot.muts$Variant_Classification %in% c("Silent","Intron","RNA","3'UTR"),]  
  }
  
  #convert columns into factors
  tot.muts$Variant_Type = as.factor(as.character(tot.muts$Variant_Type))
  tot.muts$Variant_Classification = as.factor(as.character(tot.muts$Variant_Classification))

  #split maf file into according to tumor sample barcode
  pid.split = split(tot.muts,f = as.factor(as.character(tot.muts$Tumor_Sample_Barcode))) 
  variants.df = data.frame(variants = sapply(pid.split,nrow))
  
  #create dataframes to store output
  variant.type.df = data.frame()
  variant.classification.df = data.frame()
  
  #get variant classes
  variant.classes = levels(tot.muts$Variant_Classification)
  
  #Make colors for variants classes.
  if(length(variant.classes) <= 9 ){
    type_col = structure(brewer.pal(length(variant.classes),name = "Set1"), names = variant.classes)
  } else{
    type_col = structure(rainbow(n = length(variant.classes),alpha = 1), names = variant.classes)
  }
  
  
  
  for(i in 1:length(pid.split)){
    
    mut = pid.split[[i]] #each patient
    variant.type.df = rbind(variant.type.df,sapply(split(mut,as.factor(mut$Variant_Type)),nrow))
    rownames(variant.type.df)[nrow(variant.type.df)] = names(pid.split[i])
    variant.classification.df = rbind(variant.classification.df,sapply(split(mut,as.factor(mut$Variant_Classification)),nrow))
    rownames(variant.classification.df)[nrow(variant.classification.df)] = names(pid.split[i])
  }
  
  colnames(variant.type.df) = levels(tot.muts$Variant_Type)
  colnames(variant.classification.df) = levels(tot.muts$Variant_Classification)
  
  #Reorder data acoording to number of mutations.
  variant.classification.df = variant.classification.df[order(rowSums(variant.classification.df),decreasing = T),]
  variant.classification.df = variant.classification.df[,order(colSums(variant.classification.df),decreasing = T)]
  
  #Melt data for ggplot
  vc = melt(as.matrix(variant.classification.df))
  vc$X1 = factor(vc$X1,levels = rownames(variant.classification.df)) #reorder patient levels
  vc$X2 = factor(vc$X2,levels = colnames(variant.classification.df)) #reorder variant class levels
  p = ggplot()+geom_bar(data = vc, aes(x = X1,y = value,fill = X2),stat = "identity")+theme(legend.position = "none")
  p.bar = p + scale_fill_manual(values = type_col, name = "Variant Class: ") + xlab("") + ylab("Number of Mutations") + theme(axis.text.x = element_blank())
  
  p.box = ggplot(data = vc)+geom_boxplot(aes(x = X2,y = value,fill = X2))+theme_bw()+theme(legend.position = "bottom")+scale_fill_manual(values = type_col, name = "Variant Class: ") + xlab("") + ylab("Number of Mutations")+ theme(axis.text.x = element_blank())
  
  grid.arrange(p.bar,p.box,nrow = 2, ncol=1)
  
  #Return results as a list of data frames. Self explainatory. 
  return(list(variants.per.sample = variants.df,variant.type.summary = variant.type.df,variant.classification.summary = variant.classification.df))
}