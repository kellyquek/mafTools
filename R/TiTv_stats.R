#' Classify SNVs into Transitions and Tranversions.
#' 
#'  Takes MAF file or simple six column file with chr, start, end, ref_allele, alt_allele and sample id as input
#'  and classifies them into six classes of converions.
#'  
#'  @param maf_file input file either in MAF format or six column format.
#'  @param is_maf logical. Whether input is in MAF or not. Default is TRUE.
#'  @param plot logical. Plots distribution of various conversion events as a boxplot.
#'  @return returns a list of dataframes containing raw counts classified into six different classes and fraction of conversion for each sample.
#'  @export   

TiTv_stats <- function(maf_file,is_maf = TRUE,plot = TRUE){
  
  
  if(is_maf){
    tot.muts = read.delim(file = maf_file,header = T,sep = "\t",stringsAsFactors = F,comment.char = "#",)
    tot.muts = split(tot.muts,as.factor(as.character(tot.muts$Variant_Type)))$'SNP'
    tot.muts = tot.muts[,c("Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
  } else{
    #Read the file
    tot.muts = read.delim(maf_file,header=F,sep = "\t",stringsAsFactors = F,comment.char = "#")
    tot.muts = tot.muts[,1:6]
    tot.muts = tot.muts[tot.muts[,4] %in% c("A","T","G","C"),] #select only snv's
    tot.muts = tot.muts[tot.muts[,5] %in% c("A","T","G","C"),]
  }
  
  class.summary = c()
  class.summary.df = c()
  TiTv.fractions = data.frame()
  
  mut.raw.counts = data.frame(row.names = c("A-G","T-C","C-T","G-A","A-T","T-A","A-C","T-G","C-A","G-T","C-G","G-C"))
  
  pid.split = split(tot.muts,as.factor(as.character(tot.muts[,6])))
  
  for(i in 1:length(pid.split)){
    
    mut = pid.split[[i]]
    
    mut$con = paste(mut[,4],mut[,5],sep="-")
    mut$con = as.factor(mut$con)
    mut$con = factor(mut$con,levels = c("A-G","T-C","C-T","G-A","A-T","T-A","A-C","T-G","C-A","G-T","C-G","G-C"))
    
    mut.df = data.frame(absolute_counts = summary(mut$con))
    mut.raw.counts = cbind(mut.raw.counts,mut.df)
    colnames(mut.raw.counts)[ncol(mut.raw.counts)] = names(pid.split[i])
    mut.df$class = rep(rownames(mut.df)[seq(1,12,by = 2)],each = 2)
    mut.df$TiTv = c(rep(x = "Ti",4),rep("Tv",8))
    mut.df$percent_contribution = mut.df$absolute_counts/sum(mut.df$absolute_counts)*100 #calculate contribution of each conversion type - in %)
    mut.df$conversion = row.names(mut.df)
    mut.df = mut.df[,c(5,1:4)] #Reorder the columns
    
    mut.summarized = ddply(.data = mut.df,.variables = ~class,summarise, sum(percent_contribution))
    colnames(mut.summarized) = c("class","percent_contribution")
    mut.summarized2 = ddply(.data = mut.df,.variables = ~TiTv,summarise, sum(percent_contribution))
    colnames(mut.summarized2) = c("TiTv","percent_contribution")
    TiTv.fractions = rbind(TiTv.fractions,mut.summarized2$percent_contribution)
    rownames(TiTv.fractions)[i] = names(pid.split[i])
    
    mut.summary = list(mut.df,mut.summarized,mut.summarized2) #make a summarized list of all three data frames.
    class.summary.df = rbind(class.summary.df,c(mut.summary[[2]][,2],names(pid.split)[i])) #make a table
    class.summary = rbind(class.summary,cbind(mut.summary[[2]],pid=names(pid.split)[i])) #make table for ggplot
  }
  
  colnames(TiTv.fractions) = c("Ti","Tv")
  if(plot){
    p = ggplot(data = class.summary,aes(x = class,y = percent_contribution))+geom_boxplot()+theme_minimal()+ylim(0,100)+xlab("Mutation Class")+theme(axis.line = element_line(colour = "black"))+ylab("Fraction of Mutations")
    print(p)
  }
  
  class.summary.df = as.data.frame(class.summary.df[,c(7,1:6)])
  colnames(class.summary.df) = c("pid",mut.summary[[2]][,1])
  
  return(list(fraction.contribution = class.summary.df,raw.counts = mut.raw.counts,TiTv.fractions = TiTv.fractions))
}

