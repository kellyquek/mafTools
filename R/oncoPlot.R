#' Draws oncoplot from the input maf file.
#' 
#' Takes maf file as input and plots it as a matrix similar to the output generated by oncoprint.
#' This function heavily depends on oncoplot script by Zuguang Gu for plotting (https://github.com/jokergoo/ComplexHeatmap/blob/master/vignettes/oncoprint.R)
#' while taking care of format conversion.
#' 
#' @param maf_file input file in MAF format.
#' @param removeSilent logical. Whether to discard silent mutations ("Silent","Intron","RNA","3'UTR","3'Flank","5'UTR","5'Flank","IGR") before forming matrix. Default is TRUE.
#' @param oncoPlot logical. Whether to plot mutation matrix (Upto top 20 genes). Default is TRUE. 
#' @param writeMatrix Writes the transformed maf file as a tab delimited text file, to be opened by spreadsheet applications.
#' @param top How many top n genes to plot. Defaults to 20.
#' @return Plots and writes formed matrix to a file.
#' @export


oncoPlot = function(maf_file, removeSilent = TRUE,oncoPlot = TRUE,writeMatrix = FALSE,top = 20){
  
  require(package = "ComplexHeatmap",quietly = T,warn.conflicts = F)
  
  tot.muts = read.delim(file = maf_file,header = T,sep = "\t",stringsAsFactors = F,comment.char = "#")
  tot.muts = tot.muts[,c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","Variant_Type")]
  
  if(removeSilent){
    tot.muts = tot.muts[!tot.muts$Variant_Classification %in% c("Silent","Intron","RNA","3'UTR","3'Flank","5'UTR","5'Flank","IGR"),]  
  }
  
  barcode.split = split(tot.muts,as.factor(as.character(tot.muts$Tumor_Sample_Barcode)))
  variant.classes = levels(as.factor(as.character(tot.muts$Variant_Classification)))
  mdf = data.frame(gene=as.character(unique(tot.muts$Hugo_Symbol)))
  
  for(barcode in 1:length(barcode.split))
  {
    barcode1 = barcode.split[[barcode]]
    
    if(nrow(barcode1) == 0) break
    
    barcode1$Hugo_Symbol = as.character(barcode1$Hugo_Symbol)
    
    gene.split = split(barcode1,as.factor(as.character(barcode1$Hugo_Symbol)))
    
    p = data.frame()
    for(i in 1:length(gene.split))
    {
      gene1 = gene.split[[i]]
      
      if(length(unique(gene1$Variant_Classification))=="1"){
        p=rbind(p,data.frame(gene=unique(gene1$Hugo_Symbol),mut=unique(as.character(gene1$Variant_Classification))))
      } else{
        p=rbind(p,data.frame(gene = unique(as.character(gene1$Hugo_Symbol)), mut="two_hit"))
      }
    }
    colnames(p)[2]=as.character(unique(barcode1$Tumor_Sample_Barcode))
    mdf = merge(mdf,p,all=T)
  }
  
  genes = as.character(mdf$gene)
  mdf = mdf[,-1]
  mdf = data.frame(lapply(mdf,as.character),stringsAsFactors=F)
  
  mdf.copy = mdf 
  mdf.copy[is.na(mdf.copy)] = ""
  rownames(mdf.copy) = genes

  
  mdf[is.na(mdf)] = 0
  
  variant.classes = as.list(apply(mdf,2,unique))
  variant.classes = unique(as.vector(unlist(variant.classes)))
  variant.classes = variant.classes[!variant.classes == "0"]
  
  #for(i in 1:length(s))
  #{
  #  mdf[mdf==s[i]]=i
  #}
  
  for(i in 1:length(variant.classes))
  {
    mdf[mdf==variant.classes[i]] = i
  }
  
  
  
  mdf = as.matrix(apply(mdf,2,function(x) as.numeric(as.character(x))))
  rownames(mdf) = genes
  mdf = cbind(mdf,variants = apply(mdf,1,function(x){length(x[x!="0"])}))
  mdf=mdf[order(mdf[,ncol(mdf)],decreasing=T),]
  colnames(mdf)=gsub(pattern="^X",replacement="",colnames(mdf))
  
  nMut = mdf[,ncol(mdf)]
  mdf = mdf[,-ncol(mdf)]
  
  tmdf = t(mdf)
  mdf = t(tmdf[do.call(order,c(as.list(as.data.frame(tmdf)),decreasing=T)),])
  #mdf = mdf[1:20,]
  mdf[is.na(mdf)] = 0
  
  mdf.copy = mdf.copy[rownames(mdf),]
  mdf.copy = mdf.copy[,colnames(mdf)]
  
  mat_origin = as.matrix(mdf.copy)
  mat = mat_origin[1:top,]
  
  col = c(brewer.pal(12,name = "Paired"),brewer.pal(11,name = "Spectral")[1:3],'maroon')
  names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','Intron','Missense_Mutation','IGR','Nonsense_Mutation','RNA','Splice_Site','In_Frame_Del','Frame_Shift_Ins','Silent','In_Frame_Ins','ITD','In_Frame_Ins','Translation_Start_Site',"two_hit")
  
  #if(length(variant.classes) <= 9 ){
  #  type_col = structure(brewer.pal(length(variant.classes),name = "Set1"), names = variant.classes)
  #} else{
  #  type_col = structure(rainbow(n = length(variant.classes),alpha = 1), names = variant.classes)
  #}
  
  type_col = structure(col[variant.classes], names = names(col[variant.classes]))
  type_name = structure(variant.classes,names = variant.classes)
  
  if(writeMatrix){
    write.table(mat_origin,"onco_matrix.txt",sep = "\t",quote = F)
  }
  
  if(oncoPlot){
    
    add_oncoprint = function(type, x, y, width, height) {
      
      for(i in 1:length(variant.classes)){
        if(any(type %in% variant.classes[i])) {
          grid.rect(x, y, width - unit(0.5, "mm"), height - unit(1, "mm"), gp = gpar(col = NA, fill = type_col[variant.classes[i]]))
        }
      }
      if(any(type %in% "")) {
        grid.rect(x, y, width - unit(0.5, "mm"), height - unit(1, "mm"), gp = gpar(col = NA, fill = "#CCCCCC"))
      }
    }
    
    
    #####################################################################
    # row annotation which shows percent of mutations in all samples
    anno_pct = function(index) {
      n = length(index)
      pct = apply(mat_origin[index, ], 1, function(x) sum(!grepl("^\\s*$", x))/length(x))*100
      pct = paste0(round(pct),"%")
      pushViewport(viewport(xscale = c(0, 1), yscale = c(0.5, n + 0.5)))
      grid.text(pct, x = 1, y = seq_along(index), default.units = "native", just = "right", gp = gpar(fontsize = 10))
      upViewport()
    }
    
    ha_pct = HeatmapAnnotation(pct = anno_pct, width = grobWidth(textGrob("100%", gp = gpar(fontsize = 10))), which = "row")
  
    #####################################################################
    # row annotation which is a barplot
    anno_row_bar = function(index) {
      n = length(index)
      tb = apply(mat[index, ], 1, function(x) {
        x = unlist(strsplit(x, ";"))
        x = x[!grepl("^\\s*$", x)]
        x = sort(x)
        table(x)
      })
      max_count = max(sapply(tb, sum))
      pushViewport(viewport(xscale = c(0, max_count*1.1), yscale = c(0.5, n + 0.5)))
      for(i in seq_along(tb)) {
        if(length(tb[[i]])) {
          x = cumsum(tb[[i]])
          grid.rect(x, i, width = tb[[i]], height = 0.8, default.units = "native", just = "right", gp = gpar(col = NA, fill = type_col[names(tb[[i]])]))
        }
      }
      breaks = grid.pretty(c(0, max_count))
      grid.xaxis(at = breaks, label = breaks, main = FALSE, gp = gpar(fontsize = 10))
      upViewport()
    }
    
    ha_row_bar = HeatmapAnnotation(row_bar = anno_row_bar, width = unit(4, "cm"), which = "row")
    
    ###################################################################
    # column annotation which is also a barplot
    anno_column_bar = function(index) {
      n = length(index)
      tb = apply(mat[, index], 2, function(x) {
        x = unlist(strsplit(x, ";"))
        x = x[!grepl("^\\s*$", x)]
        x = sort(x)
        table(x)
      })
      max_count = max(sapply(tb, sum))
      pushViewport(viewport(yscale = c(0, max_count*1.1), xscale = c(0.5, n + 0.5)))
      for(i in seq_along(tb)) {
        if(length(tb[[i]])) {
          y = cumsum(tb[[i]])
          grid.rect(i, y, height = tb[[i]], width = 0.8, default.units = "native", just = "top", gp = gpar(col = NA, fill = type_col[names(tb[[i]])]))
        }
      }
      breaks = grid.pretty(c(0, max_count))
      grid.yaxis(at = breaks, label = breaks, gp = gpar(fontsize = 10))
      upViewport()
    }
    
    ha_column_bar = HeatmapAnnotation(column_bar = anno_column_bar, which = "column")
    
    #####################################################################
    # the main matrix
    ht = Heatmap(mat, rect_gp = gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
      type = mat[i,j]
      add_oncoprint(type, x, y, width, height)
    }, row_names_gp = gpar(fontsize = 10), show_column_names = FALSE, show_heatmap_legend = FALSE,
    top_annotation = ha_column_bar, top_annotation_height = unit(2, "cm"))
    
    ht_list = ha_pct + ht + ha_row_bar
    
    #########################################################
    # legend
    legend = legendGrob(labels = type_name[names(type_col)], pch = 15, gp = gpar(col = type_col), nrow = 1)
    
    #pdf("oncoprint.pdf", width = 10, height = 10,bg = "white",pointsize = 9,paper = "special")
    draw(ht_list, newpage = FALSE, annotation_legend_side = "bottom", annotation_legend_list = list(legend))
    #dev.off()
    
    print("This plot is generated using oncoprint script by Zuguang Gu (https://github.com/jokergoo/ComplexHeatmap/blob/master/vignettes/oncoprint.R)",verbose = T,quote = F)
    
  }#OncoPlot End
  
}
