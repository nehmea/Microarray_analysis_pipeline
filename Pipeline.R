

###################### setting up variables ###############################################
download_series_matrix = "TRUE" #should you start by downloading the series matrix? you may otherwise use your pre-prepared pdata file
dataset = "GSE70515" #name of desired GEO dataset to be analyzed
primary_variable = "disease" #the name of the variable. you may get this info from any GSM webpage related to the dataset
#primary_replacements = c( 
#  '.*cancer.*' = 'Tumor',
#  '.*tumor.*' = 'Tumor',
#  '.*normal.*' = 'Normal',
#  '.*Normal.*' = 'Normal')
secondary_variable = NA #either specify using a character string of length 1, or NA if not available
tertiary_variable = NA #either specify using a character string of length 1, or NA if not available
write_downloaded_dataset = "TRUE" #would yu like to save the downloaded series matrix to a local file?

pheno_file = "pGSE157784.txt" #if download_series_matrix = "FALSE", what is the name of the phenodata file?

control_condition = 'Unstimulated' #the name of the control condition for your primary variable? you may get this info from a control GSM webpage related to the dataset

annotdb = 'hgu133plus2.db' #please specify using database name from bioconductor  


###################### Data download and preparation ########################################
library(tidyverse)

if(download_series_matrix == "TRUE"){
library(GEOquery)

ini_dataset = getGEO(dataset, GSEMatrix=TRUE)[[1]]

pdata = pData(ini_dataset)
pdata = subset(pdata, select = grep('.*characteristics', colnames(pdata)))
colnames(pdata) = as.character(pdata[1,]) 
pdata = pdata %>%
  rename_with(gsub, pattern = ':.*', replacement = '') %>%
  mutate_all(gsub, pattern = '.*: ', replacement = '')
#pdata[, primary_variable] = str_replace_all(pdata[, primary_variable], primary_replacements)
if(all.equal(rownames(pdata), rownames(pData(ini_dataset)))){pData(ini_dataset) <- pdata}


fData(ini_dataset)= fData(ini_dataset) %>%
  select(all_of(c("ID", 'Gene Symbol', "ENTREZ_GENE_ID"))) %>%
  rename_with(gsub, pattern = ' ', replacement = '_') %>%
  mutate(Symbol_first = gsub(pattern = ' ///.*', replacement = '', Gene_Symbol))

  
#assign(dataset, ini_dataset)
#rm(ini_dataset)


if(write_downloaded_dataset == "TRUE"){
  write.table(exprs(ini_dataset), paste0('e', dataset, '.txt'), sep = '\t', col.names = NA)
  write.table(pData(ini_dataset), pheno_file, sep = '\t', col.names = NA)
  write.table(fData(ini_dataset), paste0('f', dataset, '.txt'), sep = '\t', col.names = NA)
  
  }

}

rm(ini_dataset)

########################## Raw Data preprocessing #############################
library(affy)
library(data.table)
library(GEOquery)
{
pdata = read.table(pheno_file, header = T, row.names = 1, sep = "\t", stringsAsFactors = T)

options(timeout = 10000)
getGEOSuppFiles(dataset, makeDirectory = F, fetch_files = T)
options(timeout = 60)

untar(paste0(dataset, '_RAW.tar'))

celfiles = list.files(pattern = '.*.CEL.gz')
eset = ReadAffy(filenames = celfiles, compress = T, phenoData = pdata, verbose = T)
saveRDS(eset, 'initial_raw_affybatch.rds')}

########################## Quality control #######################################
{
##yaqc
library(yaqcaffy)
yaqc_res <- yaqc(eset, verbose=TRUE)
exclude= unique(c(names(getOutliers(yaqc_res,"sfs")), 
           names(getOutliers(yaqc_res,"biob")),
           names(getOutliers(yaqc_res,"avbg")),
           names(getOutliers(yaqc_res, 'actin')),
           names(getOutliers(yaqc_res, 'gapdh'))))

eset_filtered_hq = eset[, !colnames(eset) %in% exclude]

png(file="yaq_res.tiff", height=7, width=7, res=300, units = "in")
plot(yaqc_res)
dev.off()

saveRDS(yaqc_res, 'yaqc_res.rds')
rm(list = c('yaqc_res', 'eset', 'pdata'))
save.image()

#file.remove(celfiles)

#RLE and NUSE
library(affyPLM)
PLMqc <- fitPLM(eset_filtered_hq, normalize=T, background=T, verbosity.level = 5)
#RLE(PLMqc)
#NUSE(PLMqc)

RLE.res = data.frame(t(RLE(PLMqc, type="stats")))
write.table(RLE.res, "raw_RLE.txt", sep = "\t", col.names = NA)

NUSE.res = data.frame(t(NUSE(PLMqc, type= "stats")))
write.table(NUSE.res, "raw_NUSE.txt", sep = "\t", col.names = NA)

if(all.equal(rownames(RLE.res), rownames(NUSE.res))){
  exclude2 = unique(c(rownames(NUSE.res)[which(NUSE.res$median>1.05 | NUSE.res$median<(-1.05)| 
                                           RLE.res$median>0.15 | RLE.res$median<(-0.15))]))}

eset_filtered_hq = eset_filtered_hq[, !colnames(eset_filtered_hq) %in% exclude2]

rm(PLMqc)
save.image()

##arrayqualitymetrics
#gcrma
library(gcrma)
include_celfiles = sampleNames(eset_filtered_hq)
celfiles2 = unique(grep(paste(include_celfiles, collapse="|"), celfiles, value=TRUE))
  #paste0(sampleNames(eset_filtered_hq), '.CEL.gz')
eset_filtered_hq_gcrma = justGCRMA(filenames = celfiles2, compress = T,
                                   phenoData = pData(eset_filtered_hq),
                                   verbose = T)

if(!all.equal(colnames(exprs(eset_filtered_hq_gcrma)), rownames(pData(eset_filtered_hq_gcrma)))){
  stop('mismatch in Expressionset')}

library(arrayQualityMetrics)
aqm.res = arrayQualityMetrics(expressionset = eset_filtered_hq_gcrma,
                              outdir = "ArrayQualityMetrics_raw_hq",
                              force = T, spatial = T, do.logtransform = F,
                              intgroup = primary_variable,
                              reporttitle=paste0(dataset, "_hq_gcrma_QCmetrics"))

exclude3 = unique(names(c(aqm.res$modules$pca@outliers@which,
                   aqm.res$modules$heatmap@outliers@which,
                   aqm.res$modules$boxplot@outliers@which,
                   aqm.res$modules$density@outliers@which,
                   aqm.res$modules$meansd@outliers@which,
                   aqm.res$modules$maplot@outliers@which)))
exclude3 = exclude3[exclude3 != '']

eset_filtered_hq = eset_filtered_hq[, !colnames(eset_filtered_hq) %in% exclude3]
eset_filtered_hq_gcrma = eset_filtered_hq_gcrma[, !colnames(eset_filtered_hq_gcrma) %in% exclude3]

saveRDS(aqm.res, 'aqm_res_rds'); saveRDS(eset_filtered_hq, 'eset_filtered_hq')
rm(aqm.res, eset_filtered_hq)
save.image()


exlude_table = rbind(if(!length(exclude)==0){data.frame(SampleID = exclude, Label = 'excluded_yaqc')},
                     if(!length(exclude2)==0){data.frame(SampleID = exclude2, Label = 'excluded_PLM')},
                     if(!length(exclude3)==0){data.frame(SampleID = exclude3, Label = 'excluded_aqm')},
                     data.frame(SampleID = colnames(eset_filtered_hq_gcrma), Label = 'retained'))
                     
write.table(exlude_table, 'exclusion_table.txt', sep = '\t', col.names= NA)

#gcrma
library(gcrma)
#untar(paste0(dataset, '_RAW.tar'))
include_celfiles = sampleNames(eset_filtered_hq_gcrma)
celfiles2 = unique(grep(paste(include_celfiles, collapse="|"), celfiles, value=TRUE))
#rm(eset_filtered_hq_gcrma)
eset_filtered_hq_gcrma = justGCRMA(filenames = celfiles2, compress = T,
                                   phenoData = pData(eset_filtered_hq_gcrma),
                                   verbose = T)

if(!all.equal(colnames(exprs(eset_filtered_hq_gcrma)), rownames(pData(eset_filtered_hq_gcrma)))){
  stop('mismatch in Expressionset')}

file.remove(celfiles)

pData(eset_filtered_hq_gcrma) = data.frame(pData(eset_filtered_hq_gcrma), stringsAsFactors = T)
pData(eset_filtered_hq_gcrma)[, primary_variable] = relevel(factor(pData(eset_filtered_hq_gcrma)[, primary_variable]),
                                                            ref = control_condition)
write.table(pData(eset_filtered_hq_gcrma), 'pData_eset_filtered_hq_gcrma.txt', col.names = NA, sep = '\t')
write.table(exprs(eset_filtered_hq_gcrma), 'exprs_eset_filtered_hq_gcrma.txt', col.names = NA, sep = '\t')
save.image()

}

################################# HQ samples visualization ##################################
{

library(RColorBrewer)
library(factoextra)
library(ggplot2)
library(ggthemes)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
df_pca <- prcomp(exprs(eset_filtered_hq_gcrma), center = F)
df_out <- as.data.frame(df_pca$rotation)[,1:3]
df_out = cbind(df_out, pData(eset_filtered_hq_gcrma))

percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste0(colnames(df_out), " (", paste(as.character(percentage), " %", ")"))


variables = colnames(pData(eset_filtered_hq_gcrma))
plots = lapply(variables, function(p)
  if(is.factor(df_out[,p])|is.character(df_out[,p])){
    ggplot(df_out, aes(PC1, PC2))+
      labs(title=p, x=percentage[1], y = percentage[2], col=p)+
      geom_point(size=2, aes(col=factor(df_out[,p])))+
      scale_color_colorblind()+
      stat_ellipse(aes(col=df_out[,p]))+
      theme_minimal()+
      theme(axis.text.x = element_text(angle=45, hjust=1))
    
  }
  else if(is.numeric(df_out[,p]) | is.integer(df_out[,p])){
    ggplot(df_out, aes(PC1, PC2))+
      labs(title=p, x=percentage[1], y = percentage[2], col=p)+
      geom_point(size=2, aes(col=as.numeric(df_out[,p])))+
      scale_colour_gradient(low=colors[200], high = colors[1])+
      stat_ellipse(aes(col=df_out[,p]))+
      theme_minimal()
  }
)

library(gridExtra)
library(grid)
png(file="PCA.tiff", height=(length(variables)/4)*7, width=7, res=300, units = "in")
print(marrangeGrob(plots, nrow=length(variables)/2, ncol=2, 
                   newpage = F))
dev.off()

}
########################### statistical analysis ##############################################
if(length(sampleNames(eset_filtered_hq_gcrma))!=0){
  
library(limma)
library(plyr)

eseti = eset_filtered_hq_gcrma

limma.table = data.frame()

design = 
if(is.na(secondary_variable) & is.na(tertiary_variable)){
model.matrix(formula(paste("~", primary_variable)), eseti)} else if(!is.na(secondary_variable) & is.na(tertiary_variable)){
    model.matrix(formula(paste("~", primary_variable, '+', secondary_variable)), eseti)} else if(!is.na(secondary_variable) & !is.na(tertiary_variable)){
        model.matrix(formula(paste("~", primary_variable, '+', secondary_variable, '+', tertiary_variable)), eseti)}
    

#colnames(design) = if(is.na(secondary_variable) & is.na(tertiary_variable)){
  #c('Intercept', primary_variable)} else if(!is.na(secondary_variable) & is.na(tertiary_variable)){
   #c('Intercept', primary_variable, secondary_variable)} else if(!is.na(secondary_variable) & !is.na(tertiary_variable)){
      #c('Intercept', primary_variable, secondary_variable, tertiary_variable)}


fit = lmFit(eseti, design)
fit = eBayes(fit)

primary_factor = pData(eseti)[, primary_variable]
for(i in 2:nlevels(primary_factor)){
deg=topTable(fit, coef=paste0(primary_variable, levels(primary_factor)[i]), 
             number=Inf, sort.by = "P", resort.by = "logFC") 
limma.table = rbind(limma.table, data.frame(Group = paste0(levels(primary_factor)[i],'_vs_', control_condition), deg))
}

options(connectionObserver = NULL)
library(annotdb, character.only= T)
library(AnnotationDbi)
limma.table = data.frame(ENTREZ = mapIds(get(annotdb), rownames(limma.table), 'ENTREZID', 'PROBEID'),
                       SYMBOL = mapIds(get(annotdb), rownames(limma.table), 'SYMBOL', 'PROBEID'),
                       limma.table)
write.table(limma.table, "limma_stats.txt", sep = "\t", col.names = NA)

deg.table = limma.table[limma.table$adj.P.Val<0.05,]
deg.table$subGroup = factor(ifelse(deg.table$logFC>0, paste0(deg.table$Group,".up"),paste0(deg.table$Group,".down")))

DEGlist = split(deg.table$SYMBOL,deg.table$subGroup)
DEGlist2 = split(deg.table$SYMBOL,deg.table$Group)

write.table(deg.table, "DEGs_p.adjust0.05.txt", sep = "\t", col.names = NA)

rm(design, fit, eseti, deg)

}
#################################### enrichment analysis ###########################################
if(length(DEGlist2)!=0){

library(msigdbr)
library(clusterProfiler)
library(ggplot2)

#msigdbr_show_species()
#data.frame(msigdbr_collections())

gene_sets = data.frame()
genelist = limma.table$logFC; names(genelist) = limma.table$SYMBOL

for(category in c('H',"C2", 'C3')){
gene_sets = rbind(gene_sets, msigdbr(species = "Homo sapiens", category = category))}
gene_sets = subset(gene_sets, gs_subcat %in% c('', 'CP:KEGG', 'MIR:MIRDB', 'TFT:GTRD'))
gene_sets$gs_subcat = ifelse(gene_sets$gs_subcat=='', gene_sets$gs_cat, gene_sets$gs_subcat)
gene_sets$gs_subcat = factor(gsub(':', '_', gene_sets$gs_subcat))
  #dplyr::select(gs_name, entrez_gene)


for(gs in levels(factor(gene_sets$gs_subcat))){
  ck = compareCluster(geneCluster = DEGlist, fun = "enricher", 
                      TERM2GENE=gene_sets[gene_sets$gs_subcat==gs, c('gs_name', 'gene_symbol')],
                      pAdjustMethod = "BH",
                      universe = limma.table$SYMBOL,
                      qvalueCutoff = 0.05)
  
  svg(file=paste0('dotplot_', gs, '.svg'), height=15, width=15)
  print(dotplot(ck, showCategory = 30)+
          labs(title=paste0(gs, ': DEG sets'))+
          scale_color_continuous(low = colors[1], high = colors[75],
                                 name = "p.adjust", guide=guide_colorbar(reverse=T)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

rm(ck)

for(gs in levels(factor(gene_sets$gs_subcat))){
  print(gs)
  for(cluster in names(DEGlist)){
    print(cluster)
  enrichments = enricher(DEGlist[[cluster]], 
                      TERM2GENE=gene_sets[gene_sets$gs_subcat==gs, c('gs_name', 'gene_symbol')],
                      pAdjustMethod = "BH",
                      universe = limma.table$SYMBOL,
                      qvalueCutoff = 0.05)
  if(!min(enrichments@result$qvalue)>0.05){
    svg(file=paste0('cnetplot_', gs, '_', cluster, '.svg'), height=15, width=15)
  print(cnetplot(enrichments, categorySize="pvalue", foldChange=genelist))
  dev.off()}
  }
}
rm(enrichments)

}

