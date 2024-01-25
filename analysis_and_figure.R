#### Loading package ####

{
  library(ggrepel)
  library(RColorBrewer)
  library(ggVennDiagram)
  library(Rtsne)
  library(CEMiTool)
  library(fgsea)
  library(KEGGREST)
  library(org.Sc.sgd.db)
  library(nortest)
  library(rrvgo)
  library(ROTS)
  library(impute)
  library(mgcv)
  library(ggtree)
  library(pegas)
  library(magick)
  library(relaimpo)
  library(miscTools)
  library(fitdistrplus)
  library(MALDIquant)
  library(multcompView)
  library(topGO)
  library(gtools)
  library(pbapply)
  library(WGCNA)
  library(Biostrings)
  library(seqinr)
  library(intrval)
  library(gridExtra)
  library(rlist)
  library(GGally)
  library(dendextend)
  library(factoextra)
  library(NbClust)
  library(clipr)
  library(dplyr)
  library(bigutilsr)
  library(performance) 
  library(igraph)
  library(network)
  library(sna)
  library(ggplot2)
  library(ggfortify)
  library(tidyverse)
  library(itertools)
  library(gplots)
  library(ape)
  library(clusterProfiler)
  library(enrichplot)
  library(data.table)
  library(ggpubr)
  library(BiocManager)
  library(readr)
  library(ggcorrplot)
  library(plotly)
  library(ggsignif)
  library(pathfindR)
  library(eulerr)

  stat_box_data <- function(y, upper_limit = max(y) +sd(y)) { #can be modified if needed#
    return( 
      data.frame(
        y = 1.15 * upper_limit,
        label = paste('count =', length(y), '\n',
                      'mean =', round((mean(y)), 6), '\n')
      )
    )
  }

 
  give.n <- function(x){
    return(c(y = max(x)+sd(x), label = length(x))) 
    # experiment with the multiplier to find the perfect position
  }

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }

  
}


####data####

gene_info = fread('data_upload/CompleteGeneAnnot_strand.csv', data.table = F)
rownames(gene_info)<-gene_info$systematic_name

{
  proteomic_data <- fread('data_upload/DIA-NN_1.8/SCmedia_MBR_CommonPeptide_Approach/211025_ProteomicsData_SCmediaMBR_DIA-NN-1.8_genes_ORF.tsv', data.table = F) # proteomic data
  rownames(proteomic_data)<- proteomic_data$Protein.Group
  proteomic_data=proteomic_data[,-1]
  colnames(proteomic_data)=gsub('\\.','_', colnames(proteomic_data))
}


RNA_data_jing = read.table('data_upload/jing_data/input_vsd.csv') # RNAseq data
d = gene_info
rownames(d) = d$Gene
RNA_data_jing=RNA_data_jing[-which(duplicated(d[rownames(RNA_data_jing),"systematic_name"])),]
RNA_data_jing=RNA_data_jing[-which(is.na(d[rownames(RNA_data_jing),"systematic_name"])),]
rownames(RNA_data_jing)= d[rownames(RNA_data_jing),"systematic_name"]

a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = as.matrix(RNA_data_jing[b,a])
d = as.matrix(proteomic_data[b,a])
d = impute.knn(d)$data
c = impute.knn(c)$data
c = (log2(c))
d= (log2(d))

temp = sample(intersect(colnames(proteomic_data),colnames(RNA_data_jing)),5)


e = cbind(c,d)
boxplot(e[,colnames(e)%in%temp])

df_rank <- apply(e,2,rank,ties.method="min")
df_sorted <- data.frame(apply(e, 2, sort))
df_mean <- apply(df_sorted, 1, mean)
index_to_mean <- function(my_index, my_mean){
  return(my_mean[my_index])
}

df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean) #final normalized data
rownames(df_final)= intersect(rownames(proteomic_data), rownames(RNA_data_jing))
boxplot(df_final[,colnames(e)%in%temp])




clades = fread('data_upload/strains.csv', sep= ';', data.table = F)
clades$new = substr(clades$Clades,1,3)
rownames(clades)= clades$`Standardized name`

#a = "org.Sc.eg.db"
#library(a, character.only = TRUE)

goterms <- Term(GOTERM)
download_Yeast_GO_mapping <-
  function(yeast.GO.url='http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz'){
    ##### download Yeast Gene Ontology mapping  http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz
    GO.assocs.all  <- read.table(textConnection(readLines(gzcon(url(yeast.GO.url))))
                                 ,header=FALSE,check.names=FALSE, stringsAsFactors=FALSE, sep="\t", quote=NULL, comment.char ='!'
    )
    my.col.names <- c("Database","SGDID","DB_Object_Symbol","NOT","GOID","DB:Ref","Evidence"
                      ,"With:From","GO_Aspect","DB_Object_Name","Synonym","Type","Taxon","Date","Asigned By","Notes 1")
    GO.assocs.all <- GO.assocs.all[,1:length(my.col.names)] ### sometimes there are empty extra columns
    colnames(GO.assocs.all) <- my.col.names
    
    ### use only most important columns "SGDID", "DB_Object_Symbol", "GOID"
    GO.assocs <- GO.assocs.all[,match(c("SGDID", "GOID"),colnames(GO.assocs.all))]
    ### one can filter for a certain evidence like IDA : Inferred from Direct Assay 
    ### in our case we will accept any of the Evidences and keep unique records of the association
    GO.assocs  <- unique(GO.assocs) #### for some reasons we have same rows duplications (for example the multiple Evidence levels: IBA, IC, IDA etc)
  }

Yeast.GO.assocs =download_Yeast_GO_mapping(  yeast.GO.url = "http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz")

Yeast.GO.assocs=Yeast.GO.assocs[Yeast.GO.assocs$SGDID%in%gene_info$SGD_name,]

a = unique(Yeast.GO.assocs$SGDID)
b = lapply(a,function(i){
  b = which(gene_info$SGD_name%in%i)
  return(b)
})
b= as.character(b)
b = rownames(gene_info)[as.numeric(b)]
a= cbind(a,b)
a= as.data.frame(a)
rownames(a)=a$a
Yeast.GO.assocs$STDID= a[Yeast.GO.assocs$SGDID,"b"]
a =  unique(Yeast.GO.assocs$STDID)
b = lapply(a, function(i){
  b = Yeast.GO.assocs[Yeast.GO.assocs$STDID==i,"GOID"]
})
names(b)= unique(Yeast.GO.assocs$STDID)
geneID2GO <- b 


a =  unique(Yeast.GO.assocs$GOID)
b = lapply(a, function(i){
  b = Yeast.GO.assocs[Yeast.GO.assocs$GOID==i,"STDID"]
})
names(b)= unique(Yeast.GO.assocs$GOID)
GO2GeneID <- b 


load('data_upload/1011_tree_assets.rdata')
cladesdf2 <- cladesdf %>% mutate(inSubset='In subset')
cols <- c(gg_color_hue(length(unique(cladesdf2$clades))), 'Not present'="transparent")
names(cols) <- c(unique(cladesdf2$clades),'Not present')


more_info_strain = fread('data_upload/FinalStrainInfo1040.csv', data.table = F)
Domestication <- more_info_strain %>% select(Standardized_name, Domestication)
Domestication=Domestication[Domestication$Standardized_name%in%colnames(proteomic_data),]




chrom_loc = fread('data_upload/chrLoc_ORFs')
control_GWAS_correlated_CNV= fread('data_upload/Results_GWAS_CNVs_Proteomics_Transcriptomics_455CorrODSamples_VL_20240117.tsv',data.table = F)

control_GWAS_correlated_CNV= cbind(control_GWAS_correlated_CNV, do.call(rbind,strsplit(control_GWAS_correlated_CNV$Pheno,'\\.')))
control_GWAS_correlated_CNV$unique = paste(control_GWAS_correlated_CNV$Chr, control_GWAS_correlated_CNV$ChrPos,control_GWAS_correlated_CNV$`1`,sep='_')
control_GWAS_correlated_CNV$type = ifelse(control_GWAS_correlated_CNV$`1`==substr(control_GWAS_correlated_CNV$SNP,3,80),'CIS','TRANS')

control_GWAS_correlated_CNV$Pheno_chr <- chrom_loc$V1[match(control_GWAS_correlated_CNV$`1`,chrom_loc$V4)]
control_GWAS_correlated_CNV$Pheno_pos <- chrom_loc$V2[match(control_GWAS_correlated_CNV$`1`,chrom_loc$V4)]
control_GWAS_correlated_CNV$Pheno_pos_end <- chrom_loc$V3[match(control_GWAS_correlated_CNV$`1`,chrom_loc$V4)]


chrom_cum_p <- chrom_loc %>% 
  group_by(V1) %>% 
  summarise(max_pos = max(V3)) %>% 
  mutate(Pheno_chr = V1, pos_add_p = lag(cumsum(max_pos), default = 0))

chrom_cum_s <- chrom_loc %>% 
  group_by(V1) %>% 
  summarise(max_pos = max(V3)) %>% 
  mutate(Chr = V1, pos_add_s = lag(cumsum(max_pos), default = 0))

control_GWAS_correlated_CNV <- control_GWAS_correlated_CNV %>% 
  inner_join(chrom_cum_p[c(3:4)], by = "Pheno_chr") %>% 
  inner_join(chrom_cum_s[c(3:4)], by = "Chr") %>% 
  mutate(pos_cum_pheno = Pheno_pos + pos_add_p) %>%
  mutate(pos_cum_snp = ChrPos + pos_add_s)

x= fread('data_upload/Results_GWAS_CNVs_Proteomics_Transcriptomics_455CorrODSamples_VL_20240117.Thresholds.tsv',data.table = F)
temp = apply(control_GWAS_correlated_CNV,1,function(i){
  i<<-i
  
  return(as.numeric(i['PValue'])< (x[x$V1==i['Pheno'],"V2"]*(631/101836)))
  
})
control_GWAS_correlated_CNV$keep = temp
control_GWAS_correlated_CNV=control_GWAS_correlated_CNV[control_GWAS_correlated_CNV$keep,]


control_GWAS_correlated_CNV$type= ifelse(control_GWAS_correlated_CNV$`1`==substr(control_GWAS_correlated_CNV$SNP,3,80),'CIS','TRANS')


control_GWAS_correlated = fread('data_upload/Results_GWAS_Proteomics_Transcriptomics_455CorrODSamples_VL_20240117.tsv',data.table = F)
control_GWAS_correlated= cbind(control_GWAS_correlated, do.call(rbind,strsplit(control_GWAS_correlated$Pheno,'\\.')))
control_GWAS_correlated$unique = paste(control_GWAS_correlated$Chr, control_GWAS_correlated$ChrPos,control_GWAS_correlated$`1`,sep='_')


length(grep('RNA',control_GWAS_correlated$Pheno))
length(grep('Prot',control_GWAS_correlated$Pheno))
table(table(control_GWAS_correlated$unique))
ld_455 = fread('data_upload/full455Matrix.SNPs.Biallelic.Var.DP10.99pNonMiss.MAF5.plink.ld')
temp = control_GWAS_correlated[control_GWAS_correlated$`2`=='Prot',]
temp <- data.table(temp)
temp$Pheno <- gsub("\\.Prot","",temp$Pheno)

temp$Pheno_chr <- chrom_loc$V1[match(temp$Pheno,chrom_loc$V4)]
temp$Pheno_pos <- chrom_loc$V2[match(temp$Pheno,chrom_loc$V4)]
temp$Pheno_pos_end <- chrom_loc$V3[match(temp$Pheno,chrom_loc$V4)]



chrom_cum_p <- chrom_loc %>% 
  group_by(V1) %>% 
  summarise(max_pos = max(V3)) %>% 
  mutate(Pheno_chr = V1, pos_add_p = lag(cumsum(max_pos), default = 0))

chrom_cum_s <- chrom_loc %>% 
  group_by(V1) %>% 
  summarise(max_pos = max(V3)) %>% 
  mutate(Chr = V1, pos_add_s = lag(cumsum(max_pos), default = 0))

temp <- temp %>% 
  inner_join(chrom_cum_p[c(3:4)], by = "Pheno_chr") %>% 
  inner_join(chrom_cum_s[c(3:4)], by = "Chr") %>% 
  mutate(pos_cum_pheno = Pheno_pos + pos_add_p) %>%
  mutate(pos_cum_snp = ChrPos + pos_add_s)

temp$pos_to_trait <- 0

temp[Chr==Pheno_chr & Pheno_pos > ChrPos]$pos_to_trait <- temp[Chr==Pheno_chr & Pheno_pos > ChrPos]$Pheno_pos-temp[Chr==Pheno_chr & Pheno_pos > ChrPos]$ChrPos

temp[Chr==Pheno_chr & Pheno_pos_end < ChrPos]$pos_to_trait <- temp[Chr==Pheno_chr & Pheno_pos_end < ChrPos]$Pheno_pos_end-temp[Chr==Pheno_chr & Pheno_pos_end < ChrPos]$ChrPos




ReduceList <- function(l){
  
  require(ComplexHeatmap)
  
  for(i in seq_along(l)[-length(l)]) {
    if(length(intersect(l[[i]], l[[i+1]])) > 0) { 
      l[[i+1]] <- unique(c(l[[i]], l[[i+1]]))
      l[[i]] <- as.list(NULL)
    }
  }
  
  l = Filter(function(x) length(x) > 0, l)
  
  
  if(length(l)>1){
    
    m <- list_to_matrix(l)
    
    if(nrow(m[rowSums(m)>1,drop=F,])>0){
      
      to_merge <- unique(lapply(1:nrow(m[rowSums(m)>1,drop=F,]),function(j){
        
        names(which(m[rowSums(m)>1,drop=F,][j,]==1))
      }))
      
      repeat {
        
        for(x in seq_along(to_merge)) {
          
          end <- to_merge[[x]][length(to_merge[[x]])]
          
          l[[end]] <- unique(unlist(l[to_merge[[x]]]))
          
          for(name in to_merge[[x]][-length(to_merge[[x]])]){
            l[[name]] <- as.list(NULL)
            
          }
          
        }
        
        m1 <- list_to_matrix(l)
        
        if(nrow(m1[rowSums(m1)>1,])==0) break
        
      }
      
    }
    Filter(function(x) length(x) > 0, l)
  }else{
    l
  }
}

a = unique(temp$Pheno)
temp= as.data.table(temp)
CI_loc_snp_455 <- pblapply(a,function(i){
  
  
  sig_loci <- temp[Pheno ==i]
  sig_loci <- sig_loci[order(sig_loci$pos_cum_snp)]
  trait_range <- c(unique(sig_loci$Pheno_pos)-25000,unique(sig_loci$Pheno_pos_end)+25000)
  trait_chr <- unique(sig_loci$Pheno_chr)
  
  sig_loci$type = ""
  
  ld_sig_list <- lapply(sig_loci$SNP,function(snp){
    
    ld <- ld_455[R2 >0.6 & (SNP_A==snp | SNP_B==snp)]
    
    
    
    if(nrow(ld)==0){
      NULL
      
    }else{
      ld_pos <- unique(c(ld$SNP_A,ld$SNP_B))
      
      ld_pos
    }
    
  })
  
  names(ld_sig_list) <- sig_loci$SNP
  ld_sig_list <- Filter(Negate(is.null), ld_sig_list)
  
  ld_sig_list_reduce <- ReduceList(ld_sig_list)
  
  
  sig_loci$ld_range <- ""
  
  for(i in seq_along(ld_sig_list_reduce)){
    
    gr_name <- names(ld_sig_list_reduce[i])
    
    chr <- as.numeric(sub("chromosome","",strsplit(gr_name,"_")[[1]][2]))
    
    if(nrow(sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]])>1){
      
      minPos <- min(sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$ChrPos)
      maxPos <- max(sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$ChrPos)
      sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$ld_range = 
        paste0(minPos,"|",maxPos)
      
      if(length(intersect(trait_range[1]:trait_range[2],minPos:maxPos))==0){
        
        sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$type <- "TRANS"
        
      }else 
        if(length(intersect(trait_range[1]:trait_range[2],minPos:maxPos))>0 & chr == trait_chr){
          
          sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$type <- "CIS"
          
        }else{
          sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$type <- "TRANS"
          
        }
      
    }
    
    
  }
  
  sig_loci[ld_range=="" & Chr==Pheno_chr & abs(pos_to_trait)<25000]$type <- "CIS"
  sig_loci[type==""]$type = "TRANS"
  
  sig_loci$ld_mask <- ""
  
  for(range in unique(sig_loci$ld_range)){
    
    if(!range==""){
      sig_loci[ld_range==range][!EffectSize==max(EffectSize)]$ld_mask <- "masked"
    }
    
  }
  
  sig_loci
  
  
})



CI_loc_snp_455_Prot = do.call(rbind, CI_loc_snp_455)

table(CI_loc_snp_455_Prot$ld_mask)
CI_loc_snp_455_Prot[CI_loc_snp_455_Prot$ld_mask=='masked',type]

table(CI_loc_snp_455_Prot[CI_loc_snp_455_Prot$ld_mask=='masked',type])
CI_loc_snp_455_Prot_filt= CI_loc_snp_455_Prot[CI_loc_snp_455_Prot$ld_mask=='',]


temp = control_GWAS_correlated[control_GWAS_correlated$`2`=='RNA',]
temp <- data.table(temp)
temp$Pheno <- gsub("\\.RNA","",temp$Pheno)
chrom_loc = fread('data_upload/chrLoc_ORFs')
temp$Pheno_chr <- chrom_loc$V1[match(temp$Pheno,chrom_loc$V4)]
temp$Pheno_pos <- chrom_loc$V2[match(temp$Pheno,chrom_loc$V4)]
temp$Pheno_pos_end <- chrom_loc$V3[match(temp$Pheno,chrom_loc$V4)]



chrom_cum_p <- chrom_loc %>% 
  group_by(V1) %>% 
  summarise(max_pos = max(V3)) %>% 
  mutate(Pheno_chr = V1, pos_add_p = lag(cumsum(max_pos), default = 0))

chrom_cum_s <- chrom_loc %>% 
  group_by(V1) %>% 
  summarise(max_pos = max(V3)) %>% 
  mutate(Chr = V1, pos_add_s = lag(cumsum(max_pos), default = 0))

temp <- temp %>% 
  inner_join(chrom_cum_p[c(3:4)], by = "Pheno_chr") %>% 
  inner_join(chrom_cum_s[c(3:4)], by = "Chr") %>% 
  mutate(pos_cum_pheno = Pheno_pos + pos_add_p) %>%
  mutate(pos_cum_snp = ChrPos + pos_add_s)

temp$pos_to_trait <- 0

temp[Chr==Pheno_chr & Pheno_pos > ChrPos]$pos_to_trait <- temp[Chr==Pheno_chr & Pheno_pos > ChrPos]$Pheno_pos-temp[Chr==Pheno_chr & Pheno_pos > ChrPos]$ChrPos

temp[Chr==Pheno_chr & Pheno_pos_end < ChrPos]$pos_to_trait <- temp[Chr==Pheno_chr & Pheno_pos_end < ChrPos]$Pheno_pos_end-temp[Chr==Pheno_chr & Pheno_pos_end < ChrPos]$ChrPos




ReduceList <- function(l){
  
  require(ComplexHeatmap)
  
  for(i in seq_along(l)[-length(l)]) {
    if(length(intersect(l[[i]], l[[i+1]])) > 0) { 
      l[[i+1]] <- unique(c(l[[i]], l[[i+1]]))
      l[[i]] <- as.list(NULL)
    }
  }
  
  l = Filter(function(x) length(x) > 0, l)
  
  
  if(length(l)>1){
    
    m <- list_to_matrix(l)
    
    if(nrow(m[rowSums(m)>1,drop=F,])>0){
      
      to_merge <- unique(lapply(1:nrow(m[rowSums(m)>1,drop=F,]),function(j){
        
        names(which(m[rowSums(m)>1,drop=F,][j,]==1))
      }))
      
      repeat {
        
        for(x in seq_along(to_merge)) {
          
          end <- to_merge[[x]][length(to_merge[[x]])]
          
          l[[end]] <- unique(unlist(l[to_merge[[x]]]))
          
          for(name in to_merge[[x]][-length(to_merge[[x]])]){
            l[[name]] <- as.list(NULL)
            
          }
          
        }
        
        m1 <- list_to_matrix(l)
        
        if(nrow(m1[rowSums(m1)>1,])==0) break
        
      }
      
    }
    Filter(function(x) length(x) > 0, l)
  }else{
    l
  }
}

a = unique(temp$Pheno)
temp= as.data.table(temp)
CI_loc_snp_455 <- pblapply(a,function(i){
  
  
  sig_loci <- temp[Pheno ==i]
  sig_loci <- sig_loci[order(sig_loci$pos_cum_snp)]
  trait_range <- c(unique(sig_loci$Pheno_pos)-25000,unique(sig_loci$Pheno_pos_end)+25000)
  trait_chr <- unique(sig_loci$Pheno_chr)
  
  sig_loci$type = ""
  
  ld_sig_list <- lapply(sig_loci$SNP,function(snp){
    
    ld <- ld_455[R2 >0.6 & (SNP_A==snp | SNP_B==snp)]
    
    
    
    if(nrow(ld)==0){
      NULL
      
    }else{
      ld_pos <- unique(c(ld$SNP_A,ld$SNP_B))
      
      ld_pos
    }
    
  })
  
  names(ld_sig_list) <- sig_loci$SNP
  ld_sig_list <- Filter(Negate(is.null), ld_sig_list)
  
  ld_sig_list_reduce <- ReduceList(ld_sig_list)
  
  
  sig_loci$ld_range <- ""
  
  for(i in seq_along(ld_sig_list_reduce)){
    
    gr_name <- names(ld_sig_list_reduce[i])
    
    chr <- as.numeric(sub("chromosome","",strsplit(gr_name,"_")[[1]][2]))
    
    if(nrow(sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]])>1){
      
      minPos <- min(sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$ChrPos)
      maxPos <- max(sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$ChrPos)
      sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$ld_range = 
        paste0(minPos,"|",maxPos)
      
      if(length(intersect(trait_range[1]:trait_range[2],minPos:maxPos))==0){
        
        sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$type <- "TRANS"
        
      }else 
        if(length(intersect(trait_range[1]:trait_range[2],minPos:maxPos))>0 & chr == trait_chr){
          
          sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$type <- "CIS"
          
        }else{
          sig_loci[Chr==chr & SNP %in% ld_sig_list_reduce[[i]]]$type <- "TRANS"
          
        }
      
    }
    
    
  }
  
  sig_loci[ld_range=="" & Chr==Pheno_chr & abs(pos_to_trait)<25000]$type <- "CIS"
  sig_loci[type==""]$type = "TRANS"
  
  sig_loci$ld_mask <- ""
  
  for(range in unique(sig_loci$ld_range)){
    
    if(!range==""){
      sig_loci[ld_range==range][!EffectSize==max(EffectSize)]$ld_mask <- "masked"
    }
    
  }
  
  sig_loci
  
  
})



CI_loc_snp_455_RNA = do.call(rbind, CI_loc_snp_455)

table(CI_loc_snp_455_RNA$ld_mask)
CI_loc_snp_455_RNA[CI_loc_snp_455_RNA$ld_mask=='masked',type]

table(CI_loc_snp_455_RNA[CI_loc_snp_455_RNA$ld_mask=='masked',type])
CI_loc_snp_455_RNA_filt= CI_loc_snp_455_RNA[CI_loc_snp_455_RNA$ld_mask=='',]

CI_loc_snp_455_Prot_filt[CI_loc_snp_455_Prot_filt$unique%in%intersect(CI_loc_snp_455_RNA_filt$unique,CI_loc_snp_455_Prot_filt$unique),type]

turnover_data2=fread('data_upload/230719_turnover_data_export.txt',data.table = F)

QC_CV_control = fread('data_upload/CV_control.csv',data.table = F)



c = proteomic_data

x = apply(c,1,function(i){sd(i,na.rm = T)/mean(na.omit(i))})

x=sort(x)
c= fgsea(stat = x,minSize  = 5,
         maxSize  = 500,pathways = GO2GeneID)
c$type = 'CV'
colnames(c)[6]='NES_CV'
c=na.omit(c)
c= c[c$padj<0.01,]
c=c[order(c$pathway),]
c$term=go2term(c$pathway)[,2]
final_CV_GSEA = c


dis_genet = read.table('data_upload/1039pairwiseDiff.tab', header = T)

####S1####
#####A#####
e= as.data.frame(colnames(proteomic_data))
rownames(e)=e$`colnames(proteomic_data)`
e$clade = clades[e$`colnames(proteomic_data)`,"new"]
e$clade[is.na(e$clade)]=''
c = as.data.frame(table(e$clade))
d = as.data.frame(table(clades$new))
a = as.data.frame(colnames(RNA_data_jing))
rownames(a)=a$`colnames(RNA_data_jing)`
a$clade = clades[a$`colnames(RNA_data_jing)`,"new"]
a$clade[is.na(a$clade)]=''
i = as.data.frame(table(a$clade))
i$dataset = 'Transcriptomic'
c$dataset = 'Proteomic'
d$dataset = 'Peter et al.'
c = rbind(rbind(c,d),i)
c$Var1=as.character(c$Var1)
c$Var1[c$Var1=='']= 'NC'
c$Var1 = str_replace((c$Var1),'\\.','')
c= c[order(c$Freq,decreasing = T),]
c$Var1 = factor(c$Var1, level = unique(c$Var1))

p=ggplot(c,aes(x= Var1, y= Freq , fill = dataset))+
  geom_histogram(stat = 'identity', position = 'dodge',alpha=0.6)+
  theme_classic()+
  scale_fill_manual("Dataset",values = c('black','#EFC000FF','#0073C2FF'))+
  xlab('Clades')+
  ylab('Number of isolates')


#####B#####


a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = as.matrix(RNA_data_jing[b,a])
d = as.matrix(proteomic_data[b,a])

d = impute.knn(d)$data
c = impute.knn(c)$data
c = (log2(c))
d= (log2(d))

temp = sample(intersect(colnames(proteomic_data),colnames(RNA_data_jing)),5)


e = cbind(c,d)
e =e[,colnames(e)%in%temp]
colnames(e)[6:10]= paste(colnames(e)[6:10],'prot',sep = '_')
e= melt(e)
e$dataset= 'Transcriptomic'

e[grep('prot',e$Var2),"dataset"]='Proteomic'
ggplot(e,aes(x= Var2,y=value, fill=dataset))+
  geom_boxplot(alpha=0.6)+
  theme_classic()+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  xlab('')+
  ylab('Expression value')+
  scale_x_discrete(labels = rep(temp,2))+
  geom_vline(xintercept = 5.5)
#####C#####
e=(df_final[,colnames(df_final)%in%temp])
colnames(e)[6:10]= paste(colnames(e)[6:10],'prot',sep = '_')
e= melt(e)
e$dataset= 'Transcriptomic'

e[grep('prot',e$Var2),"dataset"]='Proteomic'
ggplot(e,aes(x= Var2,y=value, fill=dataset))+
  geom_boxplot(alpha=0.6)+
  theme_classic()+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  xlab('')+
  ylab('Expression value')+
  scale_x_discrete(labels = rep(temp,2))+
  geom_vline(xintercept = 5.5)

#####D#####


c = fread("data_upload/RNAseq_969_log2tpm_matrix.tab", data.table = F) #raw tpm data for RNA-seq
rownames(c)=c$ORF
c=c[,-1]
a= rowMeans(c, na.rm = T)
a = as.data.frame(a)
a$prot = NA
colnames(a)[1]= 'TPM'
a[rownames(a)%in%rownames(proteomic_data),"prot"] = 'Yes'
a[is.na(a$prot),"prot"]='No'

ggplot(a,aes(y=TPM, x=prot))+
  geom_violin(fill='#B5C0C1')+
  geom_boxplot(width=0.02)+
  stat_compare_means(label.x = 1.25)+
  theme_classic()+
  xlab('Present in proteomic data')+
  theme(legend.position = 'none')+
  ylab('Transcription level')+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75))

#####E#####



ho_data = fread('data_upload/mmc5.csv',data.table = F) #data from ho et al 2018
rownames(ho_data)=ho_data$V1
ho_data= ho_data[,-1]


ho_data$Encompassed = ifelse(rownames(ho_data)%in%rownames(proteomic_data),'Yes','No')
ggplot(ho_data,aes(x=`Median molecules per cell`,fill=Encompassed))+
  geom_histogram(col='black',alpha=0.5)+
  scale_x_log10()+
  scale_fill_manual(name='Encompassed is \nthis study',values = c('grey','green'))+
  theme_classic()+
  xlab('Median molecules per cell (Ho et al., 2018)')+
  ylab('Number of proteins')



####fig1####
#####C#####
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]



e= lapply(b, function(i){
  e = cor.test(as.numeric(c[i,]),as.numeric(d[i,]), method = 's')
  return(c(as.numeric(e$estimate),as.numeric(e$p.value)))
})


names(e)=b
e=do.call(rbind,e)
e = as.data.frame(e)
e$adj = p.adjust(e$V2,method = 'b')
e$sign = ifelse(e$adj<0.05, 'Significant', 'Not significant')
e[e$V1>0.42,"sign"]= 'Strongly correlated'
ggplot(e, aes(V1, fill=sign))+
  geom_histogram(col ='black',alpha=0.7)+
  scale_fill_manual(values = c('grey', '#bae1ff','#ffb3ba') ,name='')+
  theme_classic()+
  xlab('ρ')+
  geom_vline(xintercept = median(e$V1),linetype = "dashed", linewidth=1, col= '#3399ff')+
  ylab('Number of genes')+
  theme(legend.position = 'none')
median(e$V1)
IQR(e$V1)


b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]



e= lapply(b, function(i){
  e = cor.test(as.numeric(c[i,]),as.numeric(d[i,]), method = 's')
  return(c(as.numeric(e$estimate),as.numeric(e$p.value)))
})

names(e)=b
e=do.call(rbind,e)
e = as.data.frame(e)
min(e$V1)
a = cbind(d['YPL091W',],c['YPL091W',])

b= cbind(d['YLR325C',],c['YLR325C',])
#####D#####
colnames(a)=c('RNA abundance','Protein abundance')
ggplot(a, aes(`RNA abundance`,`Protein abundance`))+
  geom_point(col='#d45f53')+
  theme_classic()+ 
  ggtitle('GLR1')+
  stat_cor(method = 's')+
  scale_y_continuous(limits=c(6.5,10))+
  scale_x_continuous(limits=c(6,10))
#####E#####
colnames(b)=c('RNA abundance','Protein abundance')
ggplot(b, aes(`RNA abundance`,`Protein abundance`))+
  geom_point(col='#1fa0f0')+
  theme_classic()+ 
  ggtitle('RPL38')+
  stat_cor(method = 's')+
  scale_y_continuous(limits=c(6.5,10))+
  scale_x_continuous(limits=c(6,10))

####S2####

ggplot(QC_CV_control,aes(x= CV,fill = group))+
  geom_density(alpha=0.7)+
  scale_x_log10()+
  theme_bw()

####S3####

b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]


b = c()
j = c()
x = c()
for(i in colnames(c)){
  if(unique(rownames(c)==rownames(d))){
    a = cor.test(c[,i],d[,i], method = 's') 
    b = c(b,(as.numeric(a$estimate)))
    
  }
}

b = as.data.frame(b)
rownames(b)= colnames(c)

ggplot(b,aes(x=b))+
  geom_histogram(col='black',bins = 50,fill='lightgrey')+
  theme_classic()+
  geom_vline(xintercept = median(b$b), linetype = "dashed", linewidth=1, col= '#3399ff')+
  ylab('Number of isolates')+
  xlab('rho')

####fig2####

#####A#####

b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]

e= t(combn(colnames(c),2))

n=pbapply(e,1, function(i){
  return(c(median(na.omit(abs(log2(c[,i[1]]/c[,i[2]])))),median(na.omit(abs(log2(d[,i[1]]/d[,i[2]]))))))
})

n = t(n)
colnames(n)= c('Proteomic','Transcriptomic')
n = melt(n)


ggplot(n, aes(y =value,x=Var2))+
  geom_violin(aes(fill = Var2), alpha = 0.6)+
  geom_boxplot(width=0.02)+
  theme_classic()+
  stat_compare_means(label.x = 1.25)+
  xlab('')+
  ylab('|Log2(FC)|')+
  scale_y_log10()+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  theme(legend.position = 'none')+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75))


#####B#####


b = rownames(df_final)
c = df_final[,duplicated(colnames(df_final))]
d = df_final[,!duplicated(colnames(df_final))]

mcor <- cor(as.matrix(c), method = 's', use="complete.obs")
corlowtri <- mcor[lower.tri(mcor)]
temp = mean(corlowtri)

a = as.vector(corlowtri)

mcor <- cor(as.matrix(d), method = 's', use="complete.obs")
corlowtri <- mcor[lower.tri(mcor)]
temp = mean(corlowtri)

b= as.vector(corlowtri)
a= cbind(a, rep('Proteomic', length(a)))
b= cbind(b, rep('Transcriptomic', length(b)))
a = rbind(b,a)
a = as.data.frame(a)
a$b=as.numeric(a$b)
colnames(a)= c('corr','V2')
ggplot(a, aes(x=corr, fill=V2))+
  geom_histogram(alpha=0.5, col= 'black')+
  xlab('Rho')+
  theme_classic()+
  geom_vline(col='#EFC000FF',xintercept = median(a[a$V2=='Proteomic',1]),linetype='dashed')+
  geom_vline(col='#0073C2FF',xintercept = median(a[a$V2=='Transcriptomic',1]), linetype='dashed')+
  scale_x_continuous(breaks = c(0,0.5,round(median(a[a$V2=='Transcriptomic',1]),2),round(median(a[a$V2=='Proteomic',1]),2),1), limits = c(0,1))+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  ylab('')+
  annotate(x = mean(c(median(a[a$V2=='Proteomic',1]),median(a[a$V2=='Transcriptomic',1]))),y = 230000, geom = 'text', label='***')+
  theme(legend.position = 'none')


#####D#####


b = rownames(df_final)
c = df_final[,duplicated(colnames(df_final))]
d = df_final[,!duplicated(colnames(df_final))]


a = dist(t(c))
a = nj(a)
a= sum(a$edge.length)
e =pblapply(1:100, function(j){
  b = dist(t(c[sample(x=1:nrow(df_final),size = nrow(df_final),replace = T), ]))
  a = nj(b)
  
  sum(a$edge.length)})

b = dist(t(d))
b = nj(b)
b= sum(b$edge.length)
f =pblapply(1:100, function(j){
  b = dist(t(d[sample(x=1:nrow(df_final),size = nrow(df_final),replace = T), ]))
  a = nj(b)
  
  sum(a$edge.length)})
x = data.frame(length=c(a,b), data= c('Proteome', 'Transcriptome'), sd= c(sd(unlist(e)),sd(unlist(f))))
ggplot(x, aes(x=data, y=length, fill=data)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),alpha=0.6) +
  geom_errorbar(aes(ymin=length, ymax=length+sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('Branch length')

e = a/b

b= names(GO2GeneID)
b=go2ont(b)
b= b[b$Ontology=='BP',"go_id"]
f <- calculateSimMatrix(b,orgdb="org.Sc.sgd.db",ont=c("BP"),method="Resnik")
reducedTerms <- reduceSimMatrix(f,
                                threshold=0.5 ,
                                orgdb="org.Sc.sgd.db")
p= unique(reducedTerms$parent)
p=GO2GeneID[p]


n=lapply(p,function(i){length(intersect(i, rownames(df_final)))>5})
x= p[unlist(n)]
temp =pblapply(x, function(i){
  i<<-i
  i = intersect(rownames(df_final),i)
  c = df_final[i,duplicated(colnames(df_final))]
  d = df_final[i,!duplicated(colnames(df_final))]
  
  
  a = dist(t(c))
  a = nj(a)
  a= sum(a$edge.length)
  e =lapply(1:10, function(j){
    b = dist(t(c[sample(x=1:nrow(c),size = nrow(c),replace = T), ]))
    a = nj(b)
    
    sum(a$edge.length)})
  b = dist(t(d))
  b = nj(b)
  b= sum(b$edge.length)
  f =lapply(1:10, function(j){
    b = dist(t(d[sample(x=1:nrow(d),size = nrow(d),replace = T), ]))
    a = nj(b)
    
    sum(a$edge.length)})
  
  
  x = a/b
  return(c(a,b,x,t.test(unlist(e),unlist(f))$p.value))
})
a= do.call(rbind,temp)
a=as.data.frame(a)
a$ind = names(x)
a$ind=as.character(a$ind)
a=a[order(a$ind),]
b = go2term(a$ind)
rownames(b)=b$go_id
a$term=b[a$ind,"Term"]
a$adj= p.adjust(a$V4,method = 'b')
b = a[a$adj<0.001,]


final_PTB_tree = a
final_PTB_tree = final_PTB_tree[order(final_PTB_tree$V3),]
final_PTB_tree$type = NA
final_PTB_tree[final_PTB_tree$V3<e&final_PTB_tree$adj<0.001,"type"]='Strongly buffered variation'
final_PTB_tree[final_PTB_tree$V3>e&final_PTB_tree$adj<0.001,"type"]='Enhanced variation'
final_PTB_tree$term = factor(final_PTB_tree$term, levels = final_PTB_tree$term)

b = na.omit(final_PTB_tree)
ggplot(b,aes(term, log2(V3), fill = type))+
  geom_bar(stat = 'identity',alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values = c('#ffb3ba','#bae1ff'))+
  scale_x_discrete(guide = guide_axis(angle = 90),label = function(x) stringr::str_trunc(x, 30))+
  xlab(label = '')+
  ylab('Log2(Branch length Ratio)')

####S4####

#####A#####

b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]

n=c()
b= c()
for(i in 1:nrow(c)){
  a=as.numeric(c[i,])
  n=c(n,setNames(var(na.omit(a)), rownames(c)[i]))
  e = which(is.na(a))
  a = as.numeric(d[i,])
  if(!is_empty(e)){a= a[-e]}
  b = c(b,setNames(var(na.omit(a)), rownames(d)[i]))
}
b = cbind(stack(b),rep('Transcriptomic',length(b)))
n = cbind(stack(n),rep('Proteomic',length(n)))
colnames(b)= c('var','genes','type')
colnames(n)= c('var','genes','type')
b = rbind(b,n)
b = as.data.frame(b)
colnames(b)=c('Var','gene', 'expre')
b$Var=as.numeric(b$Var)
ggplot(b,aes(x=expre,y=Var))+
  geom_violin(aes(fill = expre), alpha=0.6)+
  geom_boxplot(width=0.02)+
  theme_classic()+
  stat_compare_means(label.x = 1.25)+
  xlab('')+
  ylab('Variance')+
  scale_y_log10()+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  theme(legend.position = 'none')+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75))

#####B#####

b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]


a = dist(t(c))
a = a[lower.tri(a)]
b= dist(t(d))
b=b[lower.tri(b)]
b = rbind(cbind((a), rep('Proteomic', length(a))),cbind((b), rep('Transcriptomic', length(b))))
b = as.data.frame(b)

b$V1=as.numeric(as.character( b$V1))

ggplot(b, aes(y =V1,x=V2))+
  geom_violin(aes(fill = V2), alpha=0.6)+
  geom_boxplot(width=0.02)+
  theme_classic()+
  stat_compare_means(label.x = 1.25)+
  xlab('')+
  ylab('Euclidean distances')+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75))+
  theme(legend.position = 'none')

#####C#####

b = rownames(df_final)
c = df_final[,duplicated(colnames(df_final))]
d = df_final[,!duplicated(colnames(df_final))]


a = dist(t(c))
a = nj(a)
a= sum(a$edge.length)
e =pblapply(1:100, function(j){
  b = dist(t(c[sample(x=1:nrow(df_final),size = nrow(df_final),replace = T), ]))
  a = nj(b)
  
  sum(a$edge.length)})

b = dist(t(d))
b = nj(b)
b= sum(b$edge.length)
f =pblapply(1:100, function(j){
  b = dist(t(d[sample(x=1:nrow(df_final),size = nrow(df_final),replace = T), ]))
  a = nj(b)
  
  sum(a$edge.length)})
x = data.frame(length=c(a,b), data= c('Proteome', 'Transcriptome'), sd= c(sd(unlist(e)),sd(unlist(f))))
ggplot(x, aes(x=data, y=length, fill=data)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),alpha=0.6) +
  geom_errorbar(aes(ymin=length, ymax=length+sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()+
  scale_fill_manual("Dataset",values = c('#EFC000FF','#0073C2FF'))+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('Branch length')


####S5####


#####ABCDEF#####

for(j in list(c(1,2),c(3,4),c(5,6))){

  n= df_final[,!duplicated(colnames(df_final))]
  n = t(n)
  n = cbind(clades[rownames(n),"new"],n)
  n = na.omit(n)
  n= as.data.frame(n)
  for(i in 2:ncol(n)){
    n[,i]=as.numeric(n[,i])
    
  }
  x = cols
  names(x) = substr(names(x),1,3)
  x= (x)[order(names(x))]
  names(x)[1:9] = str_replace(names(x)[1:9],'0','')
  names(x)[1:9] = str_replace(names(x)[1:9],'\\.','. ')
  p= autoplot(prcomp(n[,2:ncol(n)]),data= n,colour= 'V1',x=j[1],y=j[2] ,frame = TRUE, frame.type = 'norm')+
    theme_classic()+
    scale_color_manual(values=x)+
    theme(legend.position = 'none')
  p$layers[[2]]$aes_params$alpha <- 0

  plot(p)

  n= df_final[,duplicated(colnames(df_final))]
  n = t(n)
  n = cbind(clades[rownames(n),"new"],n)
  n = na.omit(n)
  n= as.data.frame(n)
  for(i in 2:ncol(n)){
    n[,i]=as.numeric(n[,i])
    
  }
  
  p = autoplot(prcomp(n[,2:ncol(n)]),data= n,colour= 'V1',x=j[1],y=j[2] ,frame = TRUE, frame.type = 'norm')+
    theme_classic()+
    theme(legend.position = 'none')
  p$layers[[2]]$aes_params$alpha <- 0
  plot(p)

}

#####GH#####

a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]

a =cor(d, method = 's',use = "complete.obs")
a = melt(a)
c = pblapply(1:nrow(a), function(i){
  b = sort(as.character(unlist(a[i,1:2])))
  return(c(b,as.numeric(a[i,3]) ))
})
c = matrix(unlist(c), byrow = T,ncol = 3)
c= c[!duplicated(c[,1:2]),]

d = pblapply(1:nrow(dis_genet), function(i){
  b = sort(as.character(unlist(dis_genet[i,1:2])))
  return(c(b,as.numeric(dis_genet[i,4]) ))
})
d = matrix(unlist(d), byrow = T,ncol = 3)
e = apply(d[,1:2],1,function(i){paste(i , collapse = '_')})
d = d[which(e %in% apply(c[,1:2],1,function(i){paste(i , collapse = '_')})),]
c = c[which(apply(c[,1:2],1,function(i){paste(i , collapse = '_')}) %in% apply(d[,1:2],1,function(i){paste(i , collapse = '_')})),]
d = d[!duplicated(c[,1:2]), ]
d= cbind(apply(d[,1:2],1,function(i){paste(i , collapse = '_')}),d[,3])
c = cbind(apply(c[,1:2],1,function(i){paste(i , collapse = '_')}),c[,3])
c= as.data.frame(c)
d= as.data.frame(d)
e =left_join(c, d,by= 'V1')
e$V2.x= as.numeric(e$V2.x)
e$V2.y=as.numeric(e$V2.y)
ggplot(e,aes(x= V2.x, y=V2.y))+
  geom_point(alpha = 0.1)+
  stat_cor(method = 's')+
  theme_classic()+
  xlab("Transcriptome Pairwise Correlation")+
  ylab("Pairwise Genetic Distances")



a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]

a =cor(proteomic_data, method = 's',use = "complete.obs")
a = melt(a)
c = pblapply(1:nrow(a), function(i){
  b = sort(as.character(unlist(a[i,1:2])))
  return(c(b,as.numeric(a[i,3]) ))
})
c = matrix(unlist(c), byrow = T,ncol = 3)
c= c[!duplicated(c[,1:2]),]

d = pblapply(1:nrow(dis_genet), function(i){
  b = sort(as.character(unlist(dis_genet[i,1:2])))
  return(c(b,as.numeric(dis_genet[i,4]) ))
})
d = matrix(unlist(d), byrow = T,ncol = 3)
e = apply(d[,1:2],1,function(i){paste(i , collapse = '_')})
d = d[which(e %in% apply(c[,1:2],1,function(i){paste(i , collapse = '_')})),]
c = c[which(apply(c[,1:2],1,function(i){paste(i , collapse = '_')}) %in% apply(d[,1:2],1,function(i){paste(i , collapse = '_')})),]
d = d[!duplicated(c[,1:2]), ]
d= cbind(apply(d[,1:2],1,function(i){paste(i , collapse = '_')}),d[,3])
c = cbind(apply(c[,1:2],1,function(i){paste(i , collapse = '_')}),c[,3])
c= as.data.frame(c)
d= as.data.frame(d)
e =left_join(c, d,by= 'V1')
e$V2.x= as.numeric(e$V2.x)
e$V2.y=as.numeric(e$V2.y)
ggplot(e,aes(x= V2.x, y=V2.y))+
  geom_point(alpha = 0.1)+
  stat_cor(method = 's')+
  theme_classic()+
  xlab("Proteome Pairwise Correlation")+
  ylab("Pairwise Genetic Distances")


####fig3####


#####A #####



a = read.tree('data_upload/1011_matrix.tree.newick')
a = cophenetic.phylo(a)
a = melt(a)
b =unique(c(a$Var1,a$Var2))[!unique(c(a$Var1,a$Var2))%in%unique(colnames(df_final))]
a = a[!a$Var1%in% b,]
a = a[!a$Var2%in% b,]
a = dcast(a,Var1~Var2)
rownames(a)= a$Var1
a = a[,-1]
a =as.matrix(a)
a= njs(a)
e = lapply(unique(clades$new), function(i){
  i<<-i
  b= clades[clades$new==i,"Standardized name"]
  a$tip.label[a$tip.label%in%b]
})
names(e)= unique(clades$new)
names(e)[5]='No Clade'
e= e[!names(e)=='16.']
a = groupOTU(a,.node=e)
x = cols
names(x) = substr(names(x),1,3)
x= (x)[order(names(x))]
names(x)[1:9] = str_replace(names(x)[1:9],'0','')
names(x)[1:9] = str_replace(names(x)[1:9],'\\.','. ')

ggtree(a,layout="circular", size= 0.2)+
  geom_tippoint(aes(color=group), size = 0.5, alpha= 0.8)+
  scale_color_manual(values = x)+
  theme(legend.position = 'none')


#####B#####
b = df_final[,!duplicated(colnames(df_final))]
b = b[,!colnames(b)%in%"XTRA_DCZ"]

c = njs(dist(t(b)))

e = lapply(unique(clades$new), function(i){
  i<<-i
  b= clades[clades$new==i,"Standardized name"]
  c$tip.label[c$tip.label%in%b]
})
names(e)= unique(clades$new)
names(e)[5]='No Clade'
e= e[!names(e)=='16.']
c = groupOTU(c,.node=e)



ggtree(c,layout="circular", size= 0.2)+
  geom_tippoint(aes(color=group), size = 0.5, alpha= 0.8)+
  scale_color_manual(values = x)+
  theme(legend.position = 'none')

#####C#####

b = df_final[,duplicated(colnames(df_final))]
b = b[,!colnames(b)%in%"XTRA_DCZ"]

d = njs(dist(t(b)))

e = lapply(unique(clades$new), function(i){
  i<<-i
  b= clades[clades$new==i,"Standardized name"]
  d$tip.label[d$tip.label%in%b]
})
names(e)= unique(clades$new)
names(e)[5]='No Clade'
e= e[!names(e)=='16.']
d = groupOTU(d,.node=e)



ggtree(d,layout="circular", size= 0.2)+
  geom_tippoint(aes(color=group), size = 0.5, alpha= 0.8)+
  scale_color_manual(values = x)+
  theme(legend.position = 'none')

#####D#####

net = blockwiseModules(t(df_final[,duplicated(colnames(df_final))]), power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "data_upload/Prot_TOM",
                       verbose = 3)
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

load('data_upload/Prot_tom-block.1.RData')
adj <- TOM
adj[adj > 0.05] = 1
adj[adj != 1] = 0



network <- graph.adjacency(adj)
network <- igraph::simplify(network) # removes self-loops

V(network)$color <- labels2colors(net$colors)
V(network)$name <- names(net$colors)
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- igraph::delete_vertices(network, igraph::degree(network)==0)

network = intergraph::asNetwork(network)

network%v%'colour'= unlist(lapply(network$val,function(i)return(unlist(i[3]))))

remove(p)
p = ggnet2(network)
p=ggplot_build(p)
p$data[[2]]$colour <- unlist(lapply(network$val,function(i)return(unlist(i[3]))))
p$data[[2]]$alpha=0.7

a = unique(p$data[[2]]$colour)
a = setNames(a, a)
names(a)[a%in%c('black','turquoise','blue','brown','green')]=c('turquoise','black','brown','blue','green')
a  = stack(a)
rownames(a)=  a$values
p$data[[2]]$colour= a[p$data[[2]]$colour,"ind"]
p$data[[2]]$size=5
p <- ggplot_gtable(p)
plot(p)

#####E#####

net = blockwiseModules(t(df_final[,!duplicated(colnames(df_final))]), power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "data_upload/RNA_tom",
                       verbose = 3)
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

load('data_upload/RNA_tom-block.1.RData')

adj <- TOM
adj[adj > 0.05] = 1
adj[adj != 1] = 0

network <- graph.adjacency(adj)
network <- igraph::simplify(network)  # removes self-loops

V(network)$color <- labels2colors(net$colors)
V(network)$name <- names(net$colors)
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- igraph::delete_vertices(network, igraph::degree(network)==0)


network = intergraph::asNetwork(network)

network%v%'colour'= unlist(lapply(network$val,function(i)return(unlist(i[3]))))

remove(p)
p = ggnet2(network)
p=ggplot_build(p)
p$data[[2]]$colour <- unlist(lapply(network$val,function(i)return(unlist(i[3]))))
p$data[[2]]$alpha=0.7
p$data[[2]]$size=5
p <- ggplot_gtable(p)
plot(p)

####S6####


net = blockwiseModules(t(df_final[,duplicated(colnames(df_final))]), power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "data_upload/Prot_TOM",
                       verbose = 3)

# open a graphics window

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

a = stack(net$colors)

a = cbind(a,mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

cem <- cemitool(as.data.frame((df_final[,duplicated(colnames(df_final))])))

cem@selected_genes= rownames(df_final)
a= cbind(unique(mergedColors),c('M1','Not.Correlated','M2','M3','M4','M5','M6','M7'))
rownames(a)= a[,1]
b= cbind(rownames(df_final), a[mergedColors,2])
colnames(b)=c('genes','modules')
rownames(b)=1:nrow(b)
b= as.data.frame(b)
cem@module = (b)


a = GO2GeneID
a= lapply(names(GO2GeneID), function(i){
  i<<-i
  b = GO2GeneID[[i]]
  cbind(b, rep(i, length(b)))
})  

a= do.call(rbind,a)
colnames(a)= c('gene', 'term')
a= a[,2:1]
a = as.data.frame(a) 

a$term=goterms[a$term]
#a=a[a$gene%in%rownames(df_final),]
cem <- mod_ora(cem, a)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")




b=cem@ora
b = apply(b,2,function(i)str_replace(i,'/', './'))
write_clip(b)
for(i in 1:length(plots)){
  i<<-i 
  j=i
  i = plots[[i]]
  i=i$pl
  i <- ggplot_build(i)
  i$data[[1]]$fill = 'grey'
  i$data[[1]]$colour = "black"
  i$plot$theme$legend.position= 'none'
  i= ggplot_gtable(i)
  i = as_ggplot(i)
  assign(x = paste0('p_',j),value =i)
}
dev.off()
ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8, ncol = 2, nrow = 4)



####S7####

net = blockwiseModules(t(df_final[,!duplicated(colnames(df_final))]), power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "data_upload/RNA_tom",
                       verbose = 3)

# open a graphics window

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

a = stack(net$colors)

a = cbind(a,mergedColors)
temp = a
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

cem <- cemitool(as.data.frame((df_final[,!duplicated(colnames(df_final))])))

cem@selected_genes= rownames(df_final)
a= cbind(unique(mergedColors),c('M1','M2','Not.Correlated','M3','M4','M5'))
rownames(a)= a[,1]
b= cbind(rownames(df_final), a[mergedColors,2])
colnames(b)=c('genes','modules')
rownames(b)=1:nrow(b)
b= as.data.frame(b)
cem@module = (b)


a = GO2GeneID
a= lapply(names(GO2GeneID), function(i){
  i<<-i
  b = GO2GeneID[[i]]
  cbind(b, rep(i, length(b)))
})  

a= do.call(rbind,a)
colnames(a)= c('gene', 'term')
a= a[,2:1]
a = as.data.frame(a) 

a$term=goterms[a$term]
a=a[a$gene%in%rownames(df_final),]
cem <- mod_ora(cem, a)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
b=cem@ora
b = apply(b,2,function(i)str_replace(i,'/', './'))
for(i in 1:length(plots)){
  i<<-i 
  j=i
  i = plots[[i]]
  i=i$pl
  i <- ggplot_build(i)
  i$data[[1]]$fill = 'grey'
  i$data[[1]]$colour = "black"
  i$plot$theme$legend.position= 'none'
  i= ggplot_gtable(i)
  i = as_ggplot(i)
  assign(x = paste0('p_',j),value =i)
}

dev.off()
ggarrange(p_1,p_2,p_3,p_4,p_5,p_6, ncol = 2, nrow = 3)

####S8####

a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))

n= df_final[,duplicated(colnames(df_final))]
n= as.data.frame(n)
c = melt(n)
c$genes<-rep(rownames(n), ncol(n))
clades$new = substr(clades$Clades,1,3)
a= pblapply(unique(clades$new), function(i){
  if(i !=''){
    cat(c(i,' '))
    a = clades[clades$new==i,"Standardized name"]
    b = clades[clades$new!=i,"Standardized name"]
    a = intersect(a, (c$variable))
    if(length(a) >1){
      b = intersect(b, (c$variable))
      a = c[c$variable%in%a,]
      a$type = i
      b = c[c$variable%in%b,]
      b$type = 'other'
      d=lapply(rownames(n),function(j){
        x = rbind(filter(a, genes ==j), filter(b, genes == j))
        if(length(na.omit(x[x$type==i,"value"]))>0){
          return(c(log2(mean(na.omit(x[x$type==i,"value"]))/mean(na.omit(x[x$type=='other',"value"]))), wilcox.test(x[x$type==i,"value"],x[x$type=='other',"value"])$p.value ,j))}
        
      })
      d = do.call(rbind,d)
      d=as.data.frame(d)
      d$V1= as.numeric(d$V1)
      d$fdr=p.adjust(d$V2, method = 'b')
      #d$log_fdr = -log10(p.adjust(d$V2, method = 'b'))
      d$sign = NA
      d[d$fdr<0.05, 'sign']= 'Diff. expressed'
      d$stand = gene_info[d$V3,'Gene']
      
      d$log_p_vlaue = -log10(as.numeric(d$V2))
      d$clade = i
      return(d)
    }
  }
}
)
temp  =do.call(rbind,a)
temp$type = ifelse(temp$V1<0, 'under','over')
temp[is.na(temp$sign),"type"]=NA

e= lapply(unique(temp$clade),function(i){
  b= temp[temp$clade==i,]
  p = ggplot(b, aes(x=V1, y = log_p_vlaue,col=type))+
    geom_point(alpha=0.4, size = 0.5)+
    theme_bw()+
    scale_color_manual(values = c('red', 'blue'), na.value = 'black')+
    theme(legend.position = 'none')+
    ylab('')+
    facet_grid(. ~ clade)+
    xlab('')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
})
do.call("grid.arrange", c(e, ncol=4))

####fig4####

#####A#####
a = readxl::read_xlsx('data_upload/TableS5.xlsx') #table s5 from caudal et al
a = as.data.frame(a)
temp = lapply(unique(a$Module), function(i){
  i<<-i
  b = a[a$Module==i,]
  b=b$geneID
  b= paste(b, collapse = '/')
  b= strsplit(b,'/')
  b=unique(unlist(b))
})
names(temp)= unique(a$Module)

fgseaRes <- fgsea(examplePathways, exampleRanks, maxSize=500)

fgseaRes$clade = NA
annot_clade_prot_2= fgseaRes[-c(1:nrow(fgseaRes)),]

a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))

n= df_final[,duplicated(colnames(df_final))]
n= as.data.frame(n)
rownames(n) = b
c = melt(n)
c$genes<-rep(rownames(n), ncol(n))
clades$new = substr(clades$Clades,1,3)
pblapply(unique(clades$new), function(i){
  if(i !=''){
    cat(c(i,' '))
    a = clades[clades$new==i,"Standardized name"]
    b = clades[clades$new!=i,"Standardized name"]
    a = intersect(a, (c$variable))
    if(length(a) >1){
      b = intersect(b, (c$variable))
      a = c[c$variable%in%a,]
      a$type = i
      b = c[c$variable%in%b,]
      b$type = 'other'
      d = matrix(ncol = 3)
      d =d[-1,]
      for(j in rownames(n)){
        x = rbind(filter(a, genes ==j), filter(b, genes == j))
        if(length(na.omit(x[x$type==i,"value"]))>0){
          d= rbind(d, c(log2(mean(na.omit(x[x$type==i,"value"]))/mean(na.omit(x[x$type=='other',"value"]))), wilcox.test(x[x$type==i,"value"],x[x$type=='other',"value"])$p.value ,j))}
      }
      d = as.data.frame(d)
      d$V1= as.numeric(d$V1)
      d$fdr=p.adjust(d$V2, method = 'b')
      d$log_fdr = -log10(p.adjust(d$V2, method = 'b'))
      d$sign = NA
      d[d$fdr<0.05, 'sign']= 'Diff. expressed'
      d$stand = gene_info[d$V3,'Gene']
      
      d$log_p_vlaue = -log10(as.numeric(d$V2))
      d=d[d$stand%in%unlist(temp),]
      e = d$V1
      e=  setNames(e,d$stand)
      e= sort(e)
      fgseaRes <- fgsea(temp, e, maxSize=500)
      fgseaRes$clade = str_replace(i,' ','')
      annot_clade_prot_2 <<- rbind(annot_clade_prot_2, fgseaRes)
    }
  }
}
)

annot_clade_prot_2$logpval = -log10(annot_clade_prot_2$padj)

annot_clade_prot_2[annot_clade_prot_2$padj>0.05,'logpval'] = NA

ggplot(annot_clade_prot_2, aes(pathway, y = clade, size = logpval, col = NES))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_gradient(low = "#005AB5", high = "#DC3220", na.value = NA,limits=c(-4,4))+
  ylab('Subpopulation')+
  xlab('Pathway')+
  scale_size_continuous(limits=c(1,20),breaks=c(5,10,15,20))+
  labs(size='-Log10(p-value)')


#####B#####


a = readxl::read_xlsx('data_upload/TableS5.xlsx')  #table s5 from caudal et al
a = as.data.frame(a)
x = lapply(unique(a$Module), function(i){
  i<<-i
  b = a[a$Module==i,]
  b=b$geneID
  b= paste(b, collapse = '/')
  b= strsplit(b,'/')
  b=unique(unlist(b))
})
names(x)= unique(a$Module)



c = df_final[,duplicated(colnames(df_final))]

b =  pblapply(rownames(c), function(i){
  a = wilcox.test(as.numeric(c[i,intersect(colnames(c),Domestication[Domestication$Domestication=='Domesticated',1])]),as.numeric(c[i,intersect(colnames(c),Domestication[Domestication$Domestication=='Wild',1])]))
  
  a$p.value
  return(c(a$p.value, log2(mean(na.omit(as.numeric(c[i,intersect(colnames(c),Domestication[Domestication$Domestication=='Domesticated',1])])))/mean(na.omit(as.numeric(c[i,intersect(colnames(c),Domestication[Domestication$Domestication=='Wild',1])]))))))
})
b=t(as.data.frame(b))
rownames(b)= rownames(c)
# neg value mean lowly expressed in domesticated
b= as.data.frame(b)
b$log10 = -log10(b$V1)
b$sign= ifelse(p.adjust(b$V1,method = 'bonferroni')<0.05,'Sign', NA)
b$type = ifelse(b$V2>0 ,'over','under')
b[is.na(b$sign),"type"]= NA
b$names = gene_info[rownames(b),"Gene"]
b[is.na(b$sign),"names"]= NA
b$names_2 = b$names
b$names_2[b$log10<10]=NA
b$stand = gene_info[rownames(b),'Gene']
e = setNames(b$V2, nm=b$stand)
e=sort(e)
a <- fgsea(pathways = x, 
           stats    = e,
           minSize  = 5,
           maxSize  = 500)

a=a[a$padj<0.05,]
p = unlist(a[a$pathway== 'respiratory electron transport chain',leadingEdge])
b$res = b$stand%in%p
b[b$res,"res"]='Respiration'
b[b$stand%in%unlist(a[a$pathway== 'chaperone mediated protein folding',leadingEdge]),"res"]='Chaperon mediated folding'
b[b$res=='FALSE',"res"]=NA
b= b[order(b$res,na.last = F),]
ggplot(b,aes(V2, log10,col= res))+
  geom_point( )+
  scale_x_continuous(limits = c(-0.15,0.15))+
  theme_classic()+
  scale_color_manual(values = c('#FD8A8A','#9EA1D4'),na.value = 'lightgrey',name=element_blank())+
  ylab('-log10(p-value)')+
  xlab('Log2(FC)')

####S9####


b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))

n= df_final[,duplicated(colnames(df_final))]
n= as.data.frame(n)
rownames(n) = b
c = melt(n)
c$genes<-rep(rownames(n), ncol(n))
clades$new = substr(clades$Clades,1,3)
e=pblapply(unique(clades$new), function(i){
  if(i !=''){
    cat(c(i,' '))
    a = clades[clades$new==i,"Standardized name"]
    b = clades[clades$new!=i,"Standardized name"]
    a = intersect(a, (c$variable))
    if(length(a) >1){
      b = intersect(b, (c$variable))
      a = c[c$variable%in%a,]
      a$type = i
      b = c[c$variable%in%b,]
      b$type = 'other'
      d = matrix(ncol = 3)
      d =d[-1,]
      for(j in rownames(n)){
        x = rbind(filter(a, genes ==j), filter(b, genes == j))
        if(length(na.omit(x[x$type==i,"value"]))>0){
          d= rbind(d, c(log2(mean(na.omit(x[x$type==i,"value"]))/mean(na.omit(x[x$type=='other',"value"]))), wilcox.test(x[x$type==i,"value"],x[x$type=='other',"value"])$p.value ,j))}
      }
      d = as.data.frame(d)
      d$V1= as.numeric(d$V1)
      d$fdr=p.adjust(d$V2, method = 'b')
      d$log_fdr = -log10(p.adjust(d$V2, method = 'b'))
      d$sign = NA
      d[d$fdr<0.05, 'sign']= 'Diff. expressed'
      d$stand = gene_info[d$V3,'Gene']
      
      d$log_p_vlaue = -log10(as.numeric(d$V2))
      d = na.omit(d)
      return(nrow(d))
    }
  }
}
)
e=setNames(e,nm = unique(clades$new))
e=unlist(e)
f=table(clades[colnames(df_final[,!duplicated(colnames(df_final))]),"new"])
f=as.matrix(f)
f=as.data.frame(f)
e = data.frame(e)
f$signa= NA
f[intersect(rownames(e),rownames(f)),"signa"]=e[intersect(rownames(e),rownames(f)),"e"]
colnames(f)=c('Number isolate','Number of signature')
f$clade = rownames(f)
f$clade = substr(f$clade,1,2)
f = f[order(f$`Number isolate`,decreasing = T),]
f$clade= factor(f$clade,levels = f$clade)
f=melt(f)
f = na.omit(f)
f=f[!f$clade=='',]
ggplot(f, aes(x= clade, y = value, fill = variable))+
  geom_bar(stat = 'identity',color="black", position=position_dodge(), width =0.8)+
  theme_classic()+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  xlab("Subpopulation")+
  ylab("Signature genes & isolate number")+
  labs(fill='')







b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))

n= df_final[,!duplicated(colnames(df_final))]
n= as.data.frame(n)
rownames(n) = b
c = melt(n)
c$genes<-rep(rownames(n), ncol(n))
clades$new = substr(clades$Clades,1,3)
e=pblapply(unique(clades$new), function(i){
  if(i !=''){
    cat(c(i,' '))
    a = clades[clades$new==i,"Standardized name"]
    b = clades[clades$new!=i,"Standardized name"]
    a = intersect(a, (c$variable))
    if(length(a) >1){
      b = intersect(b, (c$variable))
      a = c[c$variable%in%a,]
      a$type = i
      b = c[c$variable%in%b,]
      b$type = 'other'
      d = matrix(ncol = 3)
      d =d[-1,]
      for(j in rownames(n)){
        x = rbind(filter(a, genes ==j), filter(b, genes == j))
        if(length(na.omit(x[x$type==i,"value"]))>0){
          d= rbind(d, c(log2(mean(na.omit(x[x$type==i,"value"]))/mean(na.omit(x[x$type=='other',"value"]))), wilcox.test(x[x$type==i,"value"],x[x$type=='other',"value"])$p.value ,j))}
      }
      d = as.data.frame(d)
      d$V1= as.numeric(d$V1)
      d$fdr=p.adjust(d$V2, method = 'b')
      d$log_fdr = -log10(p.adjust(d$V2, method = 'b'))
      d$sign = NA
      d[d$fdr<0.05, 'sign']= 'Diff. expressed'
      d$stand = gene_info[d$V3,'Gene']
      
      d$log_p_vlaue = -log10(as.numeric(d$V2))
      d = na.omit(d)
      return(nrow(d))
    }
  }
}
)
e=setNames(e,nm = unique(clades$new))
e=unlist(e)
f=table(clades[colnames(df_final[,!duplicated(colnames(df_final))]),"new"])
f=as.matrix(f)
f=as.data.frame(f)
e = data.frame(e)
f$signa= NA
f[intersect(rownames(e),rownames(f)),"signa"]=e[intersect(rownames(e),rownames(f)),"e"]
colnames(f)=c('Number isolate','Number of signature')
f$clade = rownames(f)
f = f[order(f$`Number isolate`,decreasing = T),]
f$clade= factor(f$clade,levels = f$clade)
f=melt(f)
f = na.omit(f)
f=f[!f$clade=='',]
ggplot(f, aes(x= clade, y = value, fill = variable))+
  geom_bar(stat = 'identity',color="black", position=position_dodge(), width =0.8)+
  theme_classic()+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  xlab("Subpopulation")+
  ylab("Signature genes & isolate number")+
  labs(fill='')



####fig5####

#####A#####

ggplot(CI_loc_snp_455_Prot_filt, aes(pos_cum_snp,pos_cum_pheno))+
  geom_point(size=0.3, alpha = 0.5)+
  theme_classic()+
  geom_vline(xintercept = position_chr$start_cumul,alpha = 0.1)+
  geom_hline(yintercept = position_chr$start_cumul,alpha = 0.1)+
  #ggtitle("RNA")+
  scale_x_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  scale_y_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  xlab('pQTL position')+
  ylab('Gene position')+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),legend.position = 'none')

#####D#####
ggplot(CI_loc_snp_455_RNA_filt, aes(pos_cum_snp,pos_cum_pheno))+
  geom_point(size=0.3, alpha = 0.5)+
  theme_classic()+
  geom_vline(xintercept = position_chr$start_cumul,alpha = 0.1)+
  geom_hline(yintercept = position_chr$start_cumul,alpha = 0.1)+
  #ggtitle("RNA")+
  scale_x_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  scale_y_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  xlab('eQTL position')+
  ylab('Gene position')+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),legend.position = 'none')

#####C#####
b= as.data.frame(CI_loc_snp_455_RNA[CI_loc_snp_455_RNA$ld_mask=='',])
a= as.data.frame(CI_loc_snp_455_Prot[CI_loc_snp_455_Prot$ld_mask=='',])



a = a[,c("pos_cum_snp","1")]
a= a[!duplicated(a),]
c = cbind(seq(0,max(a$pos_cum_snp)+20000,20000),seq(0,max(a$pos_cum_snp)+20000,20000)+20000)
d =pbapply(c,1,function(i){
  
  length(unique(a[which(between(x = a$pos_cum_snp,lower = i[1],upper = i[2])),"1"]))
})

d = cbind(c,d)
d= as.data.frame(d)
colnames(d)=c('start','stop','count')


b = b[,c("pos_cum_snp","1")]
b= b[!duplicated(b),]
c = cbind(seq(0,max(b$pos_cum_snp)+20000,20000),seq(0,max(b$pos_cum_snp)+20000,20000)+20000)

e =pbapply(c,1,function(i){
  length(unique(b[which(between(x = b$pos_cum_snp,lower = i[1],upper = i[2])),"1"]))
})

e = cbind(c,e)
e= as.data.frame(e)
colnames(e)=c('start','stop','count')

d$type= 'Prot'
e$type = 'RNA'

i = rbind(d,e)
i[i$type=='RNA',"count"]= 0-(i[i$type=='RNA',"count"])

i$per=i$count/nrow(proteomic_data)*100
i$Dataset = ifelse(i$count>0, "Proteomic",'Transcriptomic')
i[i$count<=3&i$count>=-3,"Dataset"]=NA

a = i[i$count>3|i$count<(-3),]
p = ggplot(a,aes(x= start,y=per))+
  geom_bar(stat = 'identity')+
  geom_vline(xintercept = position_chr$start_cumul,alpha = 0.1,linetype ='dashed')+
  theme_classic()+
  ylab( '% of genes')+
  #scale_fill_manual(values = c('gold','darkgrey'))+
  geom_hline(yintercept = 0)+
  xlab('Chromosome')+
  scale_y_continuous(breaks = c(-12:5), labels = abs(c(-12:5)))+
  scale_x_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  annotate(geom = 'text', x = 11023260, y = 1.4,label='pQTL', color='#EFC000FF', size = 5 )+
  annotate(geom = 'text', x = 11023260, y = -1.4,label='eQTL', color='#0073C2FF' , size= 5)

p

table(paste(a$start,a$stop))

#####B#####
a= CI_loc_snp_455_Prot_filt
a$title = 'SNP-pQTL'
p1 = ggplot(a, aes(x = type,y = EffectSize))+
  geom_violin( alpha = 0.7)+
  geom_boxplot(width=0.1)+
  stat_compare_means(label.x = 1.5,label = 'p.signif')+
  theme_bw()+
  xlab('')+
  scale_fill_manual(values = c('#EFC000FF','#FFE991'))+
  ylab('SNP-pQTL effect size')+
  theme(legend.position = 'none')+
  scale_y_log10(lim=c(0.001,1))+
  scale_x_discrete(labels=c('Local','Distant'))+
  facet_grid(. ~ title)

a= control_GWAS_correlated_CNV[control_GWAS_correlated_CNV$`2`=='Prot',]
a$title = 'CNV-pQTL'
p2= ggplot(a, aes(x = type,y = EffectSize))+
  geom_violin( alpha = 0.7)+
  geom_boxplot(width=0.1)+
  stat_compare_means(label.x = 1.5,label = 'p.signif')+
  theme_bw()+
  xlab('')+
  scale_fill_manual(values = c('#EFC000FF','#FFE991'))+
  ylab('CNV-pQTL effect size')+
  theme(legend.position = 'none')+
  scale_y_log10(lim=c(0.001,1))+
  scale_x_discrete(labels=c('Local','Distant'))+
  facet_grid(. ~ title)

grid.arrange(p1,p2, ncol=2)


####S10####

p1=ggplot(control_GWAS_correlated_CNV[control_GWAS_correlated_CNV$`2`=='RNA',], aes(pos_cum_snp,pos_cum_pheno))+
  geom_point(size=0.3, alpha = 0.5)+
  theme_classic()+
  geom_vline(xintercept = position_chr$start_cumul,alpha = 0.1)+
  geom_hline(yintercept = position_chr$start_cumul,alpha = 0.1)+
  #ggtitle("RNA")+
  scale_x_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  scale_y_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  xlab('CNV-eQTL position')+
  ylab('Gene position')+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),legend.position = 'none')+
  scale_color_manual(values = c('red','black'))

p2=ggplot(control_GWAS_correlated_CNV[control_GWAS_correlated_CNV$`2`=='Prot',], aes(pos_cum_snp,pos_cum_pheno))+
  geom_point(size=0.3, alpha = 0.5)+
  theme_classic()+
  geom_vline(xintercept = position_chr$start_cumul,alpha = 0.1)+
  geom_hline(yintercept = position_chr$start_cumul,alpha = 0.1)+
  scale_x_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  scale_y_continuous(breaks = apply(position_chr[,11:12],1,mean),labels = names(apply(position_chr[,11:12],1,mean)))+
  xlab('CNV-pQTL position')+
  ylab('Gene position')+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),legend.position = 'none')+
  scale_color_manual(values = c('red','black'))


####S11####
a =CI_loc_snp_455_Prot_filt[CI_loc_snp_455_Prot_filt$type=='CIS',]
a$test= a$pos_cum_snp

a$sens=ifelse(substr(a$`1`, 7,7)=='W',yes ='+',no='-')

a$SNP_trait_relation = 't'
a$start = gene_info[a$`1`,"start"] 
a$stop = gene_info[a$`1`,"stop"]

a[between(a$ChrPos,a$start,a$stop),"SNP_trait_relation"]='ORF'


a$dist_start = 77777777
a[a$sens=='+',"dist_start"]= a[a$sens=='+',"ChrPos"]-a[a$sens=='+',"start"]
a[a$sens=='-',"dist_start"]= 0-(a[a$sens=='-',"ChrPos"]-a[a$sens=='-',"stop"])
a$gene_length = a$stop-a$start

a[a$dist_start<(-1000),"SNP_trait_relation"]='Background'
a[a$dist_start<0&a$SNP_trait_relation!= 'Background',"SNP_trait_relation"]= 'UpstreamTSS'

a[between(a$dist_start,a$gene_length,a$gene_length+200),"SNP_trait_relation"]='DownstreamTES'
a[between(a$dist_start,a$gene_length+200, rep(max(a$"dist_start")+1,nrow(a))),"SNP_trait_relation"] = 'Background'
ggplot(a,aes(x = dist_start, fill= SNP_trait_relation))+
  geom_histogram(bins =200)+
  theme_classic()+
  scale_fill_manual(values= c("ORF"="#F9B16F","DownstreamTES"="#78acd1", "UpstreamTSS"= "#FF6633", "Background" = "#CCCCCC"))+
  xlab('Distance to the start codon')


####S12####


a = control_GWAS_correlated_CNV[control_GWAS_correlated_CNV$`2`=='Prot',]

a$aneuploidy=(a$Pheno_chr==a$Chr)&(a$Chr%in%c(1,3,8,9,11))
ggplot(a,aes(aneuploidy, EffectSize))+
  geom_violin()+
  geom_boxplot(width = 0.1)+
  stat_compare_means(label.x = 1.33,label = 'p.signif')+
  theme_classic()+
  ylab('Effect Size')+
  xlab('CNV-pQTL type')+
  scale_x_discrete(label=c('Other CNV-pQTL', 'Aneuploidy related CNV-pQTL'))




table(a$Chr==a$Pheno_chr)
a$shared=a$unique%in%names(table(control_GWAS_correlated_CNV$unique)[table(control_GWAS_correlated_CNV$unique)==2])

sum(a$shared&a$aneuploidy)


a = control_GWAS_correlated_CNV[control_GWAS_correlated_CNV$`2`=='RNA',]

a$aneuploidy=(a$Pheno_chr==a$Chr)&(a$Chr%in%c(1,3,8,9,11))
table(a$Chr==a$Pheno_chr)
a$shared=a$unique%in%names(table(control_GWAS_correlated_CNV$unique)[table(control_GWAS_correlated_CNV$unique)==2])

sum(a$shared&a$aneuploidy)

####S13####

x=rownames(df_final)

b = lapply(x, function(i){
  if(i!='YPL260W'){i<<-i
  a = df_final[i,]
  cor.test(a[duplicated(names(a))],a[!duplicated(names(a))],method = 's')$estimate}
})
names(b)= x
temp=intersect(CI_loc_snp_455_RNA_filt$unique,CI_loc_snp_455_Prot_filt$unique)
temp = CI_loc_snp_455_RNA_filt[CI_loc_snp_455_RNA_filt$unique%in%temp,`1`]
x = stack(b)
x$Overlapping_QTL = 'No'
x[x$ind%in%temp,"Overlapping_QTL"]='Yes'

ggplot(x,aes(x=values, fill = Overlapping_QTL))+
  geom_histogram(col='black')+
  theme_classic()+
  xlab('Rho')+
  scale_fill_manual(name='Overlap \nin SNP-QTL',values = c('grey','pink'))


####S14####


a = intersect(colnames(proteomic_data),colnames(RNA_data_jing))
b = intersect(rownames(proteomic_data), rownames(RNA_data_jing))
c = df_final[b,duplicated(colnames(df_final))]
d= df_final[b,!duplicated(colnames(df_final))]



e= lapply(b, function(i){
  e = cor.test(as.numeric(c[i,]),as.numeric(d[i,]), method = 's')
  return(c(as.numeric(e$estimate),as.numeric(e$p.value)))
})


names(e)=b
e=do.call(rbind,e)
e = as.data.frame(e)
e$adj = p.adjust(e$V2,method = 'b')
e$sign = ifelse(e$adj<0.05, 'Significant', 'Not significant')
e[e$V1>0.42,"sign"]= 'Strongly correlated'



x2 = aggregate(data=turnover_data2, FUN = mean , kdp_value~systematic_name)
rownames(x2)=x2$systematic_name
temp=intersect(CI_loc_snp_455_RNA_filt$unique,CI_loc_snp_455_Prot_filt$unique)
x = CI_loc_snp_455_RNA_filt[CI_loc_snp_455_RNA_filt$unique%in%temp,`1`]
e$turnover = x2[rownames(e),"kdp_value"]
e$overlap = rownames(e)%in%x
e =e[!is.na(e$turnover),]
ggplot(e, aes(overlap,turnover))+
  geom_violin()+
  geom_boxplot(width = 0.1)+
  stat_compare_means(label.x = 1.5, label = 'p.signif')+
  ylab('Mean turnover rate')+
  theme_classic()+
  xlab('Overlapping SNP-QTL')
ggplot(e, aes(overlap,turnover))+
  geom_violin()+
  geom_boxplot(width = 0.1)+
  stat_compare_means(label.x = 1.5)+
  ylab('Mean turnover rate')+
  theme_classic()+
  xlab('Overlapping SNP-QTL')




x2 = aggregate(data=turnover_data2, FUN = mean , HL~systematic_name)
rownames(x2)=x2$systematic_name
temp=intersect(CI_loc_snp_455_RNA_filt$unique,CI_loc_snp_455_Prot_filt$unique)
x = CI_loc_snp_455_RNA_filt[CI_loc_snp_455_RNA_filt$unique%in%temp,`1`]
e$turnover = x2[rownames(e),"HL"]
e$overlap = rownames(e)%in%x
e =e[!is.na(e$turnover),]
ggplot(e, aes(overlap,turnover))+
  geom_violin()+
  geom_boxplot(width = 0.1)+
  stat_compare_means(label.x = 1.5, label = 'p.signif')+
  ylab('Mean protein half-life')+
  theme_classic()+
  xlab('Overlapping SNP-QTL')
ggplot(e, aes(overlap,turnover))+
  geom_violin()+
  geom_boxplot(width = 0.1)+
  stat_compare_means(label.x = 1.5)+
  ylab('Mean protein half-life')+
  theme_classic()+
  xlab('Overlapping SNP-QTL')

table(e$overlap)
