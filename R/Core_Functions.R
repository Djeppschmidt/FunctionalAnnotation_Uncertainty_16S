
#' Use gene annotations to subset phyloseq object
#'
#' This function is not finished
#' @param ps phyloseq object
#' @param reftab reference table built from get.genes()
#' @param genes genes to use in the subset
#' @param either logical, if T then returns taxa that have at least one gene in the list; if F then returns only taxa that have both
#' @param level one of "Genus" or "Species" to select from dataset
#' @param cutoff threshold - proportion of assemblies that meet genes criteria for a taxon to keep in the reference taxon. 1 represents all assemblies meet the gene critera
#' @keywords 
#' @export
#' @examples
#' subset.phyloseq()
subset.Phyloseq<-function(ps, reftab, genes, either=T, level="Species", cutoff=1){
  tax<-as.data.frame(as.matrix(tax_table(ps)))
  # add columns of reftab to tax by shared genus and species columns
  reftab$fullname<-paste(reftab$Genus, reftab$Species)
  # calculate proportion of assemblies for each taxa that have criteria
  
  
  if(either==T){
    # define reference list of names
    if(level=="Species"){ 
      
      # summarize reftable
      keepls<-table(reftab$fullname[rowSums(reftab[,colnames(reftab)%in%genes])>0])
      
      # establish list of all taxa
      totalspp<-table(reftab$fullname)
      # subset total to spp / genera with genes
      totalsppsubset<-totalspp[names(totalspp)%in%names(keepls)] 
      print(identical(names(totalsppsubset), names(keepls))) # sanity check
      
      average<-keepls/totalspp
      pdf<-data.frame("Average"=average, "Total"=totalsppsubset)
      # calculate the proportion of each genus that has AOA
      colnames(pdf)<-c("Names", "Average", "abundance")
      pm<-reshape2::melt(pdf, id.vars=c("Names"), measure.vars=c("Average"))
      pm$Names<-as.character(pm$Names)
      print("S reflist finished")
      }
    
    if(level=="Genus"){
      # genus
      # summarize reftable
      keepls<-table(reftab$Genus[rowSums(reftab[,colnames(reftab)%in%genes])>0])
      
      # establish list of all taxa
      totalspp<-table(reftab$Genus)
      # subset total to spp / genera with genes
      totalsppsubset<-totalspp[names(totalspp)%in%names(keepls)] 
      print(identical(names(totalsppsubset), names(keepls))) # sanity check
      
      average<-keepls/totalspp
      pdf<-data.frame("Average"=average, "Total"=totalsppsubset)
      # calculate the proportion of each genus that has AOA
      colnames(pdf)<-c("Names", "Average", "abundance")
      pm<-reshape2::melt(pdf, id.vars=c("Names"), measure.vars=c("Average"))
      pm$Names<-as.character(pm$Names)
      print("G reflist finished")
      
    }
    # prune taxa by reference list of names
  }
  if(either==F){
    # define reference list of names
    if(level=="Species"){ 
      
      # summarize reftable
      keepls<-table(reftab$fullname[rowSums(reftab[,colnames(reftab)%in%genes])==length(genes)])
      
      # establish list of all taxa
      totalspp<-table(reftab$fullname)
      # subset total to spp / genera with genes
      totalsppsubset<-totalspp[names(totalspp)%in%names(keepls)] 
      print(identical(names(totalsppsubset), names(keepls))) # sanity check
      
      average<-keepls/totalspp
      pdf<-data.frame("Average"=average, "Total"=totalsppsubset)
      # calculate the proportion of each genus that has AOA
      colnames(pdf)<-c("Names", "Average", "abundance")
      pm<-reshape2::melt(pdf, id.vars=c("Names"), measure.vars=c("Average"))
      pm$Names<-as.character(pm$Names)
      print("S reflist finished")
    }
    
    if(level=="Genus"){
      # genus
      # summarize reftable
      keepls<-table(reftab$Genus[rowSums(reftab[,colnames(reftab)%in%genes])==length(genes)])
      
      # establish list of all taxa
      totalspp<-table(reftab$Genus)
      # subset total to spp / genera with genes
      totalsppsubset<-totalspp[names(totalspp)%in%names(keepls)] 
      print(identical(names(totalsppsubset), names(keepls))) # sanity check
      
      average<-keepls/totalspp
      pdf<-data.frame("Average"=average, "Total"=totalsppsubset)
      # calculate the proportion of each genus with the given genes
      colnames(pdf)<-c("Names", "Average", "abundance")
      pm<-reshape2::melt(pdf, id.vars=c("Names"), measure.vars=c("Average"))
      pm$Names<-as.character(pm$Names)
      print("G reflist finished")
      
    }
  }
  if(level=="Species"){
    prune.list<-paste(tax$Genus, tax$Species, sep=" ")%in%pm$Names[pm$Average>=cutoff]
    ps.out<-prune_taxa(prune.list, ps)
    return(ps.out)
  }
  if(level=="Genus"){
    prune.list<-paste(tax$Genus, tax$Species, sep=" ")%in%pm$Names[pm$Average>=cutoff]
    ps.out<-prune_taxa(prune.list, ps)
    return(ps.out)
    }
}

#' create taxonomy vs rank-abundance curves
#'
#' 
#' @param ps phyloseq object
#' @keywords 
#' @export
#' @examples
#' RankAbundanceBurnoulliPlot()
RankAbundanceBurnoulliPlot<-function(ps){
  
  # make df of spp ordered by rank abundace, filled with sppID
  
  otu<-as.data.frame(t(as.matrix(otu_table(ps))))
  tax<-as.data.frame(as.matrix(tax_table(ps)))
  prop<-estimate_richness(ps, measures="Observed")
  #print(max(prop$Observed))
  prop3<-c(1:max(prop$Observed))
  for(i in 1:max(prop$Observed)){
    #print(sum(prop$Observed>i)/max(prop$Observed))
    prop3[i]<-sum(prop$Observed>i)/max(prop$Observed)
  }
  
  # make a matrix of names listed most to least abundant
  Rank.Matrix<-matrix(nrow=ntaxa(ps), ncol=nsamples(ps))
  for(i in 1:nsamples(ps)){
    Rank.Matrix[,i]<-rownames(otu)[order(otu[,i], decreasing=T)] # need to check that this is ordering by number !
  }
  #print(dim(Rank.Matrix))
  # make matrix with ifelse -> is.na () for tax table 
  Rank.MatrixSpp<-matrix(nrow=ntaxa(ps), ncol=nsamples(ps))
  Rank.MatrixG<-matrix(nrow=ntaxa(ps), ncol=nsamples(ps))
  Rank.MatrixF<-matrix(nrow=ntaxa(ps), ncol=nsamples(ps))
  Rank.MatrixO<-matrix(nrow=ntaxa(ps), ncol=nsamples(ps))
  
  # make taxonomy references
  Species<-ifelse(!is.na(tax$Species), 1, 0)
  names(Species)<-rownames(tax)
  #Species<-Species[Species==1]
  Genus<-ifelse(!is.na(tax$Genus), 1, 0)
  names(Genus)<-rownames(tax)
  #Genus<-Genus[Genus==1]
  Family<-ifelse(!is.na(tax$Family), 1, 0)
  names(Family)<-rownames(tax)
  #Family<-Family[Family==1]
  Order<-ifelse(!is.na(tax$Order), 1, 0)
  names(Order)<-rownames(tax)
  #Order<-Order[Order==1]
  
  # print(View(data.frame(Species,Genus,Family,Order)))
  
  for(i in 1:nsamples(ps)){
    Rank.MatrixSpp[,i]<-Species[Rank.Matrix[,i]]
  }
  # print(View(Rank.MatrixSpp))
  for(i in 1:nsamples(ps)){
    Rank.MatrixG[,i]<-Genus[Rank.Matrix[,i]]
  }
  
  for(i in 1:nsamples(ps)){
    Rank.MatrixF[,i]<-Family[Rank.Matrix[,i]]
  }
  
  for(i in 1:nsamples(ps)){
    Rank.MatrixO[,i]<-Order[Rank.Matrix[,i]]
  }
  
  # make everything zero that was not detected in a sample
  
  for(i in 1:ncol(Rank.MatrixSpp)){
    Rank.MatrixSpp[c(prop$Observed[i]:nrow(Rank.MatrixSpp)),i]<-0
  }
  # print(View(Rank.MatrixSpp))
  for(i in 1:ncol(Rank.MatrixG)){
    Rank.MatrixG[c(prop$Observed[i]:nrow(Rank.MatrixG)),i]<-0
  }
  for(i in 1:ncol(Rank.MatrixF)){
    Rank.MatrixF[c(prop$Observed[i]:nrow(Rank.MatrixF)),i]<-0
  }
  for(i in 1:ncol(Rank.MatrixO)){
    Rank.MatrixO[c(prop$Observed[i]:nrow(Rank.MatrixO)),i]<-0
  }
  
  
  datSp<-rowSums(Rank.MatrixSpp)/ncol(Rank.MatrixSpp)
  #print(datSp)
  #print(length(datSp))
  datG<-rowSums(Rank.MatrixG)/ncol(Rank.MatrixG)
  #print(length(datG))
  datF<-rowSums(Rank.MatrixF)/ncol(Rank.MatrixF)
  datO<-rowSums(Rank.MatrixO)/ncol(Rank.MatrixO)
  #print("Spp done")
  rollingmean<-function(x){
    a<-rep(NA, length(x))
    i=1
    while((i+50) < length(x)){
      a[i]<-mean(x[c(i:(i+50))])
      i=i+1}
    a
  }
  
  rollingSD<-function(x){
    a<-rep(NA, length(x))
    i=1
    while((i+50) < length(x)){
      a[i]<-sd(x[c(i:(i+50))])
      i=i+1}
    a
  }
  
  datS2<-rollingmean(datSp)
  datS.E<-rollingSD(datSp)
  # print(datS2)
  #print(max(datS2))
  datG2<-rollingmean(datG)
  datG.E<-rollingSD(datG)
  #print(length(datG))
  #print(length(datG.E))
  datF2<-rollingmean(datF)
  datF.E<-rollingSD(datF)
  
  datO2<-rollingmean(datO)
  datO.E<-rollingSD(datO)
  
  out<-NULL
  out$df<-data.frame("ref"= c(1:ntaxa(ps)), datS2, datS.E, datG2, datG.E, datF2, datF.E, datO2, datO.E)
  out$df<-out$df[c(1:max(prop)),]
  
  
  out$p1<-ggplot(out$df, aes(x=ref, y=datS2))+geom_line(col='black')+geom_ribbon(aes(ymin=datS2-datS.E, ymax=datS2+datS.E), alpha=0.1)+ggtitle("Species")+xlab("Rank Index")+ylab("Proportion")+theme_bw()
  out$p2<-ggplot(out$df, aes(x=ref, y=datG2))+geom_line(col='black')+geom_ribbon(aes(ymin=datG2-datG.E, ymax=datG2+datG.E), alpha=0.1)+ggtitle("Genus")+xlab("Rank Index")+ylab("Proportion")+theme_bw()
  out$p3<-ggplot(out$df, aes(x=ref, y=datF2))+geom_line(col='black')+geom_ribbon(aes(ymin=datF2-datF.E, ymax=datF2+datF.E), alpha=0.1)+ggtitle("Family")+xlab("Rank Index")+ylab("Proportion")+theme_bw()
  out$p4<-ggplot(out$df, aes(x=ref, y=datO2))+geom_line(col='black')+geom_ribbon(aes(ymin=datO2-datO.E, ymax=datO2+datO.E), alpha=0.1)+ggtitle("Order")+xlab("Rank Index")+ylab("Proportion")+theme_bw()
  
  #barplot(datS2, main= "Species", ylab="Proportion Known Species", xlab="Taxon Rank")
  #lines(prop3, col = "grey")
  
  #barplot(datG, main= "Genus", ylab="Proportion Known Genera", xlab="Taxon Rank")
  #barplot(datF, main= "Family", ylab="Proportion Known Family", xlab="Taxon Rank")
  #barplot(datO, main= "Order", ylab="Proportion Known Order", xlab="Taxon Rank")
  #barplot(prop3, main= "rank by samples", ylab="number of samples with given rank", xlab="Taxon Rank")
  out
}


#' this function unzips the annotation features (gene annotations) for each assembly in refseq and genbank
#'
#' 
#' @param refdir path to table of all complete prokaryote genomes from NCBI, with ftp links. Can be downloaded from https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/. For example: "~/Documents/prokaryotes.csv"
#' @param outpath = output directory
#' @param binPath = path to bin where wget is compiled.
#' @keywords 
#' @export
#' @examples
#' download.Feature.Tables()
download.Feature.Tables<-function(refdir, outpath, binPATH="/usr/local/bin/"){
  require(stringr)
  refdir<-as.data.frame(as.matrix(read.csv(refdir, sep=",")))
  setwd(outpath)
  # needs a sanity check for: make sure path to wget is good
  pattern<-"GCA_.*"
  bin<-Sys.getenv("PATH")
  Sys.setenv("PATH" = binPATH) # necessary to direct out or R bin / access BASH functions
  dir.create(paste(outpath, "RefSeq/", sep=""))
  dir.create(paste(outpath, "GenBank/", sep=""))
  
  #if(!dir.exists(paste(outpath, "RefSeq/", sep=""))){stop("Directory not made")}
  #if(!dir.exists(paste(outpath, "GenBank/", sep=""))){stop("Directory not made")}
  
  for (i in 1:length(refdir$GenBank.FTP)){
    system(paste("wget https:", substring(refdir$GenBank.FTP[i], 5), "/", str_match(refdir$GenBank.FTP[i], pattern = "GCA_.*"), "_feature_table.txt.gz", " -P ", outpath, "GenBank/", sep=""))
  }
  
  # do same for refseq
  pattern<-"GCF_.*"
  Sys.setenv("PATH" = "/usr/local/bin/")
  for (i in 1:length(refdir$RefSeq.FTP)){
    system(paste("wget https:", substring(refdir$RefSeq.FTP[i], 5), "/", str_match(refdir$RefSeq.FTP[i], pattern = "GCF_.*"), "_feature_table.txt.gz", " -P ", outpath, "RefSeq/", sep=""))
  } 
  Sys.setenv("PATH" = bin)
}


#' this function unzips the annotation features (gene annotations) for each assembly in refseq and genbank
#'
#' This is an internal subroutine for get.genes
#' @param x list of files
#' @param a 1 or 2 for GenBank or RefSeq
#' @keywords 
#' @export
#' @examples
#' unzip()
unzip<-function(x, a, directory){
  
  bin<-Sys.getenv("PATH")
  if(a==1){Sys.setenv("PATH" = "/usr/bin/")
    system(paste("gunzip ", directory, "GenBank/", x, sep=""))} # unzip the files
  if(a==2){Sys.setenv("PATH" = "/usr/bin/")
    system(paste("gunzip ", directory, "RefSeq/", x, sep=""))}
  
  Sys.setenv("PATH" = bin)
}


#' this function makes an index for each feature that indicates whether it has a genbank, refseq, both or neither annotation file. Values: 0 = neither; 1 = GenBank only; 2 = RefSeq only; 3 = both.
#'
#' This is an internal subroutine for get.genes
#' @param refdir table of NCBI prokaryote assemblies
#' @param path1 path to GenBank folder
#' @param path2 path to RefSeq folder
#' @keywords 
#' @export
#' @examples
#' compile.Flist()
compile.Flist<-function(refdir, path1, path2){
  files1<-list.files(path1)
  files2<-list.files(path2)
  L1<-rep(NA, nrow(refdir))
  for(i in 1:nrow(refdir)){
    L1[i]<-sum(ifelse(any(files2==paste(str_match(refdir$RefSeq.FTP[i], pattern = "GCF_.*"), "_feature_table.txt", sep="")), 2,0),
               ifelse(any(files1==paste(str_match(refdir$GenBank.FTP[i], pattern = "GCA_.*"), "_feature_table.txt", sep="")), 1,0))
  }
  L1
}


#' This function is a wrapper that implements the unzip, compile.Flist and compileFunctionTable functions
#'
#' 
#' @param refdir table of NCBI prokaryote assemblies
#' @param directory directory that contains the downloaded feature annotation tables
#' @param gene list of gene names, e.g. c("nifH", "nifD")
#' @keywords 
#' @export
#' @examples
#' get.genes()
get.genes<-function(refdir, directory, genes){
  # define paths to files
  path1<-paste(outpath, "GenBank/", sep="")
  path2<-paste(outpath, "RefSeq/", sep="")
  # get list of unzipped files
  files1<-list.files(paste(directory, "GenBank/", sep="")) 
  files2<-list.files(paste(directory, "RefSeq/", sep=""))
  
  # make an index for whether reference files exist for genbank and refseq. this will be used to determine which files to use (refseq first, then genbank)
  files<-compile.Flist(refdir, path1, path2) 
  refdir$index<-files # incorporate index into reference dataset
  
  # RefSeq assemblies have different file names from GenBank. This makes a column in refdir that matches the file names for the RefSeq assemblies
  refdir$assembly<-refdir$Assembly
  for(i in 1:length(refdir$Assembly)){
    refdir$assembly[i]<-paste("GCF_", str_split(refdir$RefSeq.FTP[i], pattern="_")[[1]][2], sep="")
  }
  
  # compile table
  outtab<-plyr::ldply(refdir$Assembly, compile.Functiontable, refdir, directory)
  outtab[,5:length(genes)]<-sapply(outtab[,5:length(genes)], as.numeric)
  colnames(outtab)<-c("assembly", "Genus", "Species", "Strain", names(genes))
  outtab
}


#' subroutine for compile.Functiontable
#'
#' Used to construct the index column to indicate whether GenBank or RefSeq annotation files exist
#' @param x gene symbol to search for
#' @param type whether assemblies for one or both databases exist. 1 = only one (RefSeq or GenBank); 2=both exist
#' @keywords 
#' @export
#' @examples
#' match.gene()
match.gene<-function(x, type, file1, file2=NULL){
  if(type==1){ifelse(any(str_detect(file1$symbol, pattern=x)), 1, 0)}
  
  if(type==2){ifelse(any(str_detect(file1$symbol, pattern=x))|any(str_detect(file2$symbol, pattern=x)), 1, 0)}
 
}


#' Compiles summary table of which assemblies have each gene. 
#'
#' Output table has columns for Genus, Species, Strain, Assembly ID, and each of the genes of interest
#' @param x element from refdir$Assembly. This is implemented within ldply
#' @param y will be the reference table with all metadata
#' @param directory path to files
#' @keywords 
#' @export
#' @examples
#' compile.Functiontable()
compile.Functiontable<-function(x,y,directory){
  # if neither file exists, do:
  if(y$index[y$Assembly==x]==0){
    assembly<-y$Assembly[y$Assembly==x] #paste(str_split(x, pattern="_")[[1]][1], str_split(x, pattern="_")[[1]][2], sep="_") # name of assembly
    print(assembly)
    # paste("GCF_", str_split(y$RefSeq.FTP, pattern="_")[[1]][2], sep="")
    Genus<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][1] # Genus name
    Species<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][2] # species name
    Strain<-y[y$Assembly==assembly,3] # strain name
    # test condition
    gene<-rep(NA, length(genes))
    names(gene)<-genes
    
    out<-c(assembly, Genus, Species, Strain, gene)
  }
  
  # if Genbank file exists, do:
  if(y$index[y$Assembly==x]==2){
    #print(y$RefSeq.FTP[y$Assembly==x])
    #print(y$index[y$Assembly==x])
    file1<-as.data.frame(as.matrix(read.delim(paste(directory, "RefSeq/", paste(str_match(y$RefSeq.FTP[y$Assembly==x], pattern = "GCF_.*"), "_feature_table.txt", sep=""), sep=""))))
    assembly<-y$Assembly[y$Assembly==x] #paste(str_split(x, pattern="_")[[1]][1], str_split(x, pattern="_")[[1]][2], sep="_") # name of assembly
    #print(assembly)
    # paste("GCF_", str_split(y$RefSeq.FTP, pattern="_")[[1]][2], sep="")
    Genus<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][1] # Genus name
    Species<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][2] # species name
    Strain<-y[y$Assembly==assembly,3] # strain name
    #print("labels")
    # test condition
    
    gene<-sapply(genes, match.gene, type=1, file1)
    out<-c(assembly, Genus, Species, Strain, gene)
    
  }
  # if only Genbank exists, do:
  if(y$index[y$Assembly==x]==1){
    file1<-as.data.frame(as.matrix(read.delim(paste(directory, "GenBank/", paste(str_match(y$GenBank.FTP[y$Assembly==x], pattern = "GCA_.*"), "_feature_table.txt", sep=""), sep=""))))
    assembly<-y$Assembly[y$Assembly==x] #paste(str_split(x, pattern="_")[[1]][1], str_split(x, pattern="_")[[1]][2], sep="_") # name of assembly
    print(assembly)
    # paste("GCF_", str_split(y$RefSeq.FTP, pattern="_")[[1]][2], sep="")
    Genus<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][1] # Genus name
    Species<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][2] # species name
    Strain<-y[y$Assembly==assembly,3] # strain name
    # test condition
    
    gene<-sapply(genes, match.gene, type=1, file1)
    out<-c(assembly, Genus, Species, Strain, gene)
  }
  
  
  # if both exist:
  if(y$index[y$Assembly==x]==3){
    file1<-as.data.frame(as.matrix(read.delim(paste(directory, "GenBank/", paste(str_match(y$GenBank.FTP[y$Assembly==x], pattern = "GCA_.*"), "_feature_table.txt", sep=""), sep=""))))
    file2<-as.data.frame(as.matrix(read.delim(paste(directory, "RefSeq/", paste(str_match(y$RefSeq.FTP[y$Assembly==x], pattern = "GCF_.*"), "_feature_table.txt", sep=""), sep=""))))
    assembly<-y$Assembly[y$Assembly==x] #paste(str_split(x, pattern="_")[[1]][1], str_split(x, pattern="_")[[1]][2], sep="_") # name of assembly
    #print(assembly)
    # paste("GCF_", str_split(y$RefSeq.FTP, pattern="_")[[1]][2], sep="")
    Genus<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][1] # Genus name
    Species<-str_split(y$X.Organism.Name[y$Assembly==assembly], pattern = " ")[[1]][2] # species name
    Strain<-y[y$Assembly==assembly,3] # strain name
    #print("labels")
    # test condition
    
    gene<-sapply(genes, match.gene, type=2, file1=file1, file2=file2)
    out<-c(assembly, Genus, Species, Strain, gene)
    
    
  }
  out
  
}


