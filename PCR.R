#packages needed
library(dplyr)
library(phangorn)
library(Biostrings)
library(bio3d)
library(ape)
library(seqinr)
library(DECIPHER)
library(phytools)
library(dendextend)
library(TreeDist)


#make empty dataframes to be appended to for each of the 1000 iterations
df_ent <- data.frame(Entanglement = c(''))  
df_rf <- data.frame(RF = c(''))
rows <- nrow(df)-1
#set working directory to the parent directory, then specify full paths to go into subdirectories
setwd("/MY/DIR")
#open up the pangenome gene list
roary <- read.csv('DIR/gene_presence_absence.csv')   #change to whatever the roary output is 
colnames(roary)[colnames(roary) == 'No..isolates'] = 'isolates'
#subset it to only include genes present in all strains
pan <- subset(roary, isolates == 60)
list <- subset(pan, select=c("Gene")) 
z=1
for (row in 1:1000){
  
  chr_sub <- sample_n(list, 100, replacement = FALSE) #subsample 100 randome genes from core list  ##change to whatever replicon studying
  chr_sub2 <- sample_n(list, 100, replacement = FALSE) #subsample 100 randome genes from core list
  
  # Helper function to read gff file
  read_gff <- function(name) {
    gff_file <- paste0('DIR/',name, ".gbk.gff")
    gff_f <- read.csv(gff_file, sep = "\t", quote = "", comment.char = "#")
    gff <- setNames(gff_f, c('c1','c2','c3','c4','c5','c6','c7','c8','c9'))
  }
                    

  
  # Helper function to process gff data and generate DNA sequence
  process_rep1 <- function(gff, fasta, name) {
    gff2 <- gff %>%
      filter(c3 == "CDS") %>%
      droplevels(.) %>%
      mutate(
        RefSeq_ID = sub(".*RefSeq: *(.*?) *;.*", "\\1", c9),
      )
    chr_fasta <- fasta[1]
    Chromosome <- gff2
    head(Chromosome)
    Chromosome_core_a <- Chromosome %>% #make a df with all core gene annotations
      filter(RefSeq_ID %in% chr_sub$Gene)
    positions <- Chromosome_core_a %>% select(c4, c5)
    dna <- DNAStringSet("")
    for (i in 1:nrow(Chromosome_core_a)) {
      row <- i
      s <- positions[row, "c4"]
      e <- positions[row, "c5"]
      new_gene <- subseq(chr_fasta, start = s, end = e)
      dna <- c(dna, new_gene)
      
    }
    write.fasta(dna, name, paste0('chr2chr/chr1_genes/',name, "Chr1_catgenes.fasta"), open = "w", nbchar = 100000, as.string = FALSE)
    }
    
  #same as function above, but for the second replicon in the comparison
  process_rep2 <- function(gff, fasta, name) {
     gff2 <- gff %>%
       filter(c3 == "CDS") %>%
       droplevels(.) %>%
       mutate(
         RefSeq_ID = sub(".*RefSeq: *(.*?) *;.*", "\\1", c9),
       )
    chr_fasta <- fasta[1]
    Chromosome <- gff2
    head(Chromosome)
    Chromosome_core_b <- Chromosome %>% #make a df with all core gene annotations
      filter(RefSeq_ID %in% chr_sub2$Gene)
    write.csv(Chromosome_core_b,'chr2chr/coreb.fasta')
    positions <- Chromosome_core_b %>% select(c4, c5)
    dna <- DNAStringSet("")
    for (i in 1:nrow(Chromosome_core_b)) {
      row <- i
      s <- positions[row, "c4"]
      e <- positions[row, "c5"]
      new_gene <- subseq(chr_fasta, start = s, end = e)
      dna <- c(dna, new_gene)
      
    }
    
    write.fasta(dna, name, paste0('chr2chr/chr2_genes/',name, "Chr2_catgenes.fasta"), open = "w", nbchar = 100000, as.string = FALSE)
    }
  
    
    otus1 <- list()
    otus2 <- list()
    #creat empty Multifasta string set to be appended to for each of the strains
    chr1_mf <- DNAStringSet()
    chr2_mf <- DNAStringSet()
    
    
    listofstrains <- read.csv('chr2chr/list.csv')
    head(listofstrains)
    
    chr1_mf <- DNAStringSet()
    chr2_mf <- DNAStringSet()
    for (i in 1:nrow(listofstrains)) {
      row <- i
      name <- listofstrains$Strains[row]
      
      # Read the gff file
      gff <- read_gff(name)
      
      # Read the fasta file
      fasta_file <- paste0(name, ".fasta")
      fasta <- readDNAStringSet(paste0('whole_gen/chr/',fasta_file))
      head(fasta)
      # Process the chromosome data and generate the DNA sequence
      process_rep1(gff, fasta, name)
      process_rep2(gff, fasta, name)
      
      otus1 <- append(otus1,paste0(name,'Chr1_catgenes'))
      otus2 <- append(otus2,paste0(name,'Chr2_catgenes'))
      x <- readDNAStringSet(paste0('chr2chr/chr1_genes/',name,'Chr1_catgenes.fasta'))
      chr1_mf <- DNAStringSet(c(chr1_mf,x ))
      
      head(chr1_mf)
      
      y <- readDNAStringSet(paste0('chr2chr/chr2_genes/',name,'Chr2_catgenes.fasta'))
      chr2_mf <- DNAStringSet(c(chr2_mf,y))

      head(chr2_mf)
    }
    #create file containing multifasta information
    writeXStringSet(chr1_mf,'chr2chr/chr1mf.fasta')
    writeXStringSet(chr2_mf,'chr2chr/chr2mf.fasta')
    #align the MF and create a NJ tree from it
    alig1 <- AlignTranslation(chr1_mf)
    distane_matrix1 <- DistanceMatrix(alig1)
    tree1 <- NJ(distane_matrix1)
    #plot(tree1)
    
    alig2 <- AlignTranslation(chr2_mf)
    #head(alig2)
    distane_matrix2 <- DistanceMatrix(alig2)
    tree2 <- NJ(distane_matrix2)
    #plot(tree2)
    
    #convert phylo trees into dendrograms for the creation of tanglegrams
    tree2 <- midpoint(tree2, node.labels = "support")
    tree2_um <- force.ultrametric(tree2, method=c("nnls","extend"))
    tree2_dend <- as.dendrogram(tree2_um)
    tree2_dend <- tree2_dend %>% set('labels_cex',.65)
    tree2_dend <- tree2_dend %>% set('branches_lwd', 2)
    
    
    tree1 <- midpoint(tree1, node.labels = "support")
    tree1_um <- force.ultrametric(tree1, method=c("nnls","extend"))
    tree1_dend <- as.dendrogram(tree1_um)
    tree1_dend <- tree1_dend %>% set('labels_cex',.65)
    tree1_dend <- tree1_dend %>% set('branches_lwd', 2)
    
    
    dends_list <- dendlist(tree1_dend,tree2_dend)
    
    dends_list %>% plot(main = paste("entanglement =", round(entanglement(dends_list), 2)), highlight_distinct_edges = FALSE, common_subtrees_color_branches = FALSE)
    print(row)
    
    #continously update the entanglement and RF excel files
    df_ent[z,] = entanglement(dends_list)
    write.csv(df_ent, 'chr2chr/entanglement.csv')
    df_rf[z,] = RF.dist(tree1,tree2)
    write.csv(df_rf, 'chr2chr/RF.csv')
    z= z+1
  }


