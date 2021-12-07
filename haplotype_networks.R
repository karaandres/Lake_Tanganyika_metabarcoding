### Lake Tanganyika ND2 haplotype networks
### Written by Kara Andres (kja68), last updated 7 Dec, 2021

### Clear the working environment and load required packages
rm(list = ls())

library(pegas)
library(seqinr)
library(RColorBrewer)
library(Biostrings)
library(viridis)
library(maptools)
library(mapplots)
library(sf)
library(rgdal)

### Read in files
species_annotation <- read.csv("Annotations_taxon assignment.csv") # taxonomic assignments from Katie/Pete & BLAST
ASV_ND2_profiling_table <- read.csv("ASV_ND2_profiling_table.csv") # ASV table: # reads per site per ASV
fastaFile <- readDNAStringSet("LT_ND2_ASVs.fa") # ND2 ASV sequences
ND2_seqs <- data.frame(seq_name = names(fastaFile), sequence = paste(fastaFile))
sample_metadata <- read.csv("LT_sample_metadata.csv")
lt_shapefile <- readOGR("map_data/water/") # file path to lake shapefile

### Cut down ASV_ND2_profiling_table to just samples of interest -- 92 sites
to_keep <- sample_metadata[sample_metadata$Vasco_dataset=="X",]
ASV_ND2_profiling_table <- ASV_ND2_profiling_table[,colnames(ASV_ND2_profiling_table) %in% c("OTUID", to_keep$Sample)]
ncol(ASV_ND2_profiling_table)

### Create input file for haplotpye network: fasta file with sequences representing relative ASV abundance at each site
hap_network_fasta <- function(ASV_table, annotation, species){
  # 1) Subset ASV table to ASVs matching species of interest
  ASV_subset <- annotation[grep(species, annotation$NCBI_top_match),] # subset annotation table
  ASV_table_subset <- ASV_table[ASV_table$OTUID %in% ASV_subset$OTUID,] # subset ASV table
  ASV_table_subset_prop <- data.frame(OTUID=ASV_table_subset$OTUID,
                                      lapply(ASV_table_subset[,-1], function(x) x/sum(x)*100)) # scale each sample to 100 reads
  ASV_table_subset_prop[is.na(ASV_table_subset_prop)] <- 0 # replace NAs with 0
  # 2) Combine duplicate samples per site (normalized per sample and scaled to 100 reads)
  ASV_table_subset_prop_combined <- data.frame(OTUID=ASV_table_subset_prop$OTUID)
  haplotype_frequency <- data.frame(OTUID=ASV_table_subset_prop$OTUID)
  k <- 2
  for (i in unique(sample_metadata$SampleInfor)){
    site_data <- sample_metadata[sample_metadata$SampleInfor==i&sample_metadata$Vasco_dataset=="X",]
    if (nrow(site_data)>0){
      site_reads <-  ASV_table_subset_prop[,colnames(ASV_table_subset_prop) %in% site_data$Sample]
      if(class(site_reads)=="data.frame"){ # if more than 2 sites remaining to combine
        site_reads_combined <- rowSums(site_reads, na.rm=TRUE)/sum(rowSums(site_reads, na.rm=TRUE)) # normalize read count per sample
      } else site_reads_combined <- site_reads/sum(site_reads) # normalize read count per sample
      site_reads_combined[is.na(site_reads_combined)] <- 0
      site_reads_combined_round <- round(site_reads_combined*100, digits=0) # combine replicates and scale to 100 reads
      ASV_table_subset_prop_combined <- cbind(ASV_table_subset_prop_combined,site_reads_combined_round)
      colnames(ASV_table_subset_prop_combined)[k] <- site_data$SampleInfor[1]
      haplotype_frequency <- cbind(haplotype_frequency,site_reads_combined)
      colnames(haplotype_frequency)[k] <- site_data$SampleInfor[1]
      k <- k+1
    }
  }
  # 3) Export a haplotype frequency table for species of interest
  write.csv(haplotype_frequency, file=paste("haplotype_frequency/",gsub(" ", "_",species),"_haplotype_frequency.csv", sep=""), row.names=FALSE)
  # 4) Generate a fasta file based on the read frequencies of ASVs per sample
  fasta_dat <- data.frame(seq_name=NULL, sequence=NULL)
  for (i in 2:ncol(ASV_table_subset_prop_combined)){ # for each sample in ASV table
    for (j in 1:nrow(ASV_table_subset_prop_combined)){ # for each ASV
      reads <- ASV_table_subset_prop_combined[j,i] # proportion of reads of ASV in the sample
      if (reads>0){ # if the ASV is present in the sample
        sample_metadata_subset <- sample_metadata[sample_metadata$SampleInfor==colnames(ASV_table_subset_prop_combined[i]),]
        fasta_dat_temp <- data.frame(seq_name = rep(paste(colnames(ASV_table_subset_prop_combined[i]),"_",sample_metadata_subset$Basin[1],"_",sample_metadata_subset$Color[1],"_",ASV_table_subset_prop_combined[j,1], sep = ""), reads), # paste the sample name and ASV name
                                     sequence = rep(ND2_seqs[ND2_seqs$seq_name==ASV_table_subset[j,1],]$sequence, reads)) # paste the ASV sequence
        fasta_dat <- rbind(fasta_dat, fasta_dat_temp)
      }
    }
  }
  return(fasta_dat)
}

### Plotting function: specify species and scale ratio
plot_hap_network <- function(species, scale, size){
  x <- read.dna(paste("fasta_dat/",species,".fa", sep=""), format="fasta")
  h <- haplotype(x) # extract haplotypes from the set of DNA sequences
  R <- haploFreq(x, split="_", what=3, haplo=h) # extract haplotype frequencies by site
  R <- R[, as.character(sort(as.numeric(colnames(R))))] # reorder by latitude
  cols1 <- colorRampPalette(brewer.pal(9, "Purples")[3:9])(12) # north basin colors
  cols2 <- colorRampPalette(brewer.pal(11, "BrBG")[6:11])(34) # south basin colors
  cols3 <- c(cols1, cols2) # all colors
  cols <- cols3[as.numeric(colnames(R))]
  net <- haploNet(h) # compute the haplotype network
  pdf(paste("hap_networks/",species,"_network.pdf", sep=""), height=8, width=5)
  plot(net, size = (attr(net, "freq")/size), pie=R, bg=cols, scale.ratio=scale, 
       show.mutation=2, labels=FALSE, threshold=c(1,2)) # don't display alternate connections
  title(main=bquote(~italic(.(species))))
  dev.off()
  # plot the map
  pdf(paste("hap_networks/",species,"_map.pdf", sep=""), height=8, width=5)
  par(mar=c(1,1,1.2,1))
  plot(lt_shapefile, lwd=0.7, xlim=c(29.1, 30.3), ylim=c(-6.5, -4.8), main=species)
  to_plot <- to_keep[!duplicated(to_keep$Color),]
  to_plot <- to_plot[to_plot$Color %in% colnames(R),]
  plot_cols <- cols3[as.numeric(colnames(R))]
  points(to_plot$Longitude, to_plot$Latitude, pch=21, cex=1.5, bg=plot_cols, col="black", lwd=0.5)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], lwd=2)
  dev.off()
}

### Determine which species match criteria: multiple ASVs, multiple sites, high matching %
sp_list <- NULL
for (i in unique(species_annotation$NCBI_top_match)){
  ASV_subset <- species_annotation[grep(i, species_annotation$NCBI_top_match),] # subset annotation table
  if (nrow(ASV_subset)>2){ # if more than 2 ASVs match to species
    ASV_table_subset <- ASV_ND2_profiling_table[ASV_ND2_profiling_table$OTUID %in% ASV_subset$OTUID,] # subset ASV table
    site_IDs <- sample_metadata[,c("Sample", "SampleInfor")]
    sites_present <- ASV_table_subset[,-1][,colSums(ASV_table_subset[,-1])>0] # if ASVs are present at more than 2 sites
    sites_present_IDs <- site_IDs[site_IDs$Sample %in% colnames(sites_present),]
    if(length(unique(sites_present_IDs$SampleInfor))>2){
      sp_list <- c(sp_list, i)
    }
  }
}

# Generate files for all species matching criteria: 
# fasta files, haplotype frequencies, haplotype networks, maps 
for (i in sp_list[-length(sp_list)]){
  fasta_dat <- hap_network_fasta(ASV_ND2_profiling_table, species_annotation, i)
  i_mod <- gsub(" ", "_", i)
  write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out=paste0("fasta_dat/",i_mod,".fa"), nbchar = 600, as.string = TRUE)
  plot_hap_network(i_mod, 2, 50)
}

# Once we decide on the species we want to focus on in the manuscript,
# use interactive map function to make haplotype networks prettier: 
plot_hap_network_interactive <- function(species, scale, size){
  x <- read.dna(paste("fasta_dat/",species,".fa", sep=""), format="fasta")
  h <- haplotype(x) # extract haplotypes from the set of DNA sequences
  R <- haploFreq(x, split="_", what=3, haplo=h) # extract haplotype frequencies by site
  R <- R[, as.character(sort(as.numeric(colnames(R))))] # reorder by latitude
  cols1 <- colorRampPalette(brewer.pal(9, "Purples")[3:9])(12) # north basin colors
  cols2 <- colorRampPalette(brewer.pal(11, "BrBG")[6:11])(34) # south basin colors
  cols3 <- c(cols1, cols2) # all colors
  cols <- cols3[as.numeric(colnames(R))]
  net <- haploNet(h) # compute the haplotype network
  plot(net, size = (attr(net, "freq")/size), pie=R, bg=cols, scale.ratio=scale, 
       show.mutation=2, labels=FALSE, threshold=c(1,2)) # don't display alternate connections
  o <- replot() # interactive -- drag nodes around
  plot(net, size = (attr(net, "freq")/size), pie=R, bg=cols, scale.ratio=scale, 
      show.mutation=2, labels=FALSE, threshold=c(1,2))
  replot(o) # not interactive
  title(main=bquote(~italic(.(species))))
}

