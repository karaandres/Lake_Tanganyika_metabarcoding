### Lake Tanganyika ND2 haplotype networks
### Written by Kara Andres (kja68), last updated 21 Oct, 2021

# Clear the working environment and load required packages ####
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

# Read in files
species_annnotation <- read.csv("Annotations_taxon assignment.csv") # taxonomic assignments from Katie/Pete & BLAST
ASV_ND2_profiling_table <- read.csv("ASV_ND2_profiling_table.csv") # ASV table: # reads per site per ASV
fastaFile <- readDNAStringSet("LT_ND2_ASVs.fa") # ND2 ASV sequences
ND2_seqs <- data.frame(seq_name = names(fastaFile), sequence = paste(fastaFile))
sample_metadata <- read.csv("LT_sample_metadata.csv")
lt_shapefile <- readOGR("map_data/water/") # file path to lake shapefile

# Cut down ASV_ND2_profiling_table to just samples of interest
to_keep <- sample_metadata[sample_metadata$Vasco_dataset=="X",]
ASV_ND2_profiling_table <- ASV_ND2_profiling_table[,colnames(ASV_ND2_profiling_table) %in% c("OTUID", to_keep$Sample)]

# Combine duplicate samples per site (normalized per sample and scaled to 100 reads)
ASV_ND2_profiling_table_combined <- data.frame(OTUID=ASV_ND2_profiling_table$OTUID)
j <- 2
for (i in unique(sample_metadata$SampleInfor)){
  site_data <- sample_metadata[sample_metadata$SampleInfor==i,]
  if (sum(site_data$Vasco_dataset == "X")>0){
    site_reads <-  ASV_ND2_profiling_table[,colnames(ASV_ND2_profiling_table) %in% site_data$Sample]
    site_reads_norm <- apply(site_reads, 2, function(x) x/sum(x, na.rm = TRUE)) # normalize read count per sample
    site_reads_norm_combined <- round(rowSums(site_reads_norm, na.rm=TRUE)/sum(rowSums(site_reads_norm, na.rm=TRUE))*1000, digits=0) # combine replicates and scale to 100 reads
    ASV_ND2_profiling_table_combined <- cbind(ASV_ND2_profiling_table_combined,site_reads_norm_combined)
    colnames(ASV_ND2_profiling_table_combined)[j] <- site_data$SampleInfor[1]
    j <- j+1
  }
}

# Generate a fasta file based on the read frequencies of ASVs per sample
hap_network <- function(ASV_table, annotation, species){ # start function
  ASV_subset <- annotation[grep(species, annotation$Local_database_top_match),] # subset annotation table to ASVs matching to species of interest
  # ASV_subset <- annotation[annotation$Blast_top_match==species,] # subset annotation table to ASVs matching to species of interest
  ASV_table_subset <- ASV_table[ASV_table$OTUID %in% ASV_subset$OTUID,] # subset ASV table to ASVs of interest
  # ASV_table_subset_prop <- as.data.frame(lapply(ASV_table_subset[,-1], function(x) round(x/sum(x)*100, digits=0))) # scale each sample to 100 reads
  # ASV_table_subset_prop[is.na(ASV_table_subset_prop)] <- 0 # replace NAs with 0
  rownames(ASV_table_subset) <- ASV_table_subset$OTUID # specify rownames as ASV IDs
  fasta_dat <- data.frame(seq_name=NULL, sequence=NULL)
  for (i in 2:ncol(ASV_table_subset)){ # for each sample in ASV table
    for (j in 1:nrow(ASV_table_subset)){ # for each ASV
      reads <- ASV_table_subset[j,i] # proportion of reads of ASV in the sample
      if (reads>0){ # if the ASV is present in the sample
        sample_metadata_subset <- sample_metadata[sample_metadata$SampleInfor==colnames(ASV_table_subset[i]),]
        fasta_dat_temp <- data.frame(seq_name = rep(paste(colnames(ASV_table_subset[i]),"_",sample_metadata_subset$Basin[1],"_",sample_metadata_subset$Color[1],"_",rownames(ASV_table_subset[j,]), sep = ""), reads), # paste the sample name and ASV name
                                     sequence = rep(ND2_seqs[ND2_seqs$seq_name==rownames(ASV_table_subset[j,]),]$sequence, reads)) # paste the ASV sequence
        fasta_dat <- rbind(fasta_dat, fasta_dat_temp)
      }
    }
  }
  return(fasta_dat)
}

# Get all metadata for each sample
hap_network_metadata <- function(ASV_table, annotation, species){ # start function
  ASV_subset <- annotation[grep(species, annotation$Local_database_top_match),] # subset annotation table to ASVs matching to species of interest
  # ASV_subset <- annotation[annotation$Blast_top_match==species,] # subset annotation table to ASVs matching to species of interest
  ASV_table_subset <- ASV_table[ASV_table$OTUID %in% ASV_subset$OTUID,] # subset ASV table to ASVs of interest
  # ASV_table_subset_prop <- as.data.frame(lapply(ASV_table_subset[,-1], function(x) round(x/sum(x)*100, digits=0))) # scale each sample to 100 reads
  # ASV_table_subset_prop[is.na(ASV_table_subset_prop)] <- 0 # replace NAs with 0
  rownames(ASV_table_subset) <- ASV_table_subset$OTUID # specify rownames as ASV IDs
  sample_dat <- data.frame(seq_name=NULL, Latitude=NULL, Longitude=NULL, Basin=NULL)
  for (i in 2:ncol(ASV_table_subset)){ # for each sample in ASV table
    for (j in 1:nrow(ASV_table_subset)){ # for each ASV
      reads <- ASV_table_subset[j,i] # proportion of reads of ASV in the sample
      if (reads>0){ # if the ASV is present in the sample
        sample_metadata_subset <- sample_metadata[sample_metadata$SampleInfor==colnames(ASV_table_subset[i]),]
        sample_dat_temp <- data.frame(seq_name = rep(paste(colnames(ASV_table_subset[i]),"_",sample_metadata_subset$Basin[1],"_",sample_metadata_subset$Color[1],"_",rownames(ASV_table_subset[j,]), sep = ""), reads), # paste the sample name and ASV name
                                      Latitude = rep(sample_metadata_subset$Latitude[1], reads),
                                      Longitude = rep(sample_metadata_subset$Longitude[1], reads),
                                      Basin = rep(sample_metadata_subset$Basin[1], reads))
        sample_dat <- rbind(sample_dat_temp, sample_dat)
      }
    }
  }
  return(sample_dat)
}

# Plotting function: specify species and scale ratio
plot_hap_network <- function(species, scale, size){
  x <- read.dna(paste("/Users/kbja10/Downloads/fasta_dat_",species,".fa", sep=""), format="fasta")
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
  # legend("topleft", c("North", "South"), fill=cols, bty="n")
  title(main=bquote(~italic(.(species))))
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


### Plot haplotype networks for a few different species:
### Petrochromis famula, Oreochromis tanganicae, Bathybates fasciatus,
### Cardiopharynx schoutedeni, Grammatotria lemairii, Xenotilapia flavipinnis,
### Neolamprologus toae, Petrochromis sp. kazumbae, Simochromis diagramma,
### Cyathopharynx furcifer, Neolamprologus cancellatus, Lestradea perspicax
### Limnotilapia dardennii, Lobochilotes_labiatus, Pseudosimochromis babaulti,
### Petrochromis_sp_moshi

# Petrochromis famula: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Petrochromis_famula")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Petrochromis_famula")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Petrochromis_famula.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Petrochromis_famula", 2, 50)

# Oreochromis tanganicae: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Oreochromis_tanganicae")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Oreochromis_tanganicae")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Oreochromis_tanganicae.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Oreochromis_tanganicae", 2, 100)

# Bathybates fasciatus: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Bathybates_fasciatus")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Bathybates_fasciatus")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Bathybates_fasciatus.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Bathybates_fasciatus", 3, 10)

# Cardiopharynx schoutedeni: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Cardiopharynx_schoutedeni")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Cardiopharynx_schoutedeni")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Cardiopharynx_schoutedeni.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Cardiopharynx_schoutedeni", 0.5, 20)

# Grammatotria lemairii: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Grammatotria_lemairii")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Grammatotria_lemairii")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Grammatotria_lemairii.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Grammatotria_lemairii", 2, 100)

# Xenotilapia flavipinnis: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Xenotilapia_flavipinnis")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Xenotilapia_flavipinnis")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Xenotilapia_flavipinnis.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Xenotilapia_flavipinnis", 0.5, 1)

# Neolamprologus toae: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Neolamprologus_toae")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Neolamprologus_toae")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Neolamprologus_toae.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Neolamprologus_toae", 5, 10)

# Petrochromis_sp_kazumbae: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Petrochromis_sp_kazumbae")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Petrochromis_sp_kazumbae")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Petrochromis_sp_kazumbae.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Petrochromis_sp_kazumbae", 1, 200)

# Simochromis_diagramma: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Simochromis_diagramma")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Simochromis_diagramma")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Simochromis_diagramma.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Simochromis_diagramma", 0.5, 50)

# Cyathopharynx furcifer: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Cyathopharynx_furcifer")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Cyathopharynx_furcifer")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Cyathopharynx_furcifer.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Cyathopharynx_furcifer", 0.5, 5)

# Neolamprologus cancellatus: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Neolamprologus_cancellatus")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Neolamprologus_cancellatus")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Neolamprologus_cancellatus.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Neolamprologus_cancellatus", 2, 25)

# Lestradea perspicax: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Lestradea_perspicax")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Lestradea_perspicax")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Lestradea_perspicax.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Lestradea_perspicax", 1, 1)

# Limnotilapia_dardennii: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Limnotilapia_dardennii")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Limnotilapia_dardennii")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Limnotilapia_dardennii.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Limnotilapia_dardennii", 0.5, 150)

# Lobochilotes labiatus: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Lobochilotes_labiatus")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Lobochilotes_labiatus")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Lobochilotes_labiatus.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Lobochilotes_labiatus", 1, 50)

# Pseudosimochromis babaulti: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Pseudosimochromis_babaulti")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Pseudosimochromis_babaulti")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Pseudosimochromis_babaulti.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Pseudosimochromis_babaulti", 1, 10)

# Petrochromis_sp_moshi: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Petrochromis_sp_moshi")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Petrochromis_sp_moshi")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Petrochromis_sp_moshi.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Petrochromis_sp_moshi", 2, 100)

# Eretmodus cyanostictus: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Eretmodus_cyanostictus")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Eretmodus_cyanostictus")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Eretmodus_cyanostictus.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Eretmodus_cyanostictus", 1, 10)

# Ophthalmotilapia ventralis: plot network colored by N-S sites
fasta_dat <- hap_network(ASV_ND2_profiling_table_combined, species_annnotation, "Ophthalmotilapia_ventralis")
sample_dat <- hap_network_metadata(ASV_ND2_profiling_table_combined, species_annnotation, "Ophthalmotilapia_ventralis")
write.fasta(as.list(fasta_dat$sequence), fasta_dat$seq_name, file.out="/Users/kbja10/Downloads/fasta_dat_Ophthalmotilapia_ventralis.fa", nbchar = 600, as.string = TRUE)
plot_hap_network("Ophthalmotilapia_ventralis", 2, 50)

### Full map of Lake Tanganyika and sampling sites
# pdf("LT_sample_map.pdf", height=8, width=5)
par(mar=c(1,1,1.2,1))
plot(lt_shapefile)
to_plot <- to_keep[!duplicated(to_keep$Color),]
points(to_plot$Longitude, to_plot$Latitude, pch=21, cex=1.5, bg=cols3, col="black", lwd=0.5)
# dev.off()
