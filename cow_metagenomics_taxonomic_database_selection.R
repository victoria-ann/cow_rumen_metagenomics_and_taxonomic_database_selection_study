library(devtools)
library(dada2); packageVersion('dada2')
library(beepr)
library(ggplot2)

# 1. Specify the folder path to the raw (unfiltered) sequence files by assigning the folder path to a variable.
path_raw_R <- '/Users/toripalecek/Documents/AVS554/AVS554_Project/acidosis_ITS_raw_reverse' 

# 2. Specify what file names (fastq files) to look for in that folder based on the type of file extension. As written, this looks for only zipped (gz) files.
fnsF <- list.files(path_raw_F)
fastqsF <- fnsF[grepl('.gz$', fnsF)] 

fnsR <- list.files(path_raw_R)
fastqsR <- fnsR[grepl('.gz$', fnsR)]


# 3. Fully specify the path, and file names (fns) for the reads
fnsF <- file.path(path_raw_F, fastqsF)
fnsR <- file.path(path_raw_R, fastqsR)


# 4. Pull the sample names by reading the forward read file names, cutting it, and taking the sample name chunk.
sample.names <- sapply(strsplit(basename(fastqsF), '_for'), `[`, 1) 


# 5. Assign the sample names to the list of Read 1 forward fastq files
names(fastqsF) <- sample.names 



## Check sequence quality 

# 1.  Randomly select a few samples to look at quality profile plot. 

plotQualityProfile(fnsF[[9]]); beep() 

plotQualityProfile(fnsR[[10]]); beep() 


# 2. Look at the whole dataset to aggregate all samples onto one graph
plotQualityProfile(fnsF, aggregate=T); beep() 

plotQualityProfile(fnsR, aggregate=T); beep() 

# 3. Get a full quality report for samples

# more info on this package: https://www.rdocumentation.org/packages/Rqc/versions/1.6.2

library(devtools)
library(Rqc)


# run the command to create a quality control report, be sure to update the folder path to your raw files
rqc(path = "/Users/toripalecek/Documents/AVS554/AVS554_Project/acidosis_ITS_raw_forward",
    sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
    outdir = tempdir(), file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())

# when it is done running, a webpage will open up with your sample report

# 4. Search for strings in your data
library(Biostrings)

## specify one file name to look through
DNA_to_read <- readDNAStringSet(fnsF[[9]], format="fastq", with.qualities=TRUE)

# count the number of different nucleotides
alphabetFrequency(DNA_to_read)

# look for ambiguous bases
vmatchPattern("N", DNA_to_read)

# look for a homopolymer. You can change the AAAs to other letters or number of letters.
vmatchPattern("AAAAAAAAAAA", DNA_to_read, fixed=FALSE) 

### Summary of quality ----
# This dataset contains 34 number of samples, and 1,853,959 (forward only) reads which are 300 bases long
# The R1 quality goes below 25 at cycle 250, so I will be trimming at 250. Yes the quality is good enough to use R1.
# The R2 quality goes below 25 at cycle 175. Yes the quality is good enough to use R2. There are 12 less samples for R2, so I will not be using it.


## Filter and trim raw sequences ---------------------------
# Filer and trim sequences based on various parameters

library(dada2); packageVersion('dada2')

# 1. Create folders on your computer to place your filtered data into
dir.create('/Users/toripalecek/Documents/AVS554/AVS554_Project/filtered') 
dir.create('/Users/toripalecek/Documents/AVS554/AVS554_Project/filtered/filtered_forward')
dir.create('/Users/toripalecek/Documents/AVS554/AVS554_Project/filtered/filtered_reverse')

# 2. Specify the folder path to place filtered files for reads 1 and 2.
path_filt_F <- '/Users/toripalecek/Documents/AVS554/AVS554_Project/filtered/filtered_forward'
path_filt_R <- '/Users/toripalecek/Documents/AVS554/AVS554_Project/filtered/filtered_reverse'

# 3. Choose if you want to use only the Forward R1 reads or the Forward R1 + Reverse R2.

### Filter and trim Forward only------
# Process only the Forward R1 files if the quality of your R2 was low.
# https://rdrr.io/bioc/dada2/man/filterAndTrim.html

# run this whole chunk of code from filtoutput through verbose = TRUE)
filtoutput <- filterAndTrim(
  file.path(path_raw_F,fastqsF), # path to and file names of the raw reads
  file.path(path_filt_F, paste0(sample.names, "_F_filt.fastq.gz")), # set filtered file paths and filtered file names to create
  trimLeft=10, # cuts off the first XX bases from the F reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
  trimRight=50, # cuts of last XX bases, or hash this out and use truncLen
  #truncLen=280, # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
  #minLen=10, # optional: remove any sequences longer than max length, for use with 454 if don't want to truncate
  maxEE=2, # the maximum number of expected errors allowed in a read, always >1
  maxN= 0, # max number of ambiguous bases (N) allowed, DADA2 doesn't allow any!!
  rm.phix=TRUE, # remove any PhiX DNA (used as positive control),
  verbose=TRUE, matchIDs = TRUE); beep() # verbose = print to screen

# save the filtered variable for sample totals
saveRDS(filtoutput, "filtered_output_Cow_Rumen_ITS.rds") 

# check the dimensions of the variable created
dim(filtoutput) 

# take a look at the counts
filtoutput

# order the info by the first column (reads.in)
filtoutput[order(filtoutput[,1],decreasing=FALSE),]

# order the info by the second column (reads.out) which are filtered reads
filtoutput[order(filtoutput[,2],decreasing=FALSE),]

# get a sum total for raw reads in and filtered reads out
colSums(filtoutput)

# look at trends
library(ggplot2)
ggplot(as.data.frame(filtoutput)) + geom_point(aes(row.names(filtoutput), reads.in),color = "blue") + geom_point(aes(row.names(filtoutput), reads.out), color = "orange")

### Summary of filtering ----

# What was the largest and small number of raw reads?
##The smallest is from F29_ITS_for.fastq.gz 10545
##The largest is from E53_ITS_for.fastq.gz 99336

# What was the largest and small number of filtered reads?
##The smallest is from F29_ITS_for.fastq.gz 6911
##The largest is from E53_ITS_for.fastq.gz 82582

##There are 1853959 raw reads
##There are 1248310 filtered reads

# DADA2 learn error rates----------------------------
library(dada2); packageVersion('dada2')

# 1. Set the seed for using random number generators
set.seed(1234) 

# 2. Specify the folder path
path_filt_F <- '/Users/toripalecek/Documents/AVS554/AVS554_Project/filtered/filtered_forward'

# 3. Set the names for Forward R1 filtered files
fnsF <- list.files(path_filt_F, full.names = TRUE)
filtsF <- fnsF[grepl('.gz$', fnsF)]
sample.namesF <- sapply(strsplit(basename(filtsF), '_F_filt'), `[`, 1) 
names(filtsF) <- sample.namesF 

# 4. Learn the error rates from your sequencing run.
errF <- learnErrors(filtsF, nbases = 1e6, multithread=TRUE, randomize=TRUE); beep()

# 5. Save the error profiles to output files
saveRDS(errF, 'error_profile_F_Cow_Rumen_ITS.RDS')

## DADA2 pick sequence variants ----------------------

# 1. Decide if you are running the Forward R1 only or the Forward R1 + Reverse R2.

### Sample inference of FORWARD reads only -----------------------  

# 1. Sample inference of FORWARD reads only using previously made error profile.
dadaFs <- dada(filtsF, err=errF, multithread=FALSE); beep()

# 2. Construct a sequence table.
seqtab <- makeSequenceTable(dadaFs); beep()

# 3. Get the dimensions of your table
dim(seqtab) 
# 34 samples and 11750 SVs

# 4. Save the seqtab
saveRDS(seqtab, '/Users/toripalecek/Documents/AVS554/AVS554_Project/seqtab_Cow_Rumen_ITS.rds') 


# DADA2 Remove chimeras from seqtab ----------------------------	

library(dada2); packageVersion('dada2')

# Reload your seqtab
seqtab <- readRDS('/Users/toripalecek/Documents/AVS554/AVS554_Project/seqtab_Cow_Rumen_ITS.rds')

# 1. Remove chimeras	
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) 
#My samples have 259 bimeras out of 11750 total sequences

# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab.nochim) 
#34 samples and 11491 SVs

# Calculate the percentage of chimeras identified out of the total
sum(seqtab.nochim)/sum(seqtab)

# 3. Save your chimera-free sequence table.
saveRDS(seqtab.nochim, '/Users/toripalecek/Documents/AVS554/AVS554_Project/seqtab.nochim_Cow_Rumen_ITS.rds')

#reload as needed
seqtab.nochim <- readRDS('/Users/toripalecek/Documents/AVS554/AVS554_Project/seqtab.nochim_Cow_Rumen_ITS.rds')


## Workflow verification steps--------------------------

require(tidyr)
library(ggplot2)

# 1. Track reads through the analysis here to see how many were lost at each QC step
getN <- function(x) sum(getUniques(x))

# load the filtered output file if it isn't already in your environment
filtoutput <- readRDS("filtered_output_Cow_Rumen_ITS.rds") #CHANGE file name

# 2. Load the metadata file
meta <- read.csv('/Users/toripalecek/Documents/AVS554/AVS554_Project/Cow_ITS_metadata.csv',header = TRUE, row.names = 1)


# 2.5 check the dimensions of the three data files
dim(filtoutput) # 34 rows
dim(seqtab.nochim) # 34 rows
dim(meta) # 34 rows

meta_cleaned <- meta[row.names(meta) %in% row.names(seqtab.nochim), ]

# 3. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
track <- cbind(filtoutput, rowSums(seqtab.nochim), meta$Treatment) 
# 4. Assign column names to that new variable
colnames(track) <- c("reads.in","filtered", "nonchimeras", "Treatment") 

# 5. Assign rownames for the samples in your new variable
rownames(track) <- rownames(meta)

# 6. Look at the header of that variable to check that it looks the way you want it.
head(track)

# 7. Save the tracking variable and data as an R file
saveRDS(track, 'tracked_seqs_Cow_Lumen_ITS.RDS')

# 8. Plot all reads along the QC workflow
plotData <- as.data.frame(track) %>% gather(type,totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))

# plot with Sample_type along the X axis
ggplot(plotData,aes(x=Treatment,y=as.numeric(totals))) + geom_jitter(aes(color=type)) + 
  ylab("Sequences") + xlab("Treatment") 

# or, plot with QA stage along the X axis
ggplot(plotData,aes(x=type,y=as.numeric(totals))) + geom_jitter(aes(color=Treatment)) + 
  ylab("Sequences") + xlab("QA stage") 

# 9. Randomly select a negative control sample to see what's in it
unqs.NC <- seqtab.nochim["Mock_S13_PCRneg",] 

# Sort by # seqs and drop SVs absent in the negative control sample
unqs.NC <- sort(unqs.NC[unqs.NC>0], decreasing=TRUE) 

# Print out how many reads are inferred in the negative control
cat("DADA2 inferred", length(unqs.NC), "sequence variants present in the selected sample.\n")
# There are 81 SVs in this PCR negative control

# Plot the number of sequences in the SVs found in the negative control sample.
plot(unqs.NC, ylab="Number of seqs/SV, Neg Control", xlab="SVs", main="Number of sequences/SV in Negative Control")


# 10. Randomly select a positive control to see what's in it
unqs.PC <- seqtab.nochim["E19_ITS",]

# Sort by # seqs and drop SVs absent in the positive control sample
unqs.PC <- sort(unqs.PC[unqs.PC>0], decreasing=TRUE) 

# Print out how many reads are inferred in the positive control
cat("DADA2 inferred", length(unqs.PC), "sequence variants present in the selected sample.\n")
# There are 1496 SVs in this positive control

# Plot the number of sequences in the SVs found in the positive control
plot(unqs.PC, ylab="Number of seqs/SV, Neg Control", xlab="SVs", main="Number of sequences/SV in Positive Control")


# Phyloseq First look ------------------------------
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")


library(phyloseq)

# load seqtab as needed
seqtab.nochim <- readRDS('/Users/toripalecek/Documents/AVS554/AVS554_Project/seqtab.nochim_Cow_Rumen_ITS.rds') 


# load metadata table as needed
meta <- read.csv('/Users/toripalecek/Documents/AVS554/AVS554_Project/Cow_ITS_metadata.csv',header = TRUE, row.names = 1)


# Check the sample sames and see if they match
row.names(meta)
row.names(seqtab.nochim)

meta <- meta[order(row.names(meta)),] # this orders them alphabetically
seqtab.nochim <- seqtab.nochim[order(row.names(seqtab.nochim)),] # this orders them alphabetically
row.names(seqtab.nochim) <- row.names(meta) # this copies and pastes the names
saveRDS(seqtab.nochim, 'seqtab-nochim-Cow_Rumen_ITS.RDS') # save so you don't have to do it again

## create a phyloseq object with all samples (subset out later)
EX_ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  sample_data(meta))

EX_ps
# 34 samples, 4 variables, 11491 taxa


library(ggplot2)

# Alpha diversity peek 
plot_richness(EX_ps, x="Treatment", # Factor of your choice
              measures=c("Observed","Chao1", "Shannon"), # alpha diversity measures
              color="Treatment") + # Choose factor
  theme_bw() # a ggplot theme to make the graph look nice


# Plot the taxa sums to see how populated each taxa is
plot(sort(taxa_sums(EX_ps), TRUE), 
     type="h", 
     ylim=c(0, 20)) #limit the y-axis to better see the long tail of rare taxa


# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
EX.ord <- ordinate(EX_ps, #calculate similarities
                   method ="PCoA", #ordination type
                   "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)


plot_ordination(EX_ps, EX.ord, type="samples", color="Treatment")
# Initial clustering by extraction/sequencing batch or confounding variables implies contamination issues 
# Horse-shoe patterns indicate underlying patterns to the data

# Removing contamination using negative controls-----------------------------------
### Only if you have Negative Controls

# 1, Subset the controls that would pertain to all of the samples, like ones from reagents, PCR, or other whole-batch effects. prune to be only those taxa	
EX_ps_controls = subset_samples(EX_ps, Sample_type == "Mock" | Sample_type == "NegCon_PCR") # Column name that holds the variables associated with being a negative control

EX_ps_controls <- prune_taxa(taxa_sums(EX_ps_controls) > 0, EX_ps_controls)
#6368 taxa in the 4 mock samples

# 2, Make the taxa names into a vector so you can remove them. 
control_vec <- as.vector(taxa_names(EX_ps_controls))
vec <- as.vector(taxa_names(EX_ps))
keep <- setdiff(vec, control_vec)

# 3. Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
EX_ps_NCclean <- EX_ps

# 4. Then remove the samples which are now empty, namely the NegCon_PCR and NegCon_swab.
EX_ps_NCclean <- prune_samples(sample_sums(EX_ps_NCclean) > 0, EX_ps_NCclean)

#####Because I didn't have positive controls and couldn't do a decontamination step

## 5. subset out each dna_extraction_batch 
batch1 = subset_samples(EX_ps_NCclean, DNA_extraction_batch == "1") #16 samples
batch2 = subset_samples(EX_ps_NCclean, DNA_extraction_batch == "2") #16 samples


## 6. subset the controls and prune to be only those taxa
batch1_kit = subset_samples(batch1, Sample_type == "NegCon_kit")
batch2_kit = subset_samples(batch2, Sample_type == "NegCon_kit")

batch1_kit <- prune_taxa(taxa_sums(batch1_kit) > 0, batch1_kit) # 1 sample
batch2_kit <- prune_taxa(taxa_sums(batch2_kit) > 0, batch2_kit) # 1 sample


# 7. make the taxa names into a vector so you can remove them 
b1_control_vec <- as.vector(taxa_names(batch1_kit))
b1_vec <- as.vector(taxa_names(batch1))
b1_keep <- setdiff(b1_vec, b1_control_vec)
## then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching)
b1_clean <- prune_taxa(b1_keep, batch1)

# make the taxa names into a vector so you can remove them 
b2_control_vec <- as.vector(taxa_names(batch2_kit))
b2_vec <- as.vector(taxa_names(batch2))
b2_keep <- setdiff(b2_vec, b2_control_vec)
## then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching)
b2_clean <- prune_taxa(b2_keep, batch2)

# 8. Merge the phyloseq objects back together, then remove any blank taxa or samples
EX_ps_NC_batch_clean <- merge_phyloseq(b1_clean, b2_clean)

## 9. clean out taxa/SV columns that are no longer present
EX_ps_NC_batch_clean <- prune_taxa(taxa_sums(EX_ps_NC_batch_clean) > 0, EX_ps_NC_batch_clean)
EX_ps_NC_batch_clean <- prune_samples(sample_sums(EX_ps_NC_batch_clean) > 0, EX_ps_NC_batch_clean)
# 30 number of samples and 11793 SVs remaining 

# 10. Check your ordination again
EX_cleaner.ord <- ordinate(EX_ps_NC_batch_clean, #calculate similarities
                           method ="PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(EX_ps_NC_batch_clean, EX_cleaner.ord, type="samples", color="DNA_extraction_batch", title="After cleaning out negative controls")

## ------ Decontam method -------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(decontam) # Identifies contaminants by frequency of SVs.

# 1. Look for contaminants using DNA quantification data from your metadata file
contamdf.freq <- isContaminant(EX_ps, method="frequency", conc="quant_reading") #conc= "quant_reading" to the column heading that holds the concentration information in metadata

# 2. visualize
head(contamdf.freq) 

# 3. Make it into a table
table(contamdf.freq$contaminant)

# 4. See which SVs are categorized as contaminants based on frequency
head(which(contamdf.freq$contaminant)) 

# 5. Get rid of the contaminants 
EX_ps.noncontamfreq <- prune_taxa(!contamdf.freq$contaminant, EX_ps)

# 6. How much is left
EX_ps.noncontamfreg

# 7. Identify negative controls by indicating which column/factor in metadata and which variable indicate a negative control
sample_data(EX_ps)$is.neg <- sample_data(EX_ps)$Sample_or_Control == "Control Sample"

# 8. Calculate prevalence of SVs in samples versus controls
contamdf.prev05 <- isContaminant(EX_ps, method="prevalence", neg="is.neg", threshold=0.5)

# 9. Make a table
table(contamdf.prev05$contaminant)

# 10. Look at it
head(which(contamdf.prev05$contaminant))

# 11. get rid of the contaminants 
EX_ps_NC_decontam_clean <- prune_taxa(!contamdf.prev05$contaminant, EX_ps.noncontamfreq) # remove them from the original phyloseq object, or the one with cleaned out SVs by frequency

# 12. Check how many are left
EX_ps_NC_decontam_clean

# 13. Check your ordination again
EX_cleaner.ord <- ordinate(EX_ps_NC_decontam_clean, #calculate similarities
                           method ="PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(EX_ps_NC_decontam_clean, EX_cleaner.ord, type="samples", color="dna_extraction_batch", title="After cleaning out negative controls")

# DADA2 assigning taxonomy from a reference database file---------------------------------
library(dada2); packageVersion('dada2')
library(phyloseq)
# Have a taxonomy database file downloaded from here: https://benjjneb.github.io/dada2/training.html. 

# 1. reload as needed
meta <- read.csv('/Users/toripalecek/Documents/AVS554/AVS554_Project/Cow_ITS_metadata.csv',header = TRUE, row.names = 1)

seqtab.nochim <- readRDS('/Users/toripalecek/Documents/AVS554/AVS554_Project/seqtab-nochim-Cow_Rumen_ITS.RDS') 

all.taxa <- assignTaxonomy(seqtab.nochim, 
                           '/Users/toripalecek/Documents/AVS554/AVS554_Project/sh_general_release_dynamic_04.04.2024.fasta', 
                           tryRC = TRUE, minBoot = 75, verbose = TRUE, multithread = TRUE)
all.taxa2 <- assignTaxonomy(seqtab.nochim, 
                            '/Users/toripalecek/Documents/AVS554/AVS554_Project/RDP_LSU_fixed_train_set_v2.fa', 
                            minBoot = 75, verbose = TRUE, multithread = TRUE)

# 3. 
saveRDS(all.taxa, '/Users/toripalecek/Documents/AVS554/AVS554_Project/taxa_dynamic_04.04.2024_Cow_Rumen_ITS.rds')
write.csv(all.taxa, '/Users/toripalecek/Documents/AVS554/AVS554_Project/taxa_dynamic_04.04.2024_Cow_Rumen_ITS.csv')
saveRDS(all.taxa, '/Users/toripalecek/Documents/AVS554/AVS554_Project/RDP_LSU_fixed_train_set_v2_Cow_Rumen_ITS.rds') 
write.csv(all.taxa, '/Users/toripalecek/Documents/AVS554/AVS554_Project/RDP_LSU_fixed_train_set_v2_Cow_Rumen_ITS.csv')


# Add species designation to the table.

all.taxa.sp <- addSpecies(all.taxa, '/Users/toripalecek/Documents/AVS554/AVS554_Project/pr2_version_5.0.0_SSU_dada2.fasta', allowMultiple = FALSE, verbose = FALSE); beep()

saveRDS(all.taxa.sp, '/Users/toripalecek/Documents/AVS554/AVS554_Project/taxa_silva_Cow_Rumen_ITS_sp.rds') 
write.csv(all.taxa.sp, '/Users/toripalecek/Documents/AVS554/AVS554_Project/taxa_silva_Cow_Rumen_ITS_sp.csv') 


# 4. Remake the phyloseq object with the new taxonomy file 

# reload metadata table as needed
meta <- read.csv('/Users/toripalecek/Documents/AVS554/AVS554_Project/Cow_ITS_metadata.csv',header = TRUE, row.names = 1)

#reload taxa table as needed
all.taxa <- readRDS('/Users/toripalecek/Documents/AVS554/AVS554_Project/taxa_dynamic_04.04.2024_Cow_Rumen_ITS.rds')
all.taxa2 <- readRDS('/Users/toripalecek/Documents/AVS554/AVS554_Project/RDP_LSU_fixed_train_set_v2_Cow_Rumen_ITS.rds')
# otu_t <- otu_table(EX_ps)

## create a phyloseq object with all samples
EX_ps <- phyloseq(otu_table(otu_t, taxa_are_rows=FALSE), 
                  sample_data(meta),
                  tax_table(all.taxa))

EX_ps2 <- phyloseq(otu_table(otu_t, taxa_are_rows=FALSE), 
                   sample_data(meta),
                   tax_table(all.taxa2))

saveRDS(EX_ps, "Cow_Rumen_ITS_phyloseq_object_cleaned.rds")
saveRDS(EX_ps2, "Cow_Rumen_ITS_phyloseq_object_cleaned2.rds")


# Clean out unwanted taxa------------------------------------
require(dplyr)

# Explore your taxonomy before filtering
df <- as.data.frame(EX_ps)
table(df$Kingdom)
table(df$Phylum)
table(df$Class)
table(df$Order)

#Solution for above not working ^^^
# Extract taxonomy table from phyloseq object
tax_df <- as.data.frame(tax_table(EX_ps))
tax_df2 <- as.data.frame(tax_table(EX_ps2))

# View taxonomic counts for Kingdom, Phylum, etc.
table(tax_df$Kingdom)
table(tax_df$Phylum)
table(tax_df$Class)
table(tax_df$Order)
table(tax_df$Genus)

table(tax_df2$Kingdom)
table(tax_df2$Phylum)
table(tax_df2$Class)
table(tax_df2$Order)

###### At this point I see no difference in assigning taxonomy in one database vs the other


EX_ps_clean <- EX_ps %>% 
  EX_ps_clean2 <- EX_ps2

## clean out taxa/SV columns that are no longer present
EX_ps_clean <- prune_taxa(taxa_sums(EX_ps_clean) > 0, EX_ps_clean)
EX_ps_clean <- prune_samples(sample_sums(EX_ps_clean) > 0, EX_ps_clean)

EX_ps_clean
# 30 samples and 8091 SVs left 


# Save clean phyloseq object
saveRDS(EX_ps_clean, 'EX_ps_clean_phyloseq_object.RDS') 
saveRDS(EX_ps_clean2, 'EX_ps_clean2_phyloseq_object.RDS')
# Reload as needed
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS') 
EX_ps_clean2 <-readRDS('EX_ps_clean2_phyloseq_object.RDS') 

# Rarefaction--------------------------------------------------
# Reload as needed
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS') 
EX_ps_clean2 <-readRDS('EX_ps_clean2_phyloseq_object.RDS')

# make a rarefaction curve
library(vegan)
EX-rarec <- rarecurve(otu_table(EX_ps_clean), step = 10, cex=0.5, label = FALSE) # step is major tick marks on x axis (x 1000), cex is text size, label is sample name
EX-rarec2 <- rarecurve(otu_table(EX_ps_clean2), step = 10, cex=0.5, label = FALSE)

# If code above doesnt work use this one
tab <- otu_table(Ex_ps_clean)
class(tab) <- "matrix"
tab <- t(tab) # transpose observations to rows
rarecurve(tab, step = 10, cex = 0.5, label = FALSE)

tab2 <- otu_table(Ex_ps_clean2)
class(tab2) <- "matrix"
tab2 <- t(tab2) # transpose observations to rows
rarecurve(tab2, step = 10, cex = 0.5, label = FALSE)


# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(EX_ps_clean)))
sort(rowSums(otu_table(EX_ps_clean2)))
# smallest reads (sequences) in a sample F29 (5207)
# largest in a sample E53 (80240)
# number of samples with <5000 0

# smallest reads (sequences) in a sample F29 (5207)
# largest in a sample E53 (80240)
# number of samples with <5000 0

EX_ps_clean.rar <- rarefy_even_depth(EX_ps_clean, 
                                     sample.size=18000, #5-10k is a good amount, more is better
                                     replace=FALSE, #sampling with or without replacement
                                     trimOTUs=TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                     rngseed=711, 
                                     verbose=TRUE)
EX_ps_clean.rar2 <- rarefy_even_depth(EX_ps_clean2, 
                                      sample.size=18000, # 5-10k is a good amount, more is better
                                      replace=FALSE, #sampling with or without replacement
                                      trimOTUs=TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                      rngseed=711, 
                                      verbose=TRUE)

# 4 samples and 7430 SVs were removed. Samples E21, E31, F29, F58 in both EX_ps_clean.rar and  EX_ps_clean.rar2

# Save your cleaned, rarefied phyloseq object  
saveRDS(EX_ps_clean.rar, 'EX_ps_clean_rarefied_phyloseq_object.RDS')
saveRDS(EX_ps_clean.rar2, 'EX_ps_clean_rarefied_phyloseq_object2.RDS')

# Helpful to have an SV table from the clean, rarefied phyloseq
write.csv(otu_table(EX_ps_clean.rar), 'EX_ps_clean_rarefied_phyloseq_object.csv')
write.csv(otu_table(EX_ps_clean.rar2), 'EX_ps_clean_rarefied_phyloseq_object2.csv')


# Statistical analysis here -------------------------------------
library(phyloseq)

# Non-rarefied version
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS')
EX_ps_clean2 <-readRDS('EX_ps_clean2_phyloseq_object.RDS')

# Rarefied version
EX_ps_clean.rar <-readRDS('EX_ps_clean_rarefied_phyloseq_object.RDS')
EX_ps_clean.rar2 <-readRDS('EX_ps_clean_rarefied_phyloseq_object2.RDS')

# if you have calendar dates that you want to be treated as time:
require(lubridate)

# use as.Date to import date into a usable format for continuous variable. #"%m/%d/%y" is the format of how your date is written in the metadata, and $Sample_date is the column of your time data. Repeat as needed for any dates you have.
sample_data(EX_ps_clean.rar)$Sample_date <- as.Date(sample_data(EX_ps_clean.rar)$Sample_date, "%m.%d.%Y") 


# Alpha diversity  ---------------------------------
## Alpha diversity graphics----------------------

# load ggplot and other packages to make pretty graphs
library(ggplot2)
library(RColorBrewer)

# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

plot_richness(EX_ps_clean.rar, 
              x="Treatment", # Factor on x-axis
              measures="Observed", # choose alpha diversity measure.To have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x=element_blank()) + #example, get rid of x axis labels
  geom_violin(trim=TRUE, aes(fill=Diet)) + # Factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + # Factor to group boxplots
  facet_grid(.~Week, space="free") +
  # theme(legend.position = "none") + # use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  xlab("Week")


# Use geoms to add different shapes/styles of plots: https://ggplot2-book.org/individual-geoms.html#:~:text=3.1%20Basic%20plot%20types,plot%20has%20a%20special%20name.

# Use multiple geoms to create gaph overlays in a single plot, and rearrange their order to move to the foreground or background.

# want to map factors to geoms to get it to set colors by factor? use geom_point(aes(color=Factor)). otherwise, use geom_point(color=Black) to make them all the same.

# two different richness plots grouped
library(ggplot2)
library(ggpubr)

dev.off()

plot_richness(EX_ps_clean.rar, 
              x="Treatment", # factor you want on x-axis
              measures= c("Observed","Shannon"), # richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=All), drop = FALSE) + #factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + #factor to group boxplots
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Fungal Richness (SVs)") + ggtitle("Observed and Shannon Diversity in Cow Rumen ITS Data with Unite Taxonimic Assignment")

plot_richness(EX_ps_clean.rar2, 
              x="Treatment", 
              measures=c("Observed","Shannon"), 
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=All), drop = FALSE) + 
  geom_boxplot(width = 0.1, aes(group=Treatment)) + 
  # theme(legend.position = "none") + 
  ylab("Fungal Richness (SVs)") + ggtitle("Observed and Shannon Diversity in Cow Rumen ITS Data with RDP_LSU Taxonimic Assignment")

plot2 <- plot_richness(EX_ps_clean.rar, 
                       x="Treatment", 
                       measures="Shannon",
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=All), drop = FALSE) + 
  geom_boxplot(width = 0.1, aes(group=Treatment)) +
  # theme(legend.position = "none") +
  ylab("Fungal Richness (SVs)")

plot2_2 <- plot_richness(EX_ps_clean.rar2, 
                         x="Treatment", 
                         measures="Shannon", 
                         title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=All), drop = FALSE) + 
  geom_boxplot(width = 0.1, aes(group=Treatment)) + 
  # theme(legend.position = "none") + 
  ylab("Fungal Richness (SVs)")


# plot together
ggarrange(plot1,
          ggarrange(plot2, nrow=1, labels = c("B")),
          ncol = 2, labels ="A", common.legend = TRUE)




# richness plot with significance added
library(ggplot2)
library(ggsignif)


# richness plot observed SVs with lines to fit the view screen
plot_richness(EX_ps_clean.rar, 
              x="Treatment", 
              measures="Observed", 
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=All), drop=FALSE) + 
  geom_boxplot(width = 0.1, aes(group=Treatment)) + 
  # theme(legend.position = "none") +
  ylab("Observed Bacterial Richness (SVs)") + 
  ylim(0,1500) + #define the y axis min/max
  geom_segment(aes(x = 1, y = 1200, xend = 2, yend = 1200)) +  geom_text(x = 1.5, y = 1250, label = "***") # add a drawn in line and significance tag, adjusting the x and y coordinates till it fits where you want in the window.  Add another for each line to add. As written, this will fit the view window you have, if you adjust that your segments will not adjust with it.




## Alpha diversity plotted against other metadata -----

# use phyloseq to measure alpha diversity
EX_ps_clean.rar.rich <- estimate_richness(EX_ps_clean.rar, measure=c("Observed", "Shannon")) #change to whatever measures you want
EX_ps_clean.rar.rich2 <- estimate_richness(EX_ps_clean.rar2, measure=c("Observed", "Shannon"))

# Use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
install.packages("remotes")
remotes::install_github("twbattaglia/btools")

EX_faith <- estimate_pd(EX_ps_clean.rar)

# measure evenness for each sample
EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)

EX_ps_clean.rar.even2 <- EX_ps_clean.rar.rich2$Shannon/log(EX_ps_clean.rar.rich2$Observed)

# Coerce to data.frame and add the metadata for these samples
EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

dim(EX_ps_clean.rar.rich.df)

EX_ps_clean.rar.sd2 = as(sample_data(EX_ps_clean.rar2), "matrix")
EX_ps_clean.rar.sd2 = as.data.frame(EX_ps_clean.rar.sd2)
EX_ps_clean.rar.rich.df2 <- cbind(EX_ps_clean.rar.rich2, EX_ps_clean.rar.even2, EX_ps_clean.rar.sd2)

dim(EX_ps_clean.rar.rich.df2)
# make a graph using that dataframe.
ggplot(data=EX_ps_clean.rar.rich.df, aes(x=All, y=Observed)) + 
  theme_minimal() + 
  geom_point(aes(color=Treatment), size = 3) + # sets colors by group and a set point size 
  xlab("Sample Group") + 
  ylab("Fungal Richness (SVs)") + 
  theme(text = element_text(size = 7)) + ggtitle("Observed Diversity in Cow Rumen ITS Data with Unite Taxonimic Assignment") # increases font size

ggplot(data=EX_ps_clean.rar.rich.df2, aes(x=All, y=Observed)) + 
  theme_minimal() + 
  geom_point(aes(color=Treatment), size = 3) +
  xlab("Sample Group") + 
  ylab("Fungal Richness (SVs)") + 
  theme(text = element_text(size = 7)) + ggtitle("Observed Diversity in Cow Rumen ITS Data with RDP_LSU Taxonimic Assignment")


# make that same graph but drop any samples that lack data for that FactorA
ggplot(data=subset(EX_ps_clean.rar.rich.df, !is.na(FactorA)), aes(x=Temperature, y=Observed)) + 
  theme_minimal() + 
  geom_point(aes(color=Group), size = 3) + # sets colors by group and a set point size 
  xlab("Temperature of Ocean Water (C)") + 
  ylab("Bacterial Richness (SVs)") + 
  theme(text = element_text(size = 20)) # increases font size



# Make an alpha diversity table: https://www.geeksforgeeks.org/create-table-from-dataframe-in-r/
EX_table = as.table(table(EX_ps_clean.rar.rich.df$All, EX_ps_clean.rar.rich.df$Treatment)) 
EX_table

EX_table2 = as.table(table(EX_ps_clean.rar.rich.df2$All, EX_ps_clean.rar.rich.df2$Treatment)) 
EX_table2


# Alpha diversity metrics statistics--------------
# Phyloseq can measure and visualize alpha diversity: https://joey711.github.io/phyloseq/plot_richness-examples.html

library(phyloseq)

# use phyloseq to measure alpha diversity
EX_ps_clean.rar.rich <- estimate_richness(EX_ps_clean.rar, measure=c("Observed", "Shannon"))

# Measure evenness for each SV individually
library(asbio)
library(microbiome)
EX_ps_clean.rar.even_SV <- evenness(otu_table(EX_ps_clean.rar))

# measure evenness for each sample
EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)


# Coerce to data.frame and add the metadata for these samples
EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

# make a histogram to look at the shape of the data 
hist(EX_ps_clean.rar.rich.df$Observed)

hist(EX_ps_clean.rar.rich.df2$Observed)

# Measure Kurtosis
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
kurtosis(EX_ps_clean.rar.rich.df$Observed)
kurtosis(EX_ps_clean.rar.rich.df2$Observed)
# My output is -0.4498885 which is negative kurtosis
# +/-2 would indicate a problem and is a concerning level of kurtosis

head(EX_ps_clean.rar.rich.df)

head(EX_ps_clean.rar.rich.df2)

#check the distribution of your data, which will change your approach
shapiro.test(EX_ps_clean.rar.rich.df$Shannon)
# W = 0.89889, p-value = 0.007891. My shannon diversity metric is not normally distributed

shapiro.test(EX_ps_clean.rar.rich.df$Observed)
# W = 0.98316, p-value = 0.9018, my observed data is normal 

shapiro.test(EX_ps_clean.rar.rich.df$EX_ps_clean.rar.even)
# W = 0.79864, p-value = 6.19e-05, not normally distributed


## If your alpha diversity is normally distributed (parametric)-----------------
library(vegan)

### Linear model for numeric factors -----
# Interpret linear model in R https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
# https://boostedml.com/2019/06/linear-regression-in-r-interpreting-summarylm.html

library(lme4)
library(lmerTest)
library(emmeans)

# pairwise comparisons
emmeans(lm,pairwise ~ All) 

### Wilcoxon (a.k.a. Mann-Whitney) for pairwise numeric or categorical factors -------
# Wilcoxon in R (a.k.a. Mann-Whitney): https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/


wilcox.test(y ~ A, data=EX_ps_clean.rar.rich.df)


### ANOVA for multiple categorical factors -------

# One Way Anova (Completely Randomized Design)
fit <- aov(y ~ A, data=EX_ps_clean.rar.rich.df)

# Randomized Block Design (B is the blocking factor)
fit <- aov(y ~ A + B, data=EX_ps_clean.rar.rich.df)

# Two Way Factorial Design
fit <- aov(y ~ A*B, data=EX_ps_clean.rar.rich.df) 

# Analysis of Covariance
fit <- aov(y ~ A + x, data=EX_ps_clean.rar.rich.df) # x is a series of numerical values, like pH)


# One way anova
example1 <- lm(Observed ~ Diet, data=EX_ps_clean.rar.rich.df)
example2 <- lm(Observed ~ Week, data=EX_ps_clean.rar.rich.df)

# Find the AIC value
anova(example1,example2)




#example of two way 
summary(aov(Observed ~ Diet * Week, data=EX_ps_clean.rar.rich.df))

# One Within Factor
fit <- aov(y~A+Error(Subject/A),data=EX_ps_clean.rar.rich.df) #Subject would be a Factor in your metadata, y is the value of interest, A is the Factor of interest

# Two Within Factors W1 W2, Two Between Factors B1 B2
fit <- aov(y~(W1*W2*B1*B2)+Error(Subject/(W1*W2))+(B1*B2),
           data=EX_ps_clean.rar.rich.df)

# Tukey Honestly Significant Differences
TukeyHSD(fit)

## If alpha diversity is non-parametric -----------------
### Kruskal-Wallis Test -----
# K-W is the non-parametric version of ANOVA: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
kruskal.test(Observed ~ All, data = EX_ps_clean.rar.rich.df)
# data:  Observed by All
# Kruskal-Wallis chi-squared = 12.669, df = 7, p-value = 0.806

# Conover Test for multiple comparisons
# https://rdrr.io/cran/DescTools/man/ConoverTest.html 
install.packages("conover.test")
library(conover.test)

conover.test(EX_ps_clean.rar.rich.df$Observed, EX_ps_clean.rar.rich.df$Week,
             method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) # method of correction

# Kruskal-Wallis rank sum test
# data: x and group
# Kruskal-Wallis chi-squared = 5.2867, df = 2, p-value = 0.07
# 
# Comparison of x by group                            
# (No adjustment)                                
# Col Mean-|
# Row Mean |          0          1
---------+----------------------
  #       1 |  -1.818304
  #         |     0.0417
  #         |
  #       2 |  -2.467699  -0.726045
  #         |    0.0111*     0.2379
  # 
  # alpha = 0.05
  # Reject Ho if p <= alpha/2
  
  ### Linear model for numeric factors -----
# Interpret linear model in R https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
# https://boostedml.com/2019/06/linear-regression-in-r-interpreting-summarylm.html 



library(lme4)
library(lmerTest)
library(emmeans)

# run each factor separately to find out if they are significant on their own. 
# use as.numeric() for a numeric factor and use factor() for a categorical factor
summary(lm(Observed ~ factor(All), data= EX_ps_clean.rar.rich.df)) 


# run each factor separately on just a subset of the data.
summary(lm(Observed ~ as.numeric(PHOS), data=subset(EX_ps_clean.rar.rich.df, Week=="2"))) 

### Linear mixed effects models for complicated experimental designs -----
lm <- (glmer(Observed ~ A + C + (1|Year_collected), family=poisson, data=EX_ps_clean.rar.rich.df))

summary(lm)

# if you have multiple groups, this will give you the pairwise comparisons
emmeans(lm,pairwise ~ A) 

emmeans(lm,pairwise ~ B) 

# More ways to graph alpha diversity-------
# https://www.data-to-viz.com/
# https://r-graph-gallery.com/index.html

## Heatmaps --------------
library(phyloseq)
# base plot
plot_heatmap(EX_ps_clean.rar)


EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Phylum")

plot_heatmap(EX_ps_clean.rar, fill="Phylum") +   facet_grid(~Week, space="free", scales="free") + theme(legend.position = "bottom") #, axis.text.x = element_blank()) 


plot_heatmap(EX_ps_clean.rar.glom, fill="Phylum", taxa.label = "Phylum") +   facet_grid("All", space="free", scales="free") + theme(legend.position = "bottom") #, axis.text.x = element_blank())
# Load necessary packages
library(phyloseq)
library(ggplot2)
library(reshape2)


phyloseq_phylum <- tax_glom(EX_ps_clean.rar.glom, "Phylum")

# relative abundance so heatmap is readable
phyloseq_phylum_rel <- transform_sample_counts(phyloseq_phylum, function(x) x / sum(x))

# create data frame
phylo_df <- psmelt(phyloseq_phylum_rel)

# "All" transformation to vector 
phylo_df$All <- sample_data(EX_ps_clean.rar.glom)[[ "All"]]

ggplot(phylo_df, aes(x = All, y = Phylum, fill = Abundance)) +
  geom_tile() +
  facet_grid(. ~ All, space = "free", scales = "free") + 
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Relative Abundance by Sample Group", 
       x = "Diet,Yeast Ingestion and Sample Location Groupings", 
       y = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() + 
  theme(legend.position = "bottom") 


## Correlogram --------------
# This makes a correlation matrix plot: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

library(RColorBrewer)
library(corrplot)


# corrplot requires one dataframe of metadata and/or seqtab data.  Any columns in the dataframe will be correlated against all others.  Too many columns makes the graph unreadable.
# select to top 15 most abundant SVs from your phyloseq object and extract to a dataframe
#take out top 100 taxa based on abundance
topNOTUs = names(sort(taxa_sums(EX_ps_clean.rar), TRUE)[1:15])
top15 = prune_taxa(topNOTUs, EX_ps_clean.rar)

# combine with your metadata and create one dataframe object. 
# Coerce to data.frame and add the metadata for these samples
top15.sd = as(otu_table(top15), "matrix")
top15.sd = as.data.frame(top15.sd)

# add Genus names in place of the full sequence name that is in the SV columns
top15.tax <- as.data.frame(tax_table(top15))
## if the Genus is empty, replace with the Family
top15.tax$Genus[is.na(top15.tax$Genus)] <- top15.tax$Family[is.na(top15.tax$Genus)]
colnames(top15.sd) = top15.tax$Genus

# paste all the components together
EX_ps_clean.rar.corr.df <- cbind(top15.sd, EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# change any column factor names
names(EX_ps_clean.rar.corr.df)[names(EX_ps_clean.rar.corr.df) == "EX_ps_clean.rar.even"] <- "Evenness"

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

EX_ps_clean.rar.corr.df <- subset(EX_ps_clean.rar.corr.df, select = -c(Treatment, Diet, All, Location))

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# check that all remaining columns are numeric instead of factor or character
str(EX_ps_clean.rar.corr.df)

# clean up any columns which are not registering as numeric
EX_ps_clean.rar.corr.df <- sapply(EX_ps_clean.rar.corr.df, as.numeric)


# run correlations
ex_corr_calc <- cor(EX_ps_clean.rar.corr.df, 
                    use="complete.obs", # use=Specifies the handling of missing data
                    method="spearman") # correlation method=pearson, spearman, or kendall


### test significance of correlations
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(EX_ps_clean.rar.corr.df,0.95)
res2 <- cor.mtest(EX_ps_clean.rar.corr.df,0.99)


## plot only those correlations with a significant p-value <0.05, using hierarchical clustering
corrplot(ex_corr_calc, 
         type="lower", #shape of the plot itself: full, upper, lower
         method="circle", #shape of the data: circle, square, ellipse, number, shade, color, pie 
         order="hclust", #how to cluster samples: AOE, hclust, FPC, alphabet, or leave blank
         p.mat = res1[[1]], #which correlations to use
         sig.level=0.05, #sets significance cutoff
         insig="blank", #leaves p > 0.05 blank
         tl.col="black", # text color
         tl.cex=.7, #text size
         col=brewer.pal(n=10, name="RdYlBu")) #specify a color palette and number of shades needed
mtext("Correlation Matrix of Genus or Family Assignment of Top 15 Most Abundant SVs Using the Unite Taxonomic Database", 
      side = 3, line = -2, cex = 1)

corrplot(ex_corr_calc, 
         type="lower", #shape of the plot itself: full, upper, lower
         method="circle", #shape of the data: circle, square, ellipse, number, shade, color, pie 
         order="hclust", #how to cluster samples: AOE, hclust, FPC, alphabet, or leave blank
         p.mat = res1[[1]], #which correlations to use
         sig.level=0.05, #sets significance cutoff
         insig="blank", #leaves p > 0.05 blank
         tl.col="black", # text color
         tl.cex=.7, #text size
         col=brewer.pal(n=10, name="RdYlBu")) #specify a color palette and number of shades needed
mtext("Correlation Matrix of Genus or Family Assignment of Top 15 Most Abundant SVs Using the LSU_RDP Taxonomic Database", 
      side = 3, line = -2, cex = 1)

head(res1[[1]])
head(ex_corr_calc)
rownames(res1[[1]]) <- rownames(ex_corr_calc)
colnames(res1[[1]]) <- colnames(ex_corr_calc)
# Replace NA names with empty strings for ex_corr_calc
rownames(ex_corr_calc)[is.na(rownames(ex_corr_calc))] <- "Not Identified"
colnames(ex_corr_calc)[is.na(colnames(ex_corr_calc))] <- "Not Identified"

# Replace NA names with empty strings for res1[[1]]
rownames(res1[[1]])[is.na(rownames(res1[[1]]))] <- "Not Identified"
colnames(res1[[1]])[is.na(colnames(res1[[1]]))] <- "Not Identified"


## Barplots------------------
# can add ggplot components to make pretty
plot_bar(EX_ps_clean.rar, fill="Phylum") 
plot_bar(EX_ps_clean.rar2, fill="Phylum")

# agglomerate SVs by taxa level to get rid of all the blacklines on the graph. can be done for any level of taxonomy
EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Phylum")

EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar2, "Phylum")

plot_bar(EX_ps_clean.rar.glom, fill="Phylum")


plot_bar(EX_ps_clean.rar.glom, fill="Phylum") + 
  facet_grid(~All, space="free", scales="free") + 
  theme(
    legend.position = "bottom", 
    axis.text.x = element_blank(),        # Hides the x-axis labels
    axis.text.y = element_text(size = 8), # Change y-axis text size
    strip.text = element_text(size = 8)   # Change facet label text size
  ) + ggtitle("Phylum Abundance of Each Sample in each Sample Group for Cow Rumen ITS Unite Taxonomy Data ")

plot_bar(EX_ps_clean.rar.glom2, fill="Phylum") + 
  facet_grid(~All, space="free", scales="free") + 
  theme(
    legend.position = "bottom", 
    axis.text.x = element_blank(),        # Hides the x-axis labels
    axis.text.y = element_text(size = 8), # Change y-axis text size
    strip.text = element_text(size = 8)   # Change facet label text size
  ) + ggtitle("Phylum Abundance of Each Sample in each Sample Group for Cow Rumen ITS RDP_LSU Taxonomy Data ")

EX_ps_clean.rar.stacked = transform_sample_counts(EX_ps_clean.rar.glom, function(x) x / sum(x) )

EX_ps_clean.rar.stacked2 = transform_sample_counts(EX_ps_clean.rar.glom2, function(x) x / sum(x) )

plot_bar(EX_ps_clean.rar.stacked, fill="Phylum") + 
  facet_grid(~All, space="free", scales="free") + 
  theme(
    legend.position = "bottom", 
    axis.text.x = element_text(size = 6),        # Hides the x-axis labels
    axis.text.y = element_text(size = 8), # Change y-axis text size
    strip.text = element_text(size = 8)   # Change facet label text size
  ) + ggtitle("Relative Abundance of Each Sample in each Sample Group for Cow Rumen ITS Unite Taxonomy Data ") 

plot_bar(EX_ps_clean.rar.stacked2, fill="Phylum") + 
  facet_grid(~All, space="free", scales="free") + 
  theme(
    legend.position = "bottom", 
    axis.text.x = element_text(size = 6),        # Hides the x-axis labels
    axis.text.y = element_text(size = 8), # Change y-axis text size
    strip.text = element_text(size = 8)   # Change facet label text size
  ) + ggtitle("Relative Abundance of Each Sample in each Sample Group for Cow Rumen ITS RDP_LSU Taxonomy Data ") 

# to filter by abundance and pool low abundance groups: https://github.com/joey711/phyloseq/issues/901

## Core community members -----------
#code from this website: https://microbiome.github.io/tutorials/CoremicrobiotaAmplicon.html

#load libraries as needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

library(dplyr)
library(microbiome)


# grab the taxonomy to clean up first
tax.df <- data.frame(tax_table(EX_ps_clean.rar), stringsAsFactors=FALSE)
tax.df2 <- data.frame(tax_table(EX_ps_clean.rar2), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

tax.df2$Genus[is.na(tax.df2$Genus)] <- tax.df$Family[is.na(tax.df2$Genus)]


# put the new taxa back in
tax_table(EX_ps_clean.rar) <- as.matrix(tax.df)
tax_table(EX_ps_clean.rar2) <- as.matrix(tax.df2)

# check it worked
as.data.frame(tax_table(EX_ps_clean.rar))$Genus
as.data.frame(tax_table(EX_ps_clean.rar2))$Genus


#calculate compositional version of the data (relative abundances)
core_biomeW.rel <- microbiome::transform(EX_ps_clean.rar, "compositional") 
core_biomeW.rel2 <- microbiome::transform(EX_ps_clean.rar2, "compositional")

#This returns the taxa that exceed the given prevalence and minimum abundance detection thresholds. Set your preferred thresholds.
core_taxa_standardW <- core_members(core_biomeW.rel, detection = 1/10000, prevalence = 60/100) 
core_taxa_standardW2 <- core_members(core_biomeW.rel2, detection = 1/10000, prevalence = 60/100)

#A full phyloseq object of the core microbiota at those limits is obtained as follows
phylo.coreW <- core(core_biomeW.rel, detection = 1/10000, prevalence = .6)
phylo.coreW2 <- core(core_biomeW.rel2, detection = 1/10000, prevalence = .6)

###retrieving the associated taxa names from the phyloseq object and add it to what you just made.
core.taxaW <- taxa(phylo.coreW)
core.taxaW2 <- taxa(phylo.coreW2)
class(core.taxaW)
class(core.taxaW2)

# get the taxonomy data and assign it to what you just made.
tax.matW <- tax_table(phylo.coreW)
tax.matW2 <- tax_table(phylo.coreW2)

tax.dfW <- as.data.frame(tax.matW)
tax.dfW2 <- as.data.frame(tax.matW2)

# add the SVs to last column of what you just made.
tax.dfW$SV <- rownames(tax.dfW)
tax.dfW2$SV <- rownames(tax.dfW2)

# 67 observed of 8 variables

# select taxonomy of only those OTUs that are core members based on the thresholds that were used.
core.taxa.classW <- dplyr::filter(tax.dfW, rownames(tax.dfW) %in% core.taxaW)
knitr::kable(head(core.taxa.classW))

core.taxa.classW2 <- dplyr::filter(tax.dfW2, rownames(tax.dfW2) %in% core.taxaW2)
knitr::kable(head(core.taxa.classW2))

# save the list so you can access later, can report just the list
write.csv(core.taxa.classW, "/Users/toripalecek/Documents/AVS554/AVS554_Project/core.taxa.cow_rumen.csv")

write.csv(core.taxa.classW2, "/Users/toripalecek/Documents/AVS554/AVS554_Project/core.taxa.cow_rumen2.csv")

# graph the abundance of those shared taxa, here are some example: https://microbiome.github.io/tutorials/Core.html

#aggregate the genera 
plot.gen <- aggregate_taxa(phylo.coreW, "Genus")

plot.gen2 <- aggregate_taxa(phylo.coreW2, "Genus")

# load libraries as needed
library(ggplot2)
library(RColorBrewer)
library(viridis)

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)


plot_core(plot.gen, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = .5) + #CHANGE min prevalence
  xlab("Detection Threshold (Relative Abundance (%))") + ylab("Fungal SVs") +
  theme_minimal() + scale_fill_viridis()

plot_core(plot.gen2, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = .5) + #CHANGE min prevalence
  xlab("Detection Threshold (Relative Abundance (%))") + ylab("Fungal SVs") +
  theme_minimal() + scale_fill_viridis()


# Comparing changes in taxonomy---------------------------------
## Focusing on a single taxon -------
library(phyloseq) 
library(dplyr) 
library(tidyr)  
library(ggplot2)

# transform to relative abundance
relabun.ps <- transform_sample_counts(EX_ps_clean.rar,function(x) x / sum(x)) 

relabun.ps2 <- transform_sample_counts(EX_ps_clean.rar2,function(x) x / sum(x)) 

# glom ASVs by genus
ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE) 

ps_genus2 <- tax_glom(relabun.ps2, taxrank = "Genus", NArm = FALSE) 

#subset by the genus of choice. Taxon is not present if you get error: "length of 'dimnames' [1] not equal to array extent" 
ps_genusP <- subset_taxa(ps_genus, Genus %in% c("g__Neocallimastix")) 

ps_genusP2 <- subset_taxa(ps_genus2, Genus %in% c("g__Neocallimastix")) 

# melt the data into a different configuration
genus.df <- psmelt(ps_genusP) 

genus.df2 <- psmelt(ps_genusP2) 

# grab that abundance data. 
MySummary <- genus.df %>% group_by(All, Treatment, Genus) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 

MySummary2 <- genus.df2 %>% group_by(All, Treatment, Genus) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 

#check it
head(MySummary)

head(MySummary2)

# graph it
ggplot(data=MySummary, aes(x = Treatment, y = mean_abund)) +
  geom_point(aes(color=All)) + ylab("Mean relative abundance of reads") +ggtitle("Relative Abundance")


ggplot(data=MySummary2, aes(x = Treatment, y = mean_abund)) +
  geom_point(aes(color=All)) + ylab("Mean relative abundance of reads") +ggtitle("Relative Abundance")

table(tax_df$Species)
table(tax_df2$Species)

## Differential Abundance with DESeq ---------------------------------

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
packageVersion("DESeq2")

# subset 
EX_ps_clean_subsetA = subset_samples(EX_ps, FACTOR == "Diet")
EX_ps_clean_subsetA <- prune_samples(sample_sums(EX_ps_clean_subsetA) > 0, EX_ps_clean_subsetA)
EX_ps_clean_subsetA <- prune_taxa(taxa_sums(EX_ps_clean_subsetA) > 0, EX_ps_clean_subsetA)


# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(EX_ps_clean, ~ Diet) 

diagdds2 = phyloseq_to_deseq2(EX_ps_clean2, ~ Diet)
# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans2 = apply(counts(diagdds2), 1, gm_mean)
diagdds2 = estimateSizeFactors(diagdds2, geoMeans = geoMeans)
diagdds2 = DESeq(diagdds2, fitType="local")

# calculate significance for those abundance calculations
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(EX_ps_clean)[rownames(sigtab), ], "matrix"))

head(sigtab)

dim(sigtab)

res2 = results(diagdds2)
res2 = res2[order(res2$padj, na.last=NA), ]
alpha = 0.01
sigtab2 = res2[(res2$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(EX_ps_clean2)[rownames(sigtab2), ], "matrix"))

head(sigtab2)

dim(sigtab2)


# calculate log changes and set
sigtab = sigtab[sigtab[, "log2FoldChange"] < 0, ] 
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")] 

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# if the Genus is empty, replace with the Family
sigtab$Genus = ifelse(is.na(sigtab$Genus), paste(sigtab$Family),paste(sigtab$Genus));sigtab


library("ggplot2")

## graph differential abundance
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Phylum)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) + xlab("log2 Fold Change") + labs(size = "Mean Sequence Abundance") + theme_minimal()


## Feature Prediction (Differential Abundance) with Random Forest ---------------------------------
#  https://rpubs.com/michberr/randomforestmicrobe
# https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
# https://www.rdocumentation.org/packages/randomForest/versions/4.6-12/topics/randomForest

library(phyloseq)
library(vegan)
library(plyr)
library(dplyr)
require(magrittr)
require(scales)
require(grid)
require(reshape2)
require(knitr)
library(randomForest)
library(rfPermute)
library(ggplot2)
library(RColorBrewer)

# Make a dataframe of training data with OTUs as column and samples as rows, which is the phyloseq OTU table
predictors <- otu_table(EX_ps_clean)

predictors2 <- otu_table(EX_ps_clean2)

dim(predictors)

dim(predictors2)
#samples and SVs includes: 34 samples 11491 SVs
#samples and SVs includes: 34 samples 11491 SVs

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

tax.df2 <- data.frame(tax_table(EX_ps_clean2), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

tax.df2$Genus[is.na(tax.df2$Genus)] <- tax.df2$Family[is.na(tax.df2$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)

tax.df2$Genus.species <- paste(tax.df2$Genus, tax.df2$Species)


# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(predictors) <- tax.df$Genus.species

colnames(predictors2) <- tax.df2$Genus.species

# clean up some of the other taxa info
colnames(predictors) = gsub("_unclassified", "", colnames(predictors))
colnames(predictors) = gsub("_Intercertae_Sedis", "", colnames(predictors))

colnames(predictors2) = gsub("_unclassified", "", colnames(predictors2))
colnames(predictors2) = gsub("_Intercertae_Sedis", "", colnames(predictors2))

# Make one column for our outcome/response variable.
response <- as.factor(sample_data(EX_ps_clean)$Factorial_Data) 
response <- as.numeric(sample_data(EX_ps_clean)$Numeric_Data) 

response <- as.factor(sample_data(EX_ps_clean)$All)

response2 <- as.factor(sample_data(EX_ps_clean2)$All)

# Combine response and SVs into data frame
rf.data <- data.frame(response, predictors)

rf.data2 <- data.frame(response2, predictors2)


# set seed for random number generation reproducibility
set.seed(2)

# classify for factorial data
response.pf <- rfPermute(response ~. , data = rf.data, na.action = na.omit, ntree= 500, nrep = 100) #na.omit ignores NAs in data (not tolerated). ntrees is how many forests to build, nreps generates p-value

response.pf2 <- rfPermute(response2 ~. , data = rf.data2, na.action = na.omit, ntree= 500, nrep = 100)

# or this way for numeric data
response.pf <- rfPermute(as.numeric(response) ~. , data = rf.data, na.action = na.omit, ntree= 500, nrep = 100)

print(response.pf)

print(response.pf2)
# paste the print out here, especially the OOB error. 1-(Out-of-the-box error) = accuracy of your model



# An rfPermute model

# Type of random forest: classification 
# Number of trees: 500 
# No. of variables tried at each split: 107 
# No. of permutation replicates: 100 
# Start time: 2024-11-09 10:36:34 
# End time: 2024-11-09 10:39:20 
# Run time: 2.76 mins 

# HFControlEpimural HFControlFluid HFYeastEpimural HFYeastFluid HGControlEpimural HGControlFluid
# HFControlEpimural                 0              0               5            0                 0              0
# HFControlFluid                    1              1               2            0                 0              0
# HFYeastEpimural                   2              1               4            0                 0              0
# HFYeastFluid                      0              1               0            0                 0              0
# HGControlEpimural                 0              0               0            0                 3              1
# HGControlFluid                    0              0               0            0                 0              6
# HGYeastEpimural                   0              0               1            0                 1              0
# HGYeastFluid                      0              0               0            0                 0              4
# Overall                          NA             NA              NA           NA                NA             NA
# HGYeastEpimural HGYeastFluid pct.correct LCI_0.95 UCI_0.95
# HFControlEpimural               0            0         0.0    0.000     52.2
# HFControlFluid                  0            0        25.0    0.631     80.6
# HFYeastEpimural                 0            0        57.1   18.405     90.1
# HFYeastFluid                    0            0         0.0    0.000     97.5
# HGControlEpimural               0            0        75.0   19.412     99.4
# HGControlFluid                  0            1        85.7   42.128     99.6
# HGYeastEpimural                 0            0         0.0    0.000     84.2
# HGYeastFluid                    0            0         0.0    0.000     60.2
# Overall                        NA           NA        41.2   24.647     59.3


# An rfPermute model 2

# Type of random forest: classification 
# Number of trees: 500 
# No. of variables tried at each split: 107 
# No. of permutation replicates: 100 
# Start time: 2024-11-09 10:40:15 
# End time: 2024-11-09 10:43:03 
# Run time: 2.81 mins 

# HFControlEpimural HFControlFluid HFYeastEpimural HFYeastFluid HGControlEpimural HGControlFluid
# HFControlEpimural                 0              0               5            0                 0              0
# HFControlFluid                    0              3               1            0                 0              0
# HFYeastEpimural                   1              1               5            0                 0              0
# HFYeastFluid                      0              1               0            0                 0              0
# HGControlEpimural                 0              0               0            0                 3              1
# HGControlFluid                    0              0               0            0                 0              6
# HGYeastEpimural                   0              0               0            0                 1              1
# HGYeastFluid                      0              0               0            0                 0              4
# Overall                          NA             NA              NA           NA                NA             NA
# HGYeastEpimural HGYeastFluid pct.correct LCI_0.95 UCI_0.95
# HFControlEpimural               0            0         0.0      0.0     52.2
# HFControlFluid                  0            0        75.0     19.4     99.4
# HFYeastEpimural                 0            0        71.4     29.0     96.3
# HFYeastFluid                    0            0         0.0      0.0     97.5
# HGControlEpimural               0            0        75.0     19.4     99.4
# HGControlFluid                  0            1        85.7     42.1     99.6
# HGYeastEpimural                 0            0         0.0      0.0     84.2
# HGYeastFluid                    0            0         0.0      0.0     60.2
# Overall                        NA           NA        50.0     32.4     67.6




# grab which features were labeled "important"
imp <- importance(response.pf, scale = TRUE)

imp2 <- importance(response.pf2, scale = TRUE)

# Make a data frame with predictor names and their importance
imp.df <- data.frame(predictors = rownames(imp), imp) 

imp.df2 <- data.frame(predictors2 = rownames(imp2), imp2) 


# For factorial data, grab only those features with p-value < 0.05
imp.sig <- subset(imp.df, MeanDecreaseAccuracy.pval <= 0.05) 
print(dim(imp.sig))

imp.sig2 <- subset(imp.df2, MeanDecreaseAccuracy.pval <= 0.05) 
print(dim(imp.sig2))


# how many SVs (sig features) were identified: 214 (214, 21)
# how many SVs (sig features) were identified. 214 (214, 21)

# For numeric data, grab only those features with p-value < 0.05
imp.sig <- subset(imp.df, IncNodePurity.pval <= 0.05)
print(dim(imp.sig))

# or For factorial data, sort by importance amount
imp.sort <- imp.sig[order(imp.sig$MeanDecreaseAccuracy),]

imp.sort2 <- imp.sig2[order(imp.sig2$MeanDecreaseAccuracy),]

# For numeric data, sort by importance amount
imp.sort <- imp.sig[order(imp.sig$IncNodePurity),]

#create levels to the factor based on SV table
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

imp.sort2$predictors <- factor(imp.sort2$predictors, levels = imp.sort2$predictors)

# Select the top so many predictors (more than 50 is a crowded graph)
imp.top <- imp.sort[1:50, ]

imp.top2 <- imp.sort2[1:50, ]

# figure out what they are and name them, otherwise will just be named the full SV
otunames <- imp.top$predictors

otunames2 <- imp.top2$predictors


# need to grab the abundance data out of the otu_table and make a new data.frame, and add the taxa names back in
pred.abun.colnum <- which(colnames(rf.data) %in% otunames)

pred.abun.colnum2 <- which(colnames(rf.data2) %in% otunames2)

# when you find a match, grad the abudnance data
pred.abun <- rf.data[,sort(c(pred.abun.colnum))]

pred.abun2 <- rf.data2[,sort(c(pred.abun.colnum2))]

# make this into a dataframe for manipulation
pred.abun.df <- data.frame(pred.abun, stringsAsFactors=FALSE)

pred.abun.df2 <- data.frame(pred.abun2, stringsAsFactors=FALSE)

# use the row.names (sample names) from the phyloseq object to name the samples in your forest
row.names(pred.abun.df) <- row.names(sample_data(EX_ps_clean))

row.names(pred.abun.df2) <- row.names(sample_data(EX_ps_clean2))

# set color palette  
col.bro <- (rainbow(6))

# add white to that color list
col.bro <- append(col.bro, "#ffffff", after = 6)

# add some factors that you can use to make your graph pretty
pred.abun.df$Sample <- row.names(sample_data(EX_ps_clean)) 
pred.abun.df$All <- sample_data(EX_ps_clean)$All 
pred.abun.df$Treatment <- sample_data(EX_ps_clean)$Treatment 
pred.abun.df$Diet <- sample_data(EX_ps_clean)$Diet 

pred.abun.df2$Sample <- row.names(sample_data(EX_ps_clean2)) 
pred.abun.df2$All <- sample_data(EX_ps_clean2)$All 
pred.abun.df2$Treatment <- sample_data(EX_ps_clean2)$Treatment 
pred.abun.df2$Diet <- sample_data(EX_ps_clean2)$Diet 

# reload these packages in this order, because sometimes the ddply function breaks
library(plyr)
library(dplyr)

# melt and transform the data using ALL the factors you added
m <- melt(pred.abun.df,id.vars= c("Sample", "All", "Treatment", "Diet"))
m <- ddply(m, .(variable), transform, rescale = log(1 + value)) 
m2 <- melt(pred.abun.df2, id.vars=c("Sample", "All", "Treatment", "Diet"))
m2 <- ddply(m2, .(variable), transform, rescale = log(1 + value))

# adjusted plot for factorial data, recommend using sample ID as 'factor A'
ggplot(m, aes(as.factor(Sample), variable)) + 
  theme_minimal() + 
  facet_grid(.~Sample, space="free", scale="free") + #set up graph facets to separate out levels of FactorA
  geom_tile(aes(fill = rescale), color="gray") + #add the heatmap coloring 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + #use the preset color palette
  labs(fill="Log abundance") + #rename legend heading
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, size = 8),
        axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 6)) +
  ylab('Predictor Taxa') +
  xlab('Sample') + ggtitle("Random Forest Model for Cow Rumen ITS Unite Taxonomy Data ")  

ggplot(m2, aes(as.factor(Sample), variable)) + 
  theme_minimal() + 
  facet_grid(.~Sample, space="free", scale="free") + #set up graph facets to separate out levels of FactorA
  geom_tile(aes(fill = rescale), color="gray") + #add the heatmap coloring 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + #use the preset color palette
  labs(fill="Log abundance") + #rename legend heading
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, size = 8),
        axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 6)) +
  ylab('Predictor Taxa') +
  xlab('Sample') + ggtitle("Random Forest Model for Cow Rumen ITS RDP_LSU Taxonomy Data")  


# adjusted plot for numeric data, recommend using sample ID as 'factor A'
for_moist <-ggplot(m, aes(FactorA, variable)) + 
  theme_minimal() + 
  geom_tile(aes(fill = rescale), color="gray") + 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + 
  labs(fill="Log abundance") +
  theme(legend.position = 'bottom', # could be none, top, right, left
        axis.text.x = element_text(angle = 0, size = 10),
        #axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 10))+
  # scale_x_discrete(breaks=c(0,24.4, 50.3, 74.8, 98.3), labels=c("0" = "0", "24.4" ="25", "50.3" = "50", "74.8" = "75", "98.3" ="100")) + #optional, may want to add x-axis scale breaks
  ylab('Predictor Taxa') +
  xlab('Factor A') 



## LEFSe ---------------------------------
# uses code by seahorse001x: https://github.com/seashore001x/Rrumen/blob/master/phyloseq2lefse.R

# require(dplyr)
# require(tibble)
# 

# phyloseq_to_lefse <- function(physeq){
#   # aggregate to genus level
#   ps_lefse <- physeq %>% tax_glom(taxrank = 'Genus', NArm = F)
#   
#   # extract taxonomic data from phyloseq object and then stored in a vector called lefse_tax
#   lefse_tax <- ps_lefse %>% tax_table %>% data.frame(stringsAsFactors=FALSE)
#   lefse_tax <- replace(lefse_tax, is.na(lefse_tax), 'Unclassified')
#   lefse_tax <- lefse_tax %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% mutate(id = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "|")) %>% ungroup %>% pull(id)
#   
#   # extract otu abundance matrix from phyloseq object and annotated with tax information
#   lefse_matrix <- otu_table(ps_lefse) %>% data.frame(stringsAsFactors = F) %>% t %>% data.frame
#    
#   # bookmark this is what is throwing the error
# #  colnames(lefse_matrix) <- lefse_tax
#   
# #  row.names(lefse_matrix) <- lefse_tax
#   
#   # extract sample metadata and order sample same in lefse_matrix
#   lefse_sample <- sample_data(ps_lefse)
#   
#   
#   # convert factor in lefse_sample into character in order to combine sample and abundance data
#   lefse_sample_isfactor <- sapply(lefse_sample, is.factor)
#   lefse_sample[,lefse_sample_isfactor] <- lefse_sample[,lefse_sample_isfactor] %>% lapply(as.character)
#   lefse_sample <- lefse_sample %>% data.frame
#   
#   lefse_table <- full_join(rownames_to_column(lefse_sample), rownames_to_column(lefse_matrix), by = ("rowname" = "rowname")) %>% t
#   
#   return(lefse_table)
# }
# 
# 
# EX_clean_for_lefse <- phyloseq_to_lefse(EX_ps_clean)




predictors <- otu_table(EX_ps_clean)

dim(predictors)

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)


# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(predictors) <- tax.df$Genus.species

# clean up some of the other taxa info
colnames(predictors) = gsub("_unclassified", "", colnames(predictors))
colnames(predictors) = gsub("_Intercertae_Sedis", "", colnames(predictors))

# Make one column for outcome/response variable.
response.class <- as.factor(sample_data(EX_ps_clean)$All) 
#response.subclass <- as.factor(sample_data(EX_ps_clean)$Diet)

EX_ps_clean_samples <- as.data.frame(sample_data(EX_ps_clean))

EX_ps_clean_samples$Sample <- rownames(EX_ps_clean_samples)

subjects <- as.factor(sample_data(EX_ps_clean_samples)$Sample)

Ex_ps_for_lefse.df <- data.frame(response.class, subjects, predictors)

length(response.class)      # Number of rows in response.class
length(response.subclass)   # Number of rows in response.subclass
length(subjects)            # Number of rows in subjects
dim(predictors) 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiomeMarker")
library(microbiomeMarker)

EX_ps_genus <- tax_glom(EX_ps_clean, taxrank = "Genus")

EX_ps_genus2 <- tax_glom(EX_ps_clean2, taxrank = "Genus")

# Continue with LEfSe analysis and plotting

# Run LEfSe with a subgroup
lefse_results <- run_lefse(
  ps = EX_ps_genus,                # Phyloseq object
  group = "Diet",                  # Main grouping variable (e.g., "Diet"),
  taxa_rank = "Genus",             # Taxonomic rank (e.g., "Phylum", "Genus")
  norm = "CPM",                    # Normalization method
  kw_cutoff = 0.05,                # p-value cutoff for Kruskal-Wallis test
  lda_cutoff = 2,                  # LDA score cutoff
  bootstrap_n = 30,                # Number of bootstrap iterations
  bootstrap_fraction = 2/3         # Fraction of samples to bootstrap
)


lefse_results2 <- run_lefse(
  ps = EX_ps_genus2,                # Phyloseq object
  group = "Diet",                  # Main grouping variable (e.g., "Diet"),
  taxa_rank = "Genus",             # Taxonomic rank (e.g., "Phylum", "Genus")
  norm = "CPM",                    # Normalization method
  kw_cutoff = 0.05,                # p-value cutoff for Kruskal-Wallis test
  lda_cutoff = 2,                  # LDA score cutoff
  bootstrap_n = 30,                # Number of bootstrap iterations
  bootstrap_fraction = 2/3         # Fraction of samples to bootstrap
)



lefse_significant <- lefse_results@marker_table

lefse_significant2 <- lefse_results2@marker_table

ggplot(lefse_significant, aes(x = reorder(feature, ef_lda), y = ef_lda, fill = ef_lda > 0)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Genus",
    y = "LDA Score",
    title = "LEfSe Results for Unite Taxonoimic Database Assignment"
  ) +
  theme(axis.text.y = element_text(size = 5))


ggplot(lefse_significant, aes(x = reorder(feature, ef_lda), y = ef_lda, fill = ef_lda > 0)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Genus",
    y = "LDA Score",
    title = "LEfSe Results for RDP LSU Taxonoimic Database Assignment"
  ) +
  theme(axis.text.y = element_text(size = 5))

# Extract significant features
lefse_significant <- lefse_results@marker_table


table(sample_data(EX_ps_clean_samples)$All)


# Run LEfSe analysis
lefse_results <- run_lefse(
  ps = EX_ps_clean,  # Your phyloseq object
  group = "Diet",             # Your grouping variable (e.g., diet or treatment)
  taxa_rank = "all",         # Taxonomic rank to use for analysis, can change if needed
  transform = "identity",    # Optional: transformation of the data (e.g., "log10")
  norm = "CPM",              # Normalization method
  kw_cutoff = 0.05,          # Kruskal-Wallis p-value cutoff
  lda_cutoff = 2.0,          # LDA cutoff for selecting significant features
  wilcoxon_cutoff = 0.05,    # Wilcoxon p-value cutoff for subgroup analysis
  sample_min = 10            # Minimum number of samples per group
)

# View the results
lefse_results

lefse_results_df <- lefse_results@marker_table

# Subset the significant results based on LDA score and p-value
lefse_significant <- lefse_results_df[lefse_results_df$ef_lda > 2 & lefse_results_df$pvalue < 0.05, ]

str(lefse_results)


ggplot(lefse_significant, aes(x = reorder(feature, ef_lda), y = ef_lda, fill = ef_lda > 0)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Features (OTUs/ASVs)",
    y = "LDA Score",
    title = "LEfSe Results: Most Discriminative Features"
  ) +
  theme(
    axis.text.y = element_text(size = 5) 
  )

# save this output and upload it to the LEFSe web version on Galaxy: http://huttenhower.sph.harvard.edu/galaxy/
write.table(Ex_ps_for_lefse.df, file = "Ex_ps_for_lefse.df.txt", 
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#  Beta diversity ordinations and stats ---------------------------------
# ordinations in phyloseq: https://joey711.github.io/phyloseq/plot_ordination-examples.html

## PCA  -----------------
# currently, phyloseq doesn't run a PCA, but you can manually perform one using this tutorial: https://www.datacamp.com/community/tutorials/pca-analysis-r

install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot") #takes a few minutes
library(ggbiplot)

# calculate the components
EX_clean.rar.pca <- prcomp(otu_table(EX_ps_clean.rar), center = TRUE, scale = TRUE) 

EX_clean.rar.pca2 <- prcomp(otu_table(EX_ps_clean.rar2), center = TRUE, scale = TRUE)

# take a look
summary(EX_clean.rar.pca) 

summary(EX_clean.rar.pca2)

# graph it. add ggplot2 text to make it pretty
ggbiplot(EX_clean.rar.pca)

ggbiplot(EX_clean.rar.pca2)




## PCoA in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
EX_uJ_pcoa <- ordinate(EX_ps_clean.rar, #calculate similarities
                       method ="PCoA", #ordination type
                       "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.

EX_uJ_pcoa2 <- ordinate(EX_ps_clean.rar2, #calculate similarities
                        method ="PCoA", #ordination type
                        "jaccard", binary = TRUE)

# simple ordination
plot_ordination(EX_ps_clean.rar, EX_uJ_pcoa, type="samples", color="All")

# fancy ordination
library(ggplot2) # can add components on top of the phyloseq graph
install.packages("viridis")
library(viridis) #add a cool color palette

plot_ordination(EX_ps_clean.rar, EX_uJ_pcoa, type="samples") + 
  geom_point(aes(size = as.factor(Treatment), color=as.factor(All), shape=as.factor(Diet))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  theme_minimal()+
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=as.factor(All), linetype=as.factor(All))) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Sample Group", size = "Treatment", shape = "Diet", linetype = "Sample Group") + ggtitle("PCoA Ordination with Jaccard Similarity Unite Data")

plot_ordination(EX_ps_clean.rar2, EX_uJ_pcoa, type="samples") + 
  geom_point(aes(size = as.factor(Treatment), color=as.factor(All), shape=as.factor(Diet))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  theme_minimal()+
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=as.factor(All), linetype=as.factor(All))) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Sample Group", size = "Treatment", shape = "Diet", linetype = "Sample Group") + ggtitle("PCoA Ordination with Jaccard Similarity RDP_LSU Data")


## NMDS in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
EX_uJ_nmds <- ordinate(EX_ps_clean.rar, #calculate similarities
                       method ="NMDS", #ordination type
                       "jaccard", binary = TRUE) #similarity type

EX_uJ_nmds2 <- ordinate(EX_ps_clean.rar2, #calculate similarities
                        method ="NMDS", #ordination type
                        "jaccard", binary = TRUE)

# Run 20 stress 0.145 for Unite
# Run 20 stress 0.146 for RDP_LSU

# simple ordination
plot_ordination(EX_ps_clean.rar, EX_uJ_nmds, type="samples", color="All") + ggtitle("NMDS Unite Alignment Data") 

plot_ordination(EX_ps_clean.rar2, EX_uJ_nmds2, type="samples", color="All") + ggtitle("NMDS RDP_LSU Alignment Data") 

# fancy ordination
library(ggplot2)
install.packages("viridis")
library(viridis) #add a cool color palette

plot_ordination(EX_ps_clean.rar, EX_uJ_nmds, type="samples", shape="Diet") + 
  geom_point(aes(size = as.factor(Treatment), color=as.factor(All))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=All, linetype=All)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Sample Group", size = "Treatment", shape = "Diet", linetype="Sample Group") + ggtitle("NMDS Ordination Unite Alignment Data") 

plot_ordination(EX_ps_clean.rar2, EX_uJ_nmds2, type="samples", shape="Diet") + 
  geom_point(aes(size = as.factor(Treatment), color=as.factor(All))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=All, linetype=All)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Sample Group", size = "Treatment", shape = "Diet", linetype="Sample Group") + ggtitle("NMDS Ordination RDP_LSU Alignment Data") 




### ----- permanova stats for ordinations --------
library(vegan)
# Example basic permanova test using Adonis in vegan
adonis2(EX_uJ_nmds ~ All, as(sample_data(EX_ps_clean.rar), "data.frame"), permutations=9999, na.action=na.omit, by = "terms") 

ctrl <- with(as.data.frame(sample_data(EX_ps_clean.rar)), how(blocks = Sheep_ID, nperm = 999))

adonis2(EX_uJ_nmds ~ All, data = data.frame(sample_data(EX_ps_clean.rar)), 
        permutations = ctrl, #for repeated measures
        # strata = as.data.frame(sample_data(EX_ps_clean.rar))$Blocking_factor 
) 


# Code to use Vegan to calculation ordination data and run stats, in case phyloseq and vegan aren't sharing code well.
# function and tutorial from https://rpubs.com/DKCH2020/587758
library(vegan)

veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}
# export data from phyloseq to vegan-compatible object
EX_vegan <- veganotu(EX_ps_clean.rar)

EX_vegan2 <- veganotu(EX_ps_clean.rar2)

# make a data frame that can be used in vegan from the sample_data in the phyloseq object
sampledf <- data.frame(sample_data(EX_ps_clean.rar))

sampledf2 <- data.frame(sample_data(EX_ps_clean.rar2))

# run the ordination calculation, change the variable name from EX_uJ to reflect the calculation method you choose (EX_wBC)
EX_uJ <- vegdist(wisconsin(sqrt(EX_vegan)), method = "jaccard", binary=TRUE) 

adonis2(EX_uJ ~ All, as(sample_data(EX_ps_clean.rar), "data.frame"), permutations=9999, na.action=na.omit)


EX_uJ2 <- vegdist(wisconsin(sqrt(EX_vegan2)), method = "jaccard", binary=TRUE) 

adonis2(EX_uJ2 ~ All, as(sample_data(EX_ps_clean.rar2), "data.frame"), permutations=9999, na.action=na.omit)

# Permutation test for adonis under reduced model Treatment Only 
# Permutation: free
# Number of permutations: 9999

# adonis2(formula = EX_uJ ~ Treatment, data = as(sample_data(EX_ps_clean.rar), "data.frame"), permutations = 9999, na.action = na.omit)
#          Df SumOfSqs      R2    F Pr(>F)
# Model     1   0.3896 0.03382 0.98   0.42
# Residual 28  11.1308 0.96618            
# Total    29  11.5204 1.00000 

# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 9999

# adonis2(formula = EX_uJ ~ All, data = as(sample_data(EX_ps_clean.rar), "data.frame"), permutations = 9999, na.action = na.omit)
#          Df SumOfSqs      R2      F Pr(>F)    
# Model     7   4.1789 0.36274 1.7889  1e-04 ***
# Residual 22   7.3416 0.63726                  
# Total    29  11.5204 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Permutation test for adonis under reduced model 2
# Permutation: free
# Number of permutations: 9999

# adonis2(formula = EX_uJ2 ~ Treatment, data = as(sample_data(EX_ps_clean.rar2), "data.frame"), permutations = 9999, na.action = na.omit)
#          Df SumOfSqs      R2    F Pr(>F)
# Model     1   0.3896 0.03382 0.98 0.4175
# Residual 28  11.1308 0.96618            
# Total    29  11.5204 1.00000 


# Permutation test for adonis under reduced model 2
# Permutation: free
# Number of permutations: 9999

# adonis2(formula = EX_uJ2 ~ All, data = as(sample_data(EX_ps_clean.rar2), "data.frame"), permutations = 9999, na.action = na.omit)
#          Df SumOfSqs      R2      F Pr(>F)    
# Model     7   4.1789 0.36274 1.7889  1e-04 ***
# Residual 22   7.3416 0.63726                  
# Total    29  11.5204 1.00000                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#### run betadispersion test to see how tight/loose the clusters are ----
betadisp_EX <- betadisper(EX_uJ, sampledf$All) 

betadisp_EX

betadisp_EX2 <- betadisper(EX_uJ2, sampledf2$All) 


# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = EX_uJ, group = sampledf$All)
# 
# No. of Positive Eigenvalues: 29
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   HFControlEpimural    HFControlFluid   HFYeastEpimural      HFYeastFluid HGControlEpimural    HGControlFluid   HGYeastEpimural      HGYeastFluid 
# 0.5149            0.4830            0.5408            0.0000            0.4724            0.5126            0.4413            0.4889 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 29 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 1.5959 0.7949 0.5162 0.4906 0.4090 0.3999 0.3839 0.3771 



#Homogeneity of multivariate dispersions

# Call: betadisper(d = EX_uJ2, group = sampledf2$All)
# 
# No. of Positive Eigenvalues: 29
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   HFControlEpimural    HFControlFluid   HFYeastEpimural      HFYeastFluid HGControlEpimural    HGControlFluid   HGYeastEpimural      HGYeastFluid 
# 0.5149            0.4830            0.5408            0.0000            0.4724            0.5126            0.4413            0.4889 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 29 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 1.5959 0.7949 0.5162 0.4906 0.4090 0.3999 0.3839 0.3771 

# run ANOVA to see if clusters overlap or not
anova(betadisp_EX)
anova(betadisp_EX2)

# Analysis of Variance Table
# 
# Response: Distances
# Df   Sum Sq  Mean Sq F value   Pr(>F)    
# Groups     7 0.266561 0.038080  48.918 5.86e-12 ***
#   Residuals 22 0.017126 0.000778                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Analysis of Variance Table 2
# 
# Response: Distances
# Df   Sum Sq  Mean Sq F value   Pr(>F)    
# Groups     7 0.266561 0.038080  48.918 5.86e-12 ***
#   Residuals 22 0.017126 0.000778                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                  

citation("vegan")

permutest(betadisp_EX)
permutest(betadisp_EX2)
# paste output here

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     7 0.266561 0.038080 48.918    999  0.001 ***
#   Residuals 22 0.017126 0.000778                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Permutation test for homogeneity of multivariate dispersions 2
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     7 0.266561 0.038080 48.918    999  0.001 ***
#   Residuals 22 0.017126 0.000778                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#  correct for multiple comparisons
TukeyHSD(betadisp_EX)
TukeyHSD(betadisp_EX2)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
# diff         lwr          upr     p adj
# HFControlFluid-HFControlEpimural    -0.031912225 -0.10305875  0.039234299 0.8007046
# HFYeastEpimural-HFControlEpimural    0.025859246 -0.03427054  0.085989033 0.8311845
# HFYeastFluid-HFControlEpimural      -0.514903767 -0.61905161 -0.410755921 0.0000000
# HGControlEpimural-HFControlEpimural -0.042521088 -0.10838997  0.023347793 0.4133202
# HGControlFluid-HFControlEpimural    -0.002270582 -0.06240037  0.057859204 1.0000000
# HGYeastEpimural-HFControlEpimural   -0.073560725 -0.15423330  0.007111849 0.0912239
# HGYeastFluid-HFControlEpimural      -0.025976538 -0.09184542  0.039892343 0.8828953
# HFYeastEpimural-HFControlFluid       0.057771471 -0.00809741  0.123640352 0.1143935
# HFYeastFluid-HFControlFluid         -0.482991543 -0.59055498 -0.375428110 0.0000000
# HGControlEpimural-HFControlFluid    -0.010608864 -0.08175539  0.060537659 0.9995548
# HGControlFluid-HFControlFluid        0.029641642 -0.03622724  0.095510523 0.7981876
# HGYeastEpimural-HFControlFluid      -0.041648500 -0.12668486  0.043387860 0.7253734
# HGYeastFluid-HFControlFluid          0.005935687 -0.06521084  0.077082210 0.9999910
# HFYeastFluid-HFYeastEpimural        -0.540763014 -0.64137939 -0.440146635 0.0000000
# HGControlEpimural-HFYeastEpimural   -0.068380335 -0.12851012 -0.008250548 0.0184833
# HGControlFluid-HFYeastEpimural      -0.028129829 -0.08191155  0.025651888 0.6593684
# HGYeastEpimural-HFYeastEpimural     -0.099419971 -0.17547880 -0.023361139 0.0050729
# HGYeastFluid-HFYeastEpimural        -0.051835784 -0.11196557  0.008294002 0.1258519
# HGControlEpimural-HFYeastFluid       0.472382679  0.36823483  0.576530525 0.0000000
# HGControlFluid-HFYeastFluid          0.512633185  0.41201681  0.613249563 0.0000000
# HGYeastEpimural-HFYeastFluid         0.441343042  0.32725479  0.555431291 0.0000000
# HGYeastFluid-HFYeastFluid            0.488927229  0.38477938  0.593075075 0.0000000
# HGControlFluid-HGControlEpimural     0.040250506 -0.01987928  0.100380293 0.3700881
# HGYeastEpimural-HGControlEpimural   -0.031039637 -0.11171221  0.049632938 0.8950668
# HGYeastFluid-HGControlEpimural       0.016544550 -0.04932433  0.082413432 0.9885966
# HGYeastEpimural-HGControlFluid      -0.071290143 -0.14734898  0.004768690 0.0769739
# HGYeastFluid-HGControlFluid         -0.023705956 -0.08383574  0.036423831 0.8830533
# HGYeastFluid-HGYeastEpimural         0.047584187 -0.03308839  0.128256762 0.5223260


# Tukey multiple comparisons of means 2
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
# diff         lwr          upr     p adj
# HFControlFluid-HFControlEpimural    -0.031912225 -0.10305875  0.039234299 0.8007046
# HFYeastEpimural-HFControlEpimural    0.025859246 -0.03427054  0.085989033 0.8311845
# HFYeastFluid-HFControlEpimural      -0.514903767 -0.61905161 -0.410755921 0.0000000
# HGControlEpimural-HFControlEpimural -0.042521088 -0.10838997  0.023347793 0.4133202
# HGControlFluid-HFControlEpimural    -0.002270582 -0.06240037  0.057859204 1.0000000
# HGYeastEpimural-HFControlEpimural   -0.073560725 -0.15423330  0.007111849 0.0912239
# HGYeastFluid-HFControlEpimural      -0.025976538 -0.09184542  0.039892343 0.8828953
# HFYeastEpimural-HFControlFluid       0.057771471 -0.00809741  0.123640352 0.1143935
# HFYeastFluid-HFControlFluid         -0.482991543 -0.59055498 -0.375428110 0.0000000
# HGControlEpimural-HFControlFluid    -0.010608864 -0.08175539  0.060537659 0.9995548
# HGControlFluid-HFControlFluid        0.029641642 -0.03622724  0.095510523 0.7981876
# HGYeastEpimural-HFControlFluid      -0.041648500 -0.12668486  0.043387860 0.7253734
# HGYeastFluid-HFControlFluid          0.005935687 -0.06521084  0.077082210 0.9999910
# HFYeastFluid-HFYeastEpimural        -0.540763014 -0.64137939 -0.440146635 0.0000000
# HGControlEpimural-HFYeastEpimural   -0.068380335 -0.12851012 -0.008250548 0.0184833
# HGControlFluid-HFYeastEpimural      -0.028129829 -0.08191155  0.025651888 0.6593684
# HGYeastEpimural-HFYeastEpimural     -0.099419971 -0.17547880 -0.023361139 0.0050729
# HGYeastFluid-HFYeastEpimural        -0.051835784 -0.11196557  0.008294002 0.1258519
# HGControlEpimural-HFYeastFluid       0.472382679  0.36823483  0.576530525 0.0000000
# HGControlFluid-HFYeastFluid          0.512633185  0.41201681  0.613249563 0.0000000
# HGYeastEpimural-HFYeastFluid         0.441343042  0.32725479  0.555431291 0.0000000
# HGYeastFluid-HFYeastFluid            0.488927229  0.38477938  0.593075075 0.0000000
# HGControlFluid-HGControlEpimural     0.040250506 -0.01987928  0.100380293 0.3700881
# HGYeastEpimural-HGControlEpimural   -0.031039637 -0.11171221  0.049632938 0.8950668
# HGYeastFluid-HGControlEpimural       0.016544550 -0.04932433  0.082413432 0.9885966
# HGYeastEpimural-HGControlFluid      -0.071290143 -0.14734898  0.004768690 0.0769739
# HGYeastFluid-HGControlFluid         -0.023705956 -0.08383574  0.036423831 0.8830533
# HGYeastFluid-HGYeastEpimural         0.047584187 -0.03308839  0.128256762 0.5223260

## 3D scatterplot --------

# Use the vegan3d package and any distance data to graph it in 3D
# https://rdrr.io/cran/vegan3d/man/ordiplot3d.html
install.packages('vegan3d')
library(vegan3d)
install.packages('scatterplot3d')
library(scatterplot3d)
library(vegan)


# basic version of the plot, using the object made in the Betadispersion section above
ordiplot3d(EX_uJ)


# Beta diversity components ---------------------------------

# if you need to drop samples before creating the plots below, use this:

# optional, use this if you want to convert N/A in your metadata to something else
sample_data(EX_ps_clean.rar)$Factor[is.na(sample_data(EX_ps_clean.rar)$Factor)] <- "no_data_available" 

## CCA in phyloseq -------------------------------------------

# create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data

# reload as needed
library(phyloseq)
library(vegan)
library(ggplot2)

bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray") 

bray_not_na2 <- phyloseq::distance(physeq = EX_ps_clean.rar2, method = "bray")

cca_ord <- ordinate(
  physeq = EX_ps_clean.rar,
  method = "CCA",
  distance = bray_not_na,
  formula = ~ All) 

cca_ord_diet <- ordinate(
  physeq = EX_ps_clean.rar, 
  distance = bray_not_na,
  formula = ~ Diet
)

cca_ord2 <- ordinate(
  method = "CCA",
  distance = bray_not_na2,
  formula = ~ All 
)

cca_ord

cca_ord2

cca_ord_diet


# Call: cca(formula = OTU ~ All, data = data)

# -- Model Summary --

#    Inertia Proportion Rank
# Total          7.7042     1.0000     
# Constrained    2.4588     0.3191    7
# Unconstrained  5.2454     0.6809   22

# Inertia is scaled Chi-square

# -- Eigenvalues --

# Eigenvalues for constrained axes:
# CCA1   CCA2   CCA3   CCA4   CCA5   CCA6   CCA7 
#  0.7094 0.4715 0.3409 0.2783 0.2405 0.2368 0.1813 

# Eigenvalues for unconstrained axes:
#   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
# 0.3794 0.3662 0.3270 0.3128 0.2934 0.2825 0.2717 0.2598 
# (Showing 8 of 22 unconstrained eigenvalues)


# Call: cca(formula = OTU ~ All, data = data)

# -- Model Summary --
# Inertia Proportion Rank
# Total          7.7042     1.0000     
# Constrained    2.4588     0.3191    7
# Unconstrained  5.2454     0.6809   22

# Inertia is scaled Chi-square

# -- Eigenvalues --

#  Eigenvalues for constrained axes:
#  CCA1   CCA2   CCA3   CCA4   CCA5   CCA6   CCA7 
# 0.7094 0.4715 0.3409 0.2783 0.2405 0.2368 0.1813 

# Eigenvalues for unconstrained axes:
#   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
# 0.3794 0.3662 0.3270 0.3128 0.2934 0.2825 0.2717 0.2598 
# (Showing 8 of 22 unconstrained eigenvalues)

# anova of whole model
anova(cca_ord, permu=1000)

anova(cca_ord2, permu=1000)

cca_ord_diet

# Permutation test for cca under reduced model
# Permutation: free
# Number of permutations: 999

# Model: cca(formula = OTU ~ All, data = data)
# Df ChiSquare      F Pr(>F)    
# Model     7    2.4588 1.4732  0.001 ***
# Residual 22    5.2454                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permutation test for cca under reduced model 2
# Permutation: free
# Number of permutations: 999

# Model: cca(formula = OTU ~ All, data = data)
# Df ChiSquare      F Pr(>F)    
# Model     7    2.4588 1.4732  0.001 ***
# Residual 22    5.2454                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# anova of the factors (terms) you specified
anova(cca_ord, by="terms", permu=1000)

anova(cca_ord2, by="terms", permu=1000)

anova(cca_ord_diet, by="terms", permu=1000)

# Permutation test for cca under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: cca(formula = OTU ~ Diet, data = data)
# Df ChiSquare      F Pr(>F)    
# Diet      1    0.6861 2.7372  0.001 ***
#  Residual 28    7.0181                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Permutation test for cca under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: cca(formula = OTU ~ All, data = data)
#          Df ChiSquare      F Pr(>F)    
# All       7    2.4588 1.4732  0.001 ***
# Residual 22    5.2454                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Permutation test for cca under reduced model 2
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: cca(formula = OTU ~ All, data = data)
#          Df ChiSquare      F Pr(>F)    
# All       7    2.4588 1.4732  0.001 ***
# Residual 22    5.2454                  
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# cca plot with FACTORA in the model
cca_plot <- plot_ordination(
  ordination = cca_ord, 
  color = "All", 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Treatment)) + 
  geom_point(aes(colour = All), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) +
  labs(color = "Sample Group", shape = "Treatment") + 
  ggtitle("CCA Plot Unite Alignment Data") 

cca_plot2 <- plot_ordination(
  physeq = EX_ps_clean.rar2, 
  ordination = cca_ord2, 
  color = "All", 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Treatment)) + 
  geom_point(aes(colour = All), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) +
  labs(color = "Sample Group", shape = "Treatment") + 
  ggtitle("CCA Plot RPD_LSU Alignment Data")


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cca_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CCA1, 
                 yend = CCA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CCA1, 
                 y = 1.3 * CCA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cca_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 

## RDA in phyloseq -------------------------------------------
# Redundancy analysis for linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data
bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

bray_not_na2 <- phyloseq::distance(physeq = EX_ps_clean.rar2, method = "bray")

rda_ord <- ordinate(
  physeq = EX_ps_clean.rar, 
  method = "RDA",
  distance = bray_not_na,
  formula = ~ All)

rda_ord2 <- ordinate(
  physeq = EX_ps_clean.rar2,
  method = "RDA",
  distance = bray_not_na2,
  formula = ~ All)

rda_ord
rda_ord2

# Call: rda(formula = OTU ~ All, data = data)

#  -- Model Summary --

#  Inertia Proportion Rank
# Total         5.124e+06  1.000e+00     
# Constrained   2.597e+06  5.069e-01    7
# Unconstrained 2.527e+06  4.931e-01   22

# Inertia is variance

# -- Eigenvalues --

# Eigenvalues for constrained axes:
#   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7 
# 1495333  671373  167614  123330   82895   41653   15281 

# Eigenvalues for unconstrained axes:
#  PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
# 550384 485034 340521 198452 162006 141779 120396  89623 
# (Showing 8 of 22 unconstrained eigenvalues)

# Call: rda(formula = OTU ~ All, data = data) 2

# -- Model Summary --

#  Inertia Proportion Rank
#  Total         5.124e+06  1.000e+00     
# Constrained   2.597e+06  5.069e-01    7
# Unconstrained 2.527e+06  4.931e-01   22

# Inertia is variance

# -- Eigenvalues --

#  Eigenvalues for constrained axes:
#   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7 
# 1495333  671373  167614  123330   82895   41653   15281 

# Eigenvalues for unconstrained axes:
#  PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
# 550384 485034 340521 198452 162006 141779 120396  89623 
# (Showing 8 of 22 unconstrained eigenvalues)

# anova of whole model
anova(rda_ord, permu=1000)
anova(rda_ord2 permu=1000)

# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999

# Model: rda(formula = OTU ~ All, data = data) 2
#          Df Variance      F Pr(>F)    
# Model     7  2597479 3.2307  0.001 ***
# Residual 22  2526871                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permutation: free
# Number of permutations: 999

# Model: rda(formula = OTU ~ All, data = data)
#          Df Variance      F Pr(>F)    
# Model     7  2597479 3.2307  0.001 ***
# Residual 22  2526871                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# anova of the factors (terms) you specified
anova(rda_ord, by="terms", permu=1000)

anova(rda_ord2, by="terms", permu=1000)

# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: rda(formula = OTU ~ All, data = data)
#        Df Variance      F Pr(>F)    
# All       7  2597479 3.2307  0.001 ***
# Residual 22  2526871                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permutation test for rda under reduced model 2
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: rda(formula = OTU ~ All, data = data)
#        Df Variance      F Pr(>F)    
# All       7  2597479 3.2307  0.001 ***
# Residual 22  2526871                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# rda plot with FACTORA in the model
rda_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, 
  ordination = rda_ord, 
  color = "All", 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Treatment)) + 
  geom_point(aes(colour = All), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) +
  labs(color = "Sample Group", shape = "Treatment") + 
  ggtitle("RDA Plot Unite Alignment Data")

rda_plot2 <- plot_ordination(
  physeq = EX_ps_clean.rar2,
  ordination = rda_ord2, 
  color = "All", 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Treatment)) + 
  geom_point(aes(colour = All), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) +
  labs(color = "Sample Group", shape = "Treatment") + 
  ggtitle("RDA Plot RPD_LSU Alignment Data")

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(rda_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1, 
                 yend = RDA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * RDA1, 
                 #y = 1.3 * RDA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
rda_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 


## dbRDA in phyloseq -------------------------------------------
# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data

bray_not_na2 <- phyloseq::distance(physeq = EX_ps_clean.rar2, method = "bray")

cap_ord <- ordinate(
  physeq = EX_ps_clean.rar, 
  method = "CAP", #borrows capscale from the vegan package
  distance = bray_not_na,
  formula = ~ All)

cap_ord2 <- ordinate(
  physeq = EX_ps_clean.rar2, 
  method = "CAP", #borrows capscale from the vegan package
  distance = bray_not_na2,
  formula = ~ All)

cap_ord

cap_ord2

# Call: capscale(formula = distance ~ All, data = data)

# -- Model Summary --

#  Inertia Proportion Rank
# Total          9.8964     1.0000     
# Constrained    4.8684     0.4919    7
# Unconstrained  5.0280     0.5081   22

# Inertia is squared Bray distance

# -- Eigenvalues --

# Eigenvalues for constrained axes:
# CAP1   CAP2   CAP3   CAP4   CAP5   CAP6   CAP7 
# 2.5146 0.9560 0.5351 0.3228 0.2257 0.1817 0.1324 

# Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#  0.5714 0.4561 0.3876 0.3375 0.3071 0.2838 0.2561 0.2375 
# (Showing 8 of 22 unconstrained eigenvalues)

# Call: capscale(formula = distance ~ All, data = data) 2

# -- Model Summary --

#   Inertia Proportion Rank
# Total          9.8964     1.0000     
# Constrained    4.8684     0.4919    7
# Unconstrained  5.0280     0.5081   22

# Inertia is squared Bray distance

# -- Eigenvalues --

#  Eigenvalues for constrained axes:
#   CAP1   CAP2   CAP3   CAP4   CAP5   CAP6   CAP7 
# 2.5146 0.9560 0.5351 0.3228 0.2257 0.1817 0.1324 

# Eigenvalues for unconstrained axes:
#   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
# 0.5714 0.4561 0.3876 0.3375 0.3071 0.2838 0.2561 0.2375 
# (Showing 8 of 22 unconstrained eigenvalues)

# anova of whole model
anova(cap_ord, permu=1000)

anova(cap_ord2, permu=1000)


# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999

# Model: capscale(formula = distance ~ All, data = data)
#          Df SumOfSqs      F Pr(>F)    
# Model     7   4.8684 3.0431  0.001 ***
# Residual 22   5.0280                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permutation test for capscale under reduced model 2
# Permutation: free
# Number of permutations: 999

# Model: capscale(formula = distance ~ All, data = data)
#         Df SumOfSqs      F Pr(>F)    
# Model     7   4.8684 3.0431  0.001 ***
# Residual 22   5.0280                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# anova of the factors (terms) you specified
anova(cap_ord, by="terms", permu=1000)

anova(cap_ord2, by="terms", permu=1000)


# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: capscale(formula = distance ~ All, data = data)
#        Df SumOfSqs      F Pr(>F)    
# All       7   4.8684 3.0431  0.001 ***
# Residual 22   5.0280                  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permutation test for capscale under reduced model 2
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: capscale(formula = distance ~ All, data = data)
#          Df SumOfSqs      F Pr(>F)    
# All       7   4.8684 3.0431  0.001 ***
#  Residual 22   5.0280                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# CAP plot with FACTORA in the model

cap_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, 
  ordination = cap_ord, 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Treatment)) + 
  geom_point(aes(colour = All), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) +
  labs(color = "Sample Group", shape = "Treatment") + 
  ggtitle("CAP Plot Unite Alignment Data")

cap_plot2 <- plot_ordination(
  physeq = EX_ps_clean.rar2, 
  ordination = cap_ord2, 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Treatment)) + 
  geom_point(aes(colour = All), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) +
  labs(color = "Sample Group", shape = "Treatment") + 
  ggtitle("CAP Plot RPD_LSU Alignment Data")

# add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 

## variance partitioning --------------
# http://www.hiercourse.com/docs/variationPartitioning.html 

library(vegan)

# create a dataframe from your SV table
EX.df <- as.data.frame(otu_table(EX_ps_clean.rar))

EX.df2 <- as.data.frame(otu_table(EX_ps_clean.rar2))

# extract your sample data from the phyloseq object
env.df <- as.data.frame(sample_data(EX_ps_clean.rar))

env.df2 <- as.data.frame(sample_data(EX_ps_clean.rar2))


# calculate Principal Coordinates Of Neighborhood Matrix, diversity distance data transformed into rectangular format
EX.pcnm <- pcnm(dist(bray_not_na))

EX.pcnm2 <- pcnm(dist(bray_not_na2))


# environmental variables as predictors of community similarity
cap.env <- capscale(EX.df ~ ., data=env.df, distance='bray')

# calculate CCA 
cap.pcnm <- capscale(EX.df ~ ., data=as.data.frame(scores(EX.pcnm)), distance='bray')


# make a model using SV ordination and sample data
mod0.env <- capscale(EX.df ~ 1, data=env.df, distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# make a model using SV ordination and the component scores
mod0.pcnm <- capscale(EX.df ~ 1, data=as.data.frame(scores(EX.pcnm)), distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# select variables in each predictor table
step.env <- ordistep(mod0.env, scope=formula(cap.env))

# check variance inflation factors, higher number = data are redundant to another factor, 1= data are unique. 
vif.cca(step.env)


step.pcnm <- ordistep(mod0.pcnm, scope=formula(cap.pcnm))

step.env$anova
# paste output here

step.pcnm$anova
# paste output here

EX.pcnm.sub <- scores(EX.pcnm, 
                      choices=c(1,3:16)) 
# partition variation among four predictor tables:
EX.var <- varpart(EX.df, 
                  ~ FACTORA, 
                  ~ FactorB
                  EX.pcnm.sub, data=env.df)

#plot 
par(mfrow=c(1, 2))
showvarparts(4)
plot(EX.var)

EX.var

anova(rda(EX.df  ~ FACTORA + Condition(EX.pcnm.sub), data=env.df)) # add + Condition(env.df$SAMPLE_ID) for repeated measures

## Mantel test -------

bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

# make a data frame that can be used in vegan from the sample_data in the phyloseq object
meta.df <- data.frame(sample_data(EX_ps_clean.rar))

# check the columns in your meta dataframe, need numeric data for mantel test
str(meta.df)

# make numeric any columns that are a number but as listed as "chr", for example 
meta.df$DNA_extraction_batch <- as.numeric(meta.df$DNA_extraction_batch)

# replace any binary factors with numbers, for example
meta.df["Diet"][meta.df["Diet"] == "LooseAlfalfa"] <- "1"
meta.df["Diet"][meta.df["Diet"] == "PelletedAlfalfa"] <- "2"
meta.df$Diet <- as.numeric(meta.df$Diet)

#remove any columns that you can't make numeric
meta.df <- subset(meta.df, select = -c(Sample_type, Sheep_ID, Treatment))

# Make a distance matrix for your meta data
bray.dist.meta <- vegdist(meta.df, method = "bray") 

# Mantel test to see if SV and meta distance matrixes correlate
bray.mantel <- mantel(bray_not_na, bray.dist.meta, permutations = 999)

# look at the output
bray.mantel

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = bray_not_na, ydis = bray.dist.meta, permutations = 999) 
# 
# Mantel statistic r: 0.06023 
#       Significance: 0.224 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.108 0.143 0.166 0.189 
# Permutation: free
# Number of permutations: 999


# Niche/neutral models-----
# online version: https://web.rniapps.net/webglv/

# Source tracking -----

# download the SourceTracker_modified.R

library(ape)
library(plyr)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(tidyr)
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(gridExtra)
library(phyloseq)

# need to grab the OTU_table (SV table) from your phyloseq object
ex_seqtab <- as.data.frame(otu_table(EX_ps_clean))

# take the taxa names from the tax_table and paste onto the colnames for the seqtab
# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste0(tax.df$Genus, "_", tax.df$Species)

# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(ex_seqtab) <- tax.df$Genus.species

# go from dataframe format to sourcetracker preferred format
ex_seqtab_st <- t(as.matrix(ex_seqtab))

# grab the meta
ex_meta <- as.data.frame(sample_data(EX_ps_clean))

# create a column that indicates 'source' and 'sink' if you don't already have it in there
ex_meta$SourceSink <- NA

# Designate source and sink samples based on Location
ex_meta$SourceSink <- ifelse(ex_meta$Location == "Fluid", "source", "sink")

# Subset sequence table for source and sink
train_data <- ex_seqtab_st[, train.ix]  # Source samples
test_data <- ex_seqtab_st[, test.ix]  

train.ix <- which(ex_meta$SourceSink=='source')
test.ix <- which(ex_meta$SourceSink=='sink')

# Set environmental factors (e.g., location) for the analysis
envs <- ex_meta$Location
desc <- ex_meta$Location


# Change this line to reflect wherever you downloaded the modified sourcetracker code to
source('/Users/toripalecek/Documents/AVS554/AVS554_Project/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow! like very slow!) Use 0.001 if pressed for time.
# tune.results <- tune.st(ex_seqtab_st[train.ix,], envs[train.ix])

alpha1 <- 0.05 #tune.results$best.alpha1
alpha2 <- 0.05 #tune.results$best.alpha2

# Train sourcetracker using the source samples.
st <- sourcetracker(ex_seqtab_st[train.ix,], envs[train.ix])

#Estimate source proportions in sink samples.
results <- predict(st, ex_seqtab_st[test.ix,], alpha1=alpha1, alpha2=alpha2)


# This part validates the 'source' models we made using the original training data.
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)


# plot results with legend
plot(results, type='pie', include.legend=TRUE)
plot(results.train, type='pie', include.legend=TRUE)

# Explicitly add 'Fluid' as a recognized category for plotting
results.train$Source <- factor(results.train$Source, 
                               levels = c("epimural", "fluid"),
                               labels = c("Epimural", "Fluid"))

ggplot(results.train, aes(x = Sample, y = Proportion, fill = Source)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Epimural" = "blue", "Fluid" = "green")) +
  labs(
    x = "Samples",
    y = "Proportion",
    fill = "Source"
  ) +
  theme_minimal()



# Convert proportions to a data frame
proportions <- as.data.frame(results.train$proportions)

# Rename "Unknown" to "Epimural"
colnames(proportions) <- gsub("Unknown", "Epimural", colnames(proportions))

# Add sample names
proportions$Sample <- rownames(proportions)

# Reshape to long format for plotting
library(tidyr)
results.train.df <- pivot_longer(
  proportions,
  cols = -Sample,  # All columns except Sample
  names_to = "Source",
  values_to = "Proportion"
)

# Plot
library(ggplot2)
ggplot(results.train.df, aes(x = Sample, y = Proportion, fill = Source)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Epimural" = "blue", "Fluid" = "green")) +
  labs(
    x = "Most Significant SVs (p < 0.05)",
    y = "Proportion",
    fill = "Source"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Source Tracking of Fungi with Epimural and Rumen Fluid as Source and Sink Locations")


# Convert proportions to a data frame
proportions <- as.data.frame(results.train$proportions)

# Rename "Unknown" to "Epimural"
colnames(proportions) <- gsub("Unknown", "Epimural", colnames(proportions))

# Add sample names as a column
proportions$Sample <- rownames(proportions)

# Reshape the data into long format for easier viewing and table creation
library(tidyr)
results.train.df <- pivot_longer(
  proportions,
  cols = -Sample,  # All columns except Sample
  names_to = "Source",  # Create a new column "Source" from the column names
  values_to = "Proportion"  # Store the values into a column named "Proportion"
)

# Create the table using the reshaped data
library(knitr)
kable(results.train.df, caption = "Source Tracking Proportions for Epimural and Rumen Fluid Locations")

# Save the reshaped data to a CSV file
write.csv(results.train.df, "source_tracking_results.csv", row.names = FALSE)



# Load necessary package
library(tidyr)

# Reshape the data from long to wide format
proportions_wide <- results.train.df %>%
  pivot_wider(
    names_from = "Source",  # The column for the new column names
    values_from = "Proportion"  # The column with the values to fill the new columns
  )

# View the reshaped data
head(proportions_wide)

# Save the reshaped data to a CSV file
write.csv(proportions_wide, "source_tracking_proportions_wide.csv", row.names = FALSE)


downstream <- data.frame(results$proportions)
downstream$id <- row.names(downstream)
meltdown <- melt(downstream, id.vars=c("id"))

ex_meta

colnames(meltdown) <- c("Bacterial_SV", "Diet", "value")
meltdown$Diet <- ex_meta[as.vector(meltdown$id), "Diet"]
#meltdown$Week <- ex_meta[as.vector(meltdown$id), "Week"]

meltdown$Month <- as.factor(sapply(meltdown$Date, function(x) strsplit(x, "/")[[1]][1]))
levels(meltdown$Month) <- c("June", "July", "August")
meltdown$Day<- as.factor(sapply(meltdown$Date, function(x) strsplit(x, "/")[[1]][2]))

meltdown$Rep <- as.factor(sapply(meltdown$id, function(x) strsplit(x, "_")[[1]][3]))
levels(meltdown$Rep) <- c("A", "AA", "B")

meltdown$Category <- as.factor(meltdown$Category)
levels(meltdown$Category) <- c("B", "E")

meltdown$SID <- as.factor(sapply(meltdown$id, function(x) strsplit(x, "_")[[1]][5]))

meltdown$TypeRep <- paste(meltdown$Category, meltdown$Rep, meltdown$Day, sep=".")
theme_set(theme_classic(base_size=10, base_family="Avenir"))

ggplot(meltdown, aes(x=Week, y=value, fill=Diet)) + geom_bar(stat="identity") + facet_wrap(~Downstream.Site+Month, scales="free_x") + xlab("")


# Dendogram (optional) ---------------------------------------
## Note: I typically use the UPGMA algorothm for Dendograms when I make them as stand-alone.  Some other graph functions will intergrate dendograms themselves.

#calculate dissimilarity based on your sequence data in base R
dist.object <- dist(as.data.frame(otu_table(EX_ps_clean.rar)), method = "euclidean", diag = FALSE, upper = FALSE, p = 2) # method = "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"

# cluster the dissimiliary info with base R: https://www.rdocumentation.org/packages/stats/versions/3.2.1/topics/hclust
EX_ps.clust <- hclust(dist.object, method = "average", members = NULL) #method = average for UPGMA

# visualization in phyloseq: https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html#hierarchical-clustering
plot(EX_ps.clust)



