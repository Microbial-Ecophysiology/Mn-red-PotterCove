# This script was modified to use for the analysis of bacterial amplicon data
# See: https://benjjneb.github.io/dada2/tutorial.html
# It is not optimized for:
#   amplicons which are sequenced into the rev primer at the 5' end of the fwd read, and vice versa (e.g. ITS)
#   faster computing times at slightly lower output quality (e.g. big data)

# Christiane Hassenruecks own optimization was implemented, snakemake code can be found here
#https://git.io-warnemuende.de/bio_inf/workflow_templates/src/branch/master/Amplicon_dada2_MiSeq
 # differences at merging step to not loose organisms (e.g. large SOX bacteria) with large 16S gene which cannot be merged
 # need reference for mapping for this!


## run on server

# load module on server:  

# open in separate screen, especially for large jobs:
#screen -S R
# activate environment needed for further mapping below before opening R
# module load bbmap/38.86 # important to load this before R or wrong R version will open
# conda activate metawrap-env
# module load R/4.3.1 # version number might be different, use autocomplete
# open R while in the correct working directory: 


# open R by typing R

## load packages ####
require(dada2)
require(ShortRead)
require(gridExtra)
require(tidyverse)

packageVersion("dada2")
# '1.28.0'

# save and load workspace
## path to working directory
wd <- "WORKINGDIRECTORY/"
setwd(wd)

## workspace file
tf <- "TEMPORARYFILE.Rdata"
# save.image(tf)
# load(tf) # to load workspace again

## path to scripts to load
sloc <- "SCRIPTLOCATION/Scripts_to_source/"

## path to database for mapping
dblm <- "DATABASE_LOCATION_MAP/"

## path to database for classification
dblc <- "DATABASE_LOCATION_CLAS/"

## project name to name files
prj <- "PROJECT"


# specify path to input fastq files
## this file has the names of the flowcells with the lane afterwards in each line, no header
 # created automatically in script dada2_seqprep.bash, can also be easily made by hand
seq_lanes <- read.table(paste0(wd, "flowcell_IDs.txt"),
                        col.names = "flowcell", sep = "\t")
rownames(seq_lanes) <- seq_lanes$flowcell

## here replace for Arc if analysing archaea libraries
TARGET <- "Bac"

extract_sample_paths <- function(flowcell_ID) {
path <- paste0(wd, "seq_by_flowcell/", flowcell_ID, "/Clipped_", TARGET)
fnFR_R1 <- sort(list.files(path, pattern="clip_fr_R1.fastq", full.names = TRUE))
fnFR_R2 <- sort(list.files(path, pattern="clip_fr_R2.fastq", full.names = TRUE))
fnRF_R1 <- sort(list.files(path, pattern="clip_rf_R1.fastq", full.names = TRUE))
fnRF_R2 <- sort(list.files(path, pattern="clip_rf_R2.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFR_R1), "_"), `[`, 1)
flowcell_ID <- setNames(list(fnFR_R1, fnFR_R2, fnRF_R1, fnRF_R2, sample.names), 
                        c("fnFR_R1", "fnFR_R2", "fnRF_R1", "fnRF_R2", "sample.names"))
return(flowcell_ID)
}

### create list with each seq_lane as one element
sample.paths <- apply(seq_lanes, 1, extract_sample_paths)

## quality check and trimming ####
# mkdir QualityProfiles
source(paste0(sloc, "dada2_quality_check.R"))
mapply(function(x, i) {
  quality_check(
    c(x[['fnFR_R1']], x[['fnRF_R1']]),              # need '' inside [[]] to access list elements by name
    c(x[['fnFR_R2']], x[['fnRF_R2']]),
    file_base = paste0("QualityProfiles/", i,"_QualityProfile_separate") # need to create folder before!
  )},
  x = sample.paths, i = names(sample.paths),        # need mapply instead of lapply to access the names of the list
  SIMPLIFY = F
)


### check pdfs:
# Considerations for trimming:
# expected max length: 252bp (?)
# min overlap: 30bp  --> in total R1 + R2 seq should have 285-290 bp total length
# reads should be truncated so that rev primer is not included at end of fwd reads
# It is recommended to trim to just enough for the required length for sufficient overlap
# Caution: don't remove too much


## Run parameter optimization
# Define ranges for maxEE (normally same for all)
range_maxEE <- matrix(
  c(1, 1,
    2, 1,
    2, 2,
    3, 2,
    3, 3),
  nrow = 5,
  ncol = 2,
  byrow = T
)

# add to list according parameters for testing truncation length and maxEE
# check quality profile pdfs for this:
### seq lane 1/1 HKMLNDRXY_L1
  ## forward (R1)
  # also drops at beginning of read, 
  # substantial bad peaks, cut as short as possible
  ## reverse (R2)
  # fr really good
   #rf bad peaks 140 bp, 175-180 bp
  # potential cut-off: 140, 175, 200 bp
sample.paths[['HKMLNDRXY_L1']][['range_truncLen']] <- matrix(
    c(145, 145, # first number for R1 read (= forward), second number for R2 read (=reverse)
      115, 175,
      120, 170,
      100, 190,
      90, 200),
  nrow = 5, ncol = 2, byrow = T)
sample.paths[['HKMLNDRXY_L1']][['range_maxEE']] <- range_maxEE


# Prepare directories to place filtered files in Filtered_TARGET/ subdirectory in folder seq_by_flowcell individual flowcell folders
sample.paths <- mapply(function(x, i) {
  x[['filtFR_R1']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_FR_R1_filt.fastq"))
  x[['filtFR_R2']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_FR_R2_filt.fastq"))
  x[['filtRF_R1']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_RF_R1_filt.fastq"))
  x[['filtRF_R2']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_RF_R2_filt.fastq"))
  names(x[['filtFR_R1']]) <- x[['sample.names']]
  names(x[['filtFR_R2']]) <- x[['sample.names']]
  names(x[['filtRF_R1']]) <- x[['sample.names']]
  names(x[['filtRF_R2']]) <- x[['sample.names']]
  x
},
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# Run parameter optimization
# This is quite time consuming and should only be attempted on a server with as many cores as you have samples (or at least 20)
source(paste0(sloc, "dada2_screen_settings.R"))
threads <- 60
sample.paths <- mapply(function(x, i) {
  x[['screen_filt_settings_fr']] <- 
    screen_settings(sample.names = x[['sample.names']], fnFs = x[['fnFR_R1']], fnRs = x[['fnFR_R2']], 
                    range_maxEE = x[['range_maxEE']], range_truncLen = x[['range_truncLen']], cpus = threads)
  x[['screen_filt_settings_rf']] <- 
    screen_settings(sample.names = x[['sample.names']], fnFs = x[['fnRF_R1']], fnRs = x[['fnRF_R2']], 
                    range_maxEE = x[['range_maxEE']], range_truncLen = x[['range_truncLen']], cpus = threads)
  x
},
  x = sample.paths, i = names(sample.paths),
  SIMPLIFY = F
)

save.image(tf)



## make plots
# mkdir QualityProfiles
mapply(function(x, i) {
pdf(paste0("QualityProfiles/Parameter_screening_", i, ".pdf"), width = 7, height = 7)
plot(
  x[['screen_filt_settings_fr']][, "prop.total"],
  x[['screen_filt_settings_fr']][, "q90"] - x[['screen_filt_settings_fr']][, "q10"],
  col = rep(rainbow(nrow(x[['range_maxEE']])), nrow(x[['range_truncLen']])),
  pch = 16
)
text(
  x[['screen_filt_settings_fr']][, "prop.total"],
  x[['screen_filt_settings_fr']][, "q90"] - x[['screen_filt_settings_fr']][, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
plot(
  x[['screen_filt_settings_rf']][, "prop.total"],
  x[['screen_filt_settings_rf']][, "q90"] - x[['screen_filt_settings_rf']][, "q10"],
  col = rep(rainbow(nrow(x[['range_maxEE']])), nrow(x[['range_truncLen']])),
  pch = 16
)
text(
  x[['screen_filt_settings_rf']][, "prop.total"],
  x[['screen_filt_settings_rf']][, "q90"] - x[['screen_filt_settings_rf']][, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
dev.off()
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


# This is just a gut feeling, but I would optimize for the following criteria:
#   small difference between 10 and 90 percentile of retained reads
#   high total proportion of retained reads
#   most stringent maxEE that does not result in severe loss of reads

# look at numbers in table:
sample.paths[['flowcell']][['screen_filt_settings_fr']]
sample.paths[['flowcell']][['screen_filt_settings_rf']]

## Run trimming with optimal parameters

# define optimal trunc length, these parameters will be needed again further down

# seq lane HKMLNDRXY_L1 ##
 # truncLen = c(120,170), maxEE = c(2,2);
 # index 13, fr: prop.total 0.9703; q10 0.9665; q90 0.9763; rf: prop.total 0.9576; q10 0.9534; q90 0.9646;
sample.paths[['HKMLNDRXY_L1']][['truncLen_R1']] <- 120
sample.paths[['HKMLNDRXY_L1']][['truncLen_R2']] <- 170
sample.paths[['HKMLNDRXY_L1']][['error_R1']] <- 2
sample.paths[['HKMLNDRXY_L1']][['error_R2']] <- 2


# run trimming
threads <- 60

sample.paths <- mapply(function(x, i) {
  x[['filt_FR.out']] <- filterAndTrim(
  fwd = x[['fnFR_R1']], 
  filt = x[['filtFR_R1']], 
  rev = x[['fnFR_R2']], 
  filt.rev = x[['filtFR_R2']],
  truncLen = c(x[['truncLen_R1']], x[['truncLen_R2']]),
  maxN = 0,
  minQ = 2,
  maxEE = c(x[['error_R1']], x[['error_R2']]), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = threads
  )

  x[['filt_RF.out']] <- filterAndTrim(
  fwd = x[['fnRF_R1']], 
  filt = x[['filtRF_R1']], 
  rev = x[['fnRF_R2']], 
  filt.rev = x[['filtRF_R2']],
  truncLen = c(x[['truncLen_R1']], x[['truncLen_R2']]),
  maxN = 0,
  minQ = 2,
  maxEE = c(x[['error_R1']], x[['error_R2']]), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = threads
  )
  x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

## Repeat quality check after trimming
mapply(function(x, i) {
quality_check(
  c(x[['filtFR_R1']], x[['filtRF_R1']]),
  c(x[['filtFR_R2']], x[['filtRF_R2']]),
  file_base = paste0("QualityProfiles/", i, "_QualityProfileFiltered_", x[['truncLen_R1']], "_", x[['truncLen_R2']])
)
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

save.image(tf)



### Learn error rates ####
# It is generally not necessary to increase the number of nbases used for the error estimation
# It is possible that with 10 rounds (MAX_CONSIST), the algorithm for learning the errors won't converge
# Increasing MAX_CONSIST will lead to longer run times, and may only marginally improve error estimation
# I would not recommend setting MAX_CONSIST higher than 15

## correct error estimates because of binned quality scores?
# https://github.com/benjjneb/dada2/issues/791
# https://github.com/benjjneb/dada2/issues/938
# The dada2 team is working on making some good recommendations for binned quality scores
# For now, there are 2 options:
#   1) run error learning with modified loess function (maybe more elegant)
#   Hack the loessErrfun() of dada2 package (used in both learnErrors and dada): 
#   mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)
#   2) coerce any value lower than the Q40 probability to be the Q40 value in the learnErrors() output
#   We will do this here to avoid re-running the error learning


## Use alternative Loess function, this one always gave best results from binned quality scores
source(paste0(sloc, "loessErrfun2.R"))
threads <- 60

sample.paths <- mapply(function(x, i) {
  
sink(paste0("Logfiles/", i, "_", TARGET, "_log_errorEstimation.log"))
  
x[['errFR_R1']] <- learnErrors(x[['filtFR_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)
x[['errFR_R2']] <- learnErrors(x[['filtFR_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)
x[['errRF_R1']] <- learnErrors(x[['filtRF_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)
x[['errRF_R2']] <- learnErrors(x[['filtRF_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)

sink()
x

}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# all samples reached convergense after 7-8 rounds


# plot error profiles
## do plots manually
x <- "HKMLNDRXY_L1"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()



### Dereplicate and denoise samples ####
# This step takes a while...
# For large data set (e.g. full HiSeq lane), I strongly recommend pool = "pseudo"
# I would not use pool = FALSE as this will strongly impact (i.e. lower) your alpha diversity,
# which seems to be rather an artifact of the change in parameters than any true signal
threads <- 60
sample.paths <- mapply(function(x, i) {
  
sink(paste0("Logfiles/", i, "_", TARGET, "_derep_denoise.log"))

x[['dadaFR_R1']] <- dada(x[['filtFR_R1']], err = x[['errFR_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)
x[['dadaFR_R2']] <- dada(x[['filtFR_R2']], err = x[['errFR_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)
x[['dadaRF_R1']] <- dada(x[['filtRF_R1']], err = x[['errRF_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)
x[['dadaRF_R2']] <- dada(x[['filtRF_R2']], err = x[['errRF_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)

sink()

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# it is a good idea to save your workspace here
save.image(tf)



## Merge reads ####
# mismatches are 0 by default
# verbose = TRUE prints merging pair numbers to screen
source(paste0(sloc, "extract_unmerged.R"))

sample.paths <- mapply(function(x, i) {
x[['mergers_FR0']] <- mergePairs(
  x[['dadaFR_R1']],
  x[['filtFR_R1']], 
  x[['dadaFR_R2']], 
  x[['filtFR_R2']], 
  minOverlap = 10,
  verbose = TRUE,
  returnRejects = TRUE
)
x[['mergers_RF0']] <- mergePairs(
  x[['dadaRF_R1']],
  x[['filtRF_R1']], 
  x[['dadaRF_R2']], 
  x[['filtRF_R2']], 
  minOverlap = 10,
  verbose = TRUE,
  returnRejects = TRUE
)


### rescue unmerged reads ####

#### extract them and save in fasta file
x[['unmerged_FR']] <- extract_unmerged(x[['dadaFR_R1']], x[['dadaFR_R2']], x[['mergers_FR0']])
x[['unmerged_RF']] <- extract_unmerged(x[['dadaRF_R1']], x[['dadaRF_R2']], x[['mergers_RF0']])

writeFasta(x[['unmerged_FR']][[1]], file = paste0(i, "_", TARGET, "_unmerged_FR_R1.fasta"))
writeFasta(x[['unmerged_FR']][[2]], file = paste0(i, "_", TARGET, "_unmerged_FR_R2.fasta"))
writeFasta(x[['unmerged_RF']][[1]], file = paste0(i, "_", TARGET, "_unmerged_RF_R1.fasta"))
writeFasta(x[['unmerged_RF']][[2]], file = paste0(i, "_", TARGET, "_unmerged_RF_R2.fasta"))

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)



#### outside of R reads will be mapped with bbmap.sh and insert size extracted
   # can access outside with system(), conda environment has to be active when opening R!!!
threads <- 100
map_ref <- paste0(dbl, "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta")

sample.paths <- mapply(function(x, i) {
  
system(paste0(
  "bbmap.sh threads=", threads, " ref=", map_ref, " in=", i, "_", TARGET, "_unmerged_FR_R1.fasta in2=", 
  i, "_", TARGET, "_unmerged_FR_R2.fasta out=", i, "_", TARGET, "_unmerged_FR.bam"))

system(paste0(
  "bbmap.sh threads=", threads, " ref=", map_ref, " in=", i, "_", TARGET, "_unmerged_RF_R1.fasta in2=", 
  i, "_", TARGET, "_unmerged_RF_R2.fasta out=", i, "_", TARGET, "_unmerged_RF.bam"))


##### extract insert size
system(paste0(
  "samtools view -F2304 -f66 -m50 ",
  i, "_", TARGET, "_unmerged_FR.bam | cut -f1,9 > ", i, "_", TARGET, "_unmerged_is_FR.txt"))

system(paste0(
  "samtools view -F2304 -f66 -m50 ",
  i, "_", TARGET, "_unmerged_RF.bam | cut -f1,9 > ", i, "_", TARGET, "_unmerged_is_RF.txt"))


#### read insert size back into R
x[['is.fr']] <- read.table(paste0(i, "_", TARGET, "_unmerged_is_FR.txt"), h = F, sep = "\t", col.names = c("seqID", "insert"))
x[['is.rf']] <- read.table(paste0(i, "_", TARGET, "_unmerged_is_RF.txt"), h = F, sep = "\t", col.names = c("seqID", "insert"))


#### filter to insert sizes that exceed maximum length of merged sequences
x[['is_long.fr']] <- x[['is.fr']][x[['is.fr']][['insert']] > (x[['truncLen_R1']] + x[['truncLen_R2']] - 10), ] %>% 
  separate(seqID, into = c("sample_index", "row_index"), sep = "_", remove = F, convert = T)

x[['is_long.rf']] <- x[['is.rf']][x[['is.rf']][['insert']] > (x[['truncLen_R1']] + x[['truncLen_R2']] - 10), ] %>% 
  separate(seqID, into = c("sample_index", "row_index"), sep = "_", remove = F, convert = T)

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

save.image(tf)


#### retrieve and concatenate sequence
source(paste0(sloc, "retr_concat_unmerged.R"))

sample.paths <- mapply(function(x, i) {

x[['mergers_FR']] <- retr_concat_unmerged(x[['mergers_FR0']], x[['is_long.fr']], x[['unmerged_FR']])
x[['mergers_RF']] <- retr_concat_unmerged(x[['mergers_RF0']], x[['is_long.rf']], x[['unmerged_RF']])

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)



# Create sequence table ####
# Create sequence table with actual sequences as column names
sample.paths <- mapply(function(x, i) {
  
x[['seqtab_FR']] <- makeSequenceTable(x[['mergers_FR']])
x[['seqtab_RF']] <- makeSequenceTable(x[['mergers_RF']])

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


mapply(function(x, i) {
# for dimension first number samples, second number sequences (=ASVs)
dim(x[['seqtab_FR']])
dim(x[['seqtab_RF']])
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)



# As with the tryRC option of mergeSequenceTables only those sequences which are duplicated
# will be turned, manually turn sequences of RF table
# seqtab_RF_rc

# Generate reverse complement of rf
sample.paths <- mapply(function(x, i) {
  
x[['seqtab_RF_rc']] <- x[['seqtab_RF']]
colnames(x[['seqtab_RF_rc']]) <- rc(colnames(x[['seqtab_RF']]))

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# Merge sequence tables fr and reverse rf sequences, as it is only seq lane here only in one table
seqtab <- mergeSequenceTables(sample.paths[['HKMLNDRXY_L1']][['seqtab_FR']], sample.paths[['HKMLNDRXY_L1']][['seqtab_RF_rc']] ,repeats = "sum")
dim(seqtab)


# Remove chimeras ####
# This may remove quite a bit of ASVs, but only a small fraction of your total sequences
threads <- 100
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)

ncol(seqtab.nochim)/ncol(seqtab)

summary(rowSums(seqtab.nochim)/rowSums(seqtab))



# Inspect ASV length distribution
table(nchar(colnames(seqtab.nochim))) 

table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))) 



# Check unusual sequence lengths
### Check unusual sequence lengths
# subset sequence table to most common sequences of unusual length
seqtab.unus.length <- seqtab.nochim %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "seq") %>% 
  mutate(seq_length = nchar(seq)) %>% 
  filter(seq_length < 249 | seq_length > 254) %>% 
  mutate(seq_count = rowSums(.[rownames(seqtab.nochim)]),
         seq_name = paste0("sq", 1:nrow(.))) %>% 
  mutate(seq_name = paste0(seq_name, ";size=", seq_count)) %>% 
  slice_max(order_by = seq_count, n = 50)

# simple function to write fasta file from data.frame
writeFasta <- function(data, colname_seqname, colname_seq, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum, colname_seqname], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum, colname_seq]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# save fasta file
writeFasta(seqtab.unus.length, colname_seqname = "seq_name", colname_seq = "seq",
           filename = paste0("check_", prj, ".fasta"))



# Remove potential junk sequences and singletons
# dada does not generate singletons, any singletons are introduced in the merging step
# Adjust range of sequence lengths based on expected length of marker gene fragment and read distribution (here 252 bp)
# to save unmerged reads (which were merged based on alignment) estimate upper threshold of junk:
# shortest trimming - merging overlap: 
seqtab.nochim2 <- seqtab.nochim[, colSums(seqtab.nochim) > 1 & ((nchar(colnames(seqtab.nochim)) >= 249 & nchar(colnames(seqtab.nochim)) <= 254) | nchar(colnames(seqtab.nochim)) == 300)]
dim(seqtab.nochim2)
ncol(seqtab.nochim2)/ncol(seqtab)

summary(rowSums(seqtab.nochim2)/rowSums(seqtab))



# Taxonomic classification ####
# Available options:
#   DECIPHER (IdTaxa), 
#   RDP
#   Blast
#   silvangs
#   etc.

# use silvangs (according to Christiane will give best result)
# I am disabling the bootstrap filtering, but saving the bootstrap values
# so that we can manually filter by bootstrap later
threads <- 100
tax <- assignTaxonomy(
  seqtab.nochim2, 
  paste0(dblc, "silva_nr99_v138.1_train_set.fa.gz"),
  multithread = threads,
  minBoot = 0,
  outputBootstraps = T,
  tryRC = T,
  taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
)

save.image(tf)


# Determine most suitable bootstrap cut-off
apply(tax$boot, 2, summary)


### How many sequences and ASVs will be unclassified at a minboot of e.g. 80?
cutoff <- 80
tax.filt <- tax$tax
tax.filt[tax$boot < cutoff] <- NA
tax.filt <- tax.filt[!is.na(tax.filt[, 1]), ]
otu.filt <- seqtab.nochim2[, rownames(tax.filt)]
apply(tax.filt, 2, function(x) sum(is.na(x)) )



# total counts of reads which are assigned as NA on different ranks
for(i in 1:ncol(tax.filt)) {
  print(colnames(tax.filt)[i])
  print(sum(otu.filt[, is.na(tax.filt[, i])]))
}


# percent of reads which are assigned as NA on different ranks
for(i in 1:ncol(tax.filt)) {
  print(colnames(tax.filt)[i])
  print(sum(otu.filt[, is.na(tax.filt[, i])])/sum(otu.filt))
}



## decide for bootstrap value #####
otu.filt_f <- otu.filt # final one (note here): 80
otu.filt_m <- otu.filt # 90 # alternative tried bootstrap value, just for nSeq file
tax_filt_f <- tax.filt


## Remove unwanted lineages
#### skip this step here and do it instead later in metacoder, to easily count how many chloroplast and mitochondria sequences were there
# (only Bacteria are kept, Chloroplast and Mitochondria can be deleted)


# Get nSeqs summary ####
# source functions for creating nSeq files
source(paste0(sloc, "dada2_nSeq_creation_for_list.R"))

# lib_no needs to contain everything after the default file name which is library specific, so here '_libx' including the '_'
# if the file names are the original ones as only one library is analysed use lib_no = NULL instead
nSeq_all <- mapply(function(x, i) {
  i <- create_nSeq(list_input = x, nSeq_file_input = paste0("../nSeqs_", i, "_Bac.txt"))
  
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# Warning message 
#1: Expected 1 pieces. Additional pieces discarded in 99 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#2: Expected 1 pieces. Additional pieces discarded in 99 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#3: Expected 1 pieces. Additional pieces discarded in 99 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
# is ok and does not influence outcome


# feed into function
nSeq <- amend_nSeq(lib_list_input = nSeq_all)

## calculate how many reads were lost on each step
nSeq_perc <- perc_nSeq(nSeq_input = nSeq)



# Write output ####
## string that will be used in all written output, so it has not to be changed everywhere manually
output_name <- prj

## statistic
write.table(nSeq, paste0("nSeq_dada2_", output_name, ".txt"), row.names = F, quote = F, sep = "\t")
write.table(nSeq_perc, paste0("nSeq_perc_dada2_", output_name, ".txt"), row.names = F, quote = F, sep = "\t")

## ASV sequences as fasta file
uniquesToFasta(otu.filt_f, paste0("ASVs_dada2_", output_name, ".fasta")) # sequences of asv

# read written file in again, to get correct ids for sequences
seq_unique <- readDNAStringSet(paste0("ASVs_dada2_", output_name, ".fasta"))
ASV_id_seq <- data.frame(ASV_ID = gsub(";.*", "", names(seq_unique)), 
                         seq = as.character(seq_unique))

## ASV table with ASVs in rows and samples in columns, sequence and taxonomy
asv.print <- otu.filt_f %>% 
  t() %>%                                                  # transpose table
  as.data.frame() %>% 
  tibble::rownames_to_column() %>%                         # make sequences as own column
  full_join(ASV_id_seq, by = c("rowname" = "seq")) %>%     # merge ASV ID by sequence
  full_join(                                               # add taxonomy
    (
      tax_filt_f %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column()
    ), 
    by = "rowname") %>% 
  dplyr::relocate(ASV_ID, rowname, Kingdom, Phylum, Class, Order, Family, Genus) %>% # order dataframe
  dplyr::rename("seq" = "rowname")


# full table with all extra info
write.table(asv.print, paste0("asv_table_all_info_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)

# table with only ASV IDs, could be imported with ASV IDs as rownames
asv.print %>% 
  select(-c(seq, Kingdom, Phylum, Class, Order, Family, Genus)) %>% 
  write.table(paste0("asv_table_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)

# table with taxonomy
asv.print %>% 
  select(-c(seq)) %>% 
  write.table(paste0("asv_table_tax_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)

# taxonomy table with ASV IDs, could be imported with ASV IDs as rownames
asv.print %>% 
  select(-c(seq, rownames(otu.filt_f))) %>% 
  write.table(paste0("tax_table_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)


# Write output for metacoder ####
# need to transform sequences to ASV IDs or loading into R will take super long
# seq_table with ASVs in columns and samples in rows
asv.print %>% 
  select(-c(seq, Kingdom, Phylum, Class, Order, Family, Genus)) %>% 
  column_to_rownames("ASV_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  write.table(paste0("seq_table_for_metacoder_", output_name, ".txt"), quote = F, sep = "\t", row.names = T)

# tax_table with ASVs in rows and ranks of taxonomy in separate columns, similar to above but with ASV IDs as rownames
# could also use tax_table from above and import with one column used as rownames
asv.print %>% 
  select(-c(seq, rownames(otu.filt_f))) %>% 
  column_to_rownames("ASV_ID") %>% 
  write.table(paste0("tax_table_for_metacoder_", output_name, ".txt"), quote = F, sep = "\t", row.names = T)

# end of script ####
save.image(tf)

# proceed to metacoder import R script
