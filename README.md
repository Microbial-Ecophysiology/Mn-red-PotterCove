# Manganese reduction in Potter Cove sediments
## Scripts used to analyse sequence raw data, subsequent analyses and plotting
As supplemental to manuscript Wunder, L. et al. 2024: Manganese reduction in microbial communities involved in Antarctic surface sediments from Potter Cove (submitted)

## Amplicon sequence data analysis - workflow
Involves all scripts starting with `S_` and some scripts which are sourced during this analysis in folder [Scripts to source](Scripts_to_source)


Two datasets were analysed: 
1. Potter Cove in situ, shown as `insitu` in script names
2. Slurry incubation experiments, shown as `inc` in script names.

In principle, the scripts for both datasets are the same. The scripts provided here contain the used parameters, so when the rawdata is downloaded from ENA it can be directly used to replicate the analysis process.\
[Data in situ sample](https://www.ebi.ac.uk/ena/browser/view/PRJEB72873) \
[Data incubation samples](https://www.ebi.ac.uk/ena/browser/view/PRJEB72882)

1. bash script [S_01_dada2_seqprep](S_01_dada2_seqprep.bash) for quality control, demultiplexing and primer clipping \
This is the same for in situ and incubation data analysis, as no parameters are entered here.\
Mapping files for \
[in situ samples analysis](small_data/Insitu)\
[incubation samples analysis](small_data/Inc)
3. R script `S_02_dada2`([in situ](S_insitu_02_dada2.R), [incubation](S_inc_02_dada2.R)) for sequence trimming, error correction, dereplication, denoising,  read merging, chimera removal and taxonomic classification
4. R script `S_03_metacoder_import_data` ([in situ](S_insitu_03_metacoder_import_data.R), [incubation](S_inc_03_metacoder_import_data.R)) for importing ASV table created in previous script into a `taxmap` format, remove mitochondrial and chloroplast sequences and doubletons and singletons
5. R script `S_04_rarefaction-curves` ([in situ](S_insitu_04_rarefaction-curves.R), [incubation](S_inc_04_rarefaction-curves.R)) for calculating and plotting rarefaction curves and removing samples with too low coverage if necessary
6. R script `S_05_taxonomy_plot` ([in situ](S_insitu_05_taxonomy_plot.R), [incubation](S_inc_05_taxonomy_plot.R)) for plotting the taxonomy in bar and line plots

Files with some metadata needed to run the scripts can be found here for [in situ samples](small_data/Insitu/Insitu_mdata.txt) and [incubation samples](small_data/Inc/Inc_mdata.txt).


## Geochemical data
Calculate solid phase Mn wt% data into µmol/cm^3 and plot it in this [script](Monien-2014-data_plot.R). 

Plot geochemical in situ data with this [script](insitu_geochemistry_plot.R). The data is also uploaded to PANGAEA [here](https://doi.org/10.1594/PANGAEA.941109). 

Plot and do statistical analysis on geochemical data measured during incubation experiment with this [script](inc_geochemistry_plot.R) uploaded [here](small_data/Inc/Inc_geochem_data.txt).
