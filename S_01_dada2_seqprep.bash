## Script for raw sequence data quality check, demultiplexing and primer clipping in preparation to rest of dada2 pipeline

# 0. download data
# done previousely
# make sure data is in correct directory
WDIR=""
cd $WDIR
ls $WDIR/raw_reads/


# 1. Quality control
WDIR=""
cd $WDIR

## run fastqc (this was actually done previousely with same fastqc version, copy results into folder)
mkdir quality-control_fastqc

conda activate metawrap-env
# version fastqc v0.11.8

fastqc $WDIR/raw_reads/*fq.gz -o $WDIR/quality-control_fastqc/

# looks fine, lots of duplicated Sequences


# 2. Prepare data and directories
WDIR=""
cd $WDIR

# list target kingdom in file for present libraries
nano target.txt
Bac
Arc

# look into directory 'raw_reads' and extract information which files are present
# the code is explained in detail in the explanation_dada2_seqprep.md file in the Github repository

## list all library IDs
ls -1 $WDIR/raw_reads/*fq.gz | xargs -n1 basename | while read line
 do
  res=${line//[^_]}
  no=${#res}
  lib=$(echo ${line} | cut -f1-$((no-3)) -d '_')
  echo $lib
 done | sort | uniq > lib_IDs.txt


## list all library IDs and library IDs from novogene with flowcell and lane ID
ls -1 $WDIR/raw_reads/*.fq.gz | xargs -n1 basename | while read line
 do
  res=${line//[^_]}
  no=${#res}
  
  lib=$(echo ${line} | cut -f1-$((no-3)) -d '_')
  libnovo=$(echo ${line} | cut -f$((no-2)) -d '_')
  flowcell=$(echo ${line} | cut -f$((no-1)) -d '_')
  lane=$(echo ${line} | cut -f$no -d '_')
  echo -e "${lib}\t${libnovo}\t${flowcell}_${lane}"
 done | sort | uniq > lib_flowcell_IDs.txt


## Identify which flowcell + lines are identical among libraries:
ls -1 $WDIR/raw_reads/*fq.gz | xargs -n1 basename | while read line
 do
  res=${line//[^_]}
  no=${#res}
  flowcell=$(echo ${line} | cut -f$((no-1)) -d '_')
  lane=$(echo ${line} | cut -f$no -d '_')
  echo $flowcell"_"$lane
 done | sort | uniq  > $WDIR/flowcell_IDs.txt

cat flowcell_IDs.txt


## make directories for each flowcell
while read flowcell
do
	mkdir -p seq_by_flowcell/${flowcell}/raw_reads seq_by_flowcell/${flowcell}/Logfiles seq_by_flowcell/${flowcell}/Demux_Bac seq_by_flowcell/${flowcell}/Demux_Arc seq_by_flowcell/${flowcell}/Clipped_Bac seq_by_flowcell/${flowcell}/Clipped_Arc
done < flowcell_IDs.txt


## move rawreads to according directories
while read flowcell
do
 mv $WDIR/raw_reads/*${flowcell}*".fq.gz" $WDIR/seq_by_flowcell/${flowcell}/raw_reads/ 
done < flowcell_IDs.txt


## need mappingfiles from each library
mkdir mappingfiles
# structure name: Bac_Libxx_mappingfile.txt or Arc_Libxx_mappingfile.txt
 ## fix line endings of files on windows PC first
#dos2unix *_mappingfile.txt


# 3. Demultiplexing

conda activate trimming
conda list # check if cutadapt is loaded, check version
# cutadapt version 3.1
WDIR=""
cd $WDIR
THREADS=80

while read lib libnovo flowcell
do
  
  # for bacteria
  TARGET="Bac"

  sed '1d' mappingfiles/$TARGET"_"$lib"_mappingfile.txt" | while read line
    do
    SID=$(echo "${line}" | cut -f1)
    BCD=$(echo "${line}" | cut -f2 | sed 's/^/\^/')
    OLP=$(expr ${#BCD} - 1)
  
    cutadapt -j $THREADS -O ${OLP} --no-indels -e 0 -g ${BCD} -G ${BCD} --discard-untrimmed -o seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R1.fastq.gz" -p seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R2.fastq.gz" seq_by_flowcell/${flowcell}/raw_reads/${lib}"_"${libnovo}"_"${flowcell}"_1.fq.gz" seq_by_flowcell/${flowcell}/raw_reads/${lib}"_"${libnovo}"_"${flowcell}"_2.fq.gz" > seq_by_flowcell/${flowcell}/Logfiles/${SID}".demux.log" 2>&1
 
    done
  
  # for archaea
  TARGET="Arc"

  sed '1d' mappingfiles/$TARGET"_"$lib"_mappingfile.txt" | while read line
    do
    SID=$(echo "${line}" | cut -f1)
    BCD=$(echo "${line}" | cut -f2 | sed 's/^/\^/')
    OLP=$(expr ${#BCD} - 1)
  
    cutadapt -j $THREADS -O ${OLP} --no-indels -e 0 -g ${BCD} -G ${BCD} --discard-untrimmed -o seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R1.fastq.gz" -p seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R2.fastq.gz" seq_by_flowcell/${flowcell}/raw_reads/${lib}"_"${libnovo}"_"${flowcell}"_1.fq.gz" seq_by_flowcell/${flowcell}/raw_reads/${lib}"_"${libnovo}"_"${flowcell}"_2.fq.gz" > seq_by_flowcell/${flowcell}/Logfiles/${SID}".demux.log" 2>&1
 
    done
	
done < lib_flowcell_IDs.txt


## unzip files, needed for steps after clipping anyways
gunzip seq_by_flowcell/*/Demux_Bac/*.gz
gunzip seq_by_flowcell/*/Demux_Arc/*.gz


# 4. Primer clipping
conda activate trimming
WDIR=""
cd $WDIR
THREADS=40

while read lib libnovo flowcell
do
  while read TARGET
  do

  sed '1d' mappingfiles/$TARGET"_"$lib"_mappingfile.txt" | while read line
   do
   SID=$(echo "${line}" | cut -f1)
   FWD=$(echo "${line}" | cut -f3 | sed 's/^/\^/')
   REV=$(echo "${line}" | cut -f4 | sed 's/^/\^/')
   OFWD=$(expr ${#FWD} - 2)
   OREV=$(expr ${#REV} - 2)
   ERROR=0.16
      
      # process fwd-rev orientation
      # m = minimum length, to ensure no empty sequences produced in output
      cutadapt -j $THREADS --no-indels -e ${ERROR} -g "${FWD};o=${OFWD}" -G "${REV};o=${OREV}" -m 50 --discard-untrimmed -o seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/${SID}"_clip_fr_R1.fastq" -p seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/${SID}"_clip_fr_R2.fastq" seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R1.fastq" seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R2.fastq" > seq_by_flowcell/${flowcell}/Logfiles/${SID}".clip_fr.log" 2>&1
	  
	  # process rev-fwd orientation
	  cutadapt -j $THREADS --no-indels -e ${ERROR} -g "${REV};o=${OREV}" -G "${FWD};o=${OFWD}" -m 50 --discard-untrimmed -o seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/${SID}"_clip_rf_R1.fastq" -p seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/${SID}"_clip_rf_R2.fastq" seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R1.fastq" seq_by_flowcell/${flowcell}/"Demux_"$TARGET/${SID}"_R2.fastq" > seq_by_flowcell/${flowcell}/Logfiles/${SID}".clip_rf.log" 2>&1
    
   done
   
  done < target.txt
	
done < lib_flowcell_IDs.txt


# 5. Count sequences for each step and switch to dada2 in R
WDIR=""
cd $WDIR

while read lib libnovo flowcell
do
  while read TARGET
  do
	ls -1v seq_by_flowcell/${flowcell}/"Demux_"$TARGET/*_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste <(ls -1v seq_by_flowcell/${flowcell}/"Demux_"$TARGET/*_R1.fastq | xargs -n1 basename | sed 's/_R1\.fastq//') - > tmp1
	ls -1v seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/*_clip_fr_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste tmp1 - > tmp2
	ls -1v seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/*_clip_rf_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste tmp2 - > tmp3
	echo -e 'SID\tDemux\tClipped_fr\tClipped_rf' | cat - tmp3 > "nSeqs_"${flowcell}"_"${TARGET}".txt"
	rm tmp*
	
  done < target.txt

### extract information about demultiplexing efficiency from logfiles
	ls -1v seq_by_flowcell/${flowcell}/Logfiles/*.demux*.log | xargs awk '/^Total read pairs processed/{print $NF}' | paste <(ls -1v seq_by_flowcell/${flowcell}/Logfiles/*.demux*.log | xargs -n1 basename | sed 's/.demux\.log//') - > tmp1
	ls -1v seq_by_flowcell/${flowcell}/Logfiles/*.demux*.log | xargs awk '/Read 1 with adapter:/{print $5}' | paste tmp1 - > tmp2
	ls -1v seq_by_flowcell/${flowcell}/Logfiles/*.demux*.log | xargs awk '/Read 2 with adapter:/{print $5}' | paste tmp2 - > tmp3
	ls -1v seq_by_flowcell/${flowcell}/Logfiles/*.demux*.log | xargs awk '/Pairs written/{print $5}' | paste tmp3 - > tmp4
	sed -i 's/,//g' tmp4
	echo -e 'SID\tTotal_read_pairs\tRead1_with_adapter\tRead2_with_adapter\tPairs_written' | cat - tmp4 > "reads_demux_"${flowcell}".txt"
	rm tmp*

done < lib_flowcell_IDs.txt


# 6. Check for files which have less then 100 reads and move them into other directory
WDIR=""
cd $WDIR

while read lib libnovo flowcell
do
  while read TARGET
  do
  sed '1d' "nSeqs_"$flowcell"_"$TARGET".txt" | while read line
   do
   SID=$(echo "${line}" | cut -f1)
   Demux=$(echo "${line}" | cut -f2)
     if [[ $Demux -lt 100 ]]                   # check if the sample has after demultiplexing less then 100 reads
     then
	   mkdir seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/empty/
       mv seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/${SID}*".fastq" seq_by_flowcell/${flowcell}/"Clipped_"$TARGET/empty/
	 fi
   done
  done < target.txt
done < lib_flowcell_IDs.txt

# no "empty" directory created, no sample moved

## prepare directories for subsequend analysis for Bacteria and Archaea separately
mkdir -p Bacteria/QualityProfiles Bacteria/Logfiles Bacteria/ErrorProfiles Archaea/QualityProfiles Archaea/Logfiles Archaea/ErrorProfiles


## proceed to separate dada2 R scripts for bacteria and archaea samples