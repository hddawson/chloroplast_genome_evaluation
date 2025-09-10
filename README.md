I need to analyze the chloroplast genomes available online. 

I want to assess quality through evaluation of sequence. 

I'll start with Length and GC, then - I think - taxonomy not to be silly. 

Then I think - and this will be the bioinfun part - I'll need to find the inverted repeats. 

Ideally minimal reference bias. 


First, download the genomes 

pixi run python src/1_download.py

Produces

data/genomes/
├──123.fa
├──456.fa
├──....fa

and, importantly, taxonomy_info, which contains the organism labels according to NCBI.

Second, annotate the genomes 

ls data/genomes >> data/listGenomes.txt
src/2_annotate.sh data/listGenomes.txt

This failed for three genomes:
MK560340.1
ON550389.1
ON550390.1

third, get occurrences by querying the organism name into GBIF, just get the big data

XXX_processGBIFData.py --> find overlap between taxonomy_info and gbif search
XXX_CleanOccurrences.R --> use r packages to clean occurrences


At some point, process data from TRY 

XXX_processTRYData.r

walk throughd dataProcessing.ipynb until you get to alignmet step, then,
src/aligner.sh


Then, we will also need to query GBIF
pixi run python src/3_getOccurrences.py data/taxonomy_info.csv

Fourth, process occurences into environmental data

pixi run python src/3_getOccurrences.py data/taxonomy_info.csv data/occurrences/

(head -n 1 $(find data/occurrences -name '*_clean.csv' | head -n 1) && find data/occurrences -name '*_clean.csv' | xargs -I {} tail -n +2 {}) > data/combinedOccurrences_BIEN.csv

XXX_EnvDataFromTifs.r --> 
XXX_processEnvData.r --> 




#find the species with few occurrences
pixi run python src/3_point_5_reQuery.py 

 pixi run python src/3_getOccurrences.py data/taxonomy_missing_or_under10.csv data/occurrences_query2/

#I will also look again at the species I found, just in case
pixi run python src/3_getOccurrences.py data/taxonomy_found_try1.csv data/occurrences_query3/

(head -n 1 $(find data/occurrences -name '*_clean.csv' | head -n 1) && find data/occurrences -name '*_clean.csv' | xargs -I {} tail -n +2 {}) > data/combinedOccurrences.csv

Querying GBIF and BIEN is kind of sad, because sometimes the API gives out and you recieve no results. So we check combinedOccurrences.csv occurence against taxonomy_info.csv and make a save a subset of taxonomy data for all the species (Organism in taxonomy, queryTerm in Occurrences) that have fewer than 10 occurrences 

(head -n 1 $(find data/occurrences_query2 -name '*_clean.csv' | head -n 1) && find data/occurrences_query2 -name '*_clean.csv' | xargs -I {} tail -n +2 {}) > data/combinedOccurrences_2.csv


then I get the envData from the tifs 
pulling_envData.r


make tree

/programs/FastTree-2.1.11/FastTree -nt -gtr -gamma -log tree.log < data/tmp/cds_supermatrix.fasta > tree.nwk

/programs/raxml-ng_v1.2.0/raxml-ng --all --msa ../data/tmp/cds_supermatrix.fasta --model GTR+G --bs-trees 10 --threads auto
