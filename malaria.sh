# Plasmodium data
nohup gmes_petap.pl --ES --sequence Plasmodium_berghei.genome --min_contig 10000 &
# Specify a minimum continuous sequence length of 10,000, which is used to filter out shorter sequences in the input genome sequence to avoid excessive noise or false positive results when predicting genes.
# nohup: command is used to run a command in the background without interrupting the execution of the command even if the terminal is closed. It ignores the hang signal, allowing the command to continue running.
# &: symbol puts the previous command in the background, allowing you to continue executing other commands in the terminal session without waiting for the current command to finish


# Processing of Haemoproteus tartakovskyi data
# Clean genome sequence
# the input reads derive both from the bird and the parasite, bird scaffolds should be removed
python removeScaffolds.py Haemoproteus_tartakovskyi.genome 35 Ht.genome 3000
# Use the script provided, removeScaffold.py, provided
# Conditions are: 
# the length of the sequence must be greater than or equal to 3000
# the GC content must be less than or equal to 35%

jobs # check whether it's running

# gffPaese.pl
cp /resources/binp28/Data/gffParse.pl .
# need to delete the old gffParse.pl
cp /resources/binp28/Data/gffParse.pl .
chmod +x gffParse.pl 
./gffParse.pl
ls -lh

# Still have some scaffolds that derive from the bird. These should be short
# First to do is fix the headers format in the file, otherwise the gffParse.pl can't be used
cat Haemoproteus.gff | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gff
./gffParse.pl -i Ht.genome -g ~/malaria/Ht_Prediction/Ht2.gff -f CDS -p -a gene_id -c
# output files: gffParse.fna, gffParse.log, gffParse.faa

# INFO: The gff or gtf file /home/inf-39-2023/malaria/Ht_Prediction/Ht2.gff has successfully been parsed.
#       There were 2230 scaffolds containing genes.
#       The scaffolds contained 3763 genes.
#       The genes contained the feature CDS 8539 times.

# INFO: The scaffold/genome file Ht.genome was successfully parsed.
#       There were 2343 scaffolds. This number may be higher than the one above.

# INFO: Following files were output in the directory
#       /home/inf-39-2023/malaria:
#       gffParse.fna (fasta file containing the genes)
#       gffParse.log (log file)
#       gffParse.faa (fasta file containing translated genes)

# NOTICE: There are 1296 warnings in the log file.
#         Reading frames have been adjusted for 1296 genes (with the use of -c option).


# Blast
conda create -n blast_env
conda activate blast_env
conda install -c bioconda blast

# binp28 Annotation:
# Look what NCBI offers when it comes to sequence databases (BLAST format). This takes 30 seconds or so to load:
update_blastdb.pl --showall

# Check the path of blastp
whereis blastp
# blastp: /usr/local/bin/blastp /home/inf-39-2023/anaconda3/envs/blast_env/bin/blastp
# Because have two blastp, so there is a problem in runnning

# In top
# chaeck ID: 
u inf-39-2023
c
# Blastp
nohup /usr/local/bin/blastp -query ~/malaria/Ht_Prediction/gffParse.faa -db SwissProt -out Ht.blastp -num_descriptions 10 -num_alignments 5 -num_threads 25 &


# Retrieve the host scaffolds, redirect to the output file
python datParser.py blastp/Ht.blastp Ht_Prediction/gffParse.fna taxonomy.dat uniprot_sprot.dat > scaffolds.txt
wc -l scaffolds.txt 
# 23 scaffolds.txt

# Removal: clean the fasta file from bird scaffolds
python scripts/ids_removal.py Ht.genome Ht_Prediction/scaffolds.txt Ht_removed.fna
cat Ht.genome | grep -c '^>'
# 2343
cat Ht_removed.fna | grep -c '^>'
# 2320

# Rerun a gene prediction after cleaning the birds scaffolds
$ ./gffParse.pl -i Ht_removed.fna -g ~/malaria/Ht_Prediction/Ht2.gff -f CDS -p -a gene_id -c

# INFO: The gff or gtf file /home/inf-39-2023/malaria/Ht_Prediction/Ht2.gff has successfully been parsed.
#       There were 2230 scaffolds containing genes.
#       The scaffolds contained 3763 genes.
#       The genes contained the feature CDS 8539 times.

# INFO: The scaffold/genome file Ht_removed.fna was successfully parsed.
#       There were 2320 scaffolds. This number may be higher than the one above.

# INFO: Following files were output in the directory
#       /home/inf-39-2023/malaria:
#       gffParse.fna (fasta file containing the genes)
#       gffParse.log (log file)
#       gffParse.faa (fasta file containing translated genes)

# NOTICE: There are 1281 warnings in the log file.
#         Reading frames have been adjusted for 1281 genes (with the use of -c option).

# Change the format of Pk file
# >chromosome10
# chromosome10
# CCTTTCTCCTTAACTTCCCCTCTATGAAGAGATTATTTGTCCTCACTACTCCTTACCTTCTTTCCCCTGT
# AAGAAGAGATTATTTTGCCCCTTACTACTTTTTACCTTCCTGTAAGAAGAGATTATTTTGTCCCTTACTA
cat Plasmodium_knowlesi.genome | grep -v '^chromosome' > Plasmodium_knowlesi2.genome

# gffParse.pl: Ensure that the gffParse.pl executable can be executed directly from any directory
# mkdir bin
# cd bin
# cd ..
# vi .bashrc
# source .bashrc
# echo $PATH
# /home/inf-39-2023/anaconda3/envs/GeneMark/bin:/home/inf-39-2023/anaconda3/condabin:/home/inf-39-2023/.local/bin:/home/inf-39-2023/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:~/bin
# cd bin
# chmod +x gffParse.pl 
# ls
# gffParse.pl


mkdir fasta
cd fasta/
gffParse.pl -i ~/malaria/plasmodiumGenome/Plasmodium_berghei.genome -g ~/malaria/plasmodiumGenome/Prediction/genemark.Pb.gtf -p -a gene_id -c
gffParse.pl -i ~/malaria/plasmodiumGenome/Plasmodium_cynomolgi.genome -g ~/malaria/plasmodiumGenome/Prediction/genemark.Pc.gtf -p -a gene_id -c
gffParse.pl -i ~/malaria/plasmodiumGenome/Plasmodium_faciparum.genome -g ~/malaria/plasmodiumGenome/Prediction/genemark.Pf.gtf -p -a gene_id -c
gffParse.pl -i ~/malaria/plasmodiumGenome/Plasmodium_knowlesi2.genome -g ~/malaria/plasmodiumGenome/Prediction/fixed_Pk.gtf -p -a gene_id -c
gffParse.pl -i ~/malaria/plasmodiumGenome/Plasmodium_vivax.genome -g ~/malaria/plasmodiumGenome/Prediction/genemark.Pv.gtf -p -a gene_id -c
gffParse.pl -i ~/malaria/plasmodiumGenome/Plasmodium_yoelii.genome -g ~/malaria/plasmodiumGenome/Prediction/genemark.Py.gtf -p -a gene_id -c
gffParse.pl -i ~/malaria/plasmodiumGenome/Toxoplasma_gondii.genome -g ~/malaria/plasmodiumGenome/Prediction/Tg.gff -p -a gene_id -c
# -F overwrite


### Busco
# Install Busco on the server
conda create -n Busco
conda activate Busco
conda install -c bioconda busco=5.6.1
busco --version
conda search busco
conda config -h
conda config --add channels bioconda
conda config --add channels conda-forge
conda install busco # install the latest version

# Install proteinortho on the server
conda create -n proteinortho
conda activate proteinortho
conda install -c bioconda proteinortho


# Change the format
awk '/^>/ {print $1; next} { print $0}' ~/malaria/plasmodiumGenome/fasta/Pf_fasta/gffParse.faa > Pf.faa

# Proteinortho 
nohup proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &

$ ls
# Ht.faa                 myproject.info                        myproject.proteinortho.html  Pb.faa               Pc.faa.diamond.dmnd  Pk.faa               Pv.faa.diamond.dmnd  Tg.faa
# Ht.faa.diamond.dmnd    myproject.proteinortho-graph          myproject.proteinortho.tsv   Pb.faa.diamond.dmnd  Pf.faa               Pk.faa.diamond.dmnd  Py.faa               Tg.faa.diamond.dmnd
# myproject.blast-graph  myproject.proteinortho-graph.summary  nohup.out                    Pc.faa               Pf.faa.diamond.dmnd  Pv.faa               Py.faa.diamond.dmnd

# Ht.faa.diamond.dmnd, Pb.faa.diamond.dmnd, Pc.faa.diamond.dmnd, Pf.faa.diamond.dmnd, Pk.faa.diamond.dmnd,  Pv.faa.diamond.dmnd, Py.faa.diamond.dmnd, Tg.faa.diamond.dmnd: A protein database file created by the diamond tool.
# myproject.info: A file that contains information about the project, possibly including the number of input files, parameter Settings, and so on.
# myproject. Proteinortho. HTML: HTML proteinortho6. Pl output file, may contain information about visualization of clustering results.
# myproject.proteinortho-graph: proteinortho6.pl output graph file that may contain graphical information about clustering results.
# myproject.proteinortho-graph.summary: Summary information about the proteinortho6.pl graphics file.
# myproject. Proteinortho. TSV: proteinortho6. Pl output TSV format files, may contain detailed information about the clustering results.
# nohup.out: output file of the nohup command, which contains the output information of the proteinortho6.pl command.

# BUSCO
busco -i ../faa_formated/Pb.faa -o Pb -m prot -l apicomplexa
busco -i ../../faa_formated/Pc.faa -o Pc -m prot -l apicomplexa
busco -i ../../faa_formated/Pf.faa -o Pf -m prot -l apicomplexa
busco -i ../../faa_formated/Pk.faa -o Pk -m prot -l apicomplexa
busco -i ../../faa_formated/Pv.faa -o Pv -m prot -l apicomplexa
busco -i ../../faa_formated/Py.faa -o Py -m prot -l apicomplexa
busco -i ../../faa_formated/Tg.faa -o Tg -m prot -l apicomplexa
busco -i ../../faa_formated/Ht.faa -o Ht -m prot -l apicomplexa


# protein fasta file
cat Pb_busco/Pb/run_apicomplexa_odb10/full_table.tsv | awk '/^[^#]/ && ($2 == "Complete" || $2 == "Duplicated") {print $1, $3}' > Pb_IdSeq.txt
cat *_busco/*/run_apicomplexa_odb10/full_table.tsv | awk '/^[^#]/ && ($2 == "Complete" || $2 == "Duplicated") {print $1, $3}' > *_IdSeq.txt

awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Pb_busco/Pb/run_apicomplexa_odb10/full_table.tsv > Pb_IdSeq2.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Pc_busco/Pc/run_apicomplexa_odb10/full_table.tsv > Pc_IdSeq.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Pf_busco/Pf/run_apicomplexa_odb10/full_table.tsv > Pf_IdSeq.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Pk_busco/Pk/run_apicomplexa_odb10/full_table.tsv > Pk_IdSeq.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Pv_busco/Pv/run_apicomplexa_odb10/full_table.tsv > Pv_IdSeq.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Py_busco/Py/run_apicomplexa_odb10/full_table.tsv > Py_IdSeq.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Tg_busco/Tg/run_apicomplexa_odb10/full_table.tsv > Tg_IdSeq.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1, $3}' Ht_busco/Ht/run_apicomplexa_odb10/full_table.tsv > Ht_IdSeq.txt

python ~/malaria/scripts/extract.py all_IDs.txt Ht_IdSeq.txt Pb_IdSeq.txt Pc_IdSeq.txt Pf_IdSeq.txt Pk_IdSeq.txt Pv_IdSeq.txt Py_IdSeq.txt Tg_IdSeq.txt output.txt

# Extract protein fasta
# Install the pandas library in your current conda environment
conda create -n myenv python=3.8
conda activate myenv
conda install -c conda-forge biopython
conda install pandas

# out is the directory that stores the sequence file corresponding to the busco id
ls -l out | wc -l
# 224
# 224-1=223
ls -l out
# total 1732
# -rw-r--r--. 1 inf-39-2023 students  4875 Feb 29 10:42 10090at5794.faa
# -rw-r--r--. 1 inf-39-2023 students  4673 Feb 29 10:43 10141at5794.faa
# -rw-r--r--. 1 inf-39-2023 students  4737 Feb 29 10:43 10373at5794.faa
# ...

### Alignments and trees
conda create -n clustalo
conda activate clustalo
conda install -c bioconda clustalo raxml
conda install clustalo
conda install raxml

# Multiple sequence alignment:
clustalo -i input.fasta -o output.aln --outfmt=clustal
# Read the sequence from the input.fasta file and compare them, then save the result to the output.aln file.

# Evolutionary analysis（RAxML）：
raxmlHPC-PTHREADS -f a -s output.aln -n output_tree -m GTRGAMMA -T 4
# Perform RAxML analysis using the alignment results in the output.aln file and generate an evolutionary tree.

# Clustalo
for file in *.faa; do clustalo -i "$file" -o ~/malaria/extract_fasta/clustalo/${file%.faa}_aligned.faa -v; done
# generate *_aligned.faa output files

# Raxml
for file in *_aligned.faa; do raxmlHPC -s "$file" -n ${file%_aligned.faa}.tre -o Tg -m PROTGAMMABLOSUM62 -p 12345; done

# Extract best_tree files to merge them into a file 
mkdir merge_tree
cp RAxML_bestTree.* > ../merge_tree .
cat RAxML_bestTree* > Allaligned_trees.tre
# Consense
consense Allaligned_trees.tre
# consense: can't find input tree file "intree"
# Please enter a new file name> Allaligned_trees2.tre       
# consense: can't find input tree file "Allaligned_trees2.tre"
# Please enter a new file name> Allaligned_trees.tre  


# outtree file: This file contains the generated consensus tree. Consensus tree is obtained by merging individual input trees, which represents the common feature or common branch pattern of input trees. In this file, each line represents a node, including the node's branch length and label information. You can use this file to view the topology and branch information of the merged consensus tree.
# outfile file: This file contains the detailed output of the consense run, usually including some statistics, prompts during the run, and so on. This file may provide additional information about generating the consensus tree, such as the number of input trees, merging methods, and so on. The outfile file is also a text file that you can open with a text editor to view its contents.

# Set Tg as outgroup
# Consensus tree program, version 3.697

# Settings for this run:
#  C         Consensus type (MRe, strict, MR, Ml):  Majority rule (extended)
#  O                                Outgroup root:  No, use as outgroup species  1
#  R                Trees to be treated as Rooted:  No
#  T           Terminal type (IBM PC, ANSI, none):  ANSI
#  1                Print out the sets of species:  Yes
#  2         Print indications of progress of run:  Yes
#  3                               Print out tree:  Yes
#  4               Write out trees onto tree file:  Yes

# Are these settings correct? (type Y or the letter for one to change)
# R

less outtree
# ((((((Pc:223.0,Pv:223.0):197.0,(Pb:223.0,Py:223.0):203.0):64.0,Pf:223.0):104.0,Ht:223.0):139.0,Pk:223.0):223.0,Tg:223.0);


### Questions:
# 1. Do you think that in a phylogenetic tree the parasites that use similar hosts will group together?
# Yes, in a phylogenetic tree, organisms that share a closer evolutionary relationship are expected to be grouped together.

# 2. Why?
# Sequences with more than 35% GC content are deleted to remove scaffolds that might come from bird DNA, which typically has a higher GC content. However, even when sequences with high GC content are removed, there may still be contamination or noise from other sources, such as bacterial DNA or sequencing errors

# Q: In SwissProt, the five characters after the underscore is a species abbreviation (what is before the underscore?)
# In the SwissProt database, the part before the underscore (_) is usually the short or abbreviation of the species, while the five characters after the underscore are the species-specific abbreviation or identifier.
# Typically, the section before the underline is the protein's name or description, which identifies the protein's function, subcellular location, and so on. The five characters following the underscore are usually abbreviations for the corresponding species of the protein and are used to identify the protein

# 3. Insert the missing data in the above table. Use bash, not internet!
### Genome size
for file in *.genome ; do echo $file ; grep -v \> "$file" | tr -d "\n" | wc -c; done
# Plasmodium_berghei.genome
# 17954629
# Plasmodium_cynomolgi.genome
# 26181343
# Plasmodium_faciparum.genome
# 23270305
# Plasmodium_knowlesi2.genome
# 23462187
# Plasmodium_knowlesi.genome
# 23462346
# Plasmodium_vivax.genome
# 27007701
# Plasmodium_yoelii.genome
# 22222369
# Toxoplasma_gondii.genome
# 128105889
# Haemoproteus_tartakovskyi.raw.genome
cat Haemoproteus_tartakovskyi.raw.genome | grep -v \> | tr -d "\n" | wc -c
# 27426784

### Genes
for dir in *_fasta/ ; do for file in "$dir"*.fna ; do echo "$file"; grep -c '^>' "$file"; done; done
# Pb_fasta/gffParse.fna
# 7235
# Pc_fasta/gffParse.fna
# 5787
# Pf_fasta/gffParse.fna
# 5207
# Pk_fasta/gffParse.fna
# 4952
# Pv_fasta/gffParse.fna
# 5682
# Py_fasta/gffParse.fna
# 4889
# Tg_fasta/gffParse.fna
# 15892
# Haemoproteus_tartakovskyi.raw.genome
cat Haemoproteus_tartakovskyi.raw.genome | grep -c '^>'
# 15048

### GC content
for file in *.genome ; do echo $file; awk '/^[^>]/ {at += gsub(/[ATat]/,"")} {gc+=gsub(/[GCgc]/,"")} END { printf "GC content: %.2f%%\n", (gc/(gc + at))*100}' $file; done
# Plasmodium_berghei.genome
# GC content: 23.75%
# Plasmodium_cynomolgi.genome
# GC content: 40.38%
# Plasmodium_faciparum.genome
# GC content: 19.37%
# Plasmodium_knowlesi2.genome
# GC content: 38.83%
# Plasmodium_knowlesi.genome
# GC content: 38.83%
# Plasmodium_vivax.genome
# GC content: 42.29%
# Plasmodium_yoelii.genome
# GC content: 21.77%
# Toxoplasma_gondii.genome
# GC content: 52.35%
# Haemoproteus_tartakovskyi.raw.genome 
awk '/^[^>]/ {at += gsub(/[ATat]/,"")} {gc+=gsub(/[GCgc]/,"")} END { printf "GC content: %.2f%%\n", (gc/(gc + at))*100}' Haemoproteus_tartakovskyi.raw.genome 
# GC content: 27.52%

# 4. Compare the genome sizes with other eukaryotes and bacteria. Discuss with your partner (that is student partner) the reason for the observed genome sizes.

# 5. What may cause the biased GC-contents in some of the species?
# Biased GC-content in species can be caused by biological selection favoring stability, mutation and genetic drift, horizontal gene transfer, recombination, environmental factors, and genomic architecture. These factors interact to shape the GC-content observed in different organisms.

# 6. What does the curly braces notation stand for?
# The curly braces notation {} in the command proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa is used for brace expansion in Bash. Brace expansion allows you to generate arbitrary strings based on comma-separated lists or ranges inside the braces.

# Q: What does the flag -l stands for? and why do we choose apicomplexa?
# The -l flag in the busco command stands for the lineage dataset to be used for assessing the completeness of the input sequences. In this case, -l apicomplexa indicates that the lineage dataset for the Apicomplexa phylum will be used.

# 7. Compare how many BUSCOs (orthologues proteins) that are found in each proteome. Do the investigated parasites have close to complete numbers of BUSCOs?
# Pb: 361    Complete BUSCOs (C)  C:80.9%  
# Pc: 429    Complete BUSCOs (C)  C:96.2%  
# Pf: 436    Complete BUSCOs (C)  C:97.8%
# Pk: 433    Complete BUSCOs (C)  C:97.1%
# Pv: 437    Complete BUSCOs (C)  C:97.9%
# Py: 435    Complete BUSCOs (C)  C:97.5%
# Tg: 380    Complete BU298  (C)  C:85.2%
# Ht: 298    Complete BUSCOs (C)  C:66.8%
# Complete BUSCOs account for 80.9% of the total BUSCOs. This means that a fairly high percentage of the genomes of these parasites are conserved genes shared with other species

# 8. Do you think that the assembly of the Haemoproteus tartakowskyi genome is a reasonable approximation of the true genome?
# According to the provided data, the proportion of BUSCOs in the Haemoproteus tartakowskyi genome is not very close to completeness. While there are some complete BUSCOs present, the proportion of missing BUSCOs is relatively high at 28.5%, indicating that there may be some incompleteness in the assembly of this genome. Therefore, it can be said that the assembly of this genome may be a reasonable approximation, 

# 9. How many of the BUSCOs are found in all eight organisms?
# Extract all busco ids in complete and duplicated shown once
# Use uniq -c to count the sorted ids and output the number of times each ID appears in all files
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1}' Pb_busco/Pb/run_apicomplexa_odb10/full_table.tsv > Pb_buscoIds.txt
awk '/^[^#]/ && ($2 == "Complete" || ($2 == "Duplicated" && !seen[$1]++)) {print $1}' *_busco/*/run_apicomplexa_odb10/full_table.tsv > *_buscoIds.txt
awk '{print $1}' *.tsv | sort | uniq -c | awk -v num=8 '$1 == num {print $2}' | wc -l
# 223

# 10. If Toxoplasma is removed, how many BUSCOs are shared among the remaining seven species. Interpret!
awk '{print $1}' Ht.tsv Pb.tsv Pc.tsv Pf.tsv Pk.tsv Pv.tsv Py.tsv | sort | uniq -c | awk -v num=7 '$1 == num {print $2}' | wc -l
# 252

# 11. Does all protein trees reflect the “true’ ’ species tree?
# No

# 12. What is the phylogenetic position of Plasmodium falciparum?
# First branch

# 13. Do you think that the GC contents have an impact on the tree topology?
# Yes.
# GC content affects the rate and type of mutations, particularly transitions (A <-> G, C <-> T) and transversions (A <-> T, G <-> C). Strong mutation biases towards GC-rich or AT-rich regions may result in differences in substitution rates, thereby influencing tree topology.
# Some organisms exhibit biased DNA mismatch repair favoring GC over AT alleles, known as biased gene conversion. This process can alter the pattern of nucleotide substitutions, thereby affecting tree topology.
# GC content may be influenced by selection pressures related to factors such as codon usage bias, gene expression levels, and environmental adaptation. Genomic regions under different selection pressures may exhibit varying GC content, impacting evolutionary rates and tree topology.
# Differences in GC content between genomic regions (e.g., coding vs. non-coding, introns vs. exons) can affect patterns of molecular evolution, introducing heterogeneity in substitution rates and potentially influencing tree topology.
# Gene flow and hybridization between species with different GC contents can lead to inconsistent gene trees, affecting overall tree topology.

# 14. Do you think that the host range has an impact on the tree topology?
# The host range influences the tree topology of phylogenetic analysis. It may lead to host switching events, co-evolutionary dynamics, horizontal gene transfer, niche conservation, and parasitic speciation, resulting in specific patterns and structures in tree topologies. Therefore, the host range and its effects on related organisms need to be considered when interpreting the results of phylogenetic analyses

# 15. Are the BUSCO proteins also found as orthologs in the proteinortho output?
# Yes

# 16. Make a script that concatenates the alignments for each organism and BUSCO into one fasta file that in the end should contain seven sequences. Alternatively, use bash.
# 17. Make a tree of this “superalignment’ ’. Does it correspond to the consense tree






