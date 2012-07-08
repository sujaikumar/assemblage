assemblage [əˈsɛmblɪdʒ]  
*noun*

1. a number of things or persons assembled together; collection; assembly
2. (Cookery) a list of dishes served at a meal or the dishes themselves
3. the act or process of assembling or the state of being assembled
4. (Fine Arts & Visual Arts / Art Terms) a three-dimensional work of art that combines various objects into an integrated whole

From Collins English Dictionary – Complete and Unabridged © HarperCollins Publishers 1991, 1994, 1998, 2000, 2003

About
=====

This is a set of scripts for working with genome assemblies, making taxon-annotated GC-coverage plots (a.k.a blob plots), extracting reads belonging to the blobs, etc

Example Data Sets
=================

As an example, we use the following study from the Short Read Archive:

* ERP001495 - De novo whole-genome sequence of the free-living nematode Caenorhabditis sp. 5 strain JU800 DRD-2008

This study has one sample:

* ERS147916 - Caenorhabditis sp. 5 JU800 DRD-2008  

Two libraries (or "Experiments" in SRA terminology) were run for this sample:

* ERX114449 - 300 bp library with Illumina HiSeq2000 101 bp PE sequencing
* ERX114450 - 600 bp library with Illumina HiSeq2000 101 bp PE sequencing

And these are the four files

* `g_ju800_110714HiSeq300_1.txt.gz - 300 bp library forward read`
* `g_ju800_110714HiSeq300_2.txt.gz - 300 bp library reverse read`
* `g_ju800_110714HiSeq600_1.txt.gz - 600 bp library forward read`
* `g_ju800_110714HiSeq600_1.txt.gz - 600 bp library reverse read`

How-Tos
=======

This list of How-Tos demonstrates how the scripts in this git can be used in conjunction with other installed software to assemble nematode (and other small metazoan genomes)

How to sample sequences at random for test read sets
----------------------------------------------------

The example read files are between 5 and 8 GB each, and processing the complete data set could take a long time. If you just want to run all the How-Tos in this README on smaller files to test how the scripts work, use the `pick_random.pl` script:

    zcat g_ju800_110714HiSeq300_1.txt.gz | pick_random.pl 4 0.1 > g_ju800_110714HiSeq300_1.txt.r0.1

This command will gunzip the gzipped fastq file (`zcat file.gz` is the same as `gunzip -c file.gz`) and pipe the output through `pick_random.pl`. The script takes two arguments - the first tells you how many lines at a time should be treated as a block, while the second tells you what proportion of the file to pick (0.1 or 10% in this case)

However, running this command independently on the forward (`_1`) and reverse (`_2`) read files will result in different reads being picked. To sample both pairs at the same time, first interleave the files, and then pick 8 lines at a time (4 for the forward read, and 4 for the reverse read), and then use `pick_random.pl` to get 10% of the reads:

    paste \
      <(zcat g_ju800_110714HiSeq300_1.txt.gz | multi2single -l 4) \
      <(zcat g_ju800_110714HiSeq300_2.txt.gz | multi2single -l 4) |
    pick_random.pl 8 0.1 > g_ju800_110714HiSeq300_interleaved.txt.r0.1

Note: The <() syntax is for bash process substitution. Anything inside <(...) will be executed and its output will be treated by the containing command as if it is a file.

I could have also used a separate script for interleaving two fastq files (e.g., `shuffleSequences_fastx.pl 4 <(zcat g_ju800_110714HiSeq300_1.txt.gz) <(zcat g_ju800_110714HiSeq300_2.txt.gz`). I wrote `shuffleSequences_fastx.pl` based on shuffleSequences_fastq that used to ship with Velvet. The difference is that it can be used to shuffle both fasta (1st argument should be 2) and fastq (1st argument should be 4) sequences.

How to adapter- and quality-trim Illumina fastq reads using sickle and scythe in one command with no intermediate files
-----------------------------------------------------------------------------------------------------------------------

(Hmm, that's a bit of a long How-To title, almost like article titles in the 1600s, such as this one: [An Observation and Experiment Concerning a Mineral Balsom, Found in a Mine of Italy by Signior Marc-Antonio Castagna; Inserted in the 7th. Giornale Veneto de Letterati of June 22. 1671, and Thence English'd as Follows](http://dx.doi.org/10.1098/rstl.1671.0068))

How to create a preliminary assembly using ABySS
------------------------------------------------

How to create a preliminary assembly using CLC
----------------------------------------------

How to map reads to an assembly to get coverage information for each contig/scaffold
------------------------------------------------------------------------------------

How to make a plot of insert sizes for each library
---------------------------------------------------


How to search for the best blast hit for a contig
-------------------------------------------------

How to make a taxon-annotated GC cov "blob" plot
------------------------------------------------

