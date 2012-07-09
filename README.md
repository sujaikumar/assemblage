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

As an example, we can use the following study from the Short Read Archive: 

* ERP001495 - De novo whole-genome sequence of the free-living nematode Caenorhabditis sp. 5 strain JU800 DRD-2008

This study has one sample:

* ERS147916 - Caenorhabditis sp. 5 JU800 DRD-2008  

Two libraries (or "Experiments" in SRA terminology) were run for this sample:

* ERX114449 - 300 bp library with Illumina HiSeq2000 101 bp PE sequencing
* ERX114450 - 600 bp library with Illumina HiSeq2000 101 bp PE sequencing

And these are the four files that can be accessed from http://www.ebi.ac.uk/ena/data/view/ERP001495

* `g_ju800_110714HiSeq300_1.txt.gz` - 300 bp library forward read
* `g_ju800_110714HiSeq300_2.txt.gz` - 300 bp library reverse read
* `g_ju800_110714HiSeq600_1.txt.gz` - 600 bp library forward read
* `g_ju800_110714HiSeq600_2.txt.gz` - 600 bp library reverse read

How-Tos
=======

This list of How-Tos demonstrates how the scripts in this repository can be used in conjunction with other installed software to assemble nematode (and other small metazoan genomes).

How to make a taxon-annotated GC cov "blob" plot
------------------------------------------------

This section is a meta section that describes the How-Tos you need to make a taxon-annotated GC-cov blob plot like the one at http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3294205/figure/Fig3/

1. adapter- and quality-trim Illumina fastq reads using sickle and scythe in one command with no intermediate files
2. create a preliminary assembly using ABySS
3. map reads to an assembly to get insert-size and coverage information
4. assign high-level taxon IDs to contigs
5. make blobology plots with R

How to sub-sample sequences at random for test read sets
--------------------------------------------------------

The example read files are between 5 and 8 GB each, and processing the complete data set could take a long time. If you just want to run all the How-Tos in this README on smaller files to test how the scripts work, use the `pick_random.pl` script:

    zcat g_ju800_110714HiSeq300_1.txt.gz | pick_random.pl 4 0.1 > g_ju800_110714HiSeq300_1.txt.r0.1

This command will gunzip the gzipped fastq file (`zcat file.gz` is the same as `gunzip -c file.gz`) and pipe the output through `pick_random.pl`. The script takes two arguments - the first tells you how many lines at a time should be treated as a block, while the second tells you what proportion of the file to pick (0.1 or 10% in this case)

However, running this command independently on the forward (`_1`) and reverse (`_2`) read files will result in different reads being picked. To sample both pairs at the same time, first interleave the files, and then pick 8 lines at a time (4 for the forward read, and 4 for the reverse read), and then use `pick_random.pl` to get 10% of the reads:

    shuffleSequences_fastx.pl 4 <(zcat g_ju800_110714HiSeq300_1.txt.gz) <(zcat g_ju800_110714HiSeq300_2.txt.gz`) |
    pick_random.pl 8 0.1 > g_ju800_110714HiSeq300_interleaved.txt.r0.1

Note: The <() syntax is for bash process substitution. Anything inside <(...) will be executed and its output will be treated by the containing command as if it is a file.

`shuffleSequences_fastx.pl` is based on shuffleSequences_fastq.pl that used to ship with Velvet. The difference is that it can be used to shuffle both fasta (set 1st argument to 2) and fastq (set 1st argument to 4) sequences.

`g_ju800_110714HiSeq300_interleaved.txt.r0.1` is an interleaved/shuffled file. If you want the _1 and the _2 reads in separate files, you will have to do one more step:

    unshuffleSequences_fastx.pl 4 g_ju800_110714HiSeq300_interleaved.txt.r0.1

which will create `g_ju800_110714HiSeq300_interleaved.txt.r0.1.1` and `g_ju800_110714HiSeq300_interleaved.txt.r0.1.2`

How to adapter- and quality-trim Illumina fastq reads using sickle and scythe in one command with no intermediate files
-----------------------------------------------------------------------------------------------------------------------

(Hmm, that's a bit of a long How-To title, almost like article titles in the 1600s, such as this one: [An Observation and Experiment Concerning a Mineral Balsom, Found in a Mine of Italy by Signior Marc-Antonio Castagna; Inserted in the 7th. Giornale Veneto de Letterati of June 22. 1671, and Thence English'd as Follows](http://dx.doi.org/10.1098/rstl.1671.0068))

Requirements:

1. Install scythe for adapter trimming from https://github.com/ucdavis-bioinformatics/scythe 
2. Install sickle for quality trimming from https://github.com/ucdavis-bioinformatics/sickle
3. Install GNU Parallel (optional, but highly recommended!) from http://www.gnu.org/software/parallel/
4. A file with Illumina adapters - you can use `adapters.fa` from this repository

Command:

    sickle pe -t sanger -n -l 50 \
        -f <(scythe -a adapters.fa <(zcat g_ju800_110714HiSeq300_1.txt.gz) -q sanger \
            2> g_ju800_110714HiSeq300_1.scythe.err | perl -plne 's/^$/A/; s/ 1.*/\/1/')  \
        -r <(scythe -a adapters.fa <(zcat g_ju800_110714HiSeq300_2.txt.gz) -q sanger \
            2> g_ju800_110714HiSeq300_2.scythe.err | perl -plne 's/^$/A/; s/ 2.*/\/2/')  \
        -o >(gzip >g_ju800_110714HiSeq300_1.clean.txt.gz) \
        -p >(gzip >g_ju800_110714HiSeq300_2.clean.txt.gz) \
        -s >(gzip >g_ju800_110714HiSeq300_s.clean.txt.gz) \
            &>g_ju800_110714HiSeq300.sickle.err
        
The command above will first run scythe to search for adapters.fa in `g_ju800_110714HiSeq300_1.txt.gz` and `g_ju800_110714HiSeq300_2.txt.gz` simultaneously. The stderr stream of scythe is stored in `g_ju800_110714HiSeq300_1.scythe.err` and the output of scythe is parsed through a perl one liner that replaces blank lines with a single A, and adds "/1" or "/2" to the read header. This perl one liner is needed because scythe screws up the read names, and because the next script sickle can't deal with blank sequence lines where the whole sequence has been adapter-trimmed.

`sickle pe -t sanger -n -l 50` then runs sickle on the files output by the scythe and perl one liner above:

* `pe` means treat corresponding reads from the `-f` and `-r` files as being parts of a pair
* `-t sanger` tells it that the quality values are in sanger fastq encoding
* `-n` discards reads/pairs with any Ns in them
* `-l 50` discards reads/pairs with fewer than 50 bp
* `-o` filename of the output forward read
* `-p` filename of the output paired reverse read
* `-s` filename of the singletons (where the other read in the pair was discarded if it had an N or was shorter than `-l 50`

Additional scythe and sickle parameters can be set as needed (minimum number of matches to adapter sequence, trim bases from start of read, etc).

If you do have GNU Parallel installed, you can run multiple libraries at the same time like this:

    parallel "sickle pe -t sanger -n -l 50 \
        -f <(scythe -a adapters.fa <(zcat {}_1.txt.gz) -q sanger \
            2> {}_1.scythe.err | perl -plne 's/^$/A/; s/ 1.*/\/1/')  \
        -r <(scythe -a adapters.fa <(zcat {}_2.txt.gz) -q sanger \
            2> {}_2.scythe.err | perl -plne 's/^$/A/; s/ 2.*/\/2/')  \
        -o >(gzip >{}_1.clean.txt.gz) \
        -p >(gzip >{}_2.clean.txt.gz) \
        -s >(gzip >{}_s.clean.txt.gz) \
            &>{}.sickle.err" ::: g_ju800_110714HiSeq300 g_ju800_110714HiSeq600 

This command will create three files for each library:
 - `g_ju800_110714HiSeq300_1.clean.txt.gz`
 - `g_ju800_110714HiSeq300_2.clean.txt.gz`
 - `g_ju800_110714HiSeq300_s.clean.txt.gz`
 - `g_ju800_110714HiSeq600_1.clean.txt.gz`
 - `g_ju800_110714HiSeq600_2.clean.txt.gz`
 - `g_ju800_110714HiSeq600_s.clean.txt.gz`

You can interleave the _1 and _2 files:

    shuffleSequences_fastx.pl 4 \
        <(zcat g_ju800_110714HiSeq300_1.clean.txt.gz) \
        <(zcat g_ju800_110714HiSeq300_2.clean.txt.gz) |
    gzip > g_ju800_110714HiSeq300.clean.txt.gz

How to create a preliminary assembly using ABySS
------------------------------------------------

Requirements:

1. Install ABySS from http://www.bcgsc.ca/platform/bioinfo/software/abyss . This README was tested with versions 1.3.3 and 1.3.4. Versions before 1.3.3 did not work as well. It is possible that newer versions will have different options and work even better.

To get a preliminary assembly (PASS) with no pairing information (i.e. treating all reads as single-end) from the cleaned reads above:

    name=C20rp.31.n10.se;
    mkdir $name;
    abyss-pe -C $name name=$name n=10 k=31 \
        se='g_ju800_110714HiSeq300.clean.txt.gz.keep.abundfilt.repaired \
            g_ju800_110714HiSeq600.clean.txt.gz.keep.abundfilt.repaired' &>$name/log

Making a new directory with the parameters used makes it easy to track the settings used if we later create a lot of assemblies with different settings.  

* `-C` puts the output in the specified directory (`$name` in this case) rather than in the current directory
* `n=10` is the minimum number of connecting pairs needed before two unitigs or contigs are joined
* `k=31` is the k-mer used for making the de Bruijn graph. 31 is a reasonable default. Longer k-mers may give more contiguous assemblies but the goal of the PASS is to get as much assembled sequence as possible rather than get the longest contigs.

The single-end PASS will be in `$name/$name-contigs.fa`

How to map reads to an assembly to get insert-size and coverage information
---------------------------------------------------------------------------

1. Install bowtie2 from http://bowtie-bio.sourceforge.net/bowtie2 . This README was tested with Version 2.0.0-beta5 and Version 2.0.0-beta6 but will probably work with newer versions as well if the syntax has not changed.
2. Install samtools from http://samtools.sourceforge.net
3. `sam_len_cov_gc_insert.pl` from this repository
4. Preliminary assembly - a fasta file with contigs `filename-contigs.fa`.
5. Interleaved read files for each library. They do not have to be interleaved to calculate coverage for each contig, but if they are interleaved, then generating insert-size plots (next How-To) is easier. In this example, I assume that we want to map the trimmed 300 bp and the 600 bp interleaved files:
 - `g_ju800_110714HiSeq300.clean.txt.gz`
 - `g_ju800_110714HiSeq600.clean.txt.gz`

Create a bowtie2 index

    assemblyfile=filename-contigs.fa
    bowtie2-build $assemblyfile $assemblyfile

Map reads for multiple libraries using a bash for loop

    for lib in g_ju800_110714HiSeq300.clean.txt.gz g_ju800_110714HiSeq600.clean.txt.gz
    do
        bowtie2 -x $assemblyfile -f -U $lib --very-fast-local -p 8 --reorder --mm | 
        tee >(samtools view -S -b - > $assemblyfile.$lib.bowtie2.bam) |
        tee >(samtools view -S -b - | samtools sort -m 2000000000 - $assemblyfile.$lib.bowtie2.sorted) |
        sam_len_cov_gc_insert.pl -i -f $assemblyfile -s - -out $assemblyfile.$lib
    end

The two `tee` commands take the SAM output generated by bowtie2 and write them out as read-sorted `$assemblyfile.$lib.bowtie2.bam` and contig-sorted `$assemblyfile.$lib.bowtie2.sorted.bam`.

`sam_len_cov_gc_insert.pl` takes the following options:

* `-i` use this switch if your read file is interleaved and you want to estimate insert sizes. Leave it out if you used bowtie2 without an interleaved file
* `-f` assembly fasta file
* `-s` alignments in SAM format
* `-out` output file prefix

and creates the following files:

* `outputprefix.lencovgc.txt` - tab delimited text file with these columns:
 1. library name (read filename by default)
 2. contig id
 3. contig length
 4. contig coverage (read coverage, not k-mer coverage)
 5. contig gc content (proportion from 0 to 1)
 
 This file is used to make taxon-annotated GC-coverage "blob" plots
 
* `outputprefix.lencovgc.fna` - a fasta file (with each sequence in a single line) where the contig id has a useful suffix so the fasta header looks like this `>contigid_len_cov_gc`, e.g., `>contig23_1221_67.623_0.3222349`. We found this file useful for quickly pulling out subsets of contigs that met a certain GC and coverage criteria, as described in **How to select preliminary assembly (PASS) contigs with given GC or coverage**.

If you used interleaved read files and the `-i` switch in `sam_len_cov_gc_insert.pl`, then you also get two additional files:

* `outputprefix.pairing.hist.txt` - Frequency distribution of each insert size and whether the reads were pointed at each other ("innies") or away from each other ("outies", as you would expect from Illumina's mate-pair protocol)
* `outputprefix.pairing.stat.txt` - Some statistics about what percentage of reads were mapped to the same contig and how many were innies or outies. etc.

How to make a plot of insert sizes for each library
---------------------------------------------------

If you used interleaved reads when mapping the read files to the preliminary assembly using bowtie2 (in the previous How-To), then you can use the following R script to create a pretty plot of the insert sizes.

Requirements:

1. R installed (Rscript in your path)
2. R packages ggplot2 and grid installed
3. `plot_insert_freq_txt_binned.R` from this repository
4. `outputprefix.pairing.hist.txt` from the output of `sam_len_cov_gc_insert.pl` with `-i` option (previous How-To)

    plot_insert_freq_txt_binned.R outputprefix.pairing.hist.txt

How to assign high-level taxon IDs to contigs
---------------------------------------------

Requirements

1. NCBI blast+ suite installed from ftp://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST
2. NCBI nt blast db downloaded and unpacked from ftp://ftp.ncbi.nih.gov/blast/db/
 * `wget ftp://ftp.ncbi.nih.gov/blast/db/nt.??.tar.gz`
3. NCBI taxonomy dumps downloaded and unpacked:
 * `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz`
 * `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
4. `fastaqual_select.pl` from this repository
5. `blast_taxonomy_report_topscore.pl` from this repository
 
I'm using bash variables to store filenames so that the command syntax is more generic:

    assemblyfile=filename-contigs.fa
    
Select some contigs at random from the preliminary assembly:

    fastaqual_select.pl -f $assemblyfile -s r -n 10000 >$assemblyfile.r10000

* `-f` is the name of the fasta file
* `-s r` sorts at random
* `-n 10000` picks 10000 contigs. you can pick more or less depending on how big your preliminary assembly is

Blast this file against the NCBI nt database. A megablast search should be quick, and enough to estimate the approximate taxonomic composition of the preliminary assembly:

    blastn -task megablast -query $assemblyfile.r10000 -db /path/to/blastdb/nt -max_target_seqs 1 -outfmt 6 > $assemblyfile.r10000.megablast.nt

Assign taxons to the contigs at the level of genus, order, family, superfamily, and kingdom (we find that "order" works well for visualising composition, but are keeping the rest just in case):

    blast_taxonomy_report.pl \
        -b $assemblyfile.r10000.megablast.nt \
        -nodes /path/to/ncbi/taxdmp/nodes.dmp \
        -names /path/to/ncbi/taxdmp/names.dmp \
        -gi_taxid_file /path/to/ncbi/taxdmp/gi_taxid_nucl.dmp.gz \
        -t genus=1 -t order=1 -t family=1 -t superfamily=1 -t kingdom=1 \
    >$assemblyfile.r10000.megablast.nt.taxon

How to make blobology plots with R
----------------------------------

Requirements:

1. `blobology.R` from this repository
2. `$assemblyfile.r10000.megablast.nt.taxon` - a file with taxons assigned for a set of contigs from the preliminary assembly 
3. `outputprefix.lencovgc.txt` - tab delimited file with len cov gc information obtained previously using `sam_len_cov_gc_insert.pl` script above

Take the file with taxon assignments and use it to append a taxon column to the tab-delimited file with len cov gc information. Where a taxon identification is not available, use "Not annotated":

    for file in *lencovgc.txt
    do
        cat $assemblyfile.r10000.megablast.nt.taxon $file | 
        perl -anF"\t" -e '
            chomp;
            /^(\S+).*\torder\t([^\t]+)/ and $o{$1} = $2 and next;
            if ($F[2] =~ /^\d+$/ ) { print "$_\t" . (exists $o{$F[1]} ? $o{$F[1]} : "Not annotated") . "\n" }
        ' >> lencovgc.taxon
    done    

Use ggplot2 in R to make a plot showing the taxon-annotated gc-cov blob plots for each library:

    Rscript --vanilla blobology.R lencovgc.taxon 0.005 1 2
    
The three arguments to blobology.R are:

1. `lencovgc.taxon` - the tab delimited file with columns: libraryname, contigid, length, coverage, gc, taxon 
2. `0.005` - a minimum presence threshold. A taxon must be present in at least this proportion of annotated contigs before it is reported. Without this cutoff, even a dozen random hits can quickly swamp the plot legend and make it impossible to see which taxa are present.
3. `1 2` - this last pair of numbers is purely aesthetic - if there are multiple libraries being plotted, this pair specifies the number of rows and columns in which the plots are laid out. In our example, we have only two libraries, so `1 2` or `2 1` etc would have worked fine. With only 1 library, you will need `1 1`

How to make taxon-specific blast databases
------------------------------------------

Once you've made a GC-cov blob plot, you might want to make taxon specific databases to do a more sensitive search against.

For example, the contaminants in figure http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3294205/figure/Fig3/ were all Proteobacteria (NCBI taxonomy ID 1224). On the NCBI website, we can easily restrict the target database to a particular taxon. However, for large-scale searches on a local computer/cluster, we will need a subset of the NCBI database specific to this taxon.

Requirements:

1. Taxonomy ID from NCBI for the subset you want to create: e.g., 1224 for Proteobacteria
2. Local copy of NCBI blast database (eg: nt or nr)
3. Install BLAST+ suite from ftp://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST

`curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=nuccore&dopt=gilist&qty=2000000&filter=all&term=txid1224\[Organism:exp\]" >txid1224.gids`
`blastdb_aliastool -dbtype nucl -gilist txid1224.gids -db /path/to/blastdb/nt -out nt_Proteobacteria`

Remember:

1. get the db=nuccore part right (change to db=protein for nr)
2. qty=NUM should be greater than total number of sequences

How to select preliminary assembly (PASS) contigs with given GC or coverage
---------------------------------------------------------------------------

How to separate contigs based on which taxon-specific blast database they hit better
------------------------------------------------------------------------------------

How to run khmer to digitally normalise reads
---------------------------------------------

Requirements:

1. Install screed from https://github.com/ctb/screed
2. Install khmer  from https://github.com/ctb/khmer - this tutorial site might help: http://ged.msu.edu/angus/diginorm-2012/tutorial.html

Typically, you might run *khmer* with a high-pass and a low-pass filter. First remove high coverage reads that don't contribute new k-mers beyond coverage 20, and then remove low coverage reads that have infrequent k-mers (e.g, if coverage <5 is likely to be sequencing errors or non-abundant contaminants). In some cases, only the high-pass filter is used.

khmer can be run on _1 forward and _2 reverse files separately, but we recommend that the files are interleaved before running, otherwise you might end up with almost all _2 reads discarded if all k-mers are found in the _1 file.

    shuffleSequences_fastx.pl 4 <(zcat g_ju800_110714HiSeq300_1.clean.txt.gz) <(zcat g_ju800_110714HiSeq300_2.clean.txt.gz`) | gzip >g_ju800_110714HiSeq300.clean.txt.gz

In this case, I am using the scythed and sickled cleaned read fastq files generated above.

Now, run the high-pass filter first:

    normalize-by-median.py g_ju800_110714HiSeq300.clean.txt.gz -C 20 -k 25 -N 4 -x 2e8 -s interleaved.C20.kh -R interleaved.C20.report

* The command above will take the gzipped cleaned interleaved fastq file specified and normalize-by-median all reads to a max-coverage of `-C 20`
* `-k 25` specifies that k-mers of size 25 are used for normalization (this is a reasonable default)
* `-N 4 -x 2e8` specifies how many blocks are used and how big each is (this settig will take up 4 x 2e8 = 3.2 GB of memory)
* `-s` specifies the name of the file where the k-mer counts are stored and `-R` specifies a report file

The output of this command will be a .keep file: `g_ju800_110714HiSeq300.clean.txt.gz.keep` which is a fasta file (khmer doesn't output fastq quality files).

To run the low-pass filter (in this case, remove reads with coverage less than `-C` 5X):

    filter-abund.py -C 5 interleaved.C20.kh g_ju800_110714HiSeq300.clean.txt.gz.keep

The command above takes the k-mer count file `interleaved.C20.kh` and whenever it sees reads in `g_ju800_110714HiSeq300.clean.txt.gz.keep` that have very infrequent k-mer coverage, it discards them. The final output file of this step is named `g_ju800_110714HiSeq300.clean.txt.gz.keep.abundfilt`.

One important thing to remember is that khmer works on individual reads and not on pairs. As a result, the pairing information of a read is lost if its pair is discarded. To get around this problem, we use the script `khmer_re_pair.pl` which pulls in the paired reads:

    khmer_re_pair.pl -o g_ju800_110714HiSeq300.clean.txt.gz -k g_ju800_110714HiSeq300.clean.txt.gz.keep.abundfilt > g_ju800_110714HiSeq300.clean.txt.gz.keep.abundfilt.repaired

Currently, `khmer_re_pair.pl` needs the complete file to be an interleaved fasta file (can be gzipped). For each read in the khmer-ed fasta file `g_ju800_110714HiSeq300.clean.txt.gz.keep.abundfilt`, it looks up the original interleaved fasta file and outputs both reads in the pair.

Note: all of this has to be done for the other library (600 bp) as well. i.e. we should be left with khmer-ed and repaired files:

* `g_ju800_110714HiSeq300.clean.txt.gz.keep.abundfilt.repaired`
* `g_ju800_110714HiSeq600.clean.txt.gz.keep.abundfilt.repaired`


How to find mate-paired reads that don't map to the ends of contigs and use those to conservatively scaffold the contigs
------------------------------------------------------------------------------------------------------------------------

How to rename and reorder scaffolds post-assembly
-------------------------------------------------



