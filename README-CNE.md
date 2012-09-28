Conserved Noncoding Elements
==========

Although the scripts referred to in this README were made specifically for my thesis on Next-generation nematode genomes, most of the multiple alignment format (MAF) file manipulation tools are generic and should be useful to anyone annotating any other genome.

Identify CNEs using whole-genome-alignments
--------------------------------------------

To get conserved non-coding elements (CNEs) across 20 nematode species, the basic steps are:

1. Get genomes
2. Perform pair-wise whole-genome alignments
3. Use a multiple whole-genome aligner to convert pair-wise alignments to multiple alignments (i.e. get conserved regions)
4. Mask coding regions in multiple-alignment files
5. Remove masked regions in multiple-alignment files (i.e. remove noncoding regions)
6. Report the remaining multiple alignments


Here are the details for each of these steps:

### 1. Get genomes

Requirements:

* one genome file for each species
* genome file name same as species name (I used species acronyms: ce, cbr, cbg, etc)
* each sequence in the genome file prefixed with "species.". (e.g, ">ce.chrI", "ce.chrII", etc)

These files are available at http://dx.doi.org/10.6084/m9.figshare.96089

### 2. Perform pair-wise whole-genome alignments

Requirements:

* LAST whole genome aligner (or BLASTZ/LASTZ - anything that outputs MAF files). I tried a few, and LAST seemed the fastest/most sensitive (I used version 193 from http://last.cbrc.jp)
* optional: gnu parallel to speed things up
* Multiz-TBA from http://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz (we used release 2009 Jan 21)
* list_pairwise.pl from this repository (takes arguments, returns tab separated lines with pairs of arguments)
* maf_insert_lc.pl - replaces substrings 
* maf_remove_lc.pl
* maf_select.pl
* species.coding.gff files

Make LAST databases for searching against

	parallel "lastdb -c {} {}" :::  ce cbg csp5 cr cbn csp11 cj ca pp sr bx mh mi mf as bm ls av di ts
	# -c masks lower case bases (some of the genome files have repeat regions lower case masked)

Run lastal, add maf header (LAST doesn't put it, and single_cov2 needs it), and run single_cov2

    list_pairwise.pl  ce cbg csp5 cr cbn csp11 cj ca pp sr bx mh mi mf as bm ls av di ts |
    awk '$1 != $2' | 
    perl -lane 'print "lastal -f 1 $F[0] $F[1] | tee >(gzip >$F[0].$F[1].maf.gz) | cat <(echo "##maf version=1") - | single_cov2 /dev/stdin >$F[0].$F[1].sing.maf " unless -f "$F[0].$F[1].maf.gz"' |
    parallel
    
    # list_pairwise.pl returns pairs of species from a list where the order of species is maintained:
    # i.e. 1-2 1-3 1-4 ... 2-3 2-4 ... 3-4, but never 2-1, 3-2, etc
    # lastal does the pairwise alignment
    # single_cov2 returns only those alignments that are single coverage (Multiz-TBA needs this, it can't work with multiple hits)
    # the final output is a .sing.maf file for each pair of species
    
### 3. Use a multiple whole-genome aligner to convert pair-wise alignments to multiple alignments

Multiz-TBA (Threaded Blockset Aligner) can be run in two ways, using "tba +" or "tba -". The former runs all the component commands right away. The latter outputs the list of commands that would be run. I prefer the latter because it is more flexible (e.g., you can replace /tmp with a different directory if you don't have much space in /tmp)

The command needs a guide tree of the species being aligned. And it needs those species names as the prefixes of each sequence (this was why we had to put species prefixes on each sequence in the original genome fasta files)

The closest neighbours are aligned first, followed by the next out group, followed by the next, etc. The guide tree is in newick format (but with commas replaced by spaces). The guide tree always uses the files it encounters first as the reference, so it is better to put the better assembled and better annotated genomes earlier in the list (e.g, ce before cr, for example), and to rearrange the guide tree so that the more contiguous genomes are more to the left of the list.

This is how I ran tba for the five species in Clade IV: M. hapla (mh), M. incognita (mi), M. floridensis (mf), B. xylophilus (bx), and S. ratti (sr):

    # set bash variable to make naming files easier
    tba=tba.cladeIV
    tba - "((sr bx) (mh (mi mf)))" *.sing.maf $tba.maf > $tba.bash
    # tba picks up the pairwise .sing.maf files that it needs if they are in the same folder
    
    bash $tba.bash &> $tba.log
    gzip $tba.maf
    # some of the files can be quite large, so gzipping them right away saves some space, but is not necessary

### 4. Mask coding regions in multiple-alignment files, and 5. remove masked regions

To mask coding regions, we have to use the species.coding.gff files (download from http://dx.doi.org/10.6084/m9.figshare.96089)

Use the script interval_mask.pl to create versions of the genome that have lower case bases for all the coding features:

    parallel "interval_mask.pl -ci U -co L -f {} -i <(zcat {}.coding.gff.gz | cut -f1,4,5) | gzip >{}.lc.gz" ::: ce cbg csp5 cr cbn csp11 cj ca pp sr bx mh mi mf as bm ls av di ts    
    # interval_mask.pl takes genome intervals (in three col format: chromosome/contig start end) as -i input, and masks the -f fasta file, all input bases are converted to uppercase (-U) and then all interval bases are masked as lower case. The output is the interval masked fasta file
    
Use maf_insert.pl to replace sequences in the alignment files with sequences from the codingmasked fasta files, and use maf_remove_lc.pl to remove those blocks from the alignment.

    maf_insert_lc.pl -m $tba.maf.gz -f <(zcat   {sr,bx,mh,mi,mf}.lc.gz) | maf_remove_lc.pl | gzip >$tba.lc.remove.maf.gz

### 6. Report the remaining multiple alignments

Sometimes multiz-tba returns alignments in which only 2 or 3 or 4 of the 5 species are present. To ensure that we only select those alignments where all species are present, and where the alignment is a certain minimum length or minimum identity, use maf_select.pl:

    maf_select.pl    -m $tba.lc.remove.maf.gz -s sr bx mh mi mf -l 30 -r 0.5 -i 0 2>$tba.l30.i0.r50.stats.txt | gzip >$tba.l30.i0.r50.maf.gz

* -m takes a maf file as input
* -s takes a list of species prefixes and only selects those alignment blocks that have all species present
* -l 30  selects blocks where output alignment has minimum 30 columns (this is an arbitrary cutoff)
* -r 0.5 selects blocks with a minimum relative identity of 50%
* -i 0   (optional) selects blocks with a minimum absolute identity of 0 (i.e. no filtering in this case)
* 2> redirects stderr to a new file which stores a 3 column file with: length, absolute identity, and relative identity (absolute identity is calculated as the number of alignment columns which are absolutely identical across all rows, divided by total number of columns; relative id is the number of columns where >=50% rows agree on a base, divided by the number of cols)

Repeat the same for other nodes in the phylogeny of the 20 species analysed, using these guide trees:

* All 20 nematode genomes (tba.all20): `(((((((ce (((cbg csp5) cr) (cbn csp11))) cj) ca) pp) ((sr bx) (mh (mi mf)))) (as ((bm (ls av)) di))) ts)`
* Clade IV and Clade V 20 (tba.cladeIVV): `(((((ce (((cbg csp5) cr) (cbn csp11))) cj) ca) pp) ((sr bx) (mh (mi mf))))`
* Clade V (tba.cladeV): `((((ce (((cbg csp5) cr) (cbn csp11))) cj) ca) pp)`
* Caenorhabditis (tba.caenorhabditis): `(((ce (((cbg csp5) cr) (cbn csp11))) cj) ca)`
* Elegans group (tba.elegans): `(ce (((cbg csp5) cr) (cbn csp11)))`
* Clade IV (tba.cladeIV): `((sr bx) (mh (mi mf)))`
* Meloidogyne (tba.meloidogyne): `(mh (mi mf))`
* Clade III (tba.cladeIII): `(as ((bm (ls av)) di))`
* Onchocercidae (tba.onchocercidae): `((bm (ls av)) di)`
    
Identify CNEs using MegaBLAST and clustering
--------------------------------------------

Requirements

1. genome files
2. megablast (from NCBI blast suite, I used version 2.2.25)
3. link_blast.pl
4. BEDtools (from http://code.google.com/p/bedtools/ - I used version 2.16.2)
5. map_data_ids.pl

To find all CNEs for any node on the nematode phylogenetic tree, first identify the species you want to use as the main reference alignment. For example, to find CNEs in all 8 *Caenorhabditis* species using ce as the reference, do pairwise megablasts between ce-cbg, ce-cbn, ce-cj, ce-cr, ce-ca, ce-csp5, ce-csp11

    # example (using the same megablast settings as Vavouri et al. 2007):
    megablast -W 30 -e 0.001 -i cbg -d ce -U T -m 8 -o cbg.ce.megablast.W30.e0.001.UT.m8

Once all the megablasts are done, cluster the hits together using single-linkage clustering. eg, if segmentA hits segmentB, and segmentB hits segmentC, then cluster all 3 together. The hits/links are fairly strict as they were formed using megablast with a word seed of 30.

    # cluster blast results (link_blast = blastclust, effectively)
    cat ce.{ca,cbg,cbn,cj,cr,csp11,csp5}.megablast.W30.e0.001.UT.m8 | link_blast.pl >caeno.megablast.cluster.out 2>caeno.megablast.cluster.err

link_blast.pl is better than blastclust for this application in two ways:

1. It can deal with many more sequences (blastclust starts to fail around ~10,000 >100 bp sequences)
2. It doesn't cluster the whole input sequence (e.g., all of chrI in C. elegans!). Instead it extracts the coordinates of each blast hit as intervals, merges overlapping intervals, and clusters the merged interval sequences. So, it is possible for one part of a sequence to go into one cluster, and another part to go into another cluster.

By default, the output of link_blast.pl is a file with clusters, one cluster per line. Each cluster is made up of sequence intervals. eg:

    cbn.Cbre_Contig33_514645_514842 cbg.chrV_11811066_11811127 csp11.Scaffold491_35449_35636 ce.CHROMOSOME_V_6022966_6023233
    
The next step is to convert these cluster intervals into a 3 column TAB separated BED format that bedtools can use

    # convert cluster intervals (_ separated) into tab sep for bedtools
    perl -lne 'while ( /\b(\S+)_(\d+)_(\d+)\b/g ) {print "$1\t" . ($2 - 1) . "\t$3"}' caeno.megablast.cluster.out >caeno.megablast.cluster.intervals
    
If any of these intervals overlap coding regions, then they must be removed. So, first find out which intervals those are:

    bedtools intersect -wa \
        -a caeno.megablast.cluster.intervals \
        -b <(zcat {ce,ca,cbg,cbn,cj,cr,csp11,csp5}.coding.name.gff.gz) |
    sort | uniq >caeno.megablast.cluster.intervals.coding
    
Then, (this is the sneaky bit), to remove these intervals from the clustering file, I created "map" file, (a two col tab delimited file stores maps from coding sequences to a " " (space) character. Eg:

    # format:  coding_interval<TAB>character_to_map_to
    ca.Can_chr102829_6197_6294<TAB><SPACE>
    
The map_data_ids.pl script then takes the cluster file and this map file, and transforms all cluster members that overlap coding intervals into spaces:

    perl -lane 'print "$F[0]_" . ($F[1]+1) . "_$F[2]\t "' caeno.megablast.cluster.intervals.coding |
    map_data_ids.pl -d " " -m - -i caeno.megablast.cluster.out | sed 's/ +/ /g'> caeno.megablast.cluster.out.cne
    
That last `sed 's/ +/ /g'` converts all multiple spaces into a single space. Effectively, the coding sequence intervals have been removed from the megablast clusters, and the only clusters left are clusters of non-coding elements

This next bit keeps only the clusters that have members from all species at that node in the nematode phylogeny:

    awk '/ca\./ && /cbg\./ && /cbn\./ && /ce\./ && /cj\./ && /cr\./ && /csp11\./ && /csp5\./ ' caeno.megablast.cluster.out.cne >caeno.megablast.cluster.out.cne.all

The "strict" version of the megablast-based CNE finding workflow removes a cluster even if only one sequence interval in that cluster overlaps a coding region.

    perl -lane 'print "$F[0]_" . ($F[1]+1) . "_$F[2]"' caeno.megablast.cluster.intervals.coding |
    fgrep.pl -v -f - -d " " caeno.megablast.cluster.out >caeno.megablast.cluster.out.cne.strict

The perl one liner converts 3 col (tab separated) list of coding CNEs into a "_" separated string, and fgrep.pl -v removes all lines which have those strings in them.

fgrep.pl is like fgrep, but faster (although more memory intensive as it loads the pattern file to be searched in memory) and only works on fields that are delimited by -d (in this case: " ").

Last step - pick only those CNE clusters where all species are present:

    awk '/ca\./ && /cbg\./ && /cbn\./ && /ce\./ && /cj\./ && /cr\./ && /csp11\./ && /csp5\./ ' caeno.megablast.cluster.out.cne.strict >caeno.megablast.cluster.out.cne.strict.all

Repeat this for other groups of species as well:

* All 20 nematode genomes: `ce cbg csp5 cr cbn csp11 cj ca pp sr bx mh mi mf as bm ls av di ts`
* Clade IV and Clade V 20: `ce cbg csp5 cr cbn csp11 cj ca pp sr bx mh mi mf`
* Clade V: `ce cbg csp5 cr cbn csp11 cj ca pp`
* Caenorhabditis: `ce cbg csp5 cr cbn csp11 cj ca`
* Elegans group: `ce cbg csp5 cr cbn csp11`
* Clade IV: `sr bx mh mi mf`
* Meloidogyne: `mh mi mf`
* Clade III: `as bm ls av di`
* Onchocercidae: `bm ls av di`
