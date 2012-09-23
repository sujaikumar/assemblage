Annotation
==========

Although the workflows in this README were made specifically for my thesis on Next-generation nematode genomes, most of them are generic and should be useful to anyone annotating any other genome.

How to predict genes using a two-pass (iterative) MAKER2 workflow
--------------------------------------------------------

This is a meta How-to which uses the next few How-Tos as described:

Requirements:

1. Genome assembly (nucleotide fasta file)
2. CDSs (ESTs or RNA-Seq assembly) from the same species, if possible
3. Protein set from a closely related species, if possible (we used C. briggsae
   for C. sp5, B. malayi for D. immitis and L. sigmodontis, and M. hapla for 
   M. floridensis. You can also use a well curated set of proteins like swissprot ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
4. MAKER2 pipeline from http://www.yandell-lab.org/software/maker.html (we used version 2.25)
5. CEGMA - Core Eukaryotic Genes Mapping Approach - from http://korflab.ucdavis.edu/datasets/cegma/ (we used version 2.1)
6. SNAP gene finder from http://korflab.ucdavis.edu/Software . We used the version dated  2010-07-28 
7. Augustus gene finder from http://bioinf.uni-greifswald.de/augustus/ (we used version 2.5.5)
8. GeneMark-ES gene finder from http://exon.gatech.edu/license_download.cgi (we used GeneMark-ES LINUX 64 version 2.3e - `gm_es_bp_linux64_v2.3e.tar.gz`)

Overview:

1. Run CEGMA on genome. Convert CEGMA output to SNAP HMMs
2. Run GeneMark on genome. Get GeneMark HMMs
3. Run MAKER2 first pass - with HMMs (CEGMA SNAP, GeneMark) and evidence sets (ESTs, Proteins), with permissive settings
4. Convert MAKER2 results to MAKER SNAP HMMs
5. Convert MAKER2 results to AUGUSTUS HMMs
6. Run MAKER2 second pass - with HMMs (MAKER SNAP, GeneMark, AUGUSTUS) and evidence sets (ESTs, Proteins) and more stringent settings

The rest of this README uses bash variables, so just replace them with your own file names on the command line.

How to run CEGMA to generate a SNAP HMM
--------------------------------------------------------

You can run CEGMA in the same directory as the genome file, but I find it easier to run it in a genome specific subdirectory (in case you have many CEGMA runs to do in the same folder)

    # set the bash variable for the genome assembly file:
    genome=assembly.fna

    mkdir       $genome.cegma
    cd          $genome.cegma
    cegma -g ../$genome.cegma

CEGMA creates a bunch of output files:

`output.completeness_report` is very useful for seeing what percentage of core eukaryotic genes are found.

For annotation, we need the `output.cegma.gff` file. The following commands create a SNAP HMM trained using this CEGMA file:

    cd $genome.cegma
    cegma2zff output.cegma.gff   ../$genome
    fathom genome.ann genome.dna -categorize 1000
    fathom -export 1000 -plus uni.ann uni.dna
    forge export.ann export.dna
    hmm-assembler.pl $genome . > ../$genome.cegmasnap.hmm

How to run GeneMark
--------------------------------------------------------

The command for running GeneMark is:

    /PATH/TO/gm_es_bp_linux64_v2.3e/gmes/gm_es.pl $genome

This should work fine in most cases. It may not work if you have very short contigs because the default minimum contig size is 20,000. Try reducing that to 10,000. From the GeneMark-ES README:

> --min_contig [number]
>
> Short contigs, which are frequently found in draft assemblies,
> may introduce significant noise in self-training procedure.
> All contigs shorter then "min_contig" are excluded from training procedure.

GeneMark will create many folders in working directory. The final file that you need to provide to MAKER is in `mod/es.mod`.

How to run MAKER2 first pass and second pass
--------------------------------------------------------

MAKER2 needs a set of control (ctl) files before it can run:

    # generate maker ctl files:
    maker -CTL

Check the `maker_bopts.ctl` and `maker_exe.ctl` files to make sure they have the right defaults and the right application paths.

Next, edit `maker_opts.ctl` and edit the following lines:

    est=/PATH/TO/ESTset.fna #non-redundant set of assembled ESTs in fasta format (classic EST analysis)
    protein=PATH/TO/Proteinset.faa #protein sequence file in fasta format
    snaphmm=/PATH/TO/genome.cegmasnap.hmm #SNAP HMM file
    gmhmm=/PATH/TO/GeneMark.mod #GeneMark HMM file
    est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
    protein2genome=1 #gene prediction from protein homology, 1 = yes, 0 = no
    keep_preds=1 #Add unsupported gene prediction to final annotation set, 1 = yes, 0 = no
    single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no

Run MAKER2 first pass from the same working directory where `maker_*.ctl` files are located:

    maker -g $genome -base maker1

By default, maker will use the $genome file name without the extension as an output prefix. `-base maker1` sets the output prefix base name

Take MAKER2 GFF predictions, and convert to SNAP HMM

    # set output filename variable:
    base=maker1
    
    # maker creates an output folder with a datastore index. use gff3_merge to combine all the GFFs
    gff3_merge -d $base.maker.output/$base\_master_datastore_index.log
    
    # convert this GFF to a snap model:
    mkdir $base.hmm
    cd    $base.hmm
    maker2zff ../$base.gff
    fathom genome.ann genome.dna -categorize 1000
    fathom -export 1000 -plus uni.ann uni.dna
    forge export.ann export.dna
    hmm-assembler.pl $base . > ../$base.snap.hmm

Convert MAKER2 GFF predictions into Augustus HMM

    # prepare cdna and species names for augustus:
    cdna=/PATH/TO/ESTset.fna
    species=SpeciesName
    
    cd $base.hmm
    zff2gff3.pl genome.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' >genome.gff3
    # the perl one liner is needed to add an extra GFF column which is missing in the output of zff2gff3.pl
    
    # run augustus using autoAug:
    autoAug.pl --genome=../$genome --species=$species --cdna=../$cdna --trainingset=genome.gff3 --singleCPU -v --useexisting

autoAug.pl will create a new species specific HMM in its config directory with the name SpeciesName. This HMM and the maker1 snap HMM need to be used for the second-pass of MAKER2

Edit the `maker_exe.ctl` file again (make a backup of the previous one for a complete record):

    est=/PATH/TO/ESTset.fna #non-redundant set of assembled ESTs in fasta format (classic EST analysis)
    protein=PATH/TO/Proteinset.faa #protein sequence file in fasta format
    snaphmm=/PATH/TO/maker1.snap.hmm #SNAP HMM file created after MAKER2 first-pass
    gmhmm=/PATH/TO/GeneMark.mod #GeneMark HMM file
    est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
    protein2genome=0 #gene prediction from protein homology, 1 = yes, 0 = no
    augustus_species=SpeciesName #Augustus gene prediction species model
    pred_stats=1 #report AED and QI statistics for all predictions as well as models
    min_protein=30 #require at least this many amino acids in predicted proteins
    alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
    keep_preds=1 #Add unsupported gene prediction to final annotation set, 1 = yes, 0 = no
    split_hit=4000 #length for the splitting of hits (expected max intron size for evidence alignments)
    single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
    single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
    evaluate=1 #run EVALUATOR on all annotations (very experimental), 1 = yes, 0 = no

Run MAKER2 again:

    maker -g $genome -base maker2
    
Collate gff3 and fasta output:

    base=maker2
    gff3_merge  -d $base.maker.output/$base\_master_datastore_index.log
    fasta_merge -d $base.maker.output/$base\_master_datastore_index.log

How to extract AED < 1 genes only from MAKER2 output
--------------------------------------------------------

Requirements:

* `fastaqual_select.pl`
* `fgrep.pl`

Get a list of mRNA IDs with AED < 1 from MAKER2 GFF3 file:

    gff=maker2.gff3
    perl -lne 'print $1 if /\tmRNA\t.+ID=([^;]+).+_AED=(.+?);/ and $1<1' $gff >$gff.aed-1.0.ids

Extract all GFF lines from GFF3 file that contain these IDs. This will get all children of these features (CDSs, exons etc) as well because the Parent= attribute will have these mRNA IDs. This won't get the embedded genome fasta sequences, however:

    fgrep.pl -f $gff.aed-1.0.ids $gff -d "=|;|,|:" >$gff.aed-1.0.gff

fgrep.pl is faster than fgrep in most cases because it loads the patterns to be searched in memory. The search is also faster because rather than look for the pattern in every substring of the file to be searched, it only searches in complete columns. The `-d "=|;|,|:"` parameter specifies how columns are delimited, so, effectively - each mRNA ID is in its own column.

MAKER2 also outputs protein and transcript files. To select only those sequences with AED < 1:

    fasta=$base.proteins.fasta # or $base.transcripts.fasta
    fastaqual_select.pl -f $fasta -inc $gff.aed-1.0.ids > $base.proteins.aed=1.0.fasta

How to download and standardise nematode genomes and annotations
--------------------------------------------------------

Requirements:

* `fastaqual_select.pl`

For each species, the following files were downloaded (from WormBase WS230 - ftp://ftp.wormbase.org/pub/wormbase/releases/WS230/species/ ; http://nematod.es; and from http://www.inra.fr/meloidogyne_incognita/genomic_resources )

* genome softmasked files (repeats in lowercase),
* protein fasta files
* annotations in GFF2 or GFF3

Files were renamed using these short species names (the order is for convenience in generating tables, as it corresponds to the order of species in the nematode phylogenetic tree used in the thesis figures):

    # ts as di bm ls av sr bx mh mi mf pp ca cj ce cbn csp11 cr cbg csp5
    # annotation gff files: $speciesid.gff.gz 
    # protein fasta files:  $speciesid.faa.gz
    # genome fasta files:   $speciesid.gz
    # for genome fasta files, prepend species shortcode as prefix, so that we can identify contigs/scaffolds on their own
    # randomise order so that longest sequences aren't first (helps with splitting files for parallelization)
    zcat c_elegans.WS230.genomic_softmasked.fa.gz    | fastaqual_select.pl -f - -s r -prefix "ce."    | gzip >ce.gz
    zcat c_briggsae.WS230.genomic_softmasked.fa.gz   | fastaqual_select.pl -f - -s r -prefix "cbg."   | gzip >cbg.gz
    zcat c_sp5.WS230.genomic_softmasked.fa.gz        | fastaqual_select.pl -f - -s r -prefix "csp5."  | gzip >csp5.gz
    zcat c_remanei.WS230.genomic_softmasked.fa.gz    | fastaqual_select.pl -f - -s r -prefix "cr."    | gzip >cr.gz
    zcat c_brenneri.WS230.genomic_softmasked.fa.gz   | fastaqual_select.pl -f - -s r -prefix "cbn."   | gzip >cbn.gz
    zcat c_sp11.WS230.genomic_softmasked.fa.gz       | fastaqual_select.pl -f - -s r -prefix "csp11." | gzip >csp11.gz
    zcat c_japonica.WS230.genomic_softmasked.fa.gz   | fastaqual_select.pl -f - -s r -prefix "cj."    | gzip >cj.gz
    zcat c_angaria.WS230.genomic_softmasked.fa.gz    | fastaqual_select.pl -f - -s r -prefix "ca."    | gzip >ca.gz
    zcat p_pacificus.WS230.genomic_softmasked.fa.gz  | fastaqual_select.pl -f - -s r -prefix "pp."    | gzip >pp.gz
    zcat s_ratti.WS230.genomic_softmasked.fa.gz      | fastaqual_select.pl -f - -s r -prefix "sr."    | gzip >sr.gz
    zcat b.xylophilus.WS230.genomic_softmasked.fa.gz | fastaqual_select.pl -f - -s r -prefix "bx."    | gzip >bx.gz
    zcat m_hapla.WS230.genomic_softmasked.fa.gz      | fastaqual_select.pl -f - -s r -prefix "mh."    | gzip >mh.gz
    zcat a_suum.WS230.genomic_softmasked.fa.gz       | fastaqual_select.pl -f - -s r -prefix "as."    | gzip >as.gz
    zcat b_malayi.WS230.genomic_softmasked.fa.gz     | fastaqual_select.pl -f - -s r -prefix "bm."    | gzip >bm.gz
    zcat t_spiralis.WS230.genomic_softmasked.fa.gz   | fastaqual_select.pl -f - -s r -prefix "ts."    | gzip >ts.gz
    zcat nMf.1.1.fna.gz                              | fastaqual_select.pl -f - -s r -prefix "mf."    | gzip >mf.gz
    zcat nDi.2.2.fna.gz                              | fastaqual_select.pl -f - -s r -prefix "di."    | gzip >di.gz
    zcat nLs.2.1.fna.gz                              | fastaqual_select.pl -f - -s r -prefix "ls."    | gzip >ls.gz
    zcat nAv.1.0.fna.gz                              | fastaqual_select.pl -f - -s r -prefix "av."    | gzip >av.gz
    zcat MiV1A1_2995.scaf.ref.gz                     | fastaqual_select.pl -f - -s r -prefix "mi."    | gzip >mi.gz

For each species, create a subset of the GFF file with only those entries (mRNA/CDS/exon) that correspond to the protein prediction sets for that species:

    zgrep -P "GLEAN\texon" as.gff.gz       | gzip >as.coding.gff.gz
    zgrep -P "AUGUSTUS\tCDS" av.gff.gz     | gzip >av.coding.gff.gz
    zgrep -P "WormBase\texon" bm.gff.gz    | gzip >bm.coding.gff.gz
    zgrep -P "WormBase\texon" bx.gff.gz    | gzip >bx.coding.gff.gz
    zgrep -P "WormBase\texon" ca.gff.gz    | gzip >ca.coding.gff.gz
    zgrep -P "curated\texon" cbg.gff.gz    | gzip >cbg.coding.gff.gz
    zgrep -P "curated\texon" cbn.gff.gz    | gzip >cbn.coding.gff.gz
    zgrep -P "curated\texon" ce.gff.gz     | gzip >ce.coding.gff.gz
    zgrep -P "curated\texon" cj.gff.gz     | gzip >cj.coding.gff.gz
    zgrep -P "curated\texon" cr.gff.gz     | gzip >cr.coding.gff.gz
    zgrep -P "maker\texon" csp5.gff.gz     | gzip >csp5.coding.gff.gz
    zgrep -P "WormBase\texon" csp11.gff.gz | gzip >csp11.coding.gff.gz
    zgrep -P "maker\texon" di.gff.gz       | gzip >di.coding.gff.gz
    zgrep -P "maker\texon" ls.gff.gz       | gzip >ls.coding.gff.gz
    zgrep -P "maker\texon" mf.gff.gz       | gzip >mf.coding.gff.gz
    zgrep -P "Freeze3\tCDS" mh.gff.gz      | gzip >mh.coding.gff.gz
    zgrep -P "EuGene\texon" mi.gff.gz      | gzip >mi.coding.gff.gz
    zgrep -P "curated\texon" pp.gff.gz     | gzip >pp.coding.gff.gz
    zgrep -P "WormBase\texon" sr.gff.gz    | gzip >sr.coding.gff.gz
    zgrep -P "WormBase\texon" ts.gff.gz    | gzip >ts.coding.gff.gz

For each species, parse the transcript name from the coding.gff.gz file and add it as an extra attribute at the end of the GFF entry:

    zcat as.coding.gff.gz | perl -lne '/transcript:(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >as.coding.name.gff.gz
    zcat av.coding.gff.gz | perl -lne '/Parent=(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >av.coding.name.gff.gz    
    zcat bm.coding.gff.gz | perl -lne '/transcript:(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >bm.coding.name.gff.gz
    zcat bx.coding.gff.gz | perl -lne '/transcript:(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >bx.coding.name.gff.gz
    zcat ca.coding.gff.gz | perl -lne '/transcript:(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >ca.coding.name.gff.gz
    zcat cbg.coding.gff.gz | perl -lne '(/CDS \"(.+?)\"/ or /Transcript \"(.+?)"/) and  print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >cbg.coding.name.gff.gz
    zcat cbn.coding.gff.gz | perl -lne '(/CDS \"(.+?)\"/ or /Transcript \"(.+?)"/) and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >cbn.coding.name.gff.gz
    zcat ce.coding.gff.gz | perl -lne '(/CDS \"([^.]+?\.[^.]+?)\"/ or /Transcript \"([^.]+?\.[^.]+?)(\.\w+)?\"/) and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >ce.coding.name.gff.gz
    zcat cj.coding.gff.gz | perl -lne '(/CDS \"(.+?)\"/ or /Transcript \"(.+?)"/) and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >cj.coding.name.gff.gz
    zcat cr.coding.gff.gz | perl -lne '(/CDS \"(.+?)\"/ or /Transcript \"(.+?)"/) and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >cr.coding.name.gff.gz
    zcat csp5.coding.gff.gz | perl -lne '/Parent=([^,;=\b]+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >csp5.coding.name.gff.gz
    zcat csp11.coding.gff.gz | perl -lne '/transcript:(\S+?)(,\S+)?$/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >csp11.coding.name.gff.gz
    zcat di.coding.gff.gz | perl -lne '/Parent=([^,;=\b]+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >di.coding.name.gff.gz
    zcat ls.coding.gff.gz | perl -lne '/Parent=([^,;=\b]+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >ls.coding.name.gff.gz
    zcat mf.coding.gff.gz | perl -lne '/Parent=([^,;=\b]+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >mf.coding.name.gff.gz
    zcat mh.coding.gff.gz | perl -lne '/Parent=(\S+?)(;\S+)?$/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >mh.coding.name.gff.gz
    zcat mi.coding.gff.gz | perl -lne '/Parent=mRNA:(.+?);/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >mi.coding.name.gff.gz
    zcat mi.coding.gff.gz | perl -lne '/Parent=mRNA:(.+?);/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >mi.coding.name.gff.gz
    zcat pp.coding.gff.gz | perl -lne '(/CDS \"(.+?)\"/ or /Transcript \"(.+?)"/) and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >pp.coding.name.gff.gz
    zcat sr.coding.gff.gz | perl -lne '/transcript:(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >sr.coding.name.gff.gz
    zcat ts.coding.gff.gz | perl -lne '/transcript:(\S+)/ and print "$_;Coding_Name=$1" and next; print "$_;Coding_Name="' | gzip >ts.coding.name.gff.gz

How to create a Metazoa subset of NCBI's nr database for command line blast
--------------------------------------------------------

This is easy to do on the NCBI website. But for command line use, this is one way of getting a Metazoa (or any other taxon) subset of any NCBI blast database:

Download and unpack blast db from ftp://ftp.ncbi.nlm.nih.gov/blast/db

Get NCBI taxon id for taxon of interest from http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
e.g., Metazoa = 33208

Get all nuc or protein IDs for that taxon (db=nuccore for nucleotides or db=protein for proteins). Ensure qty=2000000 is more than the number you want to download

    curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=protein&dopt=gilist&qty=2000000&filter=all&term=txid33208\[Organism:exp\]" >txid33208.gids
    blastdb_aliastool -dbtype prot -gilist txid33208.gids -db nr -out nr_metazoa  
    
How to add functional annotations to a genome
--------------------------------------------------------

Requirements:

* Genome assembly nucleotide fasta file $species.gz
* Protein fasta file (from WormBase or MAKER2 or other source)
* nr_Metazoa blast database (see above)
* `fastaqual_select.pl`
* LAST aligner - http://last.cbrc.jp/ (we used version 193)
* maf-convert.py - TBA-Multiz package from http://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz (also used by README-CNEs.md)
* Blast+ suite - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+ (we used version 2.2.25+)
* InterProScan - ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/index.html (we used version 4.8)
* Blast2GO command line - http://www.blast2go.com/data/blast2go/b2g4pipe_v2.5.zip
* tRNAscan-SE - http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz (we used version 1.3)
* Rfam_scan - Database and executables on ftp://ftp.sanger.ac.uk/pub/databases/Rfam/10.1
* `merge_blast2goannot_gff.pl`

Use a bash variable $species (with the species short id) so that the following commands can become part of a pipeline:

### Blast2GO and Interproscan

For blast2go, we need

* the protein fasta file $species.faa
* interproscan xml run on protein fasta file
* blastp xml from protein fasta run against nr proteins (I chose nr_metazoa)

To get interproscan xml

    iprscan -cli -i $species.faa -o $species.iprscan.raw -format raw -goterms -iprlookup
        #convert raw to xml:
    converter.pl -format xml -input $species.iprscan.raw -output $species.iprscan.xml
    gzip $species.iprscan.xml $species.iprscan.raw 
    # we could have directly got xml output from iprscan, but the raw output is more flexible and might come in handy later

To get blastp xml against nr_metazoa

    sge_blast blastp -query $species.faa -db nr_Metazoa -evalue 1e-5 -max_target_seqs 50 -outfmt 5 \
        -out $species.nr_Metazoa.blastp.1e-5.max50
        
To run blast2go

    # Create a folder to store split iprscan files (blast2go wants one xml file per protein):
    mkdir ipr
    zcat $species.iprscan.xml.gz | perl -lne '
        /<protein id="(.+?)"/ and open XML, ">ipr/$1.xml" and print XML "<EBIInterProScanResults><interpro_matches>";
        print XML;
        /<\/protein/ and print XML "</interpro_matches></EBIInterProScanResults>"'

    export B2G4PIPEPATH=/PATH/TO/b2g4pipe_v2.5

    java -Xmx8000m -cp $B2G4PIPEPATH/ext/*:$B2G4PIPEPATH/* es.blast2go.prog.B2GAnnotPipe \
        -in $species.nr_Metazoa.blastp.1e-5.max50 -out $species.blast2go \
        -prop $B2G4PIPEPATH/b2gPipe.properties -annot -dat -ips ipr -v

    # creates two files: $species.blast2go.annot (tab file) and $species.blast2go.dat
    # (which can be loaded in blast2go GUI to generate other analyses)

To add blast2go annotations to MAKER2-generated GFF:

    merge_blast2goannot_gff.pl -b $species.blast2go.annot -g $species.gff >$species.blast2go.gff

This GFF can now be loaded in GBrowse and the Blast2GO annotations will be searchable

### tRNA

    # gunzip genome on fly, and create tRNA output
    tRNAscan-SE <(zcat $species.gz) >$species.trna
    
    # convert tRNAscan output to GFF
    perl <$species.trna -lne '
        ($contig, $trnanum, $st, $en, $type, $codon, $intron_st, $intron_en, $score)=split /\s+/;
        next unless $trnanum =~ /^\d+$/; # get rid of trnascan headers
        print "$contig\ttRNAscan-SE-1.3\ttRNA\t" . ($st < $en ? "$st\t$en\t$score\t+\t" : "$en\t$st\t$score\t-\t") . ".\tID=$contig\_tRNA$trnanum;Name=tRNA_$type\_$codon"
    ' >$species.trna.gff

### Rfam

    rfam_scan.pl -blastdb -o $species.Rfam.gff Rfam.fasta Rfam.cm $species

    #gff produced by rfam_scan has no name= and no note= fields in the 9th col, so gbrowse won't like it. Add those like so:
    perl -i -pne 's/id=(.+?);/id=$1;Name=$1;/;s/rfam-id=(.+)/rfam-id=$1;Note=$1/;' $species.Rfam.gff


### LSU-SSU

    wget ftp://ftp.arb-silva.de/release_108/Exports/SSURef_108_NR_tax_silva_trunc_v2.fasta.tgz
    wget ftp://ftp.arb-silva.de/release_108/Exports/LSURef_108_tax_silva_trunc.fasta.tgz

Convert U to T, remove spaces, make single line, select Nematoda:

    for a in SSURef_108_NR_tax_silva_trunc_v2.fasta LSURef_108_tax_silva_trunc.fasta
    do
      cat $a | perl -i -plne 'if(not /^>/){s/ //g; tr/U/T/}' | fastaqual_multiline_to_singleline.pl | grep -A1 Nematoda -  | grep -v "--" - | pigz >$a.Nematoda.gz
    done


    # create lastdb for each species
    parallel "lastdb -c {} <(zcat {}.gz)" ::: ts as di bm ls av sr bx mh mi mf pp ca cj ce cbn csp11 cr cbg csp5
    
ALIGN TO LASTDB for each species, get >200 bp matches, merge overlaps, convert to GFF

    parallel "lastal -e 200 {} <(zcat LSU*Nematoda.gz SSU*Nematoda.gz) | maf-convert.py psl - | awk '\$1>200' | cut -f14,16,17 | sort -k1,1 -k2,2n | mergeBed | awk '{print \$1\"\tARB-LSU-SSU\tnucleotide_match\t\"\$2+1\"\t\"\$3\"\t.\t+\t.\tARB-LSU-SSU\"}' >{}.fna.ARB-LSU-SSU.gff" ::: ts as di bm ls av sr bx mh mi mf pp ca cj ce cbn csp11 cr cbg csp5

tRNA, Rfam, and LSU-SSU GFF features can be added to a gbrowse instance.

How to convert interproscan annotations to a heatmap
--------------------------------------------------------

Collect all interproscan annotations for each protein, and turn them into a 2 colun tab-delimited table of the format:

    proteinid1<TAB>IPR:000001 Description1;IPR:000002 Description2;IPR:000003 Description3
    proteinid2<TAB>IPR:000002 Description2;IPR:000005 Description5

    for species in ts as di bm ls av sr bx mh mi mf pp ca cj ce cbn csp11 cr cbg csp5
    do
        zcat $species.faa.iprscan.xml.gz | 
        perl -ne '/protein id="(.+?)"/ and print "\n$1\t"; /interpro id="(IPR\d+)" name="(.+?)"/ and $1 ne "noIPR" and print "$1 $2;"' |
        perl -lne 'print unless /^\s*$/ or /^\"/' >$species.faa.iprscan.txt
    done

To visualise these $species.faa.iprscan.txt files as heatmaps per genome, first count all the IPR domains in each genome, then use R:

    cat {ts,as,di,bm,ls,av,sr,bx,mh,mi,mf,pp,ca,cj,ce,cbn,csp11,cr,cbg,csp5}.faa.iprscan.txt |
    perl -ne '
        /^(.+?)\./ and $sp = $1;
        while (/(IPR.*?);/g) { $ipr{$1}{$sp}++; $sphash{$sp} = 1 }
    }{
        print join("\t","ipr",keys %sphash) . "\n";
        for my $iprid (keys %ipr) {
            print "$iprid";
            for my $sp (keys %sphash) {
                print "\t" . (exists($ipr{$iprid}{$sp}) ? $ipr{$iprid}{$sp} : 0)
            }
            print "\n"
        }
    ' >all.faa.iprscan.txt.counts

Heatmap:

    # in R
    counts="all.faa.iprscan.txt.counts"
    d=as.matrix(read.delim(counts,header=TRUE,row.names=1,sep="\t"))
    library(gplots, grDevices)

    # add up the IPR domain counts across all species, pick the top 100 IPR domains by total species count:
    topn=100
    ds=head(d[rev(order(rowSums(d))),],n=topn)
    pdf(paste0(counts,".",topn,".pdf"),8,18)
    heatmap.2(ds,scale="none",dendrogram="row",trace="none",Rowv=TRUE,Colv=FALSE,col=gray(seq(1,0,-0.01)))
    dev.off()

How to convert tRNA annotations to a heatmap
--------------------------------------------------------

Count the total number of times each tRNA feature occurs in each genome and then use R heatmaps to visualise:

    cat *.trna.gff |
    perl -ne '
        $trna{$2}{$1}++ if /^(.+?)\..*Name=tRNA_(\S+)/;
    }{
        @species = qw(ts as di bm ls av sr bx mh mi mf pp ca cj ce cbn csp11 cr cbg csp5);
        print join("\t","tRNA",@species) . "\n";
        for my $t (sort keys %trna) {
            print $t;
            for $sp (@species) {
                print "\t" . ( exists $trna{$t}{$sp} ? $trna{$t}{$sp} : "0" );  
            }
            print "\n";
        }
    ' > all.trna.counts

Heatmap:

    # in R
    counts="all.trna.counts"
    ds=as.matrix(read.delim(counts,header=TRUE,row.names=1,sep="\t"))
    library(gplots, grDevices)
    pdf(paste0(counts,".pdf"),8,16)
    heatmap.2(d,scale="none",trace="none",Rowv="FALSE",Colv="FALSE",col=gray(50:0 / 50))
    dev.off()

