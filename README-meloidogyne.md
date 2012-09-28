Meloidogyne floridensis comparative genomics
==========

Although the scripts referred to in this README were made specifically for my thesis on Next-generation nematode genomes, some of the workflows (e.g., using Exonerate to get protein and CDS fasta files) may be useful to others.

Extract protein coding CDSs from M. floridensis using exonerate
------------------------------------------------------------

Requirements:

* Genome assembly (nucleotide fasta file)
* Proteins from closely related species (protein fasta file)
* Exonerate (from http://www.ebi.ac.uk/~guy/exonerate/ - I used version 2.2.0)
* BEDtools (from http://code.google.com/p/bedtools/ - I used version 2.13.3 in November 2011)


Overview:

1. Align proteins from previously sequenced species (M. incognita and M. hapla) to new genome (M. floridensis) using exonerate:protein2genome model
2. Convert exonerate output format to CDS and amino acid fasta files

Details:

### 1. Align proteins to genome

Run exonerate:

    # bash variables
    protein=mh.faa
    genome=mf.fna
    exonerate --model protein2genome --percent 50 -q $protein -t $genome --showvulgar F --showalignment T --ryo '%qi,%ql,%qab,%qae,%ti,%tl,%tab,%tae,%et,%ei,%es,%em,%r,%pi,%ps,%C\n' > $protein.$genome.p2g

exonerate is run with --percent 50 to report only alignments where more than 50% of the query protein sequence is covered.
The rest of the options are for output formatting: --ryo stands for "roll your own". These are the columns shown:

* 1-4: query - id, length, alignment start, alignment end
* 5-8: target - id, length, alignment start, alignment end
* 9-12: equivalenced (the bits that align, not including introns) - total, identical, similar, mismatches (total - identical = mismatches)
* 13-16: rank of match, percent identical, percent similar, CIGAR block (like BAM, but with spaces)

Here is an example of one protein aligned with exonerate. The line at the bottom is the --ryo line:

    C4 Alignment:
    ------------
             Query: MhA1_Contig0.frz3.gene9
            Target: NODE_769_length_6021_cov_89.450089_6075_200.51_0.31 [revcomp]
             Model: protein2genome:local
         Raw score: 291
       Query range: 0 -> 94
      Target range: 524 -> 200
    
       1 : MetAlaArgIleGluLeuPheAlaHisLysSerLysArgValAsnThrIleLeuAsnGlnLeuS :  22
           ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
           MetAlaArgIleGluLeuPheAlaHisLysSerLysArgValAsnThrIleLeuAsnGlnLeuS
     524 : ATGGCACGAATTGAACTTTTTGCTCATAAATCCAAACGAGTAAATACTATCTTGAATCAGCTTA : 461
    
      23 : erSerAspPheGlnProGlyLeuThrProAsnSer-IleThrThrCysPheSerIleTyrGlnL :  43
           |||||||||||||||||||||||! !|||!:!|||#!!:! !||||||! !! !||||||||||
           erSerAspPheGlnProGlyLeuIleProSerSer#MetIleThrCysSerPheIleTyrGlnL
     460 : GCTCGGATTTCCAACCAGGTTTGATCCCTAGTTCTAATGATTACTTGCTCCTTTATTTATCAGC : 397
    
      44 : euSerLysMetSerLeuLysProAspLysIlePheTyrValProPheLys  >>>> Target  :  60
           ||! !|||||||||! !|||! !|||||||||||||||||||||||||||            41
           euLeuLysMetSerProLysLeuAspLysIlePheTyrValProPheLys++            
     396 : TGTTAAAAATGTCGCCCAAGCTGGACAAAATCTTCTACGTGCCCTTCAAAgt............ : 344
    
      61 : Intron 1 >>>>  IleSerThrGluThrIleIleTyrThrArgProLeuSerMetHisLeuA :  76
            bp            |||! !|||! !! !||||||! !||||||! !||||||||||||||||
                        ++IleLeuThrAlaIleIleIleCysThrArgHisLeuSerMetHisLeuA
     343 : .............agATTTTAACCGCAATCATTATCTGTACACGGCATCTCTCTATGCACTTGC : 257
    
      77 : rgPheLeuAlaLeuGlnProAsnArgMetPheValAsnAlaGlnSerAsnLeuLys :  94
           ||||||||!.!|||!:!! !|||!:!|||  !||||||!.!||||||!:!||||||
           rgPheLeuValLeuArgLeuAsnLysMetArgValAsnValGlnSerSerLeuLys
     256 : GTTTTCTAGTTCTACGGCTGAACAAAATGCGCGTGAATGTACAAAGCAGCTTGAAA : 201
    
    MhA1_Contig0.frz3.gene9,95,0,94,NODE_769_length_6021_cov_89.450089_6075_200.51_0.31,6075,524,200,94,73,78,21,0,77.66,82.98, M 99 D 1 M 78 D 41 M 105
    


--showalignment T is needed because we need the exact codons that correspond to the alignment

### 2. Convert exonerate alignments to CDS and amino acid fasta files

First, get the codons and the corresponding three letter amino acid abbreviation (eg: ATG and Met Respectively) from the alignment.

The command line used below assumed that the query sequences all started with M or W (Mhapla and Wormbase incognita proteins), so the `while ($_ !~ /^[MW]/)` line is used as a terminating condition, because every alignment ends with a line (that begins with one of these two letters)

    cat $protein.$genome.p2g | perl -ne '
      if (/Target range:/) {
        <>;<>;
        while ($_ !~ /^[MW]/) {
          <>;
          chomp($target_aa .= <>);
          chomp($target_na .= <>);
          <>;
          chomp($_=<>);
        }
        $target_na =~ s/[a-z\s\:\d\-\.]//g;
        $target_na =~ s/{.+?}//g;
        $target_aa =~ s/[\s\-\+]//g;
        $target_aa =~ s/{.+?}//g;
        print "$_\t$target_aa\t$target_na\n";
        $target_aa = $target_na = "";
      }
    ' >$protein.$genome.p2g.aaa

Next, take this line with 3 letters for each amino acid, and convert it to the amino acid IUPAC code:

    cat $protein.$genome.p2g.aaa | perl -ne '
      @F  = split /\t/;
      my ($aaa_pos, $aaa_seq, $fna_seq, $A_seq);
      while ($F[1] =~ /(.+?)(#*)/g) {
        $aaa_seq .= substr ($F[1], $aaa_pos, length($1));
        $fna_seq .= substr ($F[2], $aaa_pos, length($1));
        $aaa_pos += length($1) + length ($2);
      };
      if ($aaa_seq =~ s/\*\*\*.*//) {$fna_seq = substr($fna_seq, 0, length($aaa_seq))};
      %tla=("Ala"=>"A","Asx"=>"B","Cys"=>"C","Asp"=>"D","Glu"=>"E","Phe"=>"F","Gly"=>"G","His"=>"H","Ile"=>"I","Lys"=>"K","Leu"=>"L","Met"=>"M","Asn"=>"N","Pro"=>"P","Gln"=>"Q","Arg"=>"R","Ser"=>"S","Thr"=>"T","Val"=>"V","Trp"=>"W","Unk"=>"X","Tyr"=>"Y","Glx"=>"Z");
      while ($aaa_seq =~ /(...)/g) {die "$1 not found" if not exists $tla{$1}; $A_seq .= $tla{$1}};
      print "$F[0],$A_seq,$fna_seq\n";
    ' >$protein.$genome.p2g.A.nuc

The output above is a comma separated line for each alignment, with the --ryo format info taken from exonerate output, followed by the protein aa sequence and then by the nucleotide cds sequence. eg:

    MhA1_Contig999.frz3.gene2,316,21,316,NODE_35790_length_6901_cov_32.895813_6955_78.52_0.36,6955,1387,3398,285,217,231,68,0,76.14,81.05, M 99 D 3 M 42 D 48 M 59 I 4 M 51 D 31 M 44 I 2 M 102 D 4 M 105 D 3 M 24 D 646 M 77 D 3 M 87 D 352 M 93 D 58 M 80,KRLSKFLNNIPLEEAETPAPEAPAAPSYNAAEASSAEQTAPAEQPAAEHAASEEAAPSYNNNAVSAPAAAPAAPEPAPAGLF,AAAAGGTTATCAAAATTTCTAAACAATATCCCATTAGAAGAAGCCGAAACTCCTGCTCCTGAAGCACCAGCCGCACCGAGTTATAACGCAGCTGAAGCATCGTCTGCTGAACAAACTGCACCTGCCGAACAGCCGGCAGCTGAGCATGCCGCTTCAGAAGAAGCAGCTCCAAGCTACAATAATAATGCAGTTTCTGCTCCAGCTGCAGCTCCAGCAGCGCCCGAGCCGGCCCCAGCTGGTTTGTTT
    
Do the above alignment and initial processing steps for both protein files (M. hapla and M. incognita). The two sets of alignments may map to overlapping locations on the M. floridensis genome, and these coordinates need to be merged:

    cat *.$genome.p2g.A.nuc |
    perl -lanF',' -e '($F[6],$F[7])=($F[7],$F[6]) if $F[6]>$F[7]; print "$F[4]\t$F[6]\t$F[7]"' |
    sort -k1,1 -k2,2n |
    mergeBed  >mf-all.p2g.mergeBed

The `perl -lanF` bit prints a 3 col tab delimited BED format output with these columns:

1. Genome contig (or scaffold) name
2. Genome start location of alignment
3. Genome end   location of alignment (coordinates are swapped if the first number is larger than the second because BEDtools mergeBed expects the first number to always be smaller than the second number)

`sort -k1,1 -k2,2n` is needed because mergeBed expects the input to always be ordered by contig, then the start coordinate.

Pick the longest of the original alignments that intersects with this merged interval set:

    cat *.$genome.p2g.A.nuc | perl -lanF',' -e '
      ($F[6],$F[7])=($F[7],$F[6]) if $F[6]>$F[7];
      print "$F[4]\tmf\tregion\t".($F[6]+1)."\t$F[7]\t.\t.\t.\t$F[0],$F[1],$F[2],$F[3],$F[16],$F[17]"
    ' | intersectBed   -a mf-all.p2g.mergeBed -b - -wa -wb | awk '{print $1","$2","$3"\t"$12}' >mf-all.p2g.mergeBed.intersect
    
The perl `print` command creates a GFF format entry with all nuc and aa information for each p2g alignment.

The file mf-all.p2g.mergedBed.intersect has the same coordinates present multiple times (because each could have overlapped with multiple alignments). So we only report the longest alignment for each one:

    perl <mf-all.p2g.mergeBed.intersect -lane '
        @P=split(/,/,$F[1]);
        $F[0] =~ s/,/_/g;
        $faa{$F[0]}=$P[4] if length($P[4])>length($faa{$F[0]});
        $fna{$F[0]}=$P[5] if length($P[5])>length($fna{$F[0]});
    }{
        foreach (keys %faa) {
            print ">$_\n".$faa{$_};
            print STDERR ">$_\n".$fna{$_};
        }' >mf-all.p2g.mergeBed.protein.faa 2>mf-all.p2g.mergeBed.cds.fna

These two are the final CDS and Protein files that are used for all further analyses:

1. mf-all.p2g.mergeBed.protein.faa
2. mf-all.p2g.mergeBed.cds.fna

Calculate and plot CDS self-identity
------------------------------------------------------------

Requirements:

1. Species CDS file (coding transcripts) nucleotide fasta (renamed as mh.cds.fna.98 mi.cds.fna.98 mf.cds.fna.98)
2. Blastn (from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.25/)

The following bash command line will blast each cds file against itself, take the top 5 hits (1e-5) and keep only those hits where a sequence does not hit itself, and where the hit is at least 70% of the length of the query. The `perl -ane '$h{$F[0]}++;print if $h{$F[0]}==1'` bit ensures that only the top hit for each query is kept (not the others).

    for a in mh mi mf
    do makeblastdb -dbtype nucl -in $a.cds.fna.98;
       blastn -task blastn -query $a.cds.fna.98 -max_target_seqs 5 -outfmt '6 std qlen slen' \
         -num_threads 4 -db $a.cds.fna.98 -evalue 1e-5 -out $a.$a.cds.98.blastn.1e-5;
       cat $a.$a.cds.98.blastn.1e-5 |
       sort -t $'\t' -k 1,1 -k 12,12nr |
       awk '$4/$13>0.7 && $1!=$2' |
       perl -ane '$h{$F[0]}++;print if $h{$F[0]}==1' |
       awk '{print "'$a'\t"$3}' >>mhmimf.98.self.id
    done

The file mhmimf.98.self.id has a 2 col tab format, where col 1 is the species name, and col 2 is the percentage id of the best hit by that CDS to any CDS other than itself in that species. eg:

    mh      75.00
    mh      83.12
    ...
    mi      100.00
    mi      98.52
    ...
    mf      87.45
    mf      92.16
    ...

To plot these identities, use R:

    # R
    library(ggplot2)
    d=read.delim("mhmimf.98.self.id",header=F,col.names=c("Species","ID"))
    d$Species<-ordered(d$Species,levels=c("mh","mi","mf"))
    theme_set(theme_bw())

    qplot(data=d,x=ID, color=Species,geom="freqpoly", binwidth = 0.5) +
        xlim(75,101) + opts(legend.position = c(0,1), legend.justification=c(0,1)) + xlab("% Identity") + ylab("Count") +
        scale_color_manual(values = c("mh"="black","mi"="blue","mf"="red"),labels=c("M hapla","M incognita","M floridensis"))  

Cluster proteins using InParanoid and QuickParanoid
------------------------------------------------------------

Requirements:

* InParanoid (from http://software.sbc.su.se/cgi-bin/request.cgi?project=inparanoid - I used version 4.1)
* QuickParanoid (from http://pl.postech.ac.kr/QuickParanoid/#install - version downloaded on 2011-06-19)
* Protein files for all three species, renamed as mh, mi, mf
* Blastp (from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.25/)

I could have let InParanoid do the blasts on its own, but that's very slow, so I used our cluster to split the jobs. This is how you run InParanoid if you run the blasts separately:

    formatdb -i mh -p T
    formatdb -i mi -p T
    formatdb -i mf -p T
    blastall -F F -i mh -d mh -p blastp -v 12229 -b 12229 -M BLOSUM62 -z 5000000 -m7 | gzip >mh-mh.gz
    blastall -F F -i mh -d mi -p blastp -v 17999 -b 17999 -M BLOSUM62 -z 5000000 -m7 | gzip >mh-mi.gz
    blastall -F F -i mh -d mf -p blastp -v 15121 -b 15121 -M BLOSUM62 -z 5000000 -m7 | gzip >mh-mf.gz
    blastall -F F -i mi -d mh -p blastp -v 12229 -b 12229 -M BLOSUM62 -z 5000000 -m7 | gzip >mi-mh.gz
    blastall -F F -i mi -d mi -p blastp -v 17999 -b 17999 -M BLOSUM62 -z 5000000 -m7 | gzip >mi-mi.gz
    blastall -F F -i mi -d mf -p blastp -v 15121 -b 15121 -M BLOSUM62 -z 5000000 -m7 | gzip >mi-mf.gz
    blastall -F F -i mf -d mh -p blastp -v 12229 -b 12229 -M BLOSUM62 -z 5000000 -m7 | gzip >mf-mh.gz
    blastall -F F -i mf -d mi -p blastp -v 17999 -b 17999 -M BLOSUM62 -z 5000000 -m7 | gzip >mf-mi.gz
    blastall -F F -i mf -d mf -p blastp -v 15121 -b 15121 -M BLOSUM62 -z 5000000 -m7 | gzip >mf-mf.gz
    
    for a in {mh,mi,mf}-{mh,mi,mf}
    do
      zcat $a.gz | blast_parser.pl 40 >$a
    done

(Note: The blast_parser.pl script ships with InParanoid)

    perl inparanoid.pl mi mf
    perl inparanoid.pl mh mi
    perl inparanoid.pl mf mh
    
    #run quickparanoid
    qp; ./tests >results.txt

QuickParanoid gives results in this form (results.txt), eg:

    #clusterID      species gene    is_seed_ortholog        confidence_score        species_in_cluster      tree_conflict   
    1       mh      MhA1_Contig1687.frz3.gene1      1       1.000   mh-mi   No
    1       mi      Minc02622a      1       1.000   mh-mi   No
    1       mi      Minc02622b      1       1.000   mh-mi   No
    2       mh      MhA1_Contig1210.frz3.gene3      1       1.000   mh-mi   No
    2       mi      Minc13840       1       1.000   mh-mi   No
    3       mh      MhA1_Contig147.frz3.gene57      1       1.000   mh-mi-mf        No
    3       mi      Minc07000       1       1.000   mh-mi-mf        No
    3       mf      mf_2733 1       1.000   mh-mi-mf        No

Use this file above to get counts of clusters with these ratios of proteins: 1Mh:2Mi:1Mf, 1Mh:2Mi:2Mf etc:

    perl <result.txt -ane '
        next if /^#/;
        $cluster{$F[0]}{$F[1]}++;
    }{
        print "clusterid\tmh\tmi\tmf\n";
        for my $clusterid (sort {$a <=> $b} keys %cluster) {
            print $clusterid;
            for my $species qw(mh mi mf) {
                #count >4 converted to 4
                $cluster{$clusterid}{$species} = 4 if exists $cluster{$clusterid}{$species} and $cluster{$clusterid}{$species}>4;
                print "\t" . (exists $cluster{$clusterid}{$species} ? $cluster{$clusterid}{$species} : "0");
            }
            print "\n"
        }
    ' >result.txt.counts
    
    for a in {0..4}; do for b in {0..4}; do for c in {0..4}; do grep -c -P "^\d+\t$a\t$b\t$c" result.txt.counts; done; done; done >tally.txt
    R
    tally<-array(scan("tally.txt"),c(5,5,5))
    # 1Mh:_Mi:_Mf
    tally[,,2]
         [,1] [,2] [,3] [,4] [,5]
    [1,]    0  907  327   44   17
    [2,] 2196 2189  920  102   40
    [3,]  226  257  156   36   21
    [4,]   17   17   20    7   14
    [5,]    8   11    6    4   21

Use RAxML and APE to create and analyse multiple phylogenies
------------------------------------------------------------

Requirements

* RAxML (from http://sco.h-its.org/exelixis/countSource728.php)
* multi2single (from this repository)
* Clustal Omega (from http://www.clustal.org/omega/ - we used version 1.0.3)
* EMBOSS suite (from http://emboss.sourceforge.net/ - we used version 6.2.0)
* Protein files for all three species, renamed as mh, mi, mf
* CDS nucleotide fasta files for all three species, renamed: mh.cds mi.cds mf.cds. The CDS files and protein files must correspond exactly (same order, same sequence IDs, etc. Check this before starting)

Overview: 

To analyse all clusters with a particular ratio of proteins from each species (eg 1Mh:2Mi:1Mf) sets, we need to:

- work through all the protein and cds fasta files and load all nucleotide and amino acid sequences in memory (in two hashes %aa and %nuc, key = seqid)
- work through the results.txt file, and make a %clusters hash (that stores all the sequences belonging to each cluster)
- once all this is done (i.e. the bit upto the line with `}{`), then count the number of each type of protein in each cluster. If it matches the set we are working with (in this example, 1Mh:2Mi:1Mf), then print these 4 cds and protein sequences into a cluster nucleotide file (c$cluster.fna) and a cluster protein file (c$cluster.faa) respectively:

    mkdir raxml && cd raxml
    
    # set bash variables for the counts that we are pulling out results for:
    mh=1;mi=2;mf=1

    mkdir ${mh}mh${mi}mi${mf}mf

    cat \
      <(paste \
        <(multi2single -l 2 ../mh) \
        <(multi2single -l 2 ../mh.cds) \
       ) \
      <(paste \
        <(multi2single -l 2 ../mi) \
        <(multi2single -l 2 ../mi.cds) \
       ) \
      <(paste \
        <(multi2single -l 2 ../mf) \
        <(multi2single -l 2 ../mf.cds) \
       ) \
      ../result.txt | # pipe all the fasta files through
    perl -lane '
      /^#/ and next; # header of result.txt
      /^>(\S+)/ and $aa{$1} = $F[1] and $nuc{$1} = $F[3] and next;  #load fasta seqs in hash:
      $clusters{$F[0]}{$F[2]}=1;   # all other lines are from result.txt - so make clusters hash
    }{
      for my $cluster (keys %clusters) {
        @seqids = keys %{$clusters{$cluster}} and $_ = "@seqids";
        next unless s/MhA/MhA/g == '$mh' and s/Minc/Minc/g == '$mi' and s/mf_/mf_/g == '$mf';
        open NUC, ">'$mh'mh'$mi'mi'$mf'mf/c$cluster.fna";
        for my $seq (@seqids) {
          ($s = $seq) =~ s/MhA\S+/mh/;
          print NUC ">$s\n$nuc{$seq}"
        }
        close NUC;
        open  AA, ">'$mh'mh'$mi'mi'$mf'mf/c$cluster.faa";
        for my $seq (@seqids) {
          ($s = $seq) =~ s/MhA\S+/mh/;
          print AA ">$s\n$aa{$seq}"
        }
        close AA;
      }
    '

Now, we need to run RAxML for each CDS cluster file. To do that,
- We need to align each protein cluster file
-- We use clustalo to align the protein files.
- Then we use tranalign (from the emboss suite) to align the corresponding CDS file for that cluster using the protein alignment as a guide.
- Convert the CDS alignments to PHYLIP format (using MFAtoPHY.pl from omics.informatics.indiana.edu/~yuwwu/getfile.php?MFAtoPHY.pl)
- Run RAxML

Running the above steps:
    
    cd ${mh}mh${mi}mi${mf}mf

Align all protein files using clustalo:

    ls *.faa | parallel "clustalo -i {} -o {.}.clustalo --force --full --guidetree-out={.}.guidetree"
    
Use tranalign to align all cds files, using protein alignments as a guide:

    ls *.faa | parallel "tranalign -asequence {.}.fna -bsequence {.}.clustalo -outseq {.}.tranalign"
    
Convert all tranalign-ed files to phylip format:

    ls *.tranalign | parallel "MFAtoPHY.pl {}"
    
Run RAxML with GTRGAMMA model with mh as the outgroup:

    for a in *phy
    do
      raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -s $a -# 10 -n $a -T 2;
      rm RAxML_info* RAxML_log* RAxML_result* RAxML_parsimonyTree*;
      raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -s $a -# autoMRE -n $a.b -T 2 -b 12345;
      rm RAxML_info*;
      raxmlHPC-PTHREADS-SSE3 -m GTRCAT -f b -t RAxML_bestTree.$a -z RAxML_bootstrap.$a.b -n $a.l -T 2 -o mh;
      rm RAxML_info* RAxML_bestTree* RAxML_bootstrap*;
    done

Consolidate files (remove unnecessary files, tar/gzip others)

    cat RAxML_bipartitionsBranchLabels.c* > ${mh}mh${mi}mi${mf}mf.RAxML_bipartitionsBranchLabels; rm RAxML_bipartitionsBranchLabels.c*
    cat RAxML_bipartitions.c* > ${mh}mh${mi}mi${mf}mf.RAxML_bipartitions; rm RAxML_bipartitions.c*
    paste <(ls *.fna) ${mh}mh${mi}mi${mf}mf.RAxML_bipartitions ${mh}mh${mi}mi${mf}mf.RAxML_bipartitionsBranchLabels >${mh}mh${mi}mi${mf}mf.RAxML_bipartitions.all
    tar czf faa.tgz *faa; rm *faa
    tar czf fna.tgz *fna; rm *fna
    tar czf phy.tgz *phy; rm *phy
    tar czf reduced.tgz *reduced; rm *reduced
    tar czf guidetree.tgz *guidetree; rm *guidetree
    tar czf tranalign.tgz *tranalign; rm *tranalign
    tar czf clustalo.tgz *clustalo; rm *clustalo

Remember, this has to be done again for each clustering of interest (1Mh:2Mi:1Mf, 1Mh:1Mi:2Mf, 1Mh:1Mi:3Mf, 1Mh:32Mi:1Mf, 1Mh:2Mi:2Mf, 1Mh:2Mi:3Mf, 1Mh:3Mi:2Mf)

Now, we need to count the number of trees of each topology. We're just doing this for 1Mh:2Mi:1Mf for now i.e. only two alternatives are possible - (mi (mi mf)), or (mf (mi mi)) 

    cd 1mh2mi1mf
    # start R
    library(ape)
    dir="1mh2mi1mf"
    trees<-read.tree(sprintf("%s.RAxML_bipartitions",dir))
    
    # rename all the nodes to F (for floridensis), H (for hapla) and I (for incognita), i.e. two M incognita sequences will be called I1 and I2.
    # the tree ordering takes always names outward (i.e. innermost tree is I1, then I2... and so on)

    for (t in 1:length(trees))
    {
      tree.reordered <- reorder(trees[[t]],"p")
      seq.reordered  <- tree.reordered$tip.label[tree.reordered$edge[,2]]
      seq.reordered  <- seq.reordered [!is.na(seq.reordered)]
      mi_counter = 1
      mf_counter = 1
      seq.rename = vector()
      for (seq.name in seq.reordered)
      {
        if(substr(seq.name,1,2)=="Mi") {
          seq.rename[seq.name] = sprintf("I%d",mi_counter);
          mi_counter <- mi_counter + 1;
        } else
        if(substr(seq.name,1,2)=="mf") {
          seq.rename[seq.name] = sprintf("F%d",mf_counter);
          mf_counter <- mf_counter + 1;
        } else
        if(substr(seq.name,1,2)=="mh") {
          seq.rename[seq.name] = "H";
        }    
      }
      for (i in 1:length(seq.rename))
        tree.reordered$tip.label[i]=seq.rename[tree.reordered$tip.label[i]]
      trees[[t]]<-tree.reordered
    }

    # Get all the unique topologies in all the trees:
    u <- unique.multiPhylo(trees,use.edge.length=F)
    
    pdf(sprintf("%s-uniquetrees.pdf",dir),paper="a4",width=7,height=11)
    layout(matrix(c(1:15), 5,3,byrow=TRUE))
    
    # count how many trees match each unique topology, and print that topology with the count on top of it:
    
    count=vector(length=length(u))
    for (utree in 1:length(u))
    {
      count[utree] <- sum(unlist(lapply(X=trees,function(x) all.equal.phylo(x,u[[utree]],use.edge.length=F))))
      plot (u[[utree]], use.edge.length=F, font=2, cex=2);
      mtext(toString(count[utree]), font=2, cex=2)
    }
    dev.off()
    q(save="no")
