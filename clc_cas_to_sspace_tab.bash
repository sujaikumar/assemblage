#!/bin/bash

# sujai.kumar@ed.ac.uk 2012-07-09

CAS=$1 # clc cas file with mplib mapped to assemblt fasta file
FAS=$2 # assembly fasta file

DIST=600
if [ $3 ]
then DIST=$3
fi

# convert clc to SSPACE tab format,
# keeping only those read pairs that map to diff contigs
# Doing this separately because this is the step that takes time

assembly_table -a $CAS |
perl -lne '
  next unless /\/1 has (\d+) optimal match/ and $1 == 1;
  <>; $_=<>;
  /(\d+)\s+\S+\s+(\d+)\s+(\S+)/ and $st1 = $1 and $en1 = $2 and $contig1 = $3;
  <>;$_=<>;
  /reverse read/ and ($st1,$en1) = ($en1, $st1);
  <>;$_=<>;
  next unless /\/2 has (\d+) optimal match/ and $1 == 1;
  <>;$_=<>;
  /(\d+)\s+\S+\s+(\d+)\s+(\S+)/ and $st2 = $1 and $en2 = $2 and $contig2 = $3;
  <>;$_=<>;
  /reverse read/ and ($st2,$en2) = ($en2, $st2);
  next if $contig1 eq $contig2;
  print join("\t",$contig1,$st1,$en1,$contig2,$st2,$en2);
' | gzip >$CAS.tab.gz

# now remove all pairs where one of the reads maps within DIST bp from the end of the contig:

cat \
  <(fastaqual_multiline_to_singleline.pl $FAS) \
  <(zcat $CAS.tab.gz) |
perl -MList::Util='max,min' -lane '
  if(/>(\S+)/) { $id = $1; $_=<>; $len{$id} = length($_) - 1; next}

  ($contig1,$st1,$en1,$contig2,$st2,$en2) = @F;

  $l = '$DIST';
  next if
    min($st1,$en1) < $l or
    $len{$contig1} - max($st1,$en1) < $l or
    min($st2,$en2) < $l or
    $len{$contig2} - max($st2,$en2) < $l;
  print;
' >$CAS.dist$DIST.tab
