a=$1

#cat annotation.gff.features.coding <(zcat $a.gff.gz) | perl -lane '
#  next if /^#/;
#  $coding{$1}=1 and next if /^(\S+)$/;
#  print "'$a'.$_" if $coding{$F[2]} and $F[3] <= $F[4];
#' | gzip -c >$a.coding.gff.gz

cat $a.23.counts.span.coding <(zcat $a.gff.gz) | perl -F"\t" -lane '
  next if /^#/;
  $coding{$1}=1 and next if /^\S+\t\d+\t(\S+\t\S+)$/;
  next unless scalar @F==9 and $F[3] <= $F[4]; #ie.it should be valid GFF
  print "'$a'.$_" if $coding{"$F[1]\t$F[2]"};
' | cat - ~/media/WS230/2012-03-13-rfam.trna/$a.fna.*.gff | grep -v "^#" | gzip -c >$a.coding.gff.gz
