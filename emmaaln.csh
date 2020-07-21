#!/bin/csh

#Aligns *.fasta (translates, emma alns aa seqs, back translates).
#Assumes *.fasta are nt sequences corresponding to 0-frame open reading frames.

#In-frame stop codons are replaced with
#  sed 's/UAG/RHK/' | sed 's/UAA/YBS/' | sed 's/UGA/MVW/' 
#and changed back afterwards.

@ argc = `echo $argv | awk '{print NF}'`
if ($argc != 1) then
  echo "Usage: 'emmaaln muscle' or 'emmaaln clustalw'"
  exit 1
endif

set alnprog = $1
set alntype = `echo $alnprog | awk '{if ($1=="muscle") {print 1} else if ($1=="clustalw") {print 2} else {print 0}}'`
if (! $alntype) then
  echo "Usage: 'emmaaln muscle' or 'emmaaln clustalw'"
  exit 1
endif

ls *.fasta | sed 's/\.fasta//' > seqs.txt
foreach i (`cat seqs.txt`)
  echo ">$i" > tfh784hfgt
  tail -n +2 $i.fasta | awk -v ORS="" '{print $0}' | sed 's/.../&\n/g' | \
    sed 's/[a-z]/\U&/g' | sed 's/T/U/g' | grep "." | \
    awk '{if (length($1)==1) {printf "%sNN",$1} else {print $0}}' | \
    awk '{if (length($1)==2) {printf "%sN",$1} else {print $0}}' | \
    sed 's/UAG/RHK/' | sed 's/UAA/YBS/' | sed 's/UGA/MVW/' >> tfh784hfgt
  seqret tfh784hfgt $i.fasta -auto
  rm -f tfh784hfgt
end

rm -f nucseqs.txt pepseqs.txt
touch nucseqs.txt pepseqs.txt
foreach i (`cat seqs.txt`)
  transeq $i.fasta $i.pep -auto
  cat $i.fasta >> nucseqs.txt
  cat $i.pep >> pepseqs.txt
end
if (2 == $alntype) then
  emma pepseqs.txt clustalw.aln clustalw.dnd -auto
else if (1 == $alntype) then
  emma pepseqs.txt clustalw.aln clustalw.dnd -auto
  rm -f clustalw.aln
  muscle -in pepseqs.txt -out clustalw.aln -quiet
else
  echo "Can't get here."
endif
mkdir Tjg78HjdER
rm -f pepseqs.aln; touch pepseqs.aln
cd Tjg78HjdER
seqretsplit ../clustalw.aln temp -auto
foreach i (`cat ../seqs.txt`)
  set j = `echo $i | sed 's/[A-Z]/\L&/g'`
  cat ${j}_1.fasta >> ../pepseqs.aln
end
cd ..
rm -rf Tjg78HjdER
tranalign nucseqs.txt pepseqs.aln allseqs.txt -auto
awk -v ORS="@" '{print $0}' allseqs.txt | sed 's/@>/\n>/g' | grep "@" | \
  sed 's/@/ /' | sed 's/@//g' | \
  sed 's/RHK/UAG/g' | sed 's/YBS/UAA/g' | sed 's/MVW/UGA/g' | \
  awk '{printf "%s\n%s\n",$1,$2}' > allseqs1.txt
seqret allseqs1.txt allseqs2.txt -osformat aln -auto
embossaln2clustalw allseqs2.txt allseqs3.txt

emmadnd2pairs clustalw.dnd tfh784hfgt
awk '{printf "%s %s \n",$1,$2}' tfh784hfgt | sed 's/_1 / /g'  > allpairs.txt
rm -f tfh784hfgt

rm -f *.pep
rm -f allseqs.txt allseqs1.txt allseqs2.txt
mv allseqs3.txt allseqs.txt
rm -f nucseqs.txt pepseqs.txt pepseqs.aln
