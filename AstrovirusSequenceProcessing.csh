#!/bin/csh

#Astrovirus sequence retrieval and processing script relating to the manuscript
#  Lulla V, Firth AE (2020) A hidden gene in astroviruses encodes a viroporin.

#------------------------------------------------------------------------------

#26 July 2018 - tblastn of NC_001943 ORF1b vs taxid35278 (ssRNA+ viruses),
#  wordsize 3, expect threshold 1, max 500 sequences; entrez query
#  5000:10000[Sequence Length] 
#-> tblastn1.txt (truncated where it switched to avian astrovirus hits)
#-> 241 sequences 

awk '{if (NR%3==0) print $0}' tblastn1.txt  | sed 's/%//g' | \
  awk '{print $3}' | sort -n | uniq -c | more
#-> 95 to 100% ORF1b coverage
awk '{if (NR%3==0) print $0}' tlastn1.txt  | sed 's/%//g' | \
  awk '{print $5}' | sort -n | uniq -c | more
#-> 46 to 98% identity (had masking switched on)

awk '{if (NR%3==0) print $6}' tblastn1.txt > seqs.txt
#-> Download -> NR/GB/
rm NC_004579.gb NC_037655.gb 
#Since their "derived-from" sequences MG660832.gb AY179509.gb are present. Not
#  sure why these were the only 2 Refseqs that came down; perhaps others
#  removed as redundant whereas maybe these two sequences had some post-
#  submission author modifications(?).

#Discard sequences without "VRL" in the "LOCUS" field.
#-> 1 seq is ENV - looks OK (MG571777); indeed it is a good one to keep as it
#  is 10% aa divergent from other HAstVs.
#
#Discard sequences with >20 "ambiguous nucleotide" codes.
#-> Only two discarded as having too many ambig nts viz KP242038 and KY271945
#  (they both have quite a lot of Xs).
mkdir Ambig
mv KP242038.gb KY271945.gb Ambig/
#-> 237 seq left.

#For the remainder, find the 3 ORFs. Probably suffices to find all 
#  stop-stop ORFs > 1500 nt to start with.
rm -f log1.txt
foreach i (`ls *.gb | sed 's/\.gb//'`)
  getorf $i.gb temp1 -minsize 1500 -find 0 -auto
  set norfs = `grep ">" temp1 | wc -l`
  echo $i $norfs >> log1.txt
  rm -f temp1
end
awk '{print $2}' log1.txt | sort | uniq -c
#      7 2
#    232 3
#-> 7 sequences fail. Check these.
mkdir Fail3ORFs
foreach i (`awk '{if ($2!=3) print $1}' log1.txt`)
  mv $i.gb Fail3ORFs
end
cd Fail3ORFs
#Manually check them all to check we are correct to omit them -> confirmed.
cd ..
#-> 230 sequences left.

#Coords of the three stop-stop ORFs in genome location order.
rm -f log2.txt
foreach i (`ls *.gb | sed 's/\.gb//'`)
  getorf $i.gb temp1 -minsize 1500 -find 0 -auto
  echo -n "$i " >> log2.txt
  grep ">" temp1 | awk '{printf "%s..%s ",$2,$4}' | sed 's/\[//g' | \
    sed 's/\]//g' >> log2.txt
  echo >> log2.txt
  rm -f temp1
end

#Check for 5' and 3' coding completeness.
mkdir ORF1a CP
foreach line (`awk '{printf "%s:%s\n",$1,$2}' log2.txt`)
  set acc = `echo $line | awk -F: '{print $1}'`
  set cds = `echo $line | awk -F: '{print $2}'`
  echo ">$acc" > ORF1a/$acc.fasta
  extractseq $acc.gb temp1 -regions $cds -auto
  tail -n +2 temp1 >> ORF1a/$acc.fasta
  rm temp1
end
foreach line (`awk '{printf "%s:%s\n",$1,$4}' log2.txt`)
  set acc = `echo $line | awk -F: '{print $1}'`
  set cds = `echo $line | awk -F: '{print $2}'`
  echo ">$acc" > CP/$acc.fasta
  extractseq $acc.gb temp1 -regions $cds -auto
  tail -n +2 temp1 >> CP/$acc.fasta
  rm temp1
end
#In each of ORF1a/, CP/. Remember these are stop-stop ORFs currently.
emmaaln muscle
rm *.fasta
seqret clustalw.aln clustalw2.aln -osformat aln
emacs clustalw2.aln &

#ORF1a:
#From inspection of clustalw2.aln, the following look to be missing a large
#  part of ORF1a 5' end: GQ267696 JN420353 JN420355 KF157967 KJ571486 KY933399
grep " CDS " ../*gb | grep "<"
#We see exactly the same sequences as above coming up with significantly
#  shorter-than-expected ORF1a:
#    GQ267696 <1..2063 -> very similar to KM401565 (33..2696), GQ891990
#                         (1..2694) etc
#    JN420353 <1..1599 -> very similar to JN420351
#    JN420355 <1..1535 -> very similar to JN420351
#    KF157967 <1..1707 -> very similar to JQ403108 and DQ028633
#    KJ571486 <1..2208 -> close enough to LC201589 (3..2501) and LC201594
#                         (3..2510)
#    KY933399 <1..1615 -> highly divergent (but genogroup 2 so not expected to
#                         have ORFx)
#Note, compared the two with longest ORF1a (2063 and 2208 nt) against most-
#  closely related sequences to see expected ORF1a lengths, confirming that
#  these are considerably truncated.
#-> OK to delete all of these.
cd ..
mkdir ORF1a_truncated
mv {GQ267696,JN420353,JN420355,KF157967,KJ571486,KY933399}.gb ORF1a_truncated/
cd -
#MLB2
#    KX022687 <1..2338 -> just missing ~9 N-terminal aa
#These ones are just missing a little bit of 5' end of ORF1a sequence.
#    JN420358 <1..2553
#    KC692365 <1..2542
#    KT946735 <1..2742
#    KT946736 <1..2742
#    KY940075 <1..2602
#    LC201604 <1..2634
#    LC201607 <1..2597
#    LC201614 <1..2622
#    MG693176 <1..2677
#-> So OK to keep these ones.

#CP:
#From inspection of clustalw2.aln, the following look to be missing a large
#  part of CP 3' end: KY933399
grep " CDS " *gb | grep ">" | sed 's/\.gb://' | sed 's/\.\.>/ /' | \
  awk '{print $4-$3+1,$1,$3,$4}' | sort -nr 
#    2517 JN420354 4069 6585
#    2483 JN420351 4168 6650
#    2476 JN420355 2901 5376
#    2475 JN420357 4075 6549
#    2437 LC201602 4077 6513
#    2423 LC201604 4072 6494
#    2401 LC201606 4082 6482
#    2401 JN420358 3922 6322
#    2397 LC201603 4080 6476
#    2336 KY271946 4318 6653
#    2298 JN420352 4343 6640 <- probably essentially complete
#    2282 LC201607 4035 6316 <- missing about 60 aa and closely related seqs
#                               present so delete this
#    2148 KY933399 3134 5281 <- delete
#    1962 LC201591 3969 5930 <- delete
#    1955 LC201614 4060 6014 <- delete
#-> So remove only LC201614 LC201591 KY933399 LC201607
cd ..
mkdir CP_truncated
mv {LC201614,LC201591,KY933399,LC201607}.gb CP_truncated/
cd -

#-> 221 sequences left.


#Remake log2 (now called log3) with the above sequences deleted
rm -f log3.txt
foreach i (`ls *.gb | sed 's/\.gb//'`)
  getorf $i.gb temp1 -minsize 1500 -find 0 -auto
  echo -n "$i " >> log3.txt
  grep ">" temp1 | awk '{printf "%s..%s ",$2,$4}' | sed 's/\[//g' | \
    sed 's/\]//g' >> log3.txt
  echo >> log3.txt
  rm -f temp1
end

#Extract and align the RdRp ORF, then manually truncate to AAAAAAC start. Then
#  make a MrBayes tree so we can select representative sequences.
mkdir ORF1b
foreach line (`awk '{printf "%s:%s\n",$1,$3}' log3.txt`)
  set acc = `echo $line | awk -F: '{print $1}'`
  set cds = `echo $line | awk -F: '{print $2}'`
  echo ">$acc" > ORF1b/$acc.fasta
  extractseq $acc.gb temp1 -regions $cds -auto
  tail -n +2 temp1 >> ORF1b/$acc.fasta
  rm temp1
end
cd ORF1b
# For refseq alignment, need to truncate NC_025379 and NC_030922 before
#  aligning otherwise the alignment is quite messed up.
emmaaln muscle
cp allseqs.txt allseqs_shiftsitetruncated.txt
emacs allseqs_shiftsitetruncated.txt &
#-> Find the shift site and edit off 5' to that.

#Transeq and realign
transeq allseqs_shiftsitetruncated.txt allseqs_shiftsitetruncated.pep
muscle -in allseqs_shiftsitetruncated.pep -out muscle.aln -clwstrict
emacs muscle.aln &
#-> Inspect. -> 3' end looks OK (i.e. no partial ORF1b sequences).
#  Due to some of the Refseqs having indel problems fusing ORFs 1a and 1b, the
#  very N-term end is perhaps a bit messed up, so deleting about the first 10
#  aa. Also delete about the last 5 C-terminal aa.
seqret muscle.aln temp1
sed 's/_1$//' temp1 > temp2
seqret temp2 astro_orf1b.nex -osformat nexus

rm -f nameskey.txt
foreach i (`ls *fasta | sed 's/\.fasta//'`)
  echo -n "${i}:" >> nameskey.txt
  grep ORGANISM ../$i.gb | sed 's/ *ORGANISM *//' | \
     awk '{printf "%s - %s\n","'"$i"'",$0}' | sed 's/ /@/g' | \
     sed 's/\//£/g' >> nameskey.txt
end
emacs nameskey.txt &

mkdir R1000000
cp astro_orf1b.nex R1000000
cd R1000000
mb
execute astro_orf1b.nex
prset aamodelpr = mixed
lset rates=invgamma
showmodel
showparams
help mcmc
mcmc ngen = 1000000
sump
sumt

cp astro_orf1b.nex.con.tre temp1
egrep "@|£" temp1
#-> Nothing
foreach line (`cat ../nameskey.txt`)
  set oldname = `echo $line | awk -F: '{print $1}'`
  set newname = `echo $line | awk -F: '{print $2}'`
  echo $oldname $newname
  sed 's/'$oldname'/'$newname'/g' temp1 > temp2
  mv temp2 temp1
end
mv temp1 muscle2.nex.con.tre

java -jar ~/PACKAGES/FigTree_v1.4.2/lib/figtree.jar
#appearance: line weight 1
#tree: midpoint root
#tip labels: font size 5, liberation sans
#node labels: font size 5, 2 s.f., display prob, liberation sans
#scale bar: line weight 1, font size 5

#-> export as SVG
sed 's/@/ /g' muscle2.nex.con.tre.svg | sed 's/£/\//g' \
  > muscle3.nex.con.tre.svg

inkscape muscle3.nex.con.tre.svg
#Delete all node prob that are close to the tree leaves.
#-> save as eps. Also need to set 'output page size' to 'use exported
#  object's size' not 'use document's page size'.

#------------------------------------------------------------------------------
#CP tree using the same sequences we used for the ORF1b tree.

cd Bioinformatics_Jul2018/FullGenome/NR/GB
mkdir CP_tree_for_fullgenome_seqs
cd CP_tree_for_fullgenome_seqs

#Valid CP N-term peptides (verified by position wrt sgRNA promoter sequence)
mkdir RefSeqs
cd RefSeqs
cp ../../../VIRAD/*/*.gbf .
rm -f log1.txt
foreach acc (`ls *.gbf | sed 's/\.gbf//'`)
  @ c1 = `grep "#CDS " $acc.gbf | awk '{if (0+$7>0+6000) print $6}'`
  @ c2 = `grep "#CDS " $acc.gbf | awk '{if (0+$7>0+6000) print $6+149}'`
  grep -v "^#" $acc.gbf > temp1
  transeq temp1 temp2 -regions $c1..$c2 -auto
  cat temp2 >> log1.txt
  rm -f temp1 temp2
end
muscle -in log1.txt -out CP_Nterm.aln -clwstrict
cat CP_Nterm.aln
#  LC047798       ---MPRNRNGRQNRNTTRVVVQSTSGNQQAAQTPRTRRARRPTVNVNVNTRPQ-------
#  NC_004579      --MASANQAAKAEAKKVIEKVAKEVIKETK------NSAQRNQGPGKRWNSKKGRHMP--
#  NC_002469      ----MAEKPQQKAVASAAKQLAKEVVKLDKITKSNGKQHPQKNVPARKWRPRQA------
#  NC_013060      ---MAGRQPQQALPKAAAKQIAKEVVKQEK--KEPVVRKKKQFYPNPKFNNRFNK-----
#  NC_019494      --MAGDKLNASAKAAPFAKEVAKEVVKEEK--KTQ-ARRRKWYKPRRQQNQPQQQ-----
#  NC_011400      ---MANASKGVTVNINNAKRKPRFTNNQRA------RSTRPNFTPAPKFRKRRFIPNRN-
#  NC_026814      --MASKPGKDVTVEVKTSGTKST--SSRSK------SRGRNRNVKITVNSQPKTNRRRRN
#  NC_001943      --MASKSNKQVTVEVSNNGRSRSKSRARSQ------SRGRDKSVKITVNSRNRARRQP--
#  JN420356       --MASKSGNEVTVKVDSGRSRSK-SRGRSK------SRGRSKDVKITVNSKPKKQRRSG-
#  KT946734       MAGDNVAGSGGSGPRAAVGCNPR--RRRRR------ARQRGAKGPNGVGTVRSQPPPM--
#  NC_023631      MASRKQQNRRATRNTTNIVVRNGAAANQAG------AAGGQRRRRSRRNNKAPQVN----
#  NC_023674      --MANRQQKRGPRTTTNIVVRNGTAAPQAR------ASGSTAGSRRRRNRARRQPQVN--
#  NC_023675      --MANTKNNVQPQVVTTTTTVSN--RRGRR------RRRRANRIPPATTTTVRRTKITKP
#  NC_018702      --MAKAKQQQKNATTVTTTTVTGRSSRRSR------RRSVRRRAAGPSNPPTKTTTVR--
#  NC_023636      --MANRRNRRPPQRRRMRGPRPA--AETSTVTTTVSTNGQTK--PTTKTTTVKTTK----

#221 sequences
#Note a lot will have the incorrect CP start site annotated so will need to
#  truncate.
degapseq ../CP/clustalw2.aln t1
sed 's/_1//' t1 > t2
seqretsplit t2 t3
foreach lname (`ls *fasta | sed 's/\.fasta//'`)
  set uname = `echo $lname | sed 's/[a-z]/\U&/g'`
  mv $lname.fasta $uname.fasta
end
wc -l ../ORF1b/allpairs.txt
#-> 221
mkdir FASTA1
foreach name (`awk '{print $1}' ../ORF1b/allpairs.txt`)
  mv $name.fasta FASTA1
end
ls *fasta | wc -l
#-> 9
rm *.fasta
ls FASTA1/*fasta | wc -l
#-> 221
cat FASTA1/*fasta > t1
muscle -in t1 -out CP_Nterminally_truncated_to_correct_start.aln -clwstrict
emacs CP_Nterminally_truncated_to_correct_start.aln &
#-> Truncate to proper CP Nterm - use the Refseqs above for guidance.

#Split and redo alignment to refine it.
mkdir FASTA2
cd FASTA2
degapseq ../CP_Nterminally_truncated_to_correct_start.aln t1
seqretsplit t1 t2
foreach lname (`ls *fasta | sed 's/\.fasta//'`)
  set uname = `echo $lname | sed 's/[a-z]/\U&/g'`
  mv $lname.fasta $uname.fasta
end
ls *fasta | wc -l
#-> 221
cd ..
cat FASTA2/*fasta > t1
muscle -in t1 -out CP_refined.aln -clwstrict
emacs CP_refined.aln &

seqret CP_refined.aln astro_CP.nex -osformat nexus

ln -s ../ORF1b/nameskey.txt

mkdir R1000000
cp astro_CP.nex R1000000
cd R1000000
mb
execute astro_CP.nex
prset aamodelpr = mixed
lset rates=invgamma
showmodel
showparams
help mcmc
mcmc ngen = 1000000
sump
sumt

cp astro_CP.nex.con.tre temp1
egrep "@|£" temp1
#-> Nothing
foreach line (`cat ../nameskey.txt`)
  set oldname = `echo $line | awk -F: '{print $1}'`
  set newname = `echo $line | awk -F: '{print $2}'`
  echo $oldname $newname
  sed 's/'$oldname'/'$newname'/g' temp1 > temp2
  mv temp2 temp1
end
mv temp1 muscle2.nex.con.tre

java -jar ~/PACKAGES/FigTree_v1.4.2/lib/figtree.jar
#appearance: line weight 1
#tree: midpoint root
#tip labels: font size 5, liberation sans
#node labels: font size 5, 2 s.f., display prob, liberation sans
#scale bar: line weight 1, font size 5

#-> export as SVG
sed 's/@/ /g' muscle2.nex.con.tre.svg | sed 's/£/\//g' \
  > muscle3.nex.con.tre.svg

inkscape muscle3.nex.con.tre.svg
#Delete all node prob that are close to the tree leaves.
#-> save as eps. Also need to set 'output page size' to 'use exported
#  object's size' not 'use document's page size'.

#Colour clades in inkscape 4/3/19:
#44AA00 - green - G-I
#D400AA (crimson) G-VI
#0000D4 (dark blue) MLBs
#D4AA00 (yellow) G-II
#D40000 (red) G-III
#
#FF0000 (red) ORFx bars

#------------------------------------------------------------------------------

#Check the following refseqs have the correct CP start codon (questionable -
#  due to extended overlap with ORF1b). 
#-> Find the sgRNA promoter to prove it.

#NC_023631       GCCTAAAA--------------TTACTGATGAGCAGCTGGATCATCTTTGGAGGGGAGGA
#NC_023636       GGCTCTGATCTGCCAGTC---TTTTCAGATAGGATCTTGTCGTATCTTTGGGGGGGAGGA
#NC_019494       GGTGATGA---GCCACTCCGCTTTACTGATGAGATGCTTGACCGGCTTTGGGGGGGCGGA
#NC_001943       GACTCTGGCCTTCCAGCCAGACTCACAGAAGAGCAACTCCATCGCATTTGGAGGGGAGGA
#NC_026814       GATTCTGGCTTACCTACCAGATTCACAGAAGAGCAAATGCATCGCATATGGAGGGGAGGA
#                *                     *  * **   *    *        * *** **** ***
#                  sgRNA         CP
#                  v             v
#NC_023631       CCAAAGACAAATCCTAATGGCTAGCCGCAAACAGCAAAACAGGCGTGCCA----------
#NC_023636       CCAAAAAGAGA---CGATGGCCA------ATCGGCGTAACAG------------------
#NC_019494       CCAAAGGTCGG---ATATGGCTG---GTGATAAGCTCAACGCATCTGCCAAAGCGGCACC
#NC_001943       CCAAAGAAGTG---TGATGGCTA---GCA---AGTCCAACAAGCAAGTAA----------
#NC_026814       CCAAAAAATTG---CGATGGCTA---GCA---AGCCAGGCAAAGATGTCA----------
#                *****           *****            *     *                    

#So all these gb files need coordinates updating:
#  NC_023631.gb                  3908..6226 -> 3956..6226
#  NC_023636.gb                  4037..6382 -> 4175..6382
#  NC_019494.gb                  3895..6348 -> 4081..6348
#  NC_026814.gb                  4015..6519 -> 4195..6519

#  NC_026814.gb also needs ORFx coords updating: 4169..4492 -> 4217..4492

#Sanity check - extract 30 nt upstream of every gb annotated CP (which should
#  always be the 2nd longest ORF?)
cd NR/VIRAD
mkdir Check_CP_start_site
cd Check_CP_start_site
cp ../G*/BKP/*.gbf .
#For the "cp: will not overwrite just-created" ones, diff the two copies to
#  check they are identical.
foreach acc (`ls *.gbf | sed 's/\.gbf//'`)
  grep "#CDS " $acc.gbf | awk '{if (0+$7>0+6000) print $0}' | wc -l | \
    awk '{print $1}'
end
#-> All are 1
rm -f log1.txt
foreach acc (`ls *.gbf | sed 's/\.gbf//'`)
  @ c1 = `grep "#CDS " $acc.gbf | awk '{if (0+$7>0+6000) print $6-30}'`
  @ c2 = `grep "#CDS " $acc.gbf | awk '{if (0+$7>0+6000) print $6+2}'`
  grep -v "^#" $acc.gbf > temp1
  extractseq temp1 temp2 -regions $c1..$c2 -auto
  cat temp2 >> log1.txt
  rm -f temp1 temp2
end
muscle -in log1.txt -out sgRNA_CP.aln -clwstrict
emacs sgRNA_CP.aln &
#-> All look OK now.
#  KT946734    -CACTATGGTGGGGTGGACCAAA------AAAGTGTGATG
#  NC_019494   -GGCTTTGGGGGGGCGGACCAAA------GGTCGGATATG
#  NC_023636   -ATCTTTGGGGGGGAGGACCAAA------AAGAGACGATG
#  NC_011400   -ACATTTGGAGGGGCGGACCAAA------TGATGACTATG <- 3 CP-frame AUGs here
#  NC_026814   -GCATATGGAGGGGAGGACCAAA------AAATTGCGATG
#  JN420356    -GCATTTGGAGGGGAGGACCAAA------AGATTGCGATG
#  NC_001943   -GCATTTGGAGGGGAGGACCAAA------GAAGTGTGATG
#  NC_013060   -GGCTTTGGAGGGGAGGTCCAAA------GCAAAGTCATG
#  NC_023631   ----TTTGGAGGGGAGGACCAAAGAC---AAATCCTAATG
#  NC_023674   ----TTTGGAGGGGAGGACCAAAGAC---AGCATCTAATG
#  NC_004579   AAACTTTGGAGGGGAGGACCAAA-------GTGAGCTATG
#  NC_018702   -AGCTTTGGAGGGGAGGACCAAA------AGCTGTTCATG
#  LC047798    -------GGAGGGGCGGACCAAAGATGTCTAATCATAATG
#  NC_002469   -AGCTTTGGAGGGGCGGACCAAA------GTTTGATTATG
#  NC_023675   ----TTTGGAGGGGCGGACCAAAGCA---TAAGCCTAATG
#                     ** **** ** *****              ***

#------------------------------------------------------------------------------

#ORFx pepseqs in all 221 ORF1b-selected sequences -> colour alignment figures.
cd Bioinformatics_Jul2018/FullGenome/NR/GB/CP_tree_for_fullgenome_seqs
mkdir ORFX_PEPSEQS
cd ORFX_PEPSEQS
#Order of sequences in ORF1b tree:
cp ../../ORF1b/R1000000/muscle3.nex.con.tre5a.svg temp1.svg
#If you open as is in emacs, it will display as a figure not text. So remove
#  the first 100 lines.
tail -n +100 temp1.svg > temp2
#Now look for where the accession numbers start and end and delete everything
#  above and below.
#-> Hmm, unfortunately it doesn't follow the order on the tree image. But if we
#  can find the coords we should be able to sort by that.
grep transform temp2 | wc -l
#-> 221  
cp temp2 svg_truncated_to_accession_numbers.txt

awk -v ORS="@" '{print $0}' svg_truncated_to_accession_numbers.txt | \
  sed 's/<g@/\n/g' | grep "@" | wc -l
#-> 221
awk -v ORS="@" '{print $0}' svg_truncated_to_accession_numbers.txt | \
  sed 's/<g@/\n/g' | grep "@" | grep "£"
#-> none
awk -v ORS="@" '{print $0}' svg_truncated_to_accession_numbers.txt | \
  sed 's/<g@/\n/g' | grep "@" | sed 's/text/£/g' | sed 's/£[^£]*$//' | \
  sed 's/@.*£//' | sed 's/ *transform="matrix(1,0,0,1,//g' | sed 's/,/ /' | \
  sed 's/)/ /' | sort -n -k 2 | sed 's/<tspan[^>]*>//' | \
  sed 's/<\/tspan>//' | sed 's/<\///' | sed 's/.*>//' \
  > tree_ordered_accessions.txt
wc -l tree_ordered_accessions.txt

#Get corrected CP pepseqs (i.e. have inspected to find the correct CP start
#  codon for all these).
mkdir CP_correct_pepseqs
cd CP_correct_pepseqs
degapseq ../../CP_Nterminally_truncated_to_correct_start.aln t1
seqretsplit t1 t2
cd ..

#Check can access gb sequences for each CP peptide
foreach i (`awk '{print $1}' tree_ordered_accessions.txt`)
  if (! -r ../../GB/$i.gb) then
    echo "Can't find $i.gb"
  endif
end

#Find CP start coord in genome for each corrected CP pepseq.
rm -f CP_correct_start_coords.txt
foreach acc (`awk '{print $1}' tree_ordered_accessions.txt`)
  set lacc = `echo $acc | sed 's/[A-Z]/\L&/g'`
  set CPpepstart = `tail -n +2 CP_correct_pepseqs/$lacc.fasta | head -2 | \
    awk -v ORS="" '{print $0}'`
  seqret ../../GB/$acc.gb temp1 -auto
  transeq temp1 temp2.1 -frame 1 -auto
  transeq temp1 temp2.2 -frame 2 -auto
  transeq temp1 temp2.3 -frame 3 -auto
  tail -n +2 temp2.1 | awk -v ORS="" '{print $0}' | grep "." | \
    awk '{print "1",$0}' > temp3
  tail -n +2 temp2.2 | awk -v ORS="" '{print $0}' | grep "." | \
    awk '{print "2",$0}' >> temp3
  tail -n +2 temp2.3 | awk -v ORS="" '{print $0}' | grep "." | \
    awk '{print "3",$0}' >> temp3
  set CPstart = `grep $CPpepstart temp3 | sed 's/'$CPpepstart'.*//' | \
    awk '{print length($2)*3+$1}'`
  echo $acc $CPstart >> CP_correct_start_coords.txt
  rm -f temp1 temp2.[123] temp3
end
wc -l CP_correct_start_coords.txt

#Translate first 600 nt of CP in +1 frame. Find first AUG-initiated ORF, start
#  coord wrt CP start, pepseq. Flag up with 0 if no AUG found in the 600 nt.
rm -f orfx.txt
foreach line (`awk '{printf "%s:%s\n",$1,$2}' CP_correct_start_coords.txt`)
  set acc = `echo $line | awk -F: '{print $1}'`
  set c1 = `echo $line | awk -F: '{print $2+1}'`
  set c2 = `echo $line | awk -F: '{print $2+1200}'`
  transeq ../../GB/$acc.gb temp1 -regions $c1..$c2 -auto
  set anyM = `tail -n +2 temp1 | awk -v ORS="" '{print $0}' | grep M | wc -l`
  if ($anyM) then
    tail -n +2 temp1 | awk -v ORS="" '{print $0}' | sed 's/M/ M/' | \
      awk '{print 3*length($1),$2}' | sed 's/\*.*/*/' >> orfx.txt
  else
    echo $acc 0 >> orfx.txt
  endif
  rm -f temp1
end
wc -l orfx.txt
#-> 221
awk '{if (0==$1) print $0}' orfx.txt
#-> Every single sequence had at least one +1-frame AUG in the first 600 nt of
#  CP.

sort -n orfx.txt
#-> AUGs up to 99 nt from CP AUG look like they can be real ORFx. After 99 it
#  jumps to short ORFs at 147+ nt, and nothing over 22 codons in length till
#  306 nt.

awk '{print $1}' tree_ordered_accessions.txt > temp1
paste temp1 orfx.txt tree_ordered_accessions.txt | sed 's/\t/ /' \
  > orfx_named.txt

cp orfx_named.txt orfx_named_annotated.txt
emacs orfx_named_annotated.txt &
#Separate into phylogenetic clades, and delete off those which don't have ORFx.

#These MLBs just come up with a very short ORF, so presumably they have a
#  second AUG.
#  JX857870 6 MPIRV*	JX857870 - Astrovirus MLB3
#  KX273058 6 MPIRV*	KX273058 - Primate astrovirus
#  AB829252 6 MPIRV*	AB829252 - Astrovirus MLB2
#  KT224358 6 MPIRV*	KT224358 - Astrovirus MLB2
#  JF742759 6 MPIRV*	JF742759 - Astrovirus MLB2
#  KX022687 6 MPIRV*	KX022687 - Astrovirus MLB2
#These also have a short ORF and need to check whether or not they have a
#  second AUG
#  KP663426 324 MLPPTPNGGSSTLS*	KP663426 - Astrovirus Er/SZAL6/HUN/2011 
#  KT946731 516 MQT*	        KT946731 - Rodent astrovirus 

foreach acc (JX857870 KX273058 AB829252 KT224358 JF742759 KX022687 KP663426 \
  KT946731)
  echo $acc
  set c1 = `grep -w $acc CP_correct_start_coords.txt | awk '{print $2+1}'`
  set c2 = `echo $c1 | awk -F: '{print $1+600}'`
  transeq ../../GB/$acc.gb temp1 -regions $c1..$c2 -auto
  tail -n +2 temp1 | awk -v ORS="" '{print $0}' | grep M
  rm -f temp1
end
#-> Just copy and paste, where relevant, the longer ORFx pepseqs to 
#  orfx_named_annotated.txt
#-> KP663426 really does have no ORFx.
#-> KT946731 looks like might have ORFx but has lost initiation site
#-> The MLBs are all fine - uses the 3rd AUG on the message actually - 3
#  closely spaced AUGs. Maybe, like murine noro, a short 5'UTR facilitates
#  leaky scanning. Add +33 to the original start coord: 6 -> 39.


#Cp to word document. Replace non singleton groups with muscle alignments.
#  (But need to reorder back to the tree order - i.e. the order of the input
#  sequences - write a quick script for this).
#
#Paste a sequence block to temp1
sed 's/\*$//' temp1 | awk '{printf ">%s\n%s\n",$1,$2}' > temp2
muscle -in temp2 -out temp3 -clwstrict
seqret temp3 temp4 -auto
awk -v ORS="@" '{print $1}' temp4 | sed 's/>/\n/g' | grep "@" | \
  sed 's/@/ /' | sed 's/@//g' > temp5
rm -f temp6
foreach i (`awk '{print $1}' temp1`)
  grep -w $i temp5 >> temp6
end
awk '{printf "%-8s %s\n",$1,$2}' temp6 > temp7
emacs temp7 &
#Paste temp6 back into the word document
rm -f temp[1234567]

#Then colour code amino acids.
#Check for and annotate TMs.

#-> orfx_pepseqs_??.odt

#------------------------------------------------------------------------------

#Pepstats by group -> Table - orfx_pepstats_??.odt

mkdir PepStats
cd PepStats
#-> cp from orfx_pepseqs_??.odt file to *.txt

ls *.txt > alns1
#-> Actually reorder this so it is in the same order as in orfx_pepseqs_??.odt
#  (which is the same order as the tree).

#Rememember to degap seqs just in case

set qdir = `pwd`
cd $qdir
foreach i (`cat alns1 | sed 's/\.txt//'`)
  echo $i
  rm -f $i.log
  awk '{printf ">%s\n%s\n",$1,$2}' $i.txt > temp1
  degapseq temp1 temp2 -auto
  mkdir $qdir/T
  cd $qdir/T
  seqretsplit ../temp2 t1 -auto
  foreach j (`ls *.fasta`)
    pepstats $j temp3 -auto
    set acc = `head -1 temp3 | awk '{print $3}'`
    set mass = `grep "Molecular weight" temp3 | awk '{printf "%.2f",$4/1000}'`
    set pI = `grep "Isoelectric Point" temp3 | awk '{printf "%.2f",$4}'`
    set length = `grep "Residues" temp3 | head -1 | awk '{print $7}'`
    echo $acc $mass $pI $length >> ../$i.log
    rm -f temp3
  end
  cd $qdir
  rm -rf $qdir/T
  rm -f temp1 temp2
end
#acc mass pI length
  
rm -f ranges.txt medians.txt
foreach i (`cat alns1 | sed 's/\.txt//'`)
  @ n = `wc -l $i.log | awk '{print $1}'`
  set minmass = `awk '{print $2}' $i.log | sort -n | head -1`
  set maxmass = `awk '{print $2}' $i.log | sort -nr | head -1`
  set minpI = `awk '{print $3}' $i.log | sort -n | head -1`
  set maxpI = `awk '{print $3}' $i.log | sort -nr | head -1`
  set minlength = `awk '{print $4}' $i.log | sort -n | head -1`
  set maxlength = `awk '{print $4}' $i.log | sort -nr | head -1`
  echo $n $minmass-$maxmass $minpI-$maxpI $minlength-$maxlength $i \
    >> ranges.txt
  @ ntest = `echo $n | awk '{print int($1/2)}'`
  if ($n == 2 * $ntest) then
    @ ntest += 1
    set medianmass = `awk '{print $2}' $i.log | sort -n | head -$ntest | \
      tail -2 | awk '{s+=$1}END{print s/2}'`
    set medianpI = `awk '{print $3}' $i.log | sort -n | head -$ntest | \
      tail -2 | awk '{s+=$1}END{print s/2}'`
    set medianlength = `awk '{print $4}' $i.log | sort -n | head -$ntest | \
      tail -2 | awk '{s+=$1}END{print s/2}'`
  else
    @ ntest += 1
    set medianmass = `awk '{print $2}' $i.log | sort -n | head -$ntest | \
      tail -1`
    set medianpI = `awk '{print $3}' $i.log | sort -n | head -$ntest | \
      tail -1`
    set medianlength = `awk '{print $4}' $i.log | sort -n | head -$ntest | \
      tail -1`
  endif
  echo $n $medianmass $medianpI $medianlength $i | \
    awk '{printf "%3s%7.1f%7.1f%5.0f  %s\n",$1,$2,$3,$4,$5}' >> medians.txt
end
#-> range and median of mass, pI and length.

#  *    *    *    *    *    *    *    *    *    *    *    *    *    *    *

#Scatter box plot of individual lengths by group.

#Make a box whisker plot of ORFx lengths for the different groups.

cd PepStats
#Make a file with abbreviated names:
awk '{print $1}' alns1 > alns2
emacs alns2 &
#-> Add abbreviated name as $1.

rm -f length_by_group.txt
foreach line (`awk '{printf "%s:%s\n",$1,$2}' alns2`)
  set lname = `echo $line | awk -F: '{print $2}' | sed 's/\.txt//'`
  set sname = `echo $line | awk -F: '{print $1}'`
  awk '{printf "@%s@ %s\n","'"$sname"'",$4}' $lname.log | sed 's/@/"/g' \
    >> length_by_group.txt
end

#  boxplot(x, ..., range = 1.5, width = NULL, varwidth = FALSE,
#             notch = FALSE, outline = TRUE, names, plot = TRUE,
#             border = par("fg"), col = NULL, log = "",
#             pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
#             horizontal = FALSE, add = FALSE, at = NULL)
    
# range: this determines how far the plot whiskers extend out from the
#          box.  If ‘range’ is positive, the whiskers extend to the most
#          extreme data point which is no more than ‘range’ times the
#          interquartile range from the box. A value of zero causes the
#          whiskers to extend to the data extremes.

#box = interquartiles [confirmed by trial and error]
#center line = median not mean [confirmed by trial and error]
#circles = all outliers [confirmed by trial and error]

R
install.packages("ggplot2")
help(boxplot)  # Lots of parameters can alter here

R
#pdf("BoxPlot3.pdf",height=2.0,width=2.5,pointsize=12)
svg("BoxPlot3.svg",height=2.0,width=2.5,pointsize=12)
counts = read.table("length_by_group.txt")
boxplot(counts$V2~counts$V1,las=2,ylim=c(0,170),ylab="XP length (aa)",par(mar=c(4.8,4,0.2,0.1)),par(cex=0.7))
dev.off()

#With n numbers added.
R
#pdf("BoxPlot3b.pdf",height=2.0,width=2.5,pointsize=12)
svg("BoxPlot3b.svg",height=2.0,width=2.5,pointsize=11)
counts = read.table("length_by_group.txt")
boxplot(counts$V2~counts$V1,las=2,ylim=c(0,180),ylab="XP length (aa)",par(mar=c(4.8,4,0.2,0.1)),xlab="",par(cex=0.7))
#x=boxplot(counts$V2~counts$V1)
#x$n
#as.character(x$n)
#Be sure to update this line if input data set changes:!!!!
text(c(1:15),178,labels=c("1","4","39","10","5","6","17","9","19","24","25","13","4","4","2"),srt=90,adj=c(1,0.5))
dev.off()

#------------------------------------------------------------------------------
