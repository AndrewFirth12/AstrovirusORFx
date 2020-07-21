#!/bin/bash

#Bootstrap resampling of phasing patterns of reads in ORFx and reads in the
#  region of CP not overlapped by ORFX. 

#See the RiboSeqManual to get the input vRNA.bowtie files.

#------------------------------------------------------------------------------

#For the RiboSeq stats, let x be the number of codons used in the dual coding
#  region. Then bootstrap resample with replacement x codons per sample for the
#  dual-coding region. Also bootstrap resample with replacement x codons from
#  the single-coding region used. No point sampling en masse from host mRNAs,
#  because expression levels are different. For each resampling, we solve the
#  linear equations to get the estimated +0 and +1 frame translation levels and
#  plot two distributions - one for the dual- and one for the single-coding
#  region-derived resamplings.

#The astrovirus dual coding region is 4374..4712 = 112 aa + stop codon.
#We are looking at CP-frame so we take the regions 4376..4711 (CP-overlap) and
#  4715..6691 (CP-nonoverlap).
#The script avoids regions close to the initiation and termination codons.
#  Specifically, 5' end of read must be >= CDSstart and <= CDSend - 30

#So let's first check the native (i.e. not resampled) stats for the two
#  libraries:
cd Phasing-CP-overlap
set c1 = 4376
set c2 = 4711
foreach library (35_S30 36_S15)
  echo -n "$library "
  awk '{if ($3=="+") print $5+1}' $library/vRNA.bowtie | \
    awk '{if ($1+0>="'"$c1"'"+0&&$1+0<="'"$c2"'"-30) print $1-"'"$c1"'"+1}' | \
    awk '{print $1%3}' | sort -n | uniq -c > temp1
  set s = `awk '{s+=$1}END{print s}' temp1`
  awk '{printf "%s x%sx\n",$1,$2}' temp1 | sed 's/x0x/x3x/' | sed 's/x//g' | \
    awk '{print $2,$1/"'"$s"'"}' | sort -n | awk '{printf "%s ",$0}'
  rm -f temp1
  echo
end  
#BTW, I've checked, and this does select a region that is a multiple of 3.
#  Namely reads with 5' ends mapping from 4376 to 4681 (102 codons).
#
#  35_S30 1 0.607436 2 0.245493 3 0.147071
#  36_S15 1 0.603355 2 0.250242 3 0.146404


#  *    *    *    *    *    *    *    *    *    *  

#OK, now we need to bootstrap - resample 102 codons with replacement from the
#  buffered overlap region. Also for the buffered CP-nonoverlap region.

#Extract the nt-by-nt counts for the two buffered regions
cd Phasing-CP-overlap
set c1 = 4376
set c2 = 4711
set d1 = 4715
set d2 = 6691
foreach library (35_S30 36_S15)
  awk '{if ($3=="+") print $5+1}' $library/vRNA.bowtie | \
    awk '{if ($1+0>="'"$c1"'"+0&&$1+0<="'"$c2"'"-30) print $1-"'"$c1"'"+1}' | \
    sort -n | uniq -c | awk '{print $0,$2%3}' > $library.CP-overlap.txt
  awk '{if ($3=="+") print $5+1}' $library/vRNA.bowtie | \
    awk '{if ($1+0>="'"$d1"'"+0&&$1+0<="'"$d2"'"-30) print $1-"'"$d1"'"+1}' | \
    sort -n | uniq -c | awk '{print $0,$2%3}' > $library.CP-nonoverlap.txt
end  
#Check every nt is present (i.e. had a non-zero value); else need to pad with
#  zeros.
#c2 - c1 + 1 - 30 = 306
#d2 - d1 + 1 - 30 = 1947
tail -1 35_S30.CP-overlap.txt
tail -1 36_S15.CP-overlap.txt
#-> 306 306
tail -1 35_S30.CP-nonoverlap.txt
tail -1 36_S15.CP-nonoverlap.txt
#-> 1947 1947
wc -l *.CP-overlap.txt 
wc -l *.CP-nonoverlap.txt 
#-> 306 306 1945 1946
#-> Manually pad the latter two
awk '{print $2}' 36_S15.CP-nonoverlap.txt > t1
tail -n +2 t1 > t2
paste t1 t2 | awk '{if (1!=$2-$1) print $0}'
#-> 226	228
awk '{print $2}' 35_S30.CP-nonoverlap.txt > t1
tail -n +2 t1 > t2
paste t1 t2 | awk '{if (1!=$2-$1) print $0}'
#-> 898	900
#   1787 1789
#Check
awk '{print $3}' 36_S15.CP-overlap.txt | sort | uniq -c
awk '{print $3}' 35_S30.CP-overlap.txt | sort | uniq -c
awk '{print $3}' 36_S15.CP-nonoverlap.txt | sort | uniq -c
awk '{print $3}' 35_S30.CP-nonoverlap.txt | sort | uniq -c

#Now put each codon on one line for resampling by codon.
foreach file (36_S15.CP-overlap 35_S30.CP-overlap 36_S15.CP-nonoverlap \
  35_S30.CP-nonoverlap)
  awk '{if (0==$3) {print $1} else {printf "%s ",$1}}' $file.txt > $file.codons
end
#Check
cat 36_S15.CP-overlap.codons 35_S30.CP-overlap.codons \
  36_S15.CP-nonoverlap.codons 35_S30.CP-nonoverlap.codons | \
  awk '{print NF}' | uniq -c
wc -l *codons

#Observed values
foreach file (35_S30.CP-overlap 36_S15.CP-overlap 35_S30.CP-nonoverlap \
  36_S15.CP-nonoverlap)
  echo -n "$file "
  awk '{s1+=$1}{s2+=$2}{s3+=$3}END{print s1,s2,s3,s1+s2+s3}' $file.codons | \
    awk '{print $1/$4,$2/$4,$3/$4}'
end
#-> Again, these agree with the previous overlap values.

#Now bootstrap resample sets of 102 from each of the 4 files. For each
#  resampling, need to save the s1,s2,s3 values.
#Let's do 1000 samplings of each of the 4 files. Easiest way will be a C++
#  program that reads in the 4 arrays -> bootstrap.cxx.

#./bootstrap seed number_of_codons_to_sample_with_replacement
#  number_of_randomizations

@ seed = 351
foreach file (35_S30.CP-overlap 36_S15.CP-overlap 35_S30.CP-nonoverlap \
  36_S15.CP-nonoverlap)
  @ seed += 1
  cp $file.codons inputstats.txt
  ./bootstrap $seed 102 1000 > $file.bootstraps1
end
foreach file (35_S30.CP-overlap 36_S15.CP-overlap 35_S30.CP-nonoverlap \
  36_S15.CP-nonoverlap)
  awk '{print $1,$2,$3,$1+$2+$3}' $file.bootstraps1 | \
    awk '{print $1/$4,$2/$4,$3/$4}' > $file.bootstraps2
end

#Now we need to solve the linear equations for each.
#  35_S30 NT1     
#  36_S15 NT2     

#NT1 host                 0.7316049 0.1221675 0.1462275
#NT2 host                 0.7412605 0.1203969 0.1383426

#NT1 virus CP/XP overlap  0.6074359 0.2454927 0.1470713  observed
#NT2 virus CP/XP overlap  0.6033547 0.2502417 0.1464035  observed

#Solve:
#cp * (0.7316049, 0.1221675, 0.1462275) + xp * (0.1462275, 0.7316049, 0.1221675) = (0.6074359, 0.2454927, 0.1470713)
#cp * (0.7412605, 0.1203969, 0.1383426) + xp * (0.1383426, 0.7412605, 0.1203969) = (0.6033547, 0.2502417, 0.1464035)

#Rep 1

#a = c(0.7316049, 0.1221675, 0.1462275)
#b = c(0.1462275, 0.7316049, 0.1221675)
#c = c(0.6074359, 0.2454927, 0.1470713)

#Now enter cp = 1, 0.9, 0.8 etc and chose cp to minimize q
#q = cp * a + (1 - cp) * b - c
#sqrt(q[1]^2 + q[2]^2 + q[3]^2)

#Want to minimize wrt cp:
#sqrt(
#(cp * a[1] + (1 - cp) * b[1] - c[1])^2 +
#(cp * a[2] + (1 - cp) * b[2] - c[2])^2 +
#(cp * a[3] + (1 - cp) * b[3] - c[3])^2) 
#So need to take derivative wrt CP and set that to zero.

#Minimize the sqrt is equivalent to minimizing the non-sqrt value.

#cp    q
#1     0.175
#0.9   0.091
#0.85  0.046
#0.082 0.024
#0.081 0.016
#0.8   0.009
#0.79  0.0076  <- So best estimate is 79% CP, 21% XP
#0.78  0.013
#0.75  0.037
#0.7   0.079


#Rep 2

#a = c(0.7412605, 0.1203969, 0.1383426)
#b = c(0.1383426, 0.7412605, 0.1203969)
#c = c(0.6033547, 0.2502417, 0.1464035)

#Now enter cp = 1, 0.9, 0.8 etc and chose cp to minimize q
#q = cp * a + (1 - cp) * b - c
#sqrt(q[1]^2 + q[2]^2 + q[3]^2)

#cp    q
#1     0.1896
#0.9   0.1035
#0.85  0.0610
#0.082 0.0363
#0.8   0.0216
#0.79  0.0164
#0.78  0.0147 <-So best estimate is 78% CP, 22% XP
#0.77  0.0178
#0.75  0.0311
#0.7   0.0722

#OK. This can actually be solved analytically :-) Basically we want to minimize
#  sqrt(q[1]^2 + q[2]^2 + q[3]^2) which is the same as minimizing
#  q[1]^2 + q[2]^2 + q[3]^2. To do this we just need to solve dq/dx = 0 for
#  x (where x = cp).
#q = Sum_i=1_3 (x a_i + (1 - x) b_i - c_i)^2
#dq/dx = Sum_i=1_3 2 (x a_i + (1 - x) b_i - c_i) (a_i - b_i)
#dq/dx = 0 => Sum_i=1_3 (x a_i + (1 - x) b_i - c_i) (a_i - b_i) = 0
#          => Sum_i=1_3 (x a_i - x b_i)(a_i - b_i) + (b_i - c_i)(a_i - b_i) = 0
#          => x Sum_i=1_3 (a_i - b_i)^2 = Sum_i=1_3 (c_i - b_i)(a_i - b_i)
#          => x = Sum_i=1_3 (c_i - b_i)(a_i - b_i) / Sum_i=1_3 (a_i - b_i)^2
#If we plug these values into R, then I can find the exact solutions above
#  namely 0.7931541 for Rep 1 and 0.7816411 for Rep 2 giving minima of
#  0.007132978 and 0.01467661 respectively.

#So now we can convert *.bootstraps1 tables to the x = cp estimates.



R

#Rep 1 = 35_S30
c = read.table("35_S30.CP-overlap.bootstraps2")
a0 = c(0.7316049, 0.1221675, 0.1462275)
b0 = c(0.1462275, 0.7316049, 0.1221675)
a = as.data.frame(cbind(rep(a0[1],nrow(c)),rep(a0[2],nrow(c)),rep(a0[3],nrow(c))))
b = as.data.frame(cbind(rep(b0[1],nrow(c)),rep(b0[2],nrow(c)),rep(b0[3],nrow(c))))
num = (c - b) * (a - b) 
den = (a - b) * (a - b)
cp = (num$V1 + num$V2 + num$V3) / (den$V1 + den$V2 + den$V3)
write(cp,"35_S30.CP-overlap.bootstraps3")
#Temporarily appended 0.6074359, 0.2454927, 0.1470713 at the end of
#  35_S30.CP-overlap.bootstraps2 to check it recovered  0.7931541 -> OK

#Rep 1 = 35_S30
c = read.table("35_S30.CP-nonoverlap.bootstraps2")
a0 = c(0.7316049, 0.1221675, 0.1462275)
b0 = c(0.1462275, 0.7316049, 0.1221675)
a = as.data.frame(cbind(rep(a0[1],nrow(c)),rep(a0[2],nrow(c)),rep(a0[3],nrow(c))))
b = as.data.frame(cbind(rep(b0[1],nrow(c)),rep(b0[2],nrow(c)),rep(b0[3],nrow(c))))
num = (c - b) * (a - b) 
den = (a - b) * (a - b)
cp = (num$V1 + num$V2 + num$V3) / (den$V1 + den$V2 + den$V3)
write(cp,"35_S30.CP-nonoverlap.bootstraps3")

#Rep 2 = 36_S15
c = read.table("36_S15.CP-overlap.bootstraps2")
a0 = c(0.7412605, 0.1203969, 0.1383426)
b0 = c(0.1383426, 0.7412605, 0.1203969)
a = as.data.frame(cbind(rep(a0[1],nrow(c)),rep(a0[2],nrow(c)),rep(a0[3],nrow(c))))
b = as.data.frame(cbind(rep(b0[1],nrow(c)),rep(b0[2],nrow(c)),rep(b0[3],nrow(c))))
num = (c - b) * (a - b) 
den = (a - b) * (a - b)
cp = (num$V1 + num$V2 + num$V3) / (den$V1 + den$V2 + den$V3)
write(cp,"36_S15.CP-overlap.bootstraps3")
#Temporarily appended 0.6033547, 0.2502417, 0.1464035 at the end of
#  36_S15.CP-overlap.bootstraps2 to check it recovered 0.7816411 -> OK

#Rep 2 = 36_S15
c = read.table("36_S15.CP-nonoverlap.bootstraps2")
a0 = c(0.7412605, 0.1203969, 0.1383426)
b0 = c(0.1383426, 0.7412605, 0.1203969)
a = as.data.frame(cbind(rep(a0[1],nrow(c)),rep(a0[2],nrow(c)),rep(a0[3],nrow(c))))
b = as.data.frame(cbind(rep(b0[1],nrow(c)),rep(b0[2],nrow(c)),rep(b0[3],nrow(c))))
num = (c - b) * (a - b) 
den = (a - b) * (a - b)
cp = (num$V1 + num$V2 + num$V3) / (den$V1 + den$V2 + den$V3)
write(cp,"36_S15.CP-nonoverlap.bootstraps3")

#Tidy
foreach file (35_S30.CP-overlap 36_S15.CP-overlap 35_S30.CP-nonoverlap \
  36_S15.CP-nonoverlap)
  sed 's/ /\n/g' $file.bootstraps3 | grep "." > $file.bootstraps4
end
wc -l *bootstraps4
#-> 1000 1000 1000 1000



#Note, now showing median not mean, and from pnorm(-1) to pnorm(1) quantiles
#> pnorm(-1)
#[1] 0.1586553
#> pnorm(1)
#[1] 0.8413447


#See R100000 for 100000 bootstraps.

#Bootstrapping is to find the estimated distribution. The raison d'etre to do
#  it is so that we can use those distributions to put a p-value on the
#  observed.
#
#Where do the observed values for the XP:CP expression ratio for the overlap
# region fall in terms of percentile in the non-overlap distributions?
#Rep 1 observed = 0.260789
#Rep 2 observed = 0.2793595
#Rep 1 35_S30.CP-nonoverlap.bootstraps4
#Rep 2 36_S15.CP-nonoverlap.bootstraps4
cd R100000
awk '{print (1-$1)/$1}' 35_S30.CP-nonoverlap.bootstraps4 | awk '{if (0+$1>=0+0.260789) print $0}'
awk '{print (1-$1)/$1}' 36_S15.CP-nonoverlap.bootstraps4 | awk '{if (0+$1>=0+0.2793595) print $0}'
#-> Actually, the observed values are outside anything shown here. So
#  p < 1 in 100000 = 0.00001.


#Also add K-S test to show the distributions are different.
cd R100000
Rep1_olap = read.table("35_S30.CP-overlap.bootstraps4")
Rep1_nonolap = read.table("35_S30.CP-nonoverlap.bootstraps4")
Rep2_olap = read.table("36_S15.CP-overlap.bootstraps4")
Rep2_nonolap = read.table("36_S15.CP-nonoverlap.bootstraps4")
oRep1_olap     = (1 - Rep1_olap) / Rep1_olap
oRep1_nonolap  = (1 - Rep1_nonolap) / Rep1_nonolap
oRep2_olap     = (1 - Rep2_olap) / Rep2_olap
oRep2_nonolap  = (1 - Rep2_nonolap) / Rep2_nonolap
ks.test(oRep1_olap$V1,oRep1_nonolap$V1,alternative="greater")

#95% confidence intervals
#Rep1 overlap    0.12188 - 0.5060723
#Rep2 overlap    0.1144634 - 0.5888108
#Rep1 nonoverlap 0.0006841753 - 0.100235 
#Rep2 nonoverlap  0.00283727 - 0.1104908 
#> I.e. 95% confidence intervals do not overlap.

#------------------------------------------------------------------------------
