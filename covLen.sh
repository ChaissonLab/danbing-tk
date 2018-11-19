# arg1 = <.bam file>
# arg2 = <coverage>
# arg3 = <output file>

for r in `cat regions.orig.txt `; do 
    samtools depth $1 -r  $r |  
    awk -v cov=$2 '{ s += $3} END {print s/cov;}'; 
done | tee $3
