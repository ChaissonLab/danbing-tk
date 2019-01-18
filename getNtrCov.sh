# $1: locus number

NR=-1
locus=""
while read -r line; do
    samtools view -b -r ERR899717 HG00514.IL.bam $line | samtools depth /dev/stdin | head -n -200 | tail -n+200 | awk '{c += $3} END {print c/NR}'
    #NR=$((NR+1))
    #if [ $NR -eq $1 ]; then
    #    locus="$line"
    #    break
    #fi
done < regions.ntr2000.bed.txt

#echo "locus:" $NR $locus
