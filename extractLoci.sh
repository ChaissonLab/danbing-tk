#flist=("ERR899717_1.fastq" "ERR899717_2.fastq" "ERR903031_1.fastq" "ERR903031_2.fastq" "GM19240_NoIndex_L008_R1_001.fastq" "GM19240_NoIndex_L008_R2_001.fastq")
flist=("HG00514" "HG00733" "NA19240")
fsize=${#flist[@]}
loci=-1
# max = 5349
max=9
for ((ind=0; ind<$fsize; ind++)); do
	file=${flist[$ind]}
	if [ ! -f "$file"".L0.kmercount" ]; then
		awk -v LL="$loci" -v ff="$file" -v max="$max" '{ if (substr($1, 0, 1) == ">") { LL += 1; start = 1; if (LL == max+1) {exit} } else { if (start == 1) {print >>  ff".L"LL".kmercount"} } }' "$file"".21.kmers"
		for i in `seq 0 $max`; do
			sort "$file"".L""$i"".kmercount" | awk '{print $2}' >> "$file"".L""$i"".sort.count"
			# sort "$file"".L""$i"".kmercount" | awk '{print $1}' >> "$file"".L""$i"".sort.kmerstring"
			echo "$i done"
		done
	else
		echo "file already exists!"
	fi
done
