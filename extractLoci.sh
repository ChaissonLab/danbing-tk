# arg1: haplotype file name e.g. PanGenomeGenotyping.21.1950.kmers or ERR899717_1.fastq.21.kmers
# arg2: loci e.g."all" or "big or ""beg4870" or "itv100,500" or "end4870" or "4870"


nloci=( $(( $(awk '{ if (substr($0,0,1) == ">") {c++} } END {print c}' "$1") - 1)) )
lociList=()
if [ "$2" == "all" ]; then
    lociList=($(seq 0 1 $nloci))
    echo "extracting all loci"
elif [[ "$2" == "beg"* ]]; then
    beg=$2
    beg=${beg:3}
    lociList=($(seq $beg 1 $nloci))
    echo "extracting loci starting from: ""${lociList[0]}"
elif [[ "$2" == "end"* ]]; then
    end=$2
    end=${end:3}
    lociList=($(seq 0 1 $end))
    echo "extracting loci till: ""${lociList[0]}"
elif [[ "$2" == "itv"* ]]; then
    itv=$2
    itv=${itv:3}
    itv=(${itv//,/ })
    lociList=($(seq ${itv[0]} 1 ${itv[1]}))
    echo "extracting loci from: ""${lociList[0]}"" to ""${lociList[-1]}"
elif [[ "$2" == "big" ]]; then
    awk 'BEGIN {loci = -1} 
        {
            if (substr($0,0,1) == ">") {
                if (loci != -1) {print loci, c}
                c = 0;
                loci += 1
            } else {c += $2} 
        } 
        END {print loci, c}' "$1" > "$1"".loci.count"
    lociList=($(sort -nk2,2 "$1"".loci.count" | tail | awk '{print $1}'))
    echo "extracting loci: ""${lociList[@]}"
else
    num=$(($#-1))
    lociList=( ${@:2:$num} )
    echo "extracting loci: ""${lociList[@]}"
fi

lsize=${#lociList[@]}
out=$(basename -- "$1")
ii=-1
for ((ind=0; ind<$lsize; ind++)); do
    loci=${lociList[$ind]}
    if [ ! -f "$out"".L"$loci".sort.count" ]; then
        awk -v II="$ii" -v LL="$loci" -v oo="$out" '{ if (substr($1, 0, 1) == ">") { II += 1; if (II == LL) {start = 1}; if (II == LL+1) {exit} } else { if (start == 1) {print >> oo".L"LL".kmercount"} } }' "$1"
        if [ -f "$out"".L"$loci".kmercount" ]; then
            sort "$out"".L""$loci"".kmercount" | awk '{print $2}' >> "$out"".L""$loci"".sort.count"
            echo "$loci done"
        else
            echo "$loci is empty"
        fi
    else
        echo "file already exists!"
    fi
done

