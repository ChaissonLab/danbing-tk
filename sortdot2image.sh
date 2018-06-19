max=5349
cd ./image
hap="ERR899717_1."
for i in `seq 0 4` 
do
	#rm ../dot/loci."$i".sort.dot
	#rm loci."$i".sort.svg
	head -n 1 ../dot/"$hap"loci."$i".dot >> ../dot/"$hap"loci."$i".sort.dot
	head --lines=-1 ../dot/"$hap"loci."$i".dot| tail -n +2 | sort >> ../dot/"$hap"loci."$i".sort.dot
	tail -n 1 ../dot/"$hap"loci."$i".dot >> ../dot/"$hap"loci."$i".sort.dot
	cat ../dot/"$hap"loci."$i".sort.dot | dot -Tsvg -o"$hap"loci."$i".sort.svg
    
    echo "$i done"
done
