max=5349
cd ./image
for i in `seq 0 4` 
do
	rm ../dot/loci."$i".sort.dot
	rm loci."$i".sort.svg
	head -n 1 ../dot/loci."$i".dot >> ../dot/loci."$i".sort.dot
	head --lines=-1 ../dot/loci."$i".dot| tail -n +2 | sort >> ../dot/loci."$i".sort.dot
	tail -n 1 ../dot/loci."$i".dot >> ../dot/loci."$i".sort.dot
	cat ../dot/loci."$i".sort.dot | dot -Tsvg -oloci."$i".sort.svg
    
    echo "$i done"
done
