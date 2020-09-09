NAME=$1

hadd -k temp1.root ${NAME}{0..9}.root

for (( i = 0; i < 10; i++ )); do
	hadd -k temp2${i}.root ${NAME}${i}{0..9}.root
	
	for (( j = 0; j < 10; j++ )); do	
	    hadd -k temp2${i}${j}.root ${NAME}${i}${j}{0..9}.root		
	done
done

for (( i = 0; i < 10; i++ )); do
	hadd -k temp3${i}.root temp2${i}{0..9}.root		
done

hadd -k temp3.root temp3{0..9}.root
hadd -k temp2.root temp2{0..9}.root

hadd -k total.root temp{1..3}.root

rm -f temp*.root
