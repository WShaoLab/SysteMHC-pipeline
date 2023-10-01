
while read -r line
do
   grep $line $2 > out.txt

   if [ ! -s out.txt ]; then
       echo $line >> variant.txt 
   fi
done < $1
