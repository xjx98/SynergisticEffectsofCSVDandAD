indx=""
for ((i=1; i<=$1; i+=1)); do indx="$indx 1"; done
echo $indx > $2/index.txt
