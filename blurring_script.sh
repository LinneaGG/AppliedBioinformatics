e=$(bc -l <<< "e($1*l(10))")

x=$(awk -v v1=$e 'BEGIN {print 1-2*v1}')

y=$(awk -v v1=$x 'BEGIN {print log(v1)/log(10)}')
#loge=$(awk -v v1=$1 'BEGIN {print log(v1)/log(10)}')


sed -e 's/\tGT/\tGL/g' -e 's/0|0/'$y','$1','$1'/g' -e 's/0|1/'$1','$y','$1'/g' -e 's/1|0$
-e 's/1|1/'$1','$1','$y'/g' $2
