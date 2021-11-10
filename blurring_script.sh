# awk -v v1=$var1 -v v2=$var2 'BEGIN {print 141*(v1*0.0113/0.9)^-1.209*0.993^v2}' cre >> ggh

x=$(awk -v v1=$1 'BEGIN {print 1-2*v1}')

y=$(awk -v v1=$x 'BEGIN {print log(v1)/log(10)}')
loge=$(awk -v v1=$1 'BEGIN {print log(v1)/log(10)}')

#y=$(echo 'l($x)/l(10)' | bc -l)
#loge=$(echo 'l($1)/l(10)' | bc -l)

#echo $y
#echo $loge

sed -e 's/\tGT/\tGL/g' -e 's/0|0/'$y','$loge','$loge'/g' -e 's/0|1/'$loge','$y','$loge'/g' -e 's/1|0/'$loge','$y','$loge'/g' \
-e 's/1|1/'$loge','$loge','$y'/g' $2
