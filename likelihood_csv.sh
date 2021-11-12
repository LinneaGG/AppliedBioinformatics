# use:
# bash csv_script.sh first_likelihood_value interval last_likelihood_value > file.csv

seq $1 $2 $3 | tr '\n' ',' | sed 's/.$//'
