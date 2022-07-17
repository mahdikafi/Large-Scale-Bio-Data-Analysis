#!/bin/bash

echo "Enter the file name: "
read file_name
cat /home/mahdi/Work/code/MBDA/HW0/Q2/hw0_2/gene*/gene*.txt > $file_name
sort -Vk 3 $file_name > tmp_file
cat tmp_file > $file_name
rm -f tmp_file
while IFS=" " read -r name pos len
do
   rand_str=$(tr -dc "ACGT" </dev/urandom | head -c $len)
   printf '%s %s %s\n' "$name" "$pos" "$rand_str" >> tmp_file
done < "$file_name"
cat tmp_file > $file_name
rm -f tmp_file
