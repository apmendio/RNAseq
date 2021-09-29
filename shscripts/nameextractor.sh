#Notes - File organizer and sample task builder for sequencing data

#list command to organize files by name and numerics

ls -1v "file.name.here" > "file.output.txt"

###if there are duplicate names because of read pairs follow the next steps below

#Acquire unique name identifiers - usually the first 3 - 5 letters in sample name

cat "file.output.txt" | cut -c1-c"3"/"5" > "file.output2.txt"

#Remove duplicates from read pair names

cat "file.output2.txt" | uniq -d > "file.output2.txt"