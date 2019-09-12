#!/bin/bash


#Create empty folders for jobs
for i in {1..105}; do 
if [ -d "$i" ]; then
rm -rf $i
fi
mkdir $i

#Copy Rcode and PBS files to folders
cp iid_mcpt.R $i; 
cp run.pbs $i; 
cd $i; 

#Run simulation under each folder
qsub run.pbs;
cd ..
done;



