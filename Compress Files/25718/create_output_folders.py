import os

# this script should be run AFTER creating a samples.csv file
# and BEFORE submitting a CHTC job
# it creates the output folders that CHTC expects, with subfolders for all files from a single animal
# the folder name is read from the first item in the samples.csv file

# create dictionary to hold sample names
sn = []

# read file
# extract sample name as text in first field
with open('samples.csv') as f:
    for line in f:
        sn.append(line.strip().split(',')[0])

# create output folders
for i in sn:
    os.makedirs('out/' + str(i) + '/logs', exist_ok=True)
        