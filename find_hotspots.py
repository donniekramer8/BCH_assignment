#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def readFile(filename):
        with open(filename, 'r') as f:
                return [l.strip() for l in f.readlines()]

RF204_Alt055 = readFile('RF204_Alt055.tlx')
lines_RF204 = []
for line in RF204_Alt055:
        lines_RF204.append(line)

ChromInfo = readFile('ChromInfo.txt')
lines_chrom = []
for line in ChromInfo:
    lines_chrom.append(line)

reflist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
           'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX',
           'chrY']

# create list of lists containing number of translocations in 2 Mb regions of each chromosome
# this will create a list for each chromosome which shows how many translocations there are in regions divided into 2 Mb bins
list_chr_bins = [[] for _ in range(len(reflist))]
temp_list = []
b = 0
while b < len(reflist):
        a = 1
        k = 0
        i = 2000000
        for element in range(len(lines_RF204)-1):
                # list(lines_RF204[a].split())[1] is the data in the 'Rname' column of RF204_Alt055.tlx
                # 'a' increments to read through the lines of RF204_Alt055.tlx, 'a' is initialized to 1 to ignore header line
                if (list(lines_RF204[a].split())[1]) == reflist[b]:
                        # int(list(lines_RF204[a].split())[2]) is the data in the 'Junction' column of RF204_Alt055.tlx
                        while k <= int(list(lines_RF204[a].split())[2]) < i:                  
                                # append 'Junction' data to a temporary list
                                temp_list.append(int(list(lines_RF204[a].split())[2]))
                                a += 1
                        # when 'Junction' value is outside current bin width, append number of 'Junction' values to list_chr_bins[b] and clear temp_list
                        # increment bin width by 2000000
                        if (k <= int(list(lines_RF204[a].split())[2]) < i) == False:
                                list_chr_bins[b].append(len(temp_list))
                                temp_list[:] = []
                                k += 2000000
                                i += 2000000
                        # if at end of chromosome, reset bin width
                        elif int(list(lines_RF204[a+1].split())[2]) < int(list(lines_RF204[a].split())[2]):
                                k = 0
                                i = 2000000
                else:
                        a += 1
                if a == len(lines_RF204):
                        break
        b += 1

merged_list = []
y = 0
for list_ in list_chr_bins:
    merged_list = merged_list + list_chr_bins[y]
    y += 1

q75, q25 = np.percentile(merged_list, [75, 25])
iqr = q75 - q25

# values higher than upper_fence represent outliers within the data, which is my definition of a 'hotspot' here
upper_fence = q75+(1.5*iqr)

# create dictionary for data in ChromInfo.txt
chrom_number = [i.split('\t', 1)[0] for i in lines_chrom]
chrom_length = [i.split('\t', 1)[1] for i in lines_chrom]
chrom_num = [i.split('chr', 1)[1] for i in chrom_number]

z = 0
int_chrom_num = []
while z < len(chrom_num):
    if z < 22:
        int_chrom_num.append(int(chrom_num[z]))
        z += 1
    # else statement is for 'X', 'M', and 'Y' names
    else:
        int_chrom_num.append(chrom_num[z])
        z += 1

# replace M,X,Y with 23,24,25 in list/dictionary to simplify for later use
int_chrom_num.remove('M')
int_chrom_num.insert(22, 23)
int_chrom_num.remove('X')
int_chrom_num.insert(23, 24)
int_chrom_num.remove('Y')
int_chrom_num.insert(24, 25)

z = 0
int_chrom_length = []
while z < len(chrom_length):
    int_chrom_length.append(int(chrom_length[z]))
    z += 1

tuplelist = list(zip(int_chrom_num, int_chrom_length))
LengthDict = dict(tuplelist)

# function below creates a dataframe containing chr#, start_pos, end_pos, and # transpositions for each hotspot
chr = []
start_pos = []
end_pos = []
translocations = []
z = 0
t = 0
while t < len(list_chr_bins):
    for pos in range(len(list_chr_bins[t])):
        if int(list_chr_bins[t][z]) > upper_fence:
            # if the bin is the last bin of the chromosome, return length of chromome as 'end_pos'
            # list_chr_bins[t][z] is the number of transpositions in a bin, aka: number of 'Junction' values within a bin's boundaries
            if int(pos*2000000 + 2000000) > LengthDict[t+1]:
                chr.append(reflist[t])
                start_pos.append(pos * 2000000)
                end_pos.append(LengthDict[t+1])
                translocations.append(list_chr_bins[t][z])
                z += 1
            # if bin is not last bin of chromosome, simply add bin width to 'start_pos' to get 'end_pos'
            else:
                chr.append(reflist[t])
                start_pos.append(pos * 2000000)
                end_pos.append(pos * 2000000 + 2000000)
                translocations.append(list_chr_bins[t][z])
                z += 1
        else:
            z += 1
        if pos == int(len(list_chr_bins[t]))-1:
            break
    t += 1
    z = 0

# display data frame
data = {'Start Position':start_pos, 'End Position':end_pos, 'Translocation Count':translocations}
df = pd.DataFrame(data, index=chr)
pd.set_option('display.max_rows', None)
df.sort_values(by=['Translocation Count'], inplace=True, ascending=False)
print(df)

# create lists for data visualization
i = 0
unique_values = []
for x in sorted(merged_list):
    if value in unique_values:
        i += 1
    else:
        unique_values.append(x)
        i += 1

i = 0
number_occurences = []
for value in range(len(unique_values)):
    k = unique_values[i]
    number_occurences.append(merged_list.count(k))
    i += 1

# create plot with matplotlib
fig, ax = plt.subplots()
plt.plot(unique_values, number_occurences, marker='o', markersize='3')
plt.xlim(0,6000)
plt.ylim(0, 135)
# x axis is in logarithmic scale due to most values being < 100, though max is ~5500
plt.xscale('symlog')
plt.axvline(31.5, color='red', linestyle='dashed', label='hotspot cutoff')
plt.xlabel('Number of Translocations')
plt.ylabel('Number of 2 Mb bins')
ax.legend()
plt.title('Translocations per 2Mb bin')
plt.show()
