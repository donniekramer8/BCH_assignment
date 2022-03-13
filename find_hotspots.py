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

chrom_number = [i.split('\t', 1)[0] for i in lines_chrom]
chrom_length = [i.split('\t', 1)[1] for i in lines_chrom]
chrom_chromnum = [i.split('chr', 1)[1] for i in chrom_number]

z = 0
int_chrom_num = []
while z < len(chrom_chromnum):
    if z < 22:
        int_chrom_num.append(int(chrom_chromnum[z]))
        z += 1
    else:
        int_chrom_num.append(chrom_chromnum[z])
        z += 1

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

reflist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
           'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX',
           'chrY']

list_chr_bins = [[] for _ in range(len(reflist))]
temp_list = []
b = 0
while b < len(reflist):
        a = 1
        k = 0
        i = 2000000
        for element in range(len(lines_RF204)-1):
                if (list(lines_RF204[a].split())[1]) == reflist[b]:
                        while k <= int(list(lines_RF204[a].split())[2]) < i:
                                temp_list.append(int(list(lines_RF204[a].split())[2]))
                                a += 1
                        if (k <= int(list(lines_RF204[a].split())[2]) < i) == False:
                                list_chr_bins[b].append(len(temp_list))
                                temp_list[:] = []
                                k += 2000000
                                i += 2000000
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
for x in list_chr_bins:
    merged_list = merged_list + list_chr_bins[y]
    y += 1

q75, q25 = np.percentile(merged_list, [75, 25])
iqr = q75 - q25

upper_fence = q75+(1.5*iqr)

outliers = []
df0 = []
df1 = []
df2 = []
df3 = []
z = 0
t = 0
while t < len(list_chr_bins):
    for pos in range(len(list_chr_bins[t])):
        if int(list_chr_bins[t][z]) > 31.5:
            if int(pos*2000000 + 2000000) > LengthDict[t+1]:
                df0.append(reflist[t])
                df1.append(pos * 2000000)
                df2.append(LengthDict[t+1])
                df3.append(list_chr_bins[t][z])
                z += 1
            else:
                df0.append(reflist[t])
                df1.append(pos * 2000000)
                df2.append(pos * 2000000 + 2000000)
                df3.append(list_chr_bins[t][z])
                z += 1
        else:
            z += 1
        if pos == int(len(list_chr_bins[t]))-1:
            break
    t += 1
    z = 0

data = {'start pos':df1, 'end pos':df2, 'translocations':df3}
df = pd.DataFrame(data, index=df0)
pd.set_option('display.max_rows', None)
df.sort_values(by=['translocations'], inplace=True, ascending=False)
print(df)

unique_values = []
for x in sorted(merged_list):
    if x in unique_values:
        i += 1
    if x not in unique_values:
        unique_values.append(x)
        i += 1

i = 0
number_occurences = []
for x in merged_list:
    k = unique_values[i]
    number_occurences.append(merged_list.count(k))
    i += 1
    if i == len(unique_values):
        break

fig, ax = plt.subplots()
plt.plot(unique_values, number_occurences, marker='o', markersize='3')
plt.xlim(0,6000)
plt.ylim(0, 135)
plt.xscale('symlog')
plt.axvline(31.5, color='red', linestyle='dashed', label='hotspot cutoff')
plt.xlabel('Number of Translocations')
plt.ylabel('Number of 2 Mb bins')
ax.legend()
plt.title('Translocations per 2Mb bin')
plt.show()
