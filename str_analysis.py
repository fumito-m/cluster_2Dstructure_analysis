#!/usr/bin/env python

# str_analysis.py ver. 0.1
# Copyright (C) 2018 Fumito.M
# This script is aimed for visualization of correlation between clustering of RNA-Seq and its 2D stucture. 
# In: Paths to reference transcriptome (FASTA format), targeting clusters (FASTA format) and original cluster file (CSV) from PARalyzer.
# Out: Heatmap on clustering and its stemming probability
# Usage: Type "./str_analysis.py -h" in your terminal.

import argparse
import re
import gzip
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import collections

# ToDo: use tqdm for progress bar
# ToDo: replace print() with loggings

# Constants
ADD_N = 25
ALLOW_DIFF_N = 2
MFE_VAL = 2

def check_identity(arr1, arr2, allowance_n = 2):
	if len(arr1) != len(arr2):
		# print("Different lengths of arrays")
		return False
	else:
		count = 0
		for i in range(len(arr1)):
			if count > allowance_n:
				# print("Exceed allowing numbers of different sequence")
				return False
			if arr1[i] == arr2[i].upper() or arr1[i] == arr2[i].lower() :
				count = 0
			else:
				count += 1
		return True

def get_comp_strand(arr, reverse = True):
	out = ""
	if reverse == True:
		arr = arr[::-1]
	for i in arr:
		if i in ['a','A']:
			out += 'T'
		elif i in ['t','T']:
			out += 'A'
		elif i in ['g','G']:
			out += 'C'
		elif i in ['c','C']:
			out += 'G'
		else:
			out += 'N'
	return "".join(out)

parser = argparse.ArgumentParser(
    usage = './str_analysis.py <clust> <target> <reference>',
    description = 'This script is used to analysis correlation of clusters and RNA 2D structure. The output will be written as heatmap.png',
    epilog = 'str_analysis.py ver. 0.1')
parser.add_argument('clust', metavar = 'clust <STRING>', type = str,
                    help = 'path to the gzipped CSV file of PARalyzer cluster output')
parser.add_argument('target', metavar = 'clust <STRING>', type = str,
                    help = 'path to the target cluster FASTA file from previsou Cluster file')
parser.add_argument('ref', metavar = 'refrence <STRING>', type = str,
                    help = 'path to the gzipped refence sequence file of fasta file obtained from a database')
parser.add_argument('-n', '--num', metavar = '<INT>', type = int, help = 
	'# of times repeatedly make shuffled alignments (default <INT>) = 500', default = 500)
# Not implemented
parser.add_argument('-b', '--bin', metavar = '<INT>', type = int, help = 
	'# bins of stemming probability on heatmap (default <INT>) = 3', default = 3)
args = parser.parse_args()


# Fetch line numbers of reference chromosomes
ref_array = {}
r1 = re.compile(r">.+$")
with open(args.ref, "rt") as f:
	for n, line in enumerate(f):
		if r1.match(line):
			key = line.strip().replace('>','')
			ref_array.update({key: n})

# Get targeting arrays
target = {}
break_count = 0
r2 = re.compile(r">.+$")
with open(args.target, "rt") as f:
	for line in f:
		# For shortcut
		# if break_count > 10:
		# 	break
		if r2.match(line):
			key = line.strip().replace('>','')
			break_count += 1
		else:
			target.update({key: {"Array":line.strip()}})
	print('num of target clusters:', len(target))
	# print('Target: ',target.keys())

# Extend arrays based on reference transcriptome
# ToDo: Replace method with belows
# import linecache
# linecache.getline(args.clust, int(num))
r3 = re.compile(r"(chr\w+),([\+\-]),(\d+),(\d+),(\w+\.\d+),(\w+),(.+)$")
with gzip.open(args.clust, "rt") as f:
	for n, line in enumerate(f):
		m1 = r3.match(line)

		if n == 0 and not m1:
			header = line.strip()
		elif m1 and m1.group(5) in target.keys():
			clst = target[m1.group(5)]
			q_s, mod_s = divmod(int(m1.group(3))-ADD_N, 50)
			q_e, mod_e = divmod(int(m1.group(4))+ADD_N, 50)

			if m1.group(6) != clst['Array']:
				print("Warn: different arrays from original cluster file")

			target_line = ""
			com = 'sed -n '+str(ref_array[m1.group(1)]+q_s+2)+','+str(ref_array[m1.group(1)]+q_e+2)+'p '+args.ref
			# Only for linux
			ref = subprocess.check_output(com.split(" ")).decode('utf-8').strip().replace('\n','')
			ref = ref[mod_s-1:-50+mod_e]

			# ToDo: Check if getting from identical chromosome arrays.
			# ToDo: Set searching conditions
			if not check_identity(clst['Array'],ref[25:-25], ALLOW_DIFF_N):
				if check_identity(clst['Array'], get_comp_strand(ref[25:-25]), ALLOW_DIFF_N):
					clst['Array'] = get_comp_strand(ref[-25:])+clst['Array']+get_comp_strand(ref[:25])
				else:
					print("Warn: ", m1.group(5), "does not allign with reference array. Skipping extension")
					print('Target:', clst['Array'])
					print('Ref:', ref[25:-25])
			else:
				clst['Array'] = ref[:25]+clst['Array']+ref[-25:]
			# print('Head:',ref[0:25],'Ori:',ref[25:-25], 'Tail:',ref[-25:])

# Calculate 2D structures including shuffled arrays
for key in target:
	arr = target[key]['Array']
	com = 'ushuffle -s '+ arr +' -n ' + str(args.num)
	shuffled = subprocess.check_output(com.split(" ")).decode('utf-8').strip()
	with open("tmp.file", "w") as f:
		f.write(arr+'\n')
		f.write(shuffled)

	com = 'RNAfold --MEA tmp.file'
	folded = subprocess.check_output(com.split(" ")).decode('utf-8').strip().split('\n')

	# ToDo: replace with probability numbers and divide into bins
	t_counts = pd.DataFrame(columns=['b1','b2','b3'])
	for l in range(0, args.num+1):
		counts = collections.Counter(folded[l*6 + MFE_VAL])
		n1 = counts['('] + counts[')'] + counts['|']
		n2 = counts['{'] + counts['}'] + counts[',']
		n3 = counts['.']
		nn = pd.Series([n1,n2,n3], index=['b1', 'b2', 'b3'], name=l)
		t_counts = t_counts.append(nn/float(len(arr)))

	z_score = (t_counts[0:1] - t_counts[1:].mean())/t_counts[1:].std()
	target[key].update({'b1': z_score.loc[0,'b1'], 'b2': z_score.loc[0,'b2'], 'b3': z_score.loc[0,'b3']})


data = pd.DataFrame.from_dict(target).T.drop("Array", axis=1)
data2 = data[data.columns].astype(float)
plt.figure()
sns.heatmap(data2, cmap='Blues')
plt.savefig('./str_heatmap.png')
plt.close('all')