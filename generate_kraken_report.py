import pandas as pd
from glob import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import numpy as np
from collections import Counter, defaultdict
from subprocess import call
from math import log
from tqdm import tqdm
from copy import deepcopy
import gzip
import argparse

parser = argparse.ArgumentParser(description='Create Kraken2 style report from tsv')
parser.add_argument('-i', '--input', type=str, help='input tsv')
parser.add_argument('-o', '--output', type=str, help='output file')
parser.add_argument('-taxid_col', '--taxid_col', type=str, default='taxid', help='Column name with taxid')
parser.add_argument('-count_col', '--count_col', type=str, default='scaffold_count', help='Column name with counts')
args = parser.parse_args()

PATH_TO_TAXDB = '...'


####################################################################################################
#Import functions
print('Importing functions, nodes.dmp, names.dmp')

taxs={}
with open(PATH_TO_TAXDB + 'nodes.dmp') as f:
    file=f.read().strip().split('\n')
for j,i in enumerate(file):
    daughter=i.split('\t')[0]
    parent=i.split('\t')[2]
    taxs[daughter]=parent

taxid2rank={}
for i in file:
    taxid = i.split('\t')[0]
    rank = i.split('\t')[4]
    taxid2rank[taxid] = rank
taxid2rank['1'] = 'root'
    
names={}
with open(PATH_TO_TAXDB + 'names.dmp') as f:
    file=f.read().strip().split('\n')
for j,i in enumerate(file):
    daughter=i.split('\t')[0]
    parent=i.split('\t')[2]
    names[daughter]=parent

good_names={}
for j,i in enumerate(file):
    daughter=i.split('\t')[0]
    if daughter in good_names:
        good_names[daughter].append(i)
    else:
        good_names[daughter]=[i]

very_good_names={}
for i in good_names:
    nn=names[i]
    flag=False
    for j in good_names[i]:
        if j.split('\t')[-2]=='scientific name':
            very_good_names[i]=j.split('\t')[2]
            flag=True
    if not flag:
        very_good_names[i]=nn

ranks = {
    'root': 'R',
    'superkingdom': 'D',
    'kingdom': 'K',
    'phylum': 'P', 
    'class': 'C',
    'order': 'O',
    'family': 'F',
    'genus': 'G',
    'species': 'S'
}

def taxid2krakenrank(x):
    if (x=='UNIDENTIFIED') or (x=='0'):
        return 'U'
    if (x=='Host bowtie2') or (x=='Host'):
        return 'S'
    if taxid2rank[x] in ranks:
        return ranks[taxid2rank[x]]
        
    xz=0
    x_path_to_root=[x]
    xtemp=x
    while xz!='1':
        if xtemp in taxs:
            xz=taxs[xtemp]
        else: return 'Impossible to determine. Missing taxes'
        x_path_to_root.append(xz)
        xtemp=xz
    _c=0
    for i in x_path_to_root:
        if taxid2rank[i] not in ranks:
            _c+=1
            continue
        else:
            return ranks[taxid2rank[i]]+str(_c)
    raise Exception

def path_to_root(x):
    xz=0
    x_path_to_root=[x]
    xtemp=x
    while xz!='1':
        if xtemp in taxs:
            xz=taxs[xtemp]
        else: return 'Impossible to determine. Missing taxes'
        x_path_to_root.append(xz)
        xtemp=xz
    return x_path_to_root

def spaces_amount(node):
    return len(set(path_to_root(node))-set(['1']))


def percent(x, total):
    newstr = str(round(float(x)/total*100.0, 2))
    if newstr.split('.')[1]=='0':
        newstr = newstr.split('.')[0]+'.00'
    lenn = len(newstr)
    return ' '*(6-lenn)+newstr

def DFS(visited, graph, node):
    stack = []
    stack.append(node)
    while len(stack) != 0:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            global kraken_from_blast
            kraken_from_blast.append([
                str(taxid2abundance[current][0]), str(taxid2abundance[current][1]),
                taxid2krakenrank(current), current,' '*2*spaces_amount(current)+very_good_names[current]
            ])
            for nodes in graph[current]:
                if nodes not in visited:
                    stack.append(nodes)
                    
####################################################################################################
#Reading file

print('Reading file')
temp = pd.read_csv(args.input, usecols = [args.taxid_col, args.count_col], sep = '\t').set_index(args.taxid_col).astype(int)
assert len(temp.index) == len(set(temp.index)) # All uniques
new_res = {}
for i in temp.index:
    new_res[str(i)] = temp.loc[i, args.count_col]

# DOWN TO TOP
print('Processing...')
down2top = {}
for i in new_res:
    pathtoroot = path_to_root(i)
    assert pathtoroot[-1] == '1'
    for j in range(len(pathtoroot)-1):
        down2top[pathtoroot[j]] = pathtoroot[j+1]
down2top['1']= '1'

taxid2abundance = {}
for i in down2top:
    if i not in taxid2abundance:
        if i in new_res:
            taxid2abundance[i] = [new_res[i], new_res[i]]
        else:
            taxid2abundance[i] = [0, 0]
leaves = set(down2top.keys()) - set(down2top.values())

visited = set()
for leaf in leaves:
    pathtoroot = path_to_root(leaf)
    previous = leaf
    accum = 0
    for tax in range(len(pathtoroot)-1):
        if pathtoroot[tax] not in visited:
            accum += taxid2abundance[pathtoroot[tax]][1]
        taxid2abundance[pathtoroot[tax+1]][0] += accum
        visited.add(pathtoroot[tax])


# TOP TO DOWN
top2down = defaultdict(set)
for i in new_res:
    pathtoroot = path_to_root(i)
    for j in range(len(pathtoroot)-1):
        top2down[pathtoroot[j+1]].add(pathtoroot[j])
if '1' in top2down:
    top2down['1'] = top2down['1']-set(['1'])

kraken_from_blast = []
visited = set()
DFS(visited, top2down , '1')
kraken_from_blast = pd.DataFrame(kraken_from_blast)

total = sum(kraken_from_blast[1].astype(int))
kraken_from_blast[5] = kraken_from_blast[0].apply(lambda x: percent(x, total))
kraken_from_blast = kraken_from_blast[[5,0,1,2,3,4]]
kraken_from_blast.to_csv(args.output, index=None, header=None, sep='\t')

