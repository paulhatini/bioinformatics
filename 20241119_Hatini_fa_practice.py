import numpy as np
import itertools
from Bio import SeqIO

## reading .fasta files

### native python
def readfaFile(path):
    
    seqs = dict()
    key = ""
    value = ""
    
    with open(path) as file:
        for line in file:
            if line[0] == ">":
                if value != "": seqs[key] = value
                key = line[1:].strip()
                value = ""
            else:
                value += line.strip()
                
    seqs[key] = value ## this adds the final line to the seq dictionary
    
    return seqs

amplicons_preprocess = readfaFile("/Users/paulhatini/Downloads/amplicons.undup.fa")

### using itertools.groupby 
def iskey(line):
    return line[0] == ">"

def readfaFile2(path):
    
    seqs = dict()
    
    with open(path) as file:
        for label,group in itertools.groupby(file, iskey):
            if label:
                line = group.__next__()
                key = line[1:].strip()
            else:
                value = "".join(line.strip() for line in group)
                seqs[key] = value
                
    return seqs
                
amplicons_preprocess2 = readfaFile2("/Users/paulhatini/Downloads/amplicons.undup.fa")

### using biopython SeqIO
def readfaFile3(path):
    
    seqs = dict()
    
    for entry in SeqIO.parse(path, "fasta"):
        seqs[entry.id] = str(entry.seq)
        
    return seqs

amplicons_preprocess3 = readfaFile3("/Users/paulhatini/Downloads/amplicons.undup.fa")

print(amplicons_preprocess == amplicons_preprocess2 == amplicons_preprocess3)

## remove duplicate sequences (values) from the dictionary

### leverage property of unique keys in a dictionary (resulting dictionary contains the key of the last kvp where the value is seen)
def removeDuplicates(seqs):
    
    temp = dict()
    
    for key in seqs:
        temp[seqs[key]] = key
        
    useq = dict()
    
    for key in temp:
        useq[temp[key]] = key
        
    return useq

amplicons_postprocess = removeDuplicates(amplicons_preprocess)

### using sets (resulting dictionary contains the key of the first kvp where the value is seen)
def removeDuplicates2(seqs):
    vals = set()
    useq = dict()
    
    for key, value in seqs.items():
        if value not in vals:
            vals.add(value)
            useq[key] = value
    
    return useq

amplicons_postprocess2 = removeDuplicates2(amplicons_preprocess)


## report counts

### native python
counts = dict()
vals = set()

for key, value in amplicons_preprocess.items():
      if value not in vals:
          vals.add(value)
          counts[value] = 1
      else:
          counts[value] = counts[value] + 1
    
### np.unique
processedValues, processedCounts = np.unique(list(amplicons_preprocess.values()), return_counts=True)

counts2 = dict()
for i, item in enumerate(processedValues):
    counts2[processedValues[i]] = processedCounts[i]
    

## write to file

### native python
def writeFile(useqs, path):
    with open(path, "w") as file: 
        for key, value in useqs.items():
            file.write(f">{key}\n")
            file.write(f"{value}\n")

writeFile(amplicons_postprocess2, "/Users/paulhatini/Downloads/amplicons.fa")