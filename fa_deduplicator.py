import argparse

### Read FASTA file into a dictionary
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

### Remove duplicate values from the dictionary
def deduplicate(seqs):
    vals = set()
    useq = dict()
    
    for key, value in seqs.items():
        if value not in vals:
            vals.add(value)
            useq[key] = value
    
    return useq

### Write the dictionary to a new FASTA File
def writefaFile(useqs, path):
    with open(path, "w") as file: 
        for key, value in useqs.items():
            file.write(f">{key}\n")
            file.write(f"{value}\n")
            
### Main function to handle workflow execution
def main():
    parser = argparse.ArgumentParser(description="Deduplicate FASTA Sequences")
    parser.add_argument("input_fasta", help="Input FASTA file")  
    parser.add_argument("output_fasta", help="Output FASTA file")   
    args = parser.parse_args()
    
    initial = readfaFile(args.input_fasta)
    processed = deduplicate(initial)
    writefaFile(processed, args.output_fasta)
    print(f"Deduplicated sequences written to {args.output_fasta}")
    
if __name__ == "__main__":
    main()