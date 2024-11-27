#!/usr/bin/env python
#import mudules
import os
import subprocess

#set current directory
current_dir = "/home/s2647596/python_ICA"
os.chdir(current_dir)


#protein family
protein_name = "glucose-6-phosphatase"
#species
species_name = "Aves"

#use esearch and efetch in linux
def fasta_result(protein_name, species_name):
    esearch_query = f"esearch -db protein -query '{species_name}[organism] AND {protein_name}[Title]'"
    efetch_query = "efetch -format fasta"
    try:
        fa_result = subprocess.check_output(f"{esearch_query} | {efetch_query}", shell=True).decode("utf-8")
        #save the fasta result into a file
        with open(f"{species_name}.fasta", mode="w") as f:
            f.write(fa_result)
        return fa_result
    except Exception as e:
        print(e)
        return None

#fa_result = fasta_result(protein_name, species_name)
   

#align the protein sequences and plot the conservation
def plot_conserve(species_name):
    aln_query = f"clustalo -i {species_name}.fasta -o {species_name}.aln"
    try:
        subprocess.call(aln_query, shell=True)
        with open(f"{species_name}.aln", mode="r") as f:
            print("success alignment")
        #plot the conservation
        plot_query = f"plotcon -sequence {species_name}.aln -graph png"
        result = subprocess.check_output(plot_query, shell=True).decode("utf-8")
        print(result)
    except Exception as e:
        print(e)

#plot_conserve(species_name)

#find motif from PROSTIE database that associated with the protein sequences
#split and save individual sequence to independent FASTA files
#create a directory for fasta files
if not os.path.exists(f"{current_dir}/fasta_files"):
    os.mkdir("./fasta_files")

#split sequence
with open(f"{species_name}.fasta", mode="r") as f:
    fa_result = f.readlines()
    for i in range(len(fa_result)):
        line = fa_result[i].strip()
        if line.startswith(">"):
            #use accession number as file name
            acc_id = line.split(" ")[0]
            acc_id = acc_id[1:]
            #header line
            with open(f"./fasta_files/{acc_id}.fasta", mode="w") as inner_f:
                inner_f.write(line)
        elif fa_result[i-1].startswith(">"):
            #first sequence line
            with open(f"./fasta_files/{acc_id}.fasta", mode="a") as inner_f:
                inner_f.write(f"\n{line}")
        else:
            #other sequnece line
            with open(f"./fasta_files/{acc_id}.fasta", mode="a") as inner_f:
                inner_f.write(line)

