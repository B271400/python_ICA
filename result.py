#!/usr/bin/env python
#import mudules
import os
import subprocess
import re
import sys
import shutil

#set current directory
current_dir = "/home/s2647596/python_ICA"
os.chdir(current_dir)

#remove existing files
if os.path.exists("./fasta_files"):
    shutil.rmtree("./fasta_files")
if os.path.exists("./motif_files"):
    shutil.rmtree("./motif_files")
if os.path.exists("./conservation_files"):
    shutil.rmtree("./conservation_files")

'''#protein family
protein_name = "glucose-6-phosphatase"
#species
species_name = "Aves"'''

#create a directory for conservation files
if not os.path.exists(f"{current_dir}/conservation_files"):
    os.mkdir("./conservation_files")
os.chdir(f"{current_dir}/conservation_files")
#use esearch and efetch in linux
def fasta_result():
    #three attemps for user to enter the correct species name or protein name
    attempts=0
    while attempts < 3:
        #protein family
        protein_name = input("enter the protein family name correctly")
        #species
        species_name = input("enter the taxonomy correctly")
        
        esearch_query = f"esearch -db protein -query '{species_name}[organism] AND {protein_name}[Title]'"
        efetch_query = "efetch -format fasta"
        try:
            fa_result = subprocess.check_output(f"{esearch_query} | {efetch_query}", shell=True).decode("utf-8")
            #check if the result is empty
            if not fa_result.strip():
                print("no data found for the given queyr, please try agian")
                attempts += 1
                continue

            #save the fasta result into a file
            with open(f"{species_name}.fasta", mode="w") as f:
                f.write(fa_result)
            return fa_result, protein_name, species_name
        
        except Exception as e:
            print("unexpected error:",e)
            print("please try again")
            attempts += 1
    print(f"Maximum attempts ({attempts}) reached. Exiting the script.")
    sys.exit(1)

fa_result, protein_name, species_name = fasta_result()
   

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

plot_conserve(species_name)

#find motif from PROSTIE database that associated with the protein sequences
#split and save individual sequence to independent FASTA files
#create a directory for fasta files
os.chdir(current_dir)
if not os.path.exists(f"{current_dir}/fasta_files"):
    os.mkdir("./fasta_files")
#split sequence
with open(f"./conservation_files/{species_name}.fasta", mode="r") as f:
    fa_result = f.readlines()
    acc_list = []
    for i in range(len(fa_result)):
        line = fa_result[i].strip()
        if line.startswith(">"):
            #use accession number as file name
            acc_id = line.split(" ")[0]
            acc_id = acc_id[1:]
            if "|" in acc_id:
                acc_id = acc_id.split("|")[1]
            acc_list.append(acc_id)
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

#use for loop to scan each sequnce
#create a directory for scanning result
if not os.path.exists(f"{current_dir}/motif_files"):
    os.mkdir("./motif_files")

hit_acc_list = []
motif_dict = {}
for acc_id in acc_list:
    motif_query = f"patmatmotifs -sequence ./fasta_files/{acc_id}.fasta -outfile ./motif_files/{acc_id}.txt"
    subprocess.run(motif_query, shell=True, stdout=subprocess.DEVNULL)

    #only save the file with motif result
    with open(f"./motif_files/{acc_id}.txt", mode="r") as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if line.find("HitCount") != -1:
            #save the hit count number
            hit_num = line.split(" ")[-1]
            hit_num = int(hit_num)
            if hit_num > 0:
                #save the acc number if find motif from this sequence
                hit_acc_list.append(acc_id)
                print(acc_id, hit_num)
            else:
                os.remove(f"./motif_files/{acc_id}.txt")
        if line.find("Motif") != -1:
            motif = line.split(" = ")[-1]
            motif = motif.strip()
            if motif not in motif_dict.keys():
                motif_dict.update({motif:[acc_id]})
            else:
                motif_dict[motif].append(acc_id)
                

# print the final number of motif that found
print(f"total sequences scanned: {len(acc_list)}")
print(f"sequences with known motifs: {len(hit_acc_list)}")
print(f"Associated motifs ({len(motif_dict.keys())} found): ")
for motif in motif_dict.keys():
    seq_num = len(motif_dict[motif])
    print(f" - {motif} ({seq_num} sequences)")