'''
Author: Bharat Panera
Date: 23-07-2023

----------Dependencies----------
1. python3
2. logging
3. logging

--------------Usage-------------
python3 script.py
'''
import requests
from random import shuffle, choice, random
import logging

# create logger to print the logs on terminal
logging.basicConfig(format='%(asctime)s - %(levelname)s : %(message)s', level=logging.INFO)
# logging.disable()

#function to download  fasta file from URL
def download_url(url, file_name):
    '''Downloading the fasta file from the given URL'''
    url_data = requests.get(url, allow_redirects=True)
    open(file_name, 'wb').write(url_data.content)
    return print("Download completed")

#function to store the fasta file into python dictionary
def read_fasta_file(file_name):
    '''Read the fasta file and store it in python dictionary where fasta 
       header will be key and fasta sequence will be the value of the respective key'''
    fasta_dict = {}
    current_id = None
    current_sequence = ""

    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                # If the current_id is not None, add the sequence to the dictionary
                if current_id is not None:
                    fasta_dict[current_id] = current_sequence
                # Update the current_id and reset the current_sequence
                current_id = line   # Remove the ">" symbol
                current_sequence = ""
            else:
                current_sequence += line

        # Add the last sequence to the dictionary (if present)
        if current_id is not None:
            fasta_dict[current_id] = current_sequence
    return fasta_dict

#function to perform reverse compliment on a seq
def reverse_compliment(seq):
    '''Performing reverse compliment on a extracted fasta sequence'''

    #reverse the sequence
    reversed_fasta = seq[::-1]

    #replace the bases
    reversed_comp= reversed_fasta.replace("A","t").replace("T","a").replace("G","c").replace("C","g")
    reversed_comp = reversed_comp.upper()
    with open("reverse_compliment.txt", 'w') as output_file:
        output_file.write(reversed_comp)

    return reversed_comp.upper()

#function to replace a nucleotide from with random base in seq
def random_replacement(seq):
    '''Replacing a random base from the extracted fasta sequence'''
    #random bases which will be replaced
    bases = ["C", "T", "G"]
    # base_to_change = ["A", "T", "G", "C"]
    base_to_change = "A"
    shuffle(bases)
    base_to_change = choice(base_to_change)
    # print("before: " + str(seq))
    if base_to_change in seq:
        # index_to_replace = str(seq.index(base_to_change))
        index_to_replace = seq.index(base_to_change)
        # index_to_replace = choice(index_to_replace)
        new_base = choice(bases)
        seq = seq[:index_to_replace] + new_base + seq[index_to_replace + 1:]
    # print("after: " + str(seq))
    return seq

#this function will update the multifasta file with processed seq
def update_multifasta(file_name, fasta_file_data):
    '''updating and writing back the processed subsequence as a multi_fasta.fa file'''
    with open(file_name, 'w') as file:
        for fasta_header, sequence in fasta_file_data.items():
            file.write(f'{fasta_header}\n')
            file.write(f'{sequence}\n')
    return None


file_name = "multi_fasta.fa"
# file_name = "small_fasta_file.fn"
url = "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Vieuvirus.fn"

get_fasta_data = download_url(url, file_name)
logging.info("fasta file has been downloaded")

fasta_file = read_fasta_file(file_name)
logging.info("fasta file is parsed")

extracted_seq = fasta_file['>gi|768216015|gb|KP861230.1|_Acinetobacter_phage_YMC11/11/R3177,_complete_genome']
logging.info("One Sequence is extracted")

reverse_compliment_seq = reverse_compliment(extracted_seq)
logging.info("Reverse compliment is completed ")
# print(reverse_compliment_seq)

processed_subsequence = random_replacement(extracted_seq)
# print(processed_subsequence)

fasta_file['>gi|768216015|gb|KP861230.1|_Acinetobacter_phage_YMC11/11/R3177,_complete_genome'] = processed_subsequence
# print(fasta_file)

update_primary_fasta = update_multifasta(file_name, fasta_file)