import pandas as pd
import random

# male dna 

with open("male_dna.txt") as file:
    male_dna = file.read().strip()

# female dna

with open("female_dna.txt")as file:
    female_dna = file.read().strip()

#simulation of child dna

def simulate_child_dna(male_dna, female_dna):
    child_dna = ''
    for i in range(len(male_dna)):
        if random.random() < 0.5:
            child_dna += male_dna[i]
        else:
            child_dna += female_dna[i]

    return child_dna
    
child_dna = simulate_child_dna(male_dna, female_dna)
print("Simulated Child DNA:", child_dna)

# child dna data

def child_dna_stats(dna):

    gc = dna.count("G") + dna.count("C")
    at = dna.count("A") + dna.count("T")

    stats = {
        "A": dna.count("A"),
        "T": dna.count("T"),
        "C": dna.count("C"),
        "G": dna.count("G"),
        "GC_Content": (gc / len(dna)) * 100,
        "AT_Content": (at / len(dna)) * 100
    }

    return stats

stats = child_dna_stats(child_dna)

print(stats)

# mutation rate
def mutate_dna(dna, mutation_rate):
    mutated_dna = ''
    for base in dna:
        if random.random() < mutation_rate:
            mutated_dna += random.choice("ATCG")
        else:
            mutated_dna += base
    return mutated_dna
mutation_rate = 0.01
mutated_child_dna = mutate_dna(child_dna, mutation_rate)
print("Mutated Child DNA:", mutated_child_dna)

# comparing original and mutated dna
def compare_dna(original, mutated):
    differences = sum(1 for o, m in zip(original, mutated) if o != m)
    return differences
difference = compare_dna(child_dna, mutated_child_dna)
print("Number of differences between original and mutated DNA:", difference)

# mutation rate percentage
mutation_percentage = (difference / len(child_dna))*100
print("mutation percentage:", mutation_percentage)

# child dna to rna 
def dna_to_rna(dna):
    rna = dna.replace("T", "U")
    return rna
child_rna = dna_to_rna(child_dna)
print("Child RNA:", child_rna)

# child rna to proteins

from Bio.Seq import Seq
protein = Seq(child_rna).translate()
print("Child Protein Sequence:", protein)