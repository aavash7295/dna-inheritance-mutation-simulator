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

# compare mutated protein to orginal protein
mutated_rna = dna_to_rna(mutated_child_dna)
mutated_protein = Seq(mutated_rna).translate()
print("Mutated Protein Sequence:", mutated_protein)

def compare_proteins(original, mutated):
    differnces = sum(1 for o, m in zip(original, mutated) if o != m and m != "*")
    return differnces
protein_difference = compare_proteins(protein, mutated_protein)
print("Number of differences between original and mutated protein:", protein_difference )
plt.savefig("mutation_histogram.png", dpi=300)

# mutation positions
def mutation_positions(original, mutated):
    positions = [i for i, (o, m) in enumerate(zip(original, mutated)) if o!= m and m!= "*"]
    return positions
mutation_positions_list = mutation_positions(child_dna, mutated_child_dna)
print("Mutation Positions:", mutation_positions_list)

#mutation frequency plot
import matplotlib.pyplot as plt
plt.hist(mutation_positions_list, bins=20, edgecolor='black')
plt.title("Mutation Position Distribution") 
plt.xlabel("Position in DNA Sequence")
plt.ylabel("Frequency")
plt.show()

# original vs mutation base count plot
def base_count(dna):
    return {
        "A": dna.count("A"),
        "T": dna.count("T"),
        "C": dna.count("C"),
        "G": dna.count("G")
    }
original_base_count = base_count(child_dna)
mutated_base_count = base_count(mutated_child_dna)
labels = original_base_count.keys()
original_counts = original_base_count.values()
mutated_counts = mutated_base_count.values()
x = range(len(labels))
plt.bar(x, original_counts, width=0.4, label="Original", align='center')
plt.bar(x, mutated_counts, width=0.4, label="Mutated", align='edge')
plt.xticks(x, labels)   
plt.title("Base Count Comparison")
plt.xlabel("Base")
plt.ylabel("Count")
plt.legend()
plt.show()
plt.savefig("base_count_comparison.png", dpi=300)

#similarity score between original and mutated dna
def similarity_score(original, mutated):
    matches = sum(1 for o, m in zip(original, mutated) if o== m and m!= "*")
    score = (matches / len(original))*100
    return score
similarity_score = similarity_score(child_dna, mutated_child_dna)
print("Similarity Score between original and mutated DNA:", similarity_score)

#fasta file support 
def read_fasta(file_path):
    with open(file_path) as file:
        lines = file.readlines()
        header = lines[0].strip()
        sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence
def write_fasta(header, sequence, file_path):
    with open(file_path, 'w') as file:
        file.write(f">{header}\n")
        for i in range(0, len(sequence), 60):
            file.write(sequence[i:i+60] + '\n')






# =========================
# SAVE ALL OUTPUTS
# =========================

# Main analysis report
with open("child_dna_analysis.txt", "w") as file:

    file.write("========== CHILD DNA ANALYSIS ==========\n\n")

    file.write("SIMULATED CHILD DNA:\n")
    file.write(child_dna + "\n\n")

    file.write("CHILD DNA STATISTICS:\n")

    for key, value in stats.items():

        if isinstance(value, float):
            file.write(f"{key}: {value:.2f}\n")
        else:
            file.write(f"{key}: {value}\n")

    file.write("\n")

    file.write("MUTATED CHILD DNA:\n")
    file.write(mutated_child_dna + "\n\n")

    file.write("DNA COMPARISON:\n")
    file.write(f"Number of differences: {difference}\n")
    file.write(f"Mutation percentage: {mutation_percentage:.2f}%\n")
    file.write(f"Similarity score: {similarity_score:.2f}%\n\n")

    file.write("RNA SEQUENCE:\n")
    file.write(child_rna + "\n\n")

    file.write("PROTEIN SEQUENCE:\n")
    file.write(str(protein) + "\n\n")

    file.write("MUTATED PROTEIN SEQUENCE:\n")
    file.write(str(mutated_protein) + "\n\n")

    file.write("PROTEIN DIFFERENCES:\n")
    file.write(f"Protein differences: {protein_difference}\n\n")

    file.write("MUTATION POSITIONS:\n")
    file.write(str(mutation_positions_list) + "\n\n")


# Save RNA separately
with open("rna_sequence.txt", "w") as file:
    file.write(child_rna)


# Save protein separately
with open("protein_sequence.txt", "w") as file:
    file.write(str(protein))


# Save mutated protein separately
with open("mutated_protein_sequence.txt", "w") as file:
    file.write(str(mutated_protein))


# Save FASTA files
write_fasta("Child_DNA", child_dna, "child_dna.fasta")

write_fasta(
    "Mutated_Child_DNA",
    mutated_child_dna,
    "mutated_child_dna.fasta"
)

print("\nAll output files have been saved successfully.")
