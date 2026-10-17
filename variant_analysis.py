#Classifying KRAS SNVs (silent / missense / nonsense ...) using the CDS of NM_004985.5.
#Input CSV needs the columns: variant_name, ref, alt

import os
import re

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

VARIANTS_CSV = r"{File Path}"  #Enter file path of csv
OUTPUT_CSV = "kras_variants_classified.csv"
ACCESSION = "NM_004985.5"


## variant classification
def classify_variant(original_aa, altered_aa, codon_number):
    if original_aa == altered_aa:
        return "Silent"
    if altered_aa == "*":
        return "Nonsense"
    if original_aa == "*":
        return "Stop_loss"
    if codon_number == 0:
        return "Start_loss"
    return "Missense"


## fetching from ncbi
Entrez.email = "example@gmail.com"
with Entrez.efetch(db="nucleotide", id=ACCESSION, rettype="gb", retmode="text") as handle:
    record = SeqIO.read(handle, "genbank")

for feature in record.features:
    if feature.type == "CDS":
        cds_seq = Seq(str(feature.extract(record.seq)))
        break

df = pd.read_csv(VARIANTS_CSV)

c_positions = []
for name in df["variant_name"]:
    name_match = re.search(r"c\.(\d+)[ACGT]>[ACGT]", str(name))
    if name_match:
        c_positions.append(int(name_match.group(1)))
    else:
        c_positions.append(None)

df["c_positions"] = c_positions

## HGVS NOMENCLATURE STRAND VALIDATION
matches = []
for position, ref, alt in zip(df["c_positions"], df["ref"], df["alt"]):
    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()

    if position is None or pd.isna(position):
        matches.append("NO_POSITION_PARSED")
        continue
    position = int(position)
    if position > len(cds_seq):
        matches.append("OUT_OF_RANGE")
        continue
    if len(ref) != 1 or len(alt) != 1:
        matches.append("NOT_SNV")
        continue

    cds_base = str(cds_seq[position - 1])
    if cds_base == ref:
        matches.append("DIRECT_MATCH")
    elif cds_base == str(Seq(ref).complement()):
        matches.append("COMPLEMENT_MATCH")
    else:
        matches.append("MISMATCH")

df["match_status"] = matches

## obtaining original and altered aa
original_aas = []
altered_aas = []
codon_numbers = []
for alt, position, match in zip(df["alt"], df["c_positions"], df["match_status"]):
    if match not in ("DIRECT_MATCH", "COMPLEMENT_MATCH"):
        original_aas.append(None)
        altered_aas.append(None)
        codon_numbers.append(None)
        continue

    position = int(position)
    alt = str(alt).strip().upper()
    coding_alt = Seq(alt).complement() if match == "COMPLEMENT_MATCH" else Seq(alt)

    ## codon position
    codon_number = (position - 1) // 3
    position_codon = (position - 1) % 3
    original_codon = cds_seq[codon_number * 3:codon_number * 3 + 3]
    alt_codon = original_codon[:position_codon] + coding_alt + original_codon[position_codon + 1:]

    ## translating
    original_aas.append(str(original_codon.translate()))
    altered_aas.append(str(alt_codon.translate()))
    codon_numbers.append(codon_number)

df["original_aa"] = original_aas
df["altered_aa"] = altered_aas

classification = []
for original, altered, codon_number in zip(original_aas, altered_aas, codon_numbers):
    if original is None:
        classification.append("Unresolved")
    else:
        classification.append(classify_variant(original, altered, codon_number))

df["classification"] = classification

print(df[["variant_name", "match_status", "original_aa", "altered_aa", "classification"]])
df.to_csv(OUTPUT_CSV, index=False)