""" This module contains code for dealing with amino acid sequence and coding sequence (CDS) data. """
import math
import random
from typing import List

AA_SINGLE_LETTER = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "End": "*",  # Stop codon
}


def is_valid_aa_letter(aa):
    """Returns True if the given string is a valid amino acid letter, False otherwise."""
    return aa in AA_SINGLE_LETTER.values() and aa != '*'


class CodonFrequencyTable:
    """
    Class for storing codon frequency tables and calculating codon adaptation index (CAI)
    Assumes table is a text file in the format of 'CodonFrequency output in GCG Wisconsin PackageTM'
    """
    def __init__(self, table_path):
        file = open(table_path, 'r', encoding='utf-8')
        self._codons = set()
        self._codon_to_aa = {}
        self._aa_to_codons = {}
        self._codon_freq = {}
        self._aa_max_freq = {}
        # The format uses T instead of U
        self._nt_map = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}
        for line in file:
            tokens = line.strip(" \n").split()
            if len(tokens) < 3:
                continue
            aa = tokens[0]
            aa = AA_SINGLE_LETTER[aa]
            codon = ''.join([self._nt_map[nt] for nt in tokens[1]])
            freq = round(float(tokens[2]))
            self._codons.add(codon)
            self._codon_to_aa[codon] = aa
            if aa not in self._aa_to_codons:
                self._aa_to_codons[aa] = set()
            self._aa_to_codons[aa].add(codon)
            self._codon_freq[codon] = freq
            if aa not in self._aa_max_freq:
                self._aa_max_freq[aa] = 0
            self._aa_max_freq[aa] = max(self._aa_max_freq[aa], freq)
        file.close()

    def get_codon_freq(self, codon):
        return self._codon_freq[codon]

    def get_aa_max_freq(self, aa):
        return self._aa_max_freq[aa]

    def get_codons(self, aa) -> set[str]:
        return self._aa_to_codons[aa]

    # Maximum number of codons for a single amino acid
    def max_codons(self) -> int:
        return max(len(self.get_codons(aa)) for aa in self._aa_to_codons)

    def get_aa(self, codon):
        return self._codon_to_aa[codon]

    def codon_adaption_weight(self, codon):
        return self.get_codon_freq(codon) / self.get_aa_max_freq(self.get_aa(codon))

    def codon_adaptation_index(self, cds: List[str]) -> float:
        cai = 1
        for codon in cds:
            cai *= self.codon_adaption_weight(codon)
        return cai**(1/len(cds))

    def log_codon_adaptation_index(self, cds):
        cai = 0
        for codon in cds:
            cai += math.log(self.codon_adaption_weight(codon))
        return cai / len(cds)


def random_aa_seq(length: int) -> str:
    """Generates a random amino acid sequence of the given length."""
    assert length >= 2, "Length must be at least 2 to allow for a start and stop codon."
    internal = ''.join(random.choices(list(AA_SINGLE_LETTER.values())[:-1], k=length-2))
    return f'M{internal}*'

def random_cds(aa_seq, freq_table):
    cds = []
    for aa in aa_seq:
        cds.append(random.choice(list(freq_table.aa_to_codons[aa])))
    return cds

def rna_to_cds(rna_seq):
    return [rna_seq[i:i+3] for i in range(0, len(rna_seq), 3)]

def main():
    path = "../data/codon_tables/homosapiens.txt"
    ct = CodonFrequencyTable(path)
    print(ct.get_codon_freq("AUG"))

    cds = ["AUG", "GGC", "UUG", "UCC", "CGG", "AGC", "GAG"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)

    cds = ["AUG", "CCC", "CCC", "GGG", "GGC", "UUA", "AAA"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)

    cds = ["AUG", "GUU", "AAA", "GUA", "GUA", "GGG", "GUA"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)


if __name__ == "__main__":
    main()
