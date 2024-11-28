
import subprocess
import os
from dataclasses import dataclass
import tempfile
import time
import protein


class DernaException(Exception):
    pass

@dataclass
class DernaResult:
    rna_seq: str = ''
    db: str = ''
    mfe: float = 0.0
    cai: float = 0.0
    time_s: float = 0.0
    
    
def make_derna_cft_csv(cft: protein.CodonFrequencyTable) -> str:
    """Creates a csv file in the format required by DENRA"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        for aa in protein.AA_SINGLE_LETTER.values():
            # DERNA doesn't like stop codons
            if aa == '*':
                continue
            codons = cft.get_codons(aa)
            codon_row = [""]
            codon_row.extend(codons)
            file.write(','.join(codon_row) + '\n')
            aa_row = [aa] + [0]*6
            for i, cdn in enumerate(codons):
                aa_row[i+1] = cft.get_codon_freq(cdn)
            file.write(','.join(map(str, aa_row)) + '\n')
                
    return file.name

def make_derna_fasta(aa_seq: str) -> str:
    """Creates a fasta file in the format required by DENRA"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write('>seq\n')
        file.write(aa_seq)
    return file.name


def call_derna(cft: protein.CodonFrequencyTable, path: str, aa_seq: str, lambda_value:float | None = 0.5) -> DernaResult:
    """Calls DERNA via a subprocess"""
    csv_cft = make_derna_cft_csv(cft)
    fasta = make_derna_fasta(aa_seq)
    fname = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        ts = time.time()
        result = subprocess.run([os.path.join(path, 'build/derna'), '-c', csv_cft, '-i', fasta, '-o', file.name, '-m', '1', '-s', '2', '-l', str(lambda_value)], capture_output=True, text=True, check=False)
        te = time.time()
        if result.returncode != 0:
            raise DernaException(f'DERNA failed with return code: {result.returncode}, and stderror: {result.stderr}')
        fname = file.name
    res = DernaResult(time_s = te-ts)
    with open(fname, 'r', encoding='utf-8') as file:
        for ln in file:
            if ln.startswith('zuker cai rna:'):
                res.rna_seq = ln[len('zuker cai rna: '):].split('.size')[0].strip()
            elif ln.startswith('Codon Adaptation Index: '):
                res.cai = float(ln[len('Codon Adaptation Index: '):].strip())
            elif ln.startswith('Minimum Free Energy: '):
                res.mfe = float(ln[len('Minimum Free Energy: '):].strip())
            elif ln.startswith('zuker cai bp: '):
                res.db = ln[len('zuker cai bp: '):].split(',size')[0].strip()
    # Delete garbage DERNA creates
    os.remove("dd.txt")
    return res
            
            
    
def main():
    cft = protein.CodonFrequencyTable('../data/homosapiens.txt')
    print(call_derna(cft, "../extern/derna-main", "MLLLLV"))
    
if __name__ == "__main__":
    main()
    