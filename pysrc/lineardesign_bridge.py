
import subprocess
import os
from dataclasses import dataclass
import tempfile
import time
import protein


class LinearDesignException(Exception):
    pass

@dataclass
class LinearDesignResult:
    rna_seq: str = ''
    db: str = ''
    mfe: float = 0.0
    cai: float = 0.0
    time_s: float = 0.0
    
    
def make_linear_design_cft_csv(cft: protein.CodonFrequencyTable) -> str:
    """Creates a csv file in the format required by LinearDesign"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write('#,,\n')
        for aa in protein.AA_SINGLE_LETTER.values():
            for codon in cft.get_codons(aa):
                file.write(f'{codon},{aa},{cft.get_codon_freq(codon)}\n')
    return file.name


def call_lineardesign(cft: protein.CodonFrequencyTable, path: str, aa_seq: str, lambda_value:float = 0.0) -> LinearDesignResult:
    """Calls LinearDesign via a subprocess"""
    # Change the working directory to the LinearDesign directory since it can't load .so files otherwise
    orig_path = os.getcwd()
    csv_cft = make_linear_design_cft_csv(cft)
    os.chdir(path)
    ts = time.time()
    result = subprocess.run(['bin/LinearDesign_2D', '0', str(lambda_value), csv_cft], input=aa_seq, capture_output=True, text=True, check=False)
    te = time.time()
    os.chdir(orig_path)
    if result.returncode != 0:
        raise LinearDesignException(f'LinearDesign failed with return code: {result.returncode}, and stderror: {result.stderr}')
    
    # Parse result
    res = LinearDesignResult(time_s = te-ts)
    for ln in result.stdout.split("\n"):
        if ln.startswith('mRNA sequence:  '):
            res.rna_seq = ln[len('mRNA sequence:  '):]
        elif ln.startswith('mRNA structure: '):
            res.db = ln[len('mRNA structure: '):]
        elif ln.startswith('mRNA folding free energy: '):
            toks = ln[len('mRNA folding free energy: '):].split('; mRNA CAI: ')
            res.mfe = float(toks[0].split(' kcal/mol')[0].strip())
            res.cai = float(toks[1].strip())
    return res
            
            
    
def main():
    cft = protein.CodonFrequencyTable('../data/homosapiens.txt')
    print(call_lineardesign(cft, "../extern/LinearDesign-main/", "MLLLLV"))
    
if __name__ == "__main__":
    main()
    