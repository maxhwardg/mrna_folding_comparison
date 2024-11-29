
import subprocess
import os
from dataclasses import dataclass
import tempfile
import time
import protein


class MrnaFoldException(Exception):
    pass

@dataclass
class MrnaFoldResult:
    rna_seq: str = ''
    db: str = ''
    mfe: float = 0.0
    cai: float = 0.0
    time_s: float = 0.0
    

    
def make_mrnafold_config(aa_seq: str, parallel:bool, lambda_value: float) -> str:
    """Creates a folding configuration file in the format required by MrnaFold"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write(f'aa_seq {aa_seq}\n')
        file.write(f'parallel {str(parallel).lower()}\n')
        file.write(f'lambda {lambda_value}\n')
        file.write('num_subopt_traces 1\n')
        file.write('banned_motif_tests 0\n')
    return file.name


def call_mrnafold(path: str, aa_seq: str, parallel: bool = True, lambda_value:float = 0.0) -> MrnaFoldResult:
    """Calls mRNAFold via a subprocess"""
    fname = make_mrnafold_config(aa_seq, parallel, lambda_value)
    ts = time.time()
    result = subprocess.run([os.path.join(path, 'build/exe/fold_codon_graph'), fname], capture_output=True, text=True, check=False)
    te = time.time()
    
    # Clean up tmp file
    os.remove(fname)
    
    if result.returncode != 0:
        raise MrnaFoldException(f'mRNAFold failed with return code: {result.returncode}, and stderror: {result.stderr}')
    
    # Parse result
    res = MrnaFoldResult(time_s = te-ts)
    lns = result.stdout.split("\n")
    res.rna_seq = lns[0]
    res.db = lns[1]
    res.cai = float(lns[3].split('CAI: ')[1])
    res.mfe = float(lns[4].split('MFE: ')[1])
    
    
    return res
            
            
    
def main():
    print(call_mrnafold("../extern/mrnafold-main/", "MLLLLV"))
    
if __name__ == "__main__":
    main()
    