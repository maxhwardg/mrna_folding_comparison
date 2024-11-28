
import subprocess
import os
from dataclasses import dataclass
import tempfile
import time


class CdsFoldException(Exception):
    pass


@dataclass
class CdsFoldResult:
    rna_seq: str = ''
    db: str = ''
    mfe: float = 0.0
    time_s: float = 0.0


def make_cdsfold_fasta(aa_seq: str) -> str:
    """Creates a fasta file in the format required by DENRA"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write('>seq\n')
        file.write(aa_seq)
    return file.name


def call_cdsfold(path: str, aa_seq: str) -> CdsFoldResult:
    """Calls CdsFold via a subprocess"""
    fasta = make_cdsfold_fasta(aa_seq)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        ts = time.time()
        result = subprocess.run([os.path.join(
            path, 'src/CDSfold'), fasta], capture_output=True, text=True, check=False)
        te = time.time()
        if result.returncode != 0:
            raise CdsFoldException(f'CdsFold failed with return code: {result.returncode}, and stderror: {result.stderr}')
    ret = CdsFoldResult(time_s=te-ts)
    lns = result.stdout.split('\n')
    ret.rna_seq = lns[-5]
    ret.db = lns[-4]
    ret.mfe = float(lns[-3].split('MFE:')[1].split(' kcal/mol')[0])
    return ret


def main():
    print(call_cdsfold("../extern/CDSfold-main", "MLLLLV"))


if __name__ == "__main__":
    main()
