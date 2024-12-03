"""Bridge code for calling mrna folding algorithms"""
from dataclasses import dataclass
from typing import Sequence
import subprocess
import threading
import tempfile
import time
import os
import psutil
import protein


class FoldException(Exception):
    """Exception class for mrna folding algorithms"""


@dataclass
class FoldResult:
    """Dataclass for storing the results of an mrna folding algorithm"""
    rna_seq: str = ''
    db: str = ''
    mfe: float = float('inf')
    cai: float = -1.0
    time_s: float = -1.0
    memory_bytes: int = -1


def monitor_memory_usage(pid: int, memory_log: list[int]):
    """Monitor memory usage of the process with the given PID. Appends the memory usage to memory_log"""
    try:
        process = psutil.Process(pid)
        while process.is_running():
            memory_log.append(process.memory_info().rss)
            time.sleep(0.1)
    except psutil.NoSuchProcess:
        pass


def call_subprocess(args: Sequence[str], input_str: str = '') -> tuple[subprocess.CompletedProcess, int, float]:
    """
    Calls a subprocess with a timeout and commandline input.
    Returns a tuple of the CompletedProcess, the memory usage in bytes, and the time taken in seconds.
    """
    ts = time.time()
    p = subprocess.Popen(
        args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    # Launch background thread to monitor memory usage
    memory_log = []
    monitor_thread = threading.Thread(
        target=monitor_memory_usage, args=(p.pid, memory_log))
    monitor_thread.start()
    # Send stdin
    stdout, stderr = p.communicate(input=input_str)
    # Wait for the process to finish
    error_code = p.wait()
    te = time.time()
    return subprocess.CompletedProcess(args, error_code, stdout, stderr), max(memory_log), te-ts


def make_fasta_file(aa_seq: str) -> str:
    """Creates a fasta file in the format"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write('>seq\n')
        file.write(aa_seq)
    return file.name


def call_cdsfold(path: str, aa_seq: str) -> FoldResult:
    """Calls CdsFold via a subprocess"""
    fasta = make_fasta_file(aa_seq)
    # result = subprocess.run([os.path.join(
    #     path, 'src/CDSfold'), fasta], capture_output=True, text=True, check=False)
    result, mem_b, time_s = call_subprocess(
        [os.path.join(path, 'src/CDSfold'), fasta])

    # Clean up tmp file
    os.remove(fasta)

    if result.returncode != 0:
        raise FoldException(f"CdsFold failed with return code: "
                            f"{result.returncode}, and stderror: {result.stderr}")
    ret = FoldResult(time_s=time_s, memory_bytes=mem_b)
    lns = result.stdout.split('\n')
    ret.rna_seq = lns[-5]
    ret.db = lns[-4]
    ret.mfe = float(lns[-3].split('MFE:')[1].split(' kcal/mol')[0])
    return ret


def make_derna_cft_csv(cft: protein.CodonFrequencyTable) -> str:
    """Creates a csv file in the format required by DENRA"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        for aa in protein.AA_SINGLE_LETTER.values():
            # DERNA doesn't like stop codons
            if aa == '*':
                continue
            codons = cft.get_codons(aa)
            codon_row = ['']
            codon_row.extend(codons)
            file.write(','.join(codon_row) + '\n')
            aa_row = [aa] + [0]*6
            for i, cdn in enumerate(codons):
                aa_row[i+1] = cft.get_codon_freq(cdn)
            aa_row.append('')
            file.write(','.join(map(str, aa_row)) + '\n')

    return file.name


def call_derna(cft: protein.CodonFrequencyTable, path: str, aa_seq: str, lambda_value: float = 1.0) -> FoldResult:
    """Calls DERNA via a subprocess"""
    csv_cft = make_derna_cft_csv(cft)
    fasta = make_fasta_file(aa_seq)
    fname = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        result, mem_b, time_s = call_subprocess([os.path.join(
            path, 'build/derna'), '-c', csv_cft, '-i', fasta,
                                                 '-o', file.name,
                                                 '-m', '1',
                                                 '-s', '2',
                                                 '-l', str(lambda_value)])
        fname = file.name

    # Delete tmp input files
    os.remove(csv_cft)
    os.remove(fasta)
    with open(fname, 'r', encoding='utf-8') as file:
        lns = file.readlines()
    # Delete tmp output file
    os.remove(fname)
    # Delete garbage DERNA creates
    os.remove("dd.txt")

    if result.returncode != 0:
        os.remove(fname)
        raise FoldException(f'DERNA failed with return code: '
                            f'{result.returncode}, and stderror: {result.stderr}')

    res = FoldResult(time_s=time_s, memory_bytes=mem_b)
    for ln in lns:
        if ln.startswith('zuker cai rna:'):
            res.rna_seq = ln[len('zuker cai rna: '):].split('.size')[0].strip()
        elif ln.startswith('Codon Adaptation Index: '):
            res.cai = float(ln[len('Codon Adaptation Index: '):].strip())
        elif ln.startswith('Minimum Free Energy: '):
            res.mfe = float(ln[len('Minimum Free Energy: '):].strip())
        elif ln.startswith('zuker cai bp: '):
            res.db = ln[len('zuker cai bp: '):].split(',size')[0].strip()
    return res


def make_linear_design_cft_csv(cft: protein.CodonFrequencyTable) -> str:
    """Creates a csv file in the format required by LinearDesign"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write('#,,\n')
        for aa in protein.AA_SINGLE_LETTER.values():
            for codon in cft.get_codons(aa):
                file.write(f'{codon},{aa},{cft.get_codon_freq(codon)}\n')
    return file.name


def call_lineardesign(cft: protein.CodonFrequencyTable, path: str, aa_seq: str, lambda_value: float = 0.0) -> FoldResult:
    """Calls LinearDesign via a subprocess"""
    # Change the working directory to the LinearDesign directory since it can't load .so files otherwise
    orig_path = os.getcwd()
    csv_cft = make_linear_design_cft_csv(cft)
    os.chdir(path)
    result, mem_b, time_s = call_subprocess(
        ['bin/LinearDesign_2D', '0', str(lambda_value), csv_cft], input_str=aa_seq)
    os.chdir(orig_path)

    # Clean up tmp file
    os.remove(csv_cft)

    if result.returncode != 0:
        raise FoldException(f'LinearDesign failed with return code: '
                            f'{result.returncode}, and stderror: {result.stderr}')

    # Parse result
    res = FoldResult(time_s=time_s, memory_bytes=mem_b)
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


def make_mrnafold_config(aa_seq: str, parallel: bool, lambda_value: float) -> str:
    """Creates a folding configuration file in the format required by MrnaFold"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as file:
        file.write(f'aa_seq {aa_seq}\n')
        file.write(f'parallel {str(parallel).lower()}\n')
        file.write(f'lambda {lambda_value}\n')
        file.write('num_subopt_traces 1\n')
        file.write('banned_motif_tests 0\n')
    return file.name


def call_mrnafold(path: str, aa_seq: str, parallel: bool = True, lambda_value: float = 0.0) -> FoldResult:
    """Calls mRNAFold via a subprocess"""
    fname = make_mrnafold_config(aa_seq, parallel, lambda_value)
    result, mem_b, time_s = call_subprocess(
        [os.path.join(path, 'build/exe/fold_codon_graph'), fname])

    # Clean up tmp file
    os.remove(fname)

    if result.returncode != 0:
        raise FoldException(f'mRNAFold failed with return code: {
                            result.returncode}, and stderror: {result.stderr}')

    # Parse result
    res = FoldResult(time_s=time_s, memory_bytes=mem_b)
    lns = result.stdout.split("\n")
    res.rna_seq = lns[0]
    res.db = lns[1]
    res.cai = float(lns[3].split('CAI: ')[1])
    res.mfe = float(lns[4].split('MFE: ')[1])

    return res


def main():
    aa_str = "MLVLVLVLVL"
    print(call_cdsfold("../extern/CDSfold-main", aa_str))
    cft = protein.CodonFrequencyTable(
        "../data/homosapiens.txt")
    print(call_derna(cft, "../extern/derna-main", aa_str))
    print(call_lineardesign(cft, "../extern/LinearDesign-main/", aa_str))
    print(call_mrnafold("../extern/mrnafold-main/", aa_str))


if __name__ == "__main__":
    main()
