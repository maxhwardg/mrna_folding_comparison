"""Checks that mRNAfold is consistent with ViennaRNA and the codon frequency table."""
from bridge import call_mrnafold
import protein
import vienna

def main():
    eps = 1e-3
    aa_len = 30
    cft = protein.CodonFrequencyTable('../data/homosapiens.txt')
    for _ in range(1000):
        aa_seq = protein.random_aa_seq(aa_len)
        print(aa_seq)
        res = call_mrnafold("../extern/mrnafold-main", aa_seq, lambda_value=0.0)
        assert len(res.rna_seq) == len(aa_seq)*3, f'{len(res.rna_seq)} != {len(aa_seq)*3}'
        assert len(res.rna_seq) == len(res.db), f'{len(res.rna_seq)} != {len(res.db)}'
        cds = protein.rna_to_cds(res.rna_seq)
        for i, codon in enumerate(cds):
            assert codon in cft.get_codons(aa_seq[i]), f'{codon} not in {cft.get_codons(aa_seq[i])} for amino acid {aa_seq[i]} (index {i})'
        print(res)
        ctx = vienna.ViennaContext(res.rna_seq, dangles=0)
        vienna_fe = ctx.free_energy(res.db)
        vienna_mfe = ctx.free_energy(ctx.mfe())
        print(res.mfe, vienna_fe, vienna_mfe)
        assert abs(vienna_fe - res.mfe) < eps, f'{vienna_fe} != {res.mfe}'
        assert abs(vienna_mfe - res.mfe) < eps, f'{vienna_mfe} != {res.mfe}'
        ref_cai = cft.codon_adaptation_index(cds)
        print(res.cai, ref_cai)
        assert abs(ref_cai - res.cai) < eps, f'{ref_cai} != {res.cai}'
        

if __name__ == '__main__':
    main()
