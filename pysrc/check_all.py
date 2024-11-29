from cdsfold_bridge import call_cdsfold
from derna_bridge import call_derna
from lineardesign_bridge import call_lineardesign
from mrnafold_bridge import call_mrnafold
import protein
import vienna

def validate_res(seq, mfe, db, eps, alg):
    ctx = vienna.ViennaContext(seq, dangles=0)
    vienna_fe = ctx.free_energy(db)
    vienna_mfe = ctx.free_energy(ctx.mfe())
    assert abs(vienna_fe - mfe) < eps, f'{vienna_fe} != {mfe} for {alg}'
    assert abs(vienna_mfe - mfe) < eps, f'{vienna_mfe} != {mfe} for {alg}'
    

def main():
    eps = 1e-3
    aa_len = 30
    cft = protein.CodonFrequencyTable('../data/homosapiens.txt')
    for _ in range(1000):
        aa_seq = protein.random_aa_seq(aa_len)
        print(aa_seq)
        cds_res = call_cdsfold("../extern/CDSfold-main", aa_seq)
        linear_res = call_lineardesign(cft, "../extern/LinearDesign-main/", aa_seq)
        derna_res = call_derna(cft, "../extern/derna-main", aa_seq)
        mrna_res = call_mrnafold("../extern/mrnafold-main", aa_seq)
        
        
        print(cds_res.mfe, linear_res.mfe, derna_res.mfe, mrna_res.mfe)
        assert abs(cds_res.mfe - mrna_res.mfe) < eps, f'{cds_res.mfe} != {mrna_res.mfe}'
        assert abs(cds_res.mfe - linear_res.mfe) < eps, f'{cds_res.mfe} != {linear_res.mfe}'
        assert abs(cds_res.mfe - derna_res.mfe) < eps, f'{cds_res.mfe} != {derna_res.mfe}'
        
        validate_res(mrna_res.rna_seq, mrna_res.mfe, mrna_res.db, eps, 'mrnafold')
        validate_res(cds_res.rna_seq, cds_res.mfe, cds_res.db, eps, 'cdsfold')
        validate_res(linear_res.rna_seq, linear_res.mfe, linear_res.db, eps, 'lineardesign')
        validate_res(derna_res.rna_seq, derna_res.mfe, derna_res.db, eps, 'derna')
        
        
        

if __name__ == '__main__':
    main()
