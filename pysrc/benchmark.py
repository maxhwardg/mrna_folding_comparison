from cdsfold_bridge import call_cdsfold
from derna_bridge import call_derna
from lineardesign_bridge import call_lineardesign
from mrnafold_bridge import call_mrnafold
import protein
    

def main():
    cft = protein.CodonFrequencyTable('../data/homosapiens.txt')
    for aa_len in range(50, 2000, 50):
        aa_seq = protein.random_aa_seq(aa_len)
        print('aa_len:', aa_len)
        cds_res = call_cdsfold("../extern/CDSfold-main", aa_seq)
        linear_res = call_lineardesign(cft, "../extern/LinearDesign-main/", aa_seq)
        derna_res = call_derna(cft, "../extern/derna-main", aa_seq)
        mrna_res = call_mrnafold("../extern/mrnafold-main", aa_seq, parallel=True)
        print('cdsfold time(s):', cds_res.time_s)
        print('lineardesign time(s):', linear_res.time_s)
        print('derna time(s):', derna_res.time_s)
        print('mrnafold time(s):', mrna_res.time_s)
        # For the output to stop buffering
        print('', flush=True)
        
        
        

if __name__ == '__main__':
    main()
