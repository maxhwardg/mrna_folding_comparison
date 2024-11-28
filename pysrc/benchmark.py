from cdsfold_bridge import call_cdsfold
from derna_bridge import call_derna
from lineardesign_bridge import call_lineardesign
import protein
    

def main():
    cft = protein.CodonFrequencyTable('../data/homosapiens.txt')
    for aa_len in range(50, 1000, 50):
        aa_seq = protein.random_aa_seq(aa_len)
        print('aa_len:', aa_len)
        cds_res = call_cdsfold("../extern/CDSfold-main", aa_seq)
        linear_res = call_lineardesign(cft, "../extern/LinearDesign-main/", aa_seq)
        derna_res = call_derna(cft, "../extern/derna-main", aa_seq)
        print('cdsfold time(s):', cds_res.time_s)
        print('lineardesign time(s):', linear_res.time_s)
        print('derna time(s):', derna_res.time_s)
        
        
        

if __name__ == '__main__':
    main()
