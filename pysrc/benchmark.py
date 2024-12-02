from cdsfold_bridge import call_cdsfold
from derna_bridge import call_derna
from lineardesign_bridge import call_lineardesign, LinearDesignException
from mrnafold_bridge import call_mrnafold
import protein


def main():
    cft = protein.CodonFrequencyTable("../data/homosapiens.txt")

    derna_mx = 0
    cds_mx = 0
    threshold = 1800

    for aa_len in range(50, 2000, 50):
        print("aa_len:", aa_len)

        # LinearDesign tends to fail an assert and crash, so retry until it works
        while True:
            aa_seq = protein.random_aa_seq(aa_len)
            try:
                linear_res = call_lineardesign(
                    cft, "../extern/LinearDesign-main/", aa_seq
                )
            except LinearDesignException as e:
                print(f"LinearDesign failed with error: {e}. Trying again.")
                continue
            break

        print("lineardesign time(s):", linear_res.time_s)

        if cds_mx < threshold:
            cds_res = call_cdsfold("../extern/CDSfold-main", aa_seq)
            print("cdsfold time(s):", cds_res.time_s)
            cds_mx = max(cds_mx, cds_res.time_s)

        if derna_mx < threshold:
            derna_res = call_derna(cft, "../extern/derna-main", aa_seq)
            print("derna time(s):", derna_res.time_s)
            derna_mx = max(derna_mx, derna_res.time_s)

        mrna_res = call_mrnafold("../extern/mrnafold-main", aa_seq, parallel=True)
        print("mrnafold time(s):", mrna_res.time_s)
        # For the output to stop buffering
        print("", flush=True)


if __name__ == "__main__":
    main()
