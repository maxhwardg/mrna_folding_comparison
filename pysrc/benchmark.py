from bridge import call_cdsfold, call_mrnafold, call_lineardesign, call_derna, FoldException
import protein


def main():
    cft = protein.CodonFrequencyTable("../data/homosapiens.txt")

    derna_mx = 0
    cds_mx = 0
    timeout_s = 3600
    random_seq = False
    
    def gen_ml_seq(sz: int) -> str:
        return "M" + "L" * (sz - 1)
    
    def gen_seq(sz: int) -> str:
        return protein.random_aa_seq(sz) if random_seq else gen_ml_seq(sz)

    for aa_len in range(50, 1501, 50):
        print("aa_len:", aa_len)

        # LinearDesign tends to fail an assert and crash, so retry until it works
        while True:
            aa_seq = gen_seq(aa_len)
            try:
                linear_res = call_lineardesign(
                    cft, "../extern/LinearDesign-main/", aa_seq
                )
            except FoldException as e:
                if random_seq:
                    print(f"LinearDesign failed with error: {e}. Trying again.")
                    continue
                else:
                    # Can't retry with deterministic sequences
                    raise e
            break

        print("lineardesign time(s):", linear_res.time_s)
        print("lineardesign memory(bytes):", linear_res.memory_bytes, flush=True)

        if cds_mx < timeout_s:
            cds_res = call_cdsfold("../extern/CDSfold-main", aa_seq)
            print("cdsfold time(s):", cds_res.time_s)
            print("cdsfold memory(bytes):", cds_res.memory_bytes, flush=True)
            cds_mx = max(cds_mx, cds_res.time_s)

        if derna_mx < timeout_s:
            derna_res = call_derna(cft, "../extern/derna-main", aa_seq, lambda_value=1.0)
            print("derna time(s):", derna_res.time_s)
            print("derna memory(bytes):", derna_res.memory_bytes, flush=True)
            derna_mx = max(derna_mx, derna_res.time_s)

        mrna_res = call_mrnafold("../extern/mrnafold-main", aa_seq, parallel=True)
        print("mrnafold time(s):", mrna_res.time_s)
        print("mrnafold memory(bytes):", mrna_res.memory_bytes, flush=True)


if __name__ == "__main__":
    main()
