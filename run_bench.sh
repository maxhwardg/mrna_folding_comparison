python3 pysrc/benchmark.py --mode=random --codon_table=data/homosapiens.txt --bin_root=extern > data/bench_runs/random1.txt
python3 pysrc/benchmark.py --mode=random --codon_table=data/homosapiens.txt --bin_root=extern > data/bench_runs/random2.txt
python3 pysrc/benchmark.py --mode=mll --codon_table=data/homosapiens.txt --bin_root=extern > data/bench_runs/mll1.txt
python3 pysrc/benchmark.py --mode=mll --codon_table=data/homosapiens.txt --bin_root=extern > data/bench_runs/mll2.txt
python3 pysrc/benchmark.py --mode=random --codon_table=data/homosapiens.txt --bin_root=extern > data/bench_runs/random3.txt
python3 pysrc/benchmark.py --mode=mll --codon_table=data/homosapiens.txt --bin_root=extern > data/bench_runs/mll3.txt