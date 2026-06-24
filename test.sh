build/tests/generate_known_relationships --snps 20000 --seed 42 --prefix kit_20k
bin/fastlckin relatedness -p tests/data/kit_20k --no-ld-prune --classify -o tests/data/kit_20k_output.tsv

python3 tests/validate_ground_truth.py --fastlckin bin/fastlckin --data-dir tests/data --skip-generation
python tests/validate_ground_truth.py --fastlckin bin/fastlckin --data-dir tests/data
