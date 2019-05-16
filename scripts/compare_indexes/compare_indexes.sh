if [ $# -ne 5 ]; then
    echo "Usage: bash compare_indexes.sh <path_to_build/bin/compare_indexes> <path_to_idx_1> <path_to_kmer_prgs_of_idx_1> <path_to_idx_2> <path_to_kmer_prgs_of_idx_2>"
    exit 1
fi

set -eux
$1 $2 $4
diff -rq $3 $5
echo "If you are seeing this message, then the given indexes are identical!"

