#! /bin/bash -x

# USING THE MOST RECENT mer_counts.jf, this script generates all k-mers occurring with the specified number. Jellyfish then dumps these
# into a .fasta format, and then they are located on the reference genome using Smalt.
# This script should be run in the directory containing the reference 

NUM_OCCS=$1
NAME=$2
EXTENSION=$3

echo "Started for "$NUM_OCCS
/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish dump -L $NUM_OCCS -U $NUM_OCCS -ct "/lustre/scratch110/sanger/gh10/Data/mer_counts/"$NAME"_mer_counts.jf" > "generate_occurrence_locations_"$NUM_OCCS".tmp.dump.fasta"

cat "generate_occurrence_locations_"$NUM_OCCS".tmp.dump.fasta" | awk '{print ">try.dat " " \n"$1}' > generate_occurrence_locations.tmp.dat
/software/hpag/icas/0.61/icas/bin/rename_fastq -name kmer_reads -len 10 generate_occurrence_locations.tmp.dat generate_occurrence_locations.tmp.fastq

# Generate index if requried (but hopefully will already be there)
if [ ! -f $NAME"_hash_file"* ]; then
	echo "Creating hash file..."
	/software/hpag/bin/smalt-0.7.4 index -k 17 -s 17 $NAME"_hash_file" $NAME"."$EXTENSION
	echo "Finished creating hash file"
fi

/software/hpag/bin/smalt-0.7.4 map -m 20 -f ssaha -n 4 -O -d -0 $NAME"_hash_file" generate_occurrence_locations.tmp.fastq > "/lustre/scratch110/sanger/gh10/Data/k_mer_locations_and_words/"$NAME"_"$NUM_OCCS"_occs.tmp.ssaha" 
sed -i "s|$| $NUM_OCCS|" "/lustre/scratch110/sanger/gh10/Data/k_mer_locations_and_words/"$NAME"_"$NUM_OCCS"_occs.tmp.ssaha" 
rm generate_occurrence_locations.tmp.dat generate_occurrence_locations.tmp.fastq

echo "Finished for "$NUM_OCCS
