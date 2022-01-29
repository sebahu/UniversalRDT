#!/bin/sh
echo -n "" > all_mapping.txt
for rxn in $(ls -1 reaction_intermediates)
do
rxn_folder=reaction_intermediates/$rxn
sed 's/,/\n/g' $rxn_folder/mapping.txt | sort > $rxn_folder/mapping.sorted.txt
cat $rxn_folder/mapping.sorted.txt | sed "s/^/$rxn /"  >> all_mapping.txt
done

sort all_mapping.txt > all_mapping.sorted.txt
grep ':N#' all_mapping.sorted.txt > all_mapping.N.sorted.txt
sed 's/ .*//' all_mapping.N.sorted.txt | sort | uniq -c | sort -n > all_rxn_N_count.txt
sed 's/ *//; s/ .*//' all_rxn_N_count.txt | sort -n | uniq -c > all_rxn_N_count.histo

sed 's/.* //; s/=/\n/' all_mapping.N.sorted.txt | sort -u > all_atoms.N.sorted.txt
sed 's/.* //; s/=/\n/' all_mapping.N.sorted.txt | sort | uniq -c | sort -n > all_atoms.N.count.txt
sed 's/ *//; s/ .*//' all_atoms.N.count.txt | sort -n | uniq -c > all_atoms.N.count.histo
