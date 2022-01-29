#!/bin/sh
echo -n "" > all_metacyc_mapping.txt
echo -n "" > all_mapping.txt
for rxn_id in $(cat rxns_filtered1.txt)
do
rxn_folder=reaction_intermediates/$rxn_id
cat $rxn_folder/metacyc_mapping.sorted.txt | sed "s/^/$rxn_id /" >> all_metacyc_mapping.txt
cat $rxn_folder/mapping.sorted.txt | sed "s/^/$rxn_id /"  >> all_mapping.txt
done

sort all_metacyc_mapping.txt > all_metacyc_mapping.sorted.txt
sort all_mapping.txt > all_mapping.sorted.txt
