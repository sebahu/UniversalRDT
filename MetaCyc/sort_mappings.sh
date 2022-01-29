#!/bin/bash
for rxn_id in $(cat rxns_filtered1.txt)
do
rxn_folder=reaction_intermediates/$rxn_id

cat ${rxn_folder}/mapping.txt | sed 's/,/\n/g' | sort > ${rxn_folder}/mapping.sorted.txt
cat ${rxn_folder}/metacyc_mapping.txt | sed 's/,/\n/g' | sort > ${rxn_folder}/metacyc_mapping.sorted.txt


done
