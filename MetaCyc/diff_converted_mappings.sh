#!/bin/sh
for rxn_id in $(cat rxns_filtered1.txt)
do
echo -n ">>>>>>>>>>> from " $rxn_id " <<<<<"
rxn_folder=reaction_intermediates/$rxn_id
diff $rxn_folder/metacyc_mapping.sorted.txt $rxn_folder/mapping.sorted.txt
done
