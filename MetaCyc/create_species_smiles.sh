#!/bin/sh

while IFS= read -r line
do
rxn_id="$( echo "$line" | cut -f 1)"
mkdir -p reaction_intermediates/$rxn_id
echo "$line" | cut -f 3 | sed 's% // %\n%g' > reaction_intermediates/$rxn_id/species_ids.txt
echo "$line" | cut -f 2 | sed 's% // %\n%g' > reaction_intermediates/$rxn_id/species_smiles.txt
done < rxn_species_smiles_ids.txt
