#!/bin/bash

echo -n "" > species_id_without_cmp_inchikey.txt
while IFS= read -r line
do
rxn_id=$( echo "$line" | cut -f 1)
smiles=$( echo "$line" | cut -f 2)
inchikey=$(obabel -:"$smiles" -oinchikey)
echo "${rxn_id}	${inchikey}" >> species_id_without_cmp_inchikey.txt
done < species_id_smiles.txt
