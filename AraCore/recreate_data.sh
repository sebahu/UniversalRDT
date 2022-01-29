#!/bin/sh
while IFS= read -r line
do
species_id=$( echo "$line" | cut -f 1)
smiles=$( echo "$line" | cut -f 2)
echo -n $species_id"	"
obabel -:"$smiles" -oinchikey
done < species_id_smiles.txt | sed 's/^.*	M_/M_/' > species_id_without_cmp_inchikey.txt
