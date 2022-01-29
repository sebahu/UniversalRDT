#!/bin/bash
#species_id_smiles.txt
for rxn_id in $(cat rxns_with_one_mapping.txt)
do
rxn_folder=reaction_intermediates/$rxn_id
# first create rxn.smiles from mapping
#FROM-SIDE - (CPD-14766 0 17) (PHOSPHORYL-ETHANOLAMINE 18 25) 
#TO-SIDE - (PHTYOSPHINGOSINE-1-P 0 25)
fromside=$(grep 'FROM-SIDE' atom_mappings/$rxn_id | sed 's/^[^(]*(/(/; s/(\([^ ]* \)[^(]*/(\1/g; s/(//g' )
toside=$(grep 'TO-SIDE' atom_mappings/$rxn_id | sed 's/^[^(]*(/(/; s/(\([^ ]* \)[^(]*/(\1/g; s/(//g' )
echo $fromside' ' $toside' '
speciesn=$(echo $fromside' ' $toside' ' | sed 's/ /\n/g' | sort -u)
echo -n "" > ${rxn_folder}/species_inchikey.txt
for species in $speciesn
do
smiles=$(grep "^$species	" species_id_smiles.txt | cut -f2)
inchikey=$(obabel -:"$smiles" -oinchikey)
echo "${species}	${inchikey}" >> ${rxn_folder}/species_inchikey.txt
done
done
