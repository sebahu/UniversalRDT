#!/bin/bash
echo '#!/bin/sh' > replace_cpd_with_smiles.sh; cat species_id_smiles.txt | sed "s/\([\]\)/\1\1/g; s/^/sed 's,(/; s/\t/ ,:/; s/$/ ,g' |\\\/; " | head -c -3 | sed 's/_DASH_/-/g' >> replace_cpd_with_smiles.sh
echo -n  "" > rxn_table.flat.txt
for rxn_id in $(cat rxns_with_one_mapping.txt)
do
rxn_folder=reaction_intermediates/$rxn_id
# first create rxn.smiles from mapping
#FROM-SIDE - (CPD-14766 0 17) (PHOSPHORYL-ETHANOLAMINE 18 25) 
#TO-SIDE - (PHTYOSPHINGOSINE-1-P 0 25)
fromside=$(grep 'FROM-SIDE' atom_mappings/$rxn_id | sed 's/^[^(]*(/(/; s/(\([^ ]* \)[^(]*/(\1/g' )
toside=$(grep 'TO-SIDE' atom_mappings/$rxn_id | sed 's/^[^(]*(/(/; s/(\([^ ]* \)[^(]*/(\1/g' )
speciesn=$(echo $fromside' ' $toside' ' | sed 's/ /\n/g' | sort -u)
echo -n "" > ${rxn_folder}/species_inchikey.txt
for species in $speciesn
do
smiles=$(grep "^$species	" species_id_smiles.txt | cut -f2)
inchikey=$(obabel -:"$smiles" -oinchikey)
echo "${species}	${inchikey}" >> ${rxn_folder}/species_inchikey.txt
done
#remove (-duplications - however they were created
echo $rxn_id'	'$fromside '>>' $toside ' '  | sed 's/((/(/g ' >> rxn_table.flat.txt
done

#remove unmatched parts: they still start with (
cat rxn_table.flat.txt | sed 's/((/(/g' | sh replace_cpd_with_smiles.sh | sed 's/[\t][(][^ ]* /\t /g; s/ [(][^ ]* / /g; s/ [(][^ ]* / /g; s/ /./g; s/[.][.]/./g; s/[.]$//; s/[.]>>/>>/; s/[\t][.]/\t/; s/>>[.]/>>/; s/://g' > rxn_table.smiles.txt

while IFS= read -r line
do
rxn_id=$( echo "$line" | cut -f 1)
smiles=$( echo "$line" | cut -f 2)
mkdir -p reaction_intermediates/$rxn_id
echo $smiles > reaction_intermediates/$rxn_id/rxn.smiles
done < rxn_table.smiles.txt

echo -n "" > species_id_without_cmp_inchikey.txt
while IFS= read -r line
do
rxn_id=$( echo "$line" | cut -f 1)
smiles=$( echo "$line" | cut -f 2)
inchikey=$(obabel -:"$smiles" -oinchikey)
echo "${rxn_id}	${inchikey}" >> species_id_without_cmp_inchikey.txt
done < species_id_smiles.txt
