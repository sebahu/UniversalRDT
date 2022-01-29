#!/bin/bash
for rxn_id in $(cat rxns_filtered1.txt) # "1.1.1.152-RXN" #"2.3.1.92-RXN" #
do
rxn_folder=reaction_intermediates/$rxn_id
cd $rxn_folder
from_to=$(grep -m1 -B1 '$MOL' ECBLAST_smiles_AAM.rxn | head -n 1)
# handle the stuff before the first $MOL
# it contains in the last line the number of from-side and to-side molecules
from_num=$(echo "$from_to" | head -c3 | sed 's/ //g')

echo -n "" > mapping_lines.txt
counter=1

for fn2 in MOL_??
do

if [ -s ${fn2}.inchikey ]
then
#now we need the name of the molecule metabolite
#species_id=$(grep $(cat ${fn2}.inchikey) ../../species_id_without_cmp_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
species_id=$(grep $(cat ${fn2}.inchikey) species_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
if [ -z "$species_id" ]
then
species_id=$(grep $(head -c14 ${fn2}.inchikey) species_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
fi
#grep $(cat ${fn2}.inchikey) ../../species_id_without_cmp_inchikey.txt
#cat ${fn2}.inchikey
#echo $species_id_without_cmp
if [ $counter -le $from_num ]
then
# if we have the same inchikey from and to, then this will select the from species
species_id=$(echo $species_id | sed 's/ .*//')
mapping_side="from"
mapping_end="="
else
# if we have the same inchikey from and to, then this will select the to species
species_id=$(echo $species_id | sed 's/.* //')
mapping_side="to"
mapping_end=","
fi
echo $species_id > ${fn2}.species_id

atom_counter=0
last_element=
#from each line in InChi-Atom order we want the element and the mapping-index
#inchi-order is element-wise, so first all C, then all N, etc.
#so we use a format: <species>:<element>#<index in InChI for this atom per element - i.e. starts with 1 for every element>
if grep -q 'AuxInfo.*/N:' ${fn2}.inchi
then
	inchi_index=$(grep 'AuxInfo' ${fn2}.inchi | sed 's/^.*\/N://; s/\/.*$//; s/,/ /g')
else
	inchi_index=1
fi

#echo $inchi_index
for rdt_line in $inchi_index
do 
element_and_index=$(tail -n +$rdt_line ${fn2}.rdt_index | head -n 1)
element=$(echo "$element_and_index" | cut -f1)
mapping_index=$(echo "$element_and_index" | cut -f2)
#echo $element
#echo $mapping_index

if [ "$last_element" = "$element" ]
then
atom_counter=$(($atom_counter + 1))
else
last_element=$element
atom_counter=1
fi
echo ${mapping_index}"	${mapping_side}	"${species_id}":"${element}"#"${atom_counter}${mapping_end} >> mapping_lines.txt
done

fi #if [ -s ${fn2}.inchikey ]
counter=$(($counter + 1))
done

#now create the actual mapping, remove Hydrogen atoms (which likely will not be matched)
sort -n  mapping_lines.txt | grep -v ':H#' | cut -f3 | tr '\n' ' ' | sed 's/ //g; s/,$//' > mapping.txt
cd -
done
