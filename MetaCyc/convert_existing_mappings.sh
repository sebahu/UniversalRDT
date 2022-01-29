#!/bin/bash
for rxn_id in $(cat rxns_filtered1.txt) # "1.1.1.176-RXN" "1.3.1.63-RXN" "1.1.3.32-RXN" #
do
 rxn_folder=reaction_intermediates/$rxn_id
 #FROM-SIDE - (CPD-14766 0 17) (PHOSPHORYL-ETHANOLAMINE 18 25) 
 #TO-SIDE - (PHTYOSPHINGOSINE-1-P 0 25)
 #INDICES - 1 0 2 6 9 4 3 5 7 8 13 14 15 11 10 12 16 18 17 
 # or:
 #FROM-SIDE - (Ox-NADPH-Hemoprotein-Reductases 0 18) (S-STYLOPINE 19 42) ((WATER 1) 43 43) ((WATER 2) 44 44)
 # containing a molecule without inchi, and a molecule in 2 instances
 # we don't need the identification of WATER 1 and WATER 2
 # remove the text before the first (
 # remove the second ( and the number before the closing ) in ((NAME n)
 # separate the parts by a newline instead of space
 fromside=$(grep 'FROM-SIDE' atom_mappings/$rxn_id | sed 's/^[^(]*(/(/; s/[(][(]\([^ ]*\) [^)]*[)]/(\1/g; s/[)] [(]/\n/g; s/^(//; s/)$//' )
 toside=$(grep 'TO-SIDE' atom_mappings/$rxn_id | sed 's/^[^(]*(/(/; s/[(][(]\([^ ]*\) [^)]*[)]/(\1/g; s/[)] [(]/\n/g; s/^(//; s/)$//'  )
 mapping_indices=( $(grep 'INDICES' atom_mappings/$rxn_id | tail -c+11) )

 echo -n "" > ${rxn_folder}/metacyc_mapping_lines.txt

 
 #echo -e 1 2 "\n" 3 4 5 "\n" 678 | while IFS= read -r line; do echo $line; done
 echo "$fromside" | while IFS= read -r line
 do
  # formula can be extreme, like: C10H14N5O6PS2.Mo.2O.S or C62H91N13O17P2.C10H12N5O3.Co - but we exclude these
  species_id=$(echo $line | sed 's/ .*//');
  mapping_index=$(echo $line | sed 's/[^ ]* //; s/ .*//');
  formula=$(grep -- "^$species_id	" species_id_formula.txt | cut -f2)
  if [ -n "$formula" ]
  then
   formula_parts=$(echo $formula | sed 's/\([A-Z]\)/\n\1/g' | sort)
   for element_and_count in $formula_parts
   do
    element=$(echo $element_and_count | sed 's/[0-9]*//g')
    if [ "$element" != "H" ]
    then
     num_of_atoms=$(echo $element_and_count | sed 's/[^0-9]*//g')
     if [ -z "$num_of_atoms" ]
     then
      num_of_atoms=1
     fi

     for atom_counter in $(seq 1 $num_of_atoms)
     do
      echo ${mapping_index}"	from	"${species_id}":"${element}"#"${atom_counter}"=" >> ${rxn_folder}/metacyc_mapping_lines.txt
      mapping_index=$(($mapping_index + 1))
     done
    fi
   done
  fi
 done

 echo "$toside" | while IFS= read -r line
 do
  # formula can be extreme, like: C10H14N5O6PS2.Mo.2O.S or C62H91N13O17P2.C10H12N5O3.Co
  species_id=$(echo $line | sed 's/ .*//');
  mapping_counter=$(echo $line | sed 's/[^ ]* //; s/ .*//');
  formula=$(grep -- "^$species_id	" species_id_formula.txt | cut -f2)
  if [ -n "$formula" ]
  then
   formula_parts=$(echo $formula | sed 's/\([A-Z]\)/\n\1/g' | sort)
   for element_and_count in $formula_parts
   do
    element=$(echo $element_and_count | sed 's/[0-9]*//g')
    if [ "$element" != "H" ]
    then
     num_of_atoms=$(echo $element_and_count | sed 's/[^0-9]*//g')
     if [ -z "$num_of_atoms" ]
     then
      num_of_atoms=1
     fi

     for atom_counter in $(seq 1 $num_of_atoms)
     do
      mapping_index=${mapping_indices[${mapping_counter}]}
      echo ${mapping_index}"	to	"${species_id}":"${element}"#"${atom_counter}"," >> ${rxn_folder}/metacyc_mapping_lines.txt
      mapping_counter=$(($mapping_counter + 1))
     done
    fi
   done
  fi
 done

 #now create the actual mapping
 sort -n  ${rxn_folder}/metacyc_mapping_lines.txt | cut -f3 | tr '\n' ' ' | sed 's/ //g; s/,$//' > ${rxn_folder}/metacyc_mapping.txt

 cat ${rxn_folder}/mapping.txt | sed 's/,/\n/g; s/COA-GROUPCO-A/CO-A/g' | sort > ${rxn_folder}/mapping.sorted.txt
 cat ${rxn_folder}/metacyc_mapping.txt | sed 's/,/\n/g; s/COA-GROUPCO-A/CO-A/g' | sort > ${rxn_folder}/metacyc_mapping.sorted.txt

done
