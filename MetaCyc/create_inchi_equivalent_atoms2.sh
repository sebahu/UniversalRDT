#sort -u all_species_id_inchi_equiv.txt > species_id_inchi_equiv.txt

echo -n "" > species_inchi_number_to_element_atom_number.txt
tail -n+2 species_id_formula.txt | while IFS= read -r line
 do
  # formula can be extreme, like: C10H14N5O6PS2.Mo.2O.S or C62H91N13O17P2.C10H12N5O3.Co
  #echo "$line" | xxd
  species_id=$(echo "$line" | cut -f1);  
  formula=$(echo "$line" | cut -f2)
  #echo $species_id
  #echo $formula
  if [ -n "$formula" ]
  then
   formula_parts=$(echo $formula | sed 's/\([A-Z]\)/\n\1/g' | sort)
   total_atom_count=1
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
     #echo $element $num_of_atoms

     for atom_counter in $(seq 1 $num_of_atoms)
     do
      echo ${species_id}"	${total_atom_count}	"${element}"#"${atom_counter} >> species_inchi_number_to_element_atom_number.txt
      total_atom_count=$(($total_atom_count + 1))
     done
    fi
   done
  fi
 done
