echo '' > replace_equiv_atoms.sed
cat species_id_inchi_equiv_filtered.txt | while IFS= read -r line
 do
  # species_id  equiv groups
  species_id=$(echo "$line" | cut -f1);  
  equiv_groups=$(echo "$line" | cut -f2)
  for equiv_group in $(echo $equiv_groups | tr '()' '  ')
  do
   new_equiv_group=""
   for equiv_nr in $(echo $equiv_group | tr ',' ' ')
   do
    new_equiv_group="$new_equiv_group $(grep -- "^$species_id	$equiv_nr	" species_inchi_number_to_element_atom_number.txt | cut -f3)"
   done
   #we need a replacement of every single atom of the group with the whole group
   replacement="("$(echo $new_equiv_group | sed 's/ /;/g')")"
   for equiv_atom in $new_equiv_group
   do
    echo "s/${species_id}:${equiv_atom}"'\([=,]\)/'"${species_id}:${replacement}"'\\1/g' >> replace_equiv_atoms.sed
   done
  done
 done
