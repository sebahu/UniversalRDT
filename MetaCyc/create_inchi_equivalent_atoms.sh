echo -n "" > all_species_id_inchi_equiv.txt
for fn in reaction_intermediates/*/*.inchi
do
	if [ -s $fn ]
	then
		species_fn=$(echo "$fn" | sed 's/[.]inchi$/.species_id/');
		species_id=$(cat $species_fn);
		equiv=$(grep 'AuxInfo.*/E:' $fn | sed 's/^.*\/E://; s/\/.*$//');
		echo "$species_id	$equiv" >> all_species_id_inchi_equiv.txt
	fi
done 
sort -u all_species_id_inchi_equiv.txt > species_id_inchi_equiv.txt
grep -v '^       ' species_id_inchi_equiv.txt | grep -v '        $' > species_id_inchi_equiv_filtered.txt
