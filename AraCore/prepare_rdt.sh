#!/bin/sh
echo '#!/bin/sh' > generated/replace_cpd_with_smiles.sh; cat species_id_smiles.txt | sed "s/\([\]\)/\1\1/g; s/^/sed 's%/; s/\t/ %/; s/$/ %g' |\\\/; " | head -c -3 | sed 's/_DASH_/-/g' >> generated/replace_cpd_with_smiles.sh

cat generated/rxn_table_with_cmp_flat.txt | sed 's/\[.\]//g' > generated/rxn_table_without_cmp_flat.txt
cat generated/rxn_table_without_cmp_flat.txt | sh generated/replace_cpd_with_smiles.sh | \
sed 's/ + /./g; s/ <=> />>/; s/ -->/>>/; s/<==>/>>/; s/ //g' |sed 's/[.][.]*/./g; s/[.]>>/>>/; s/>>[.]/>>/; s/[.]$//' | grep '>>' > generated/rxn_table.smiles.txt

while IFS= read -r line
do
rxn_id=$( echo "$line" | cut -f 1)
smiles=$( echo "$line" | cut -f 2)
mkdir -p reaction_intermediates/$rxn_id
echo $smiles > reaction_intermediates/$rxn_id/rxn.smiles
grep "$rxn_id	" generated/rxn_table_with_cmp_flat.txt | cut -f2 > reaction_intermediates/$rxn_id/rxn.flat_species_with_cmp
sed 's/>.*$//' reaction_intermediates/$rxn_id/rxn.flat_species_with_cmp | tr '+' '\n' | sed 's/[<=> ]//g; s/--//' | sort -u > reaction_intermediates/$rxn_id/from_species_with_cmp
sed 's/^.*>//' reaction_intermediates/$rxn_id/rxn.flat_species_with_cmp | tr '+' '\n' | sed 's/[<=> ]//g; s/--//' | sort -u > reaction_intermediates/$rxn_id/to_species_with_cmp
cat reaction_intermediates/$rxn_id/from_species_with_cmp  reaction_intermediates/$rxn_id/to_species_with_cmp | sed 's/\[.*//; s/-/_DASH_/g; s/^/^/; s/$/	/' | sort -u > reaction_intermediates/$rxn_id/species_without_cmp.txt
grep -f reaction_intermediates/$rxn_id/species_without_cmp.txt species_id_without_cmp_inchikey.txt > reaction_intermediates/$rxn_id/species_id_inchikey.txt
done < generated/rxn_table.smiles.txt
#in the end, we need to have reaction_intermediates/$rxn_id/  with rxn.smiles, from_species_with_cmp and to_species_with_cmp
# and here: species_id_without_cmp_inchikey.txt
