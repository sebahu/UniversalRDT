#!/bin/sh
sed 's/:N#.//g' all_mapping.N.sorted.txt | sort -u > all_mapping.N.species.txt
sed 's/=.*//' all_mapping.N.species.txt | uniq -c | sort -r -n | grep '  2 ' > rxn_species_N_split.txt
sed 's/ .*=/ /' all_mapping.N.species.txt | uniq -c | sort -r -n | grep '  2 ' > rxn_species_N_join.txt
cat rxn_species_N_split.txt rxn_species_N_join.txt | sed 's/.* 2 //; s/ .*//' | sort -u > rxns_split_join.txt
sed 's/ /./; s/=/:N#.=/; s/[][]/./g' all_mapping.N.species.txt > all_mapping.N.species.grep
for combo in $(cat all_mapping.N.species.grep)
do 
trouble=$(grep "$combo" all_mapping.N.sorted.txt | sed 's/=.*#/=/; s/.*#//; s/=/\n/' | sort | uniq -c | grep -v '     2' | head -n1)
if [ -n "$trouble" ]
then
echo $combo
fi
done > trouble.txt
grep '=' trouble.txt | grep -v 'M_ATP.*M_ADP' | grep -v 'M_ADP.*M_ATP' | sed 's/[.].*//' | sort -u > trouble_rxns.txt
cat rxns_split_join.txt trouble_rxns.txt | sort -u | grep -v AMPK_ > critical_rxns.txt
