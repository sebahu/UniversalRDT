#!/bin/sh
#for metaCycRxn_Id in $(cut -f1 All_reactions_EC_MetaCyc.txt)
for metaCycRxn_Id in $(tail -n 2700 rxns_filtered1.txt)   #cat
do
wget --no-check-certificate https://metacyc.org/META/download-atom-mappings?object=${metaCycRxn_Id} -O atom_mappings2/${metaCycRxn_Id}
done
