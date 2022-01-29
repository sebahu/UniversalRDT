#:/bin/bash
echo -ne "Number of reactions from MetaCyc, which could be mapped by RDT without obvious errors:\t";
cat test1c.txt | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc, which are mapped, but not identically by RDT:\t";
cat test3c.txt | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc containing C atoms, which are mapped by RDT:\t";
cat test1c.txt | grep '[:(]C#' | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc containing C atoms, which are mapped, but not identically by RDT:\t";
cat test3c.txt | grep '[:(]C#' | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc containing O atoms, which are mapped by RDT:\t";
cat test1c.txt | grep '[:(]O#' | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc containing O atoms, which are mapped, but not identically by RDT:\t";
cat test3c.txt | grep '[:(]O#' | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc containing N atoms, which are mapped by RDT:\t";
cat test1c.txt | grep '[:(]N#' | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from MetaCyc containing N atoms, which are mapped, but not identically by RDT:\t";
cat test3c.txt | grep '[:(]N#' | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc, which are mapped by RDT:\t";
cat test1c.txt | sed 's/ .*/ /' | grep -f rxns_ModelSEED_and_MetaCyc.grep |sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in ModelSEED, which are not identically mapped by RDT and in MetaCyc:\t";
cat test3c.txt | sed 's/ .*/ /' | grep -f rxns_ModelSEED_and_MetaCyc.grep |sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc containing C atoms, which are mapped by RDT:\t";
cat test1c.txt | grep '[:(]C#' | grep -f rxns_ModelSEED_and_MetaCyc.grep | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc containing C atoms, which are mapped, but not identically by RDT:\t";
cat test3c.txt | grep '[:(]C#' | grep -f rxns_ModelSEED_and_MetaCyc.grep | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc containing O atoms, which are mapped by RDT:\t";
cat test1c.txt | grep '[:(]O#' | grep -f rxns_ModelSEED_and_MetaCyc.grep | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc containing O atoms, which are mapped, but not identically by RDT:\t";
cat test3c.txt | grep '[:(]O#' | grep -f rxns_ModelSEED_and_MetaCyc.grep | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc containing N atoms, which are mapped by RDT:\t";
cat test1c.txt | grep '[:(]N#' | grep -f rxns_ModelSEED_and_MetaCyc.grep | sed 's/ .*/ /' | sort -u | wc -l
echo -ne "Number of reactions from ModelSEED in MetaCyc containing N atoms, which are mapped, but not identically by RDT:\t";
cat test3c.txt | grep '[:(]N#' | grep -f rxns_ModelSEED_and_MetaCyc.grep | sed 's/ .*/ /' | sort -u | wc -l
cat test3c.txt | grep '[:(]N#' | sed 's/ .*//' | sort | uniq -c | sed 's/^ *//; s/ .*//' | sort -n | uniq -c > N_diff.hist
cat test1c.txt | grep '[:(]N#' | sed 's/ .*//' | sort | uniq -c | sed 's/^ *//; s/ .*//' | sort -n | uniq -c > N_count.hist
cat test3c.txt | grep '[:(]N#' | sed 's/ .*//' | sort | uniq -c | grep ' 4 ' | sed 's/.* //' > rxns_4N_diff.txt
cat test3c.txt | grep '[:(]N#' | sed 's/ .*//' | sort | uniq -c | grep ' 14 ' | sed 's/.* //' > rxns_14N_diff.txt
