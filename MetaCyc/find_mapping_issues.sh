grep -f elemental_mismatch.grep test1.txt test2.txt | sed 's/ .*/ /; s/^.*:/^/' | sort -u > rxns_with_element_mismatch.grep
grep ':.*:.*:' test1.txt test2.txt | sed 's/ .*/ /; s/^.*:/^/' | sort -u > rxns_with_3side_match.grep
grep '^[^:]*:[^:]*$' test1.txt test2.txt | sed 's/ .*/ /; s/^.*:/^/' | sort -u > rxns_with_1side_match.grep

cat rxns_with_element_mismatch.grep rxns_with_3side_match.grep rxns_with_1side_match.grep | sort -u > rxns_with_any_mismatch.grep

grep -v -f rxns_with_any_mismatch.grep test1.txt > test1c.txt
grep -v -f rxns_with_any_mismatch.grep test2.txt > test2c.txt
comm -3 -2 test1c.txt test2c.txt > test3c.txt
