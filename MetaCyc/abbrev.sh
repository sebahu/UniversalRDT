#!/bin/sh
sed 's/&alpha;/alpha/g; s/&beta;/beta/g; s/&chi;/chi/g; s/&delta;/delta/g; s/&Delta;/Delta/g; s/&epsilon;/epsilon/g;
s/&gamma;/gamma/g; s/&harr;/<-->/g; s/&iota;/iota/g; s/&kappa;/kappa/g; s/&lambda;/lambda/g; s/&mdash;/-/g;
s/&mu;/mu/g; s/&ndash;/-/g; s/&nu;/nu/g; s/&omega;/omega/g; s/&phi;/phi/g; s/&pi;/pi/g; s/&psi;/psi/g;
s/&rarr;/->/g; s/&tau;/tau/g; s/&xi;/xi/g; s/&zeta;/zeta/g; s/<[a-zA-Z]*>//g; s/<\/[a-zA-Z]*>//g;' All_reactions_of_MetaCyc.txt | sed "s/\\&prime;/'/g" > All_reactions_of_MetaCyc_plain.txt
