function [S_N, S_N_rxns, S_N_rxns_i, atom_met_inx] = create_S_N(model, ...
    atom_names, atom_map_rxns, atom_map_mapping)
atom_name_prefix_length = 2;

atom_cpds = extractBefore(atom_names, ":");
atom_met_inx = zeros(size(atom_names));
cpd_atoms = {[]};

for met_i = 1:length(model.mets)
    met = string(model.mets(met_i));
    cpd_atoms{met_i} = [];
    if ismember(met, atom_cpds)
        atom_met_inx(atom_cpds == met) = met_i;
        cpd_atoms{met_i} = find(atom_cpds == met);
    end
end

S_N = zeros(length(atom_names), length(atom_map_mapping));
S_N_rxns = string(zeros(size(S_N,2),1));
S_N_rxns_i = zeros(size(S_N,2),1);

for map_i = 1:length(atom_map_rxns) %284:284 %
    rxn = atom_map_rxns(map_i);
    sourceAtom = extractAfter(extractBefore(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    destAtom = extractAfter(extractAfter(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    source_i = (atom_names == sourceAtom);
    dest_i = (atom_names == destAtom);
    S_N(source_i,map_i) = -1;
    S_N(dest_i,map_i) = 1;
    rxn_i = find(model.rxns == rxn);
    if model.S(atom_met_inx(source_i), rxn_i) > 0
        S_N(:,current_single_transition) = -1*S_N(:,map_i);
    end
    S_N_rxns(map_i) = rxn;
    S_N_rxns_i(map_i) = rxn_i;
end

end