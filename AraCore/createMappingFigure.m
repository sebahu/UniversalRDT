if 1
    %addpath(fullfile('~/apps/gurobi903','linux64', 'matlab'));
    addpath('~/apps/cobratoolbox/');
    initCobraToolbox(false);
    %changeCobraSolver('gurobi', 'all');
end
if 1
    model = readSBML('~/Documents/Promotion/nfm/data/AraCore2/ArabidopsisCoreModel.xml', 1000);
end
relevant_rxns = ["GlnS_c", "AsnS_c", "Asnase_c"];
relevant_rxn_is = find(ismember(model.rxns,relevant_rxns));
relevant_mappings = find(ismember(S_N_rxns,relevant_rxns));
relevant_met_is = find(sum(abs(model.S(:,relevant_rxn_is)),2));
relevant_atoms = find(sum(abs(S_N(:,relevant_mappings)),2));
relevant_S_N = S_N(relevant_atoms,relevant_mappings);

diA = zeros((1+length(relevant_rxn_is))*length(relevant_atoms));
node_names = atom_names(relevant_atoms);
rxn_sizes = zeros(size(relevant_rxn_is));
curr_rxn_node = length(relevant_atoms);
node_levels = -90*ones(size(relevant_atoms));
node_x = [ 2, 2.1, 2.2, 2.3, 2.4, 2, 2.1, 2.2, 2.3, 2.4, 5, 5.1, 3, 2, 2.1, 2.2, 2.3, 2.4, 6, 6.1, 7, 4 ];
node_labels = string(atom_cpds(relevant_atoms));
for i = 1:(length(relevant_atoms)-1)
    if node_labels(i) == node_labels(i+1)
        node_labels(i) = "";
    end
end
for i=1:length(relevant_rxns)
    currMappings = find(ismember(S_N_rxns(relevant_mappings), relevant_rxns(i)));
    substrates = find(sum(relevant_S_N(:,currMappings)<0,2));
    rxn_sizes(i) = length(currMappings);
    for j = 1:length(currMappings)
        curr_rxn_node = curr_rxn_node+1;
        substrate = find(sum(relevant_S_N(:,currMappings(j))<0,2));
        node_levels(substrate) = max(node_levels(substrate), -i);
        diA(substrate,curr_rxn_node)=1;
        product=find(sum(relevant_S_N(:,currMappings(j))>0,2));
        node_levels(product) = max(node_levels(product), -i-1);
        diA(curr_rxn_node, product)=1;
        node_names(curr_rxn_node) = string(relevant_rxns(i))+j;
        node_levels(curr_rxn_node) = -i-0.5;
        node_x(curr_rxn_node) = 4 + 0.1*j;
        if j < length(currMappings)
            node_labels(curr_rxn_node) = "";
        else
            node_labels(curr_rxn_node) = string(relevant_rxns(i));
        end
    end
end
diA = diA(1:curr_rxn_node, 1:curr_rxn_node);
relevant_N_digraph=digraph(diA,node_names);
H = plot(relevant_N_digraph);layout(H,'layered');
H.YData=node_levels;
H.XData = node_x;
H.NodeLabel = node_labels;

