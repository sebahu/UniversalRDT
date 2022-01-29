function [N_graph, N_digraph, N_rel_inx] = graph_analysis(atom_names, S_N, S_N_rxns_i, fluxes, bidirectional, nh3_import_id)
diA = zeros(length(atom_names));
undiA = zeros(length(atom_names));
for i=1:size(S_N,2)
if sum(abs(S_N(:,i))) == 2
if fluxes(S_N_rxns_i(i)) < 0
    diA(find(S_N(:,i)==1),find(S_N(:,i)==-1))=1;
elseif fluxes(S_N_rxns_i(i)) > 0
    diA(find(S_N(:,i)==-1),find(S_N(:,i)==1))=1;
    if bidirectional(S_N_rxns_i(i))
        diA(find(S_N(:,i)==1),find(S_N(:,i)==-1))=1;
    end
end
if fluxes(S_N_rxns_i(i)) ~= 0
    undiA((S_N(:,i)==-1),(S_N(:,i)==1))=1;
    undiA((S_N(:,i)==1),(S_N(:,i)==-1))=1;
end
end
end
% this is depending on the flux distribution!!!
N_graph=graph(undiA,atom_names);
N_digraph=digraph(diA,atom_names);
figure
plot(N_graph)
figure
plot(N_digraph)

source_node = find(atom_names == "NO3[c]:N#1");
N_rel_inx = [];
for node = 1:length(atom_names)
    if ~isempty(source_node) && ~isempty(shortestpath(N_digraph,source_node,node))
        N_rel_inx(end+1)=node;
    end
end

