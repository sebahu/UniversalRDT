if 1
if 1
    %addpath(fullfile('~/apps/gurobi903','linux64', 'matlab'));
    addpath('~/apps/cobratoolbox/');
    initCobraToolbox(false);
    %changeCobraSolver('gurobi', 'all');
end
if 1
    model = readSBML('~/Documents/Promotion/nfm/data/AraCore2/ArabidopsisCoreModel.xml', 1000);
end

amino_acids = ["Gln"; "Asp"; "Glu"; "Asn"; "Ser"; "Cys"; "Thr"; "Gly"; "Met"; "Pro"; 
"Ala"; "Arg"; "Lys"; "His"; "Ile"; "Leu"; "Phe"; "Trp"; "Tyr"; "Val"];
aa_atoms = ["Asp:N#1"; "Glu:N#1"; "Ser:N#1"; "Cys:N#1"; "Thr:N#1"; "Gly:N#1";
    "Met:N#1"; "Pro:N#1"; "Ala:N#1"; "Ile:N#1"; "Leu:N#1"; "Phe:N#1"; "Tyr:N#1"; "Val:N#1"
    "Gln:N#1"; "Gln:N#2"; "Lys:N#1"; "Lys:N#2"; "Asn:N#1"; "Asn:N#2"; "Trp:N#1"; "Trp:N#2"; 
    "His:N#1"; "His:N#2"; "His:N#3"; "Arg:N#1"; "Arg:N#2"; "Arg:N#3"; "Arg:N#4" ];


species_lit_conc_table = readtable('species_no_cmp_literature_concentration.txt', 'ReadVariableNames', false, 'Delimiter', '\t');
species_lit_id = string(species_lit_conc_table.Var1);
species_lit_conc = species_lit_conc_table.Var2;

atom_name_prefix_length = 2;

atom_N_id_table = readtable('~/Documents/Promotion/nfm/data/AraCore2/all_atoms.N.sorted.txt', 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_N_id_table.Var1),atom_name_prefix_length);

AtomTransitionRDT_table = readtable("all_mapping.N.sorted.txt","Delimiter"," ", 'ReadVariableNames', false);
atom_map_rxns = string(AtomTransitionRDT_table.Var1);
atom_map_mapping = string(AtomTransitionRDT_table.Var2);

[S_N, S_N_rxns, S_N_rxns_i, atom_met_inx] = create_S_N(model, atom_names, atom_map_rxns, atom_map_mapping);

%nitrate and ammonium seem to be the source of choice of N in this setting ...
nitrate_feed_rxn_id = 450;
nh3_feed_rxn_id = 451;
N15_feed_rxn_ids = [ nitrate_feed_rxn_id, nh3_feed_rxn_id];

nitrate_import_id = find(atom_names == 'NO3[c]:N#1');
nh3_import_id = find(atom_names == 'NH4[c]:N#1');
N15_feed_atoms = [ nitrate_import_id, nh3_import_id ];

model2 = model;
model_bio1 = 'Bio_opt';
model_bio2 = 'Bio_NLim';
model_bio1_rxn_index = find(string(model.rxns) == model_bio1);
model_bio2_rxn_index = find(string(model.rxns) == model_bio2);

[minFlux, maxFlux] = fluxVariability(model,100);
opt_bio1 = maxFlux(model_bio1_rxn_index);
model.lb(model_bio1_rxn_index)=opt_bio1;
variable_rxns = (minFlux < -0.01 & maxFlux > 0.01);
model.c=ones(size(model.rxns));
solution = optimizeCbModel(model, 'min', 'zero');

diA = zeros(length(atom_names));
undiA = zeros(length(atom_names));
for i=1:size(S_N,2)
if sum(abs(S_N(:,i))) == 2
if solution.v(S_N_rxns_i(i)) < 0
    diA(find(S_N(:,i)==1),find(S_N(:,i)==-1))=1;
elseif solution.v(S_N_rxns_i(i)) > 0
    diA(find(S_N(:,i)==-1),find(S_N(:,i)==1))=1;
end
if solution.v(S_N_rxns_i(i)) ~= 0
    undiA((S_N(:,i)==-1),(S_N(:,i)==1))=1;
    undiA((S_N(:,i)==1),(S_N(:,i)==-1))=1;
end
end
end
% this is depending on the flux distribution!!!
N_graph=graph(undiA,atom_names);
N_digraph=digraph(diA,atom_names);
bins = conncomp(N_graph);
N_rel_inx=find(bins == bins(nh3_import_id));
N_graph_rel=graph(undiA(N_rel_inx,N_rel_inx),atom_names(N_rel_inx));
N_digraph_rel=graph(diA(N_rel_inx,N_rel_inx),atom_names(N_rel_inx));
rel_met_inx = atom_met_inx(N_rel_inx);
rel_met_no_cmp = unique(mets_no_cmp(rel_met_inx));

undiAFull = zeros(length(atom_names));
for i=1:size(S_N,2)
    undiAFull((S_N(:,i)==-1),(S_N(:,i)==1))=1;
    undiAFull((S_N(:,i)==1),(S_N(:,i)==-1))=1;
end
% this is depending on the flux distribution!!!
N_graph_full=graph(undiAFull,atom_names);
bins_full = conncomp(N_graph_full);
N_irrel_inx=find(~(bins_full == bins_full(nh3_import_id)));
N_rel_full_inx=find(bins_full == bins_full(nh3_import_id));
N_irrel_met_inx = atom_met_inx(N_irrel_inx);
N_irrel_met = unique(model.mets(N_irrel_met_inx));
N_irrel_met_no_cmp = unique(mets_no_cmp(N_irrel_met_inx));


mets_diA = zeros(length(rel_met_no_cmp));
for i=1:length(rel_met_no_cmp)
    for j=1:length(rel_met_no_cmp)
        if i ~= j
            mets_diA(i,j) = sum(sum(diA(mets_no_cmp(atom_met_inx) == rel_met_no_cmp(i), ...
                            mets_no_cmp(atom_met_inx) == rel_met_no_cmp(j))>0));
        end
    end
end
mets_diA(mets_diA>0)=1;
N_met_digraph=digraph(mets_diA,rel_met_no_cmp);

compact_diA = mets_diA;
compact_mets = rel_met_no_cmp;
i = 1;
while i < length(compact_mets)
    if sum(compact_diA(:,i)) == 1 %only one way into this node
        if sum(compact_diA(i,:)) == 1 %only one way out
            next_met = find(compact_diA(i,:));
            if sum(compact_diA(:,next_met)) == 1 % the next met has only this met as input
                % join them
                compact_mets(i) = compact_mets(i)+","+compact_mets(next_met);
                compact_mets = [ compact_mets(1:next_met -1); compact_mets(next_met+1:end)];
                compact_diA(i,:) = compact_diA(next_met,:);
                compact_diA = [compact_diA(1:next_met-1,:);compact_diA(next_met+1:end,:)];
                compact_diA = [compact_diA(:,1:next_met-1),compact_diA(:,next_met+1:end)];
                continue;
            end
        end
    end
    i = i+1;
end

N_compact_digraph=digraph(compact_diA,compact_mets);
sources = [find(compact_mets == "NH4"), find(compact_mets == "NO3,NO2")];
sinks = find(sum(compact_diA,2)==0);

target_node = find(compact_mets == "NH4");
sinks2 = [];
for node = 1:length(compact_mets)
    if isempty(shortestpath(N_compact_digraph,node,target_node))
        sinks2(end+1)=node;
    end
end
non_sinks = setdiff(1:length(compact_mets), sinks2);
compact_diA(non_sinks,sinks2) = compact_diA(non_sinks,sinks2)*10;

N_compact_digraph=digraph(compact_diA,compact_mets);
H = plot(N_compact_digraph);
layout(H,'layered', 'Sources', sources, 'Sinks', sinks);
delete(H);

last_met_throughput = abs(model.S)*abs(solution.v);

mets_no_cmp = extractBefore(string(model.mets),"[");

simFluxesPerHour = solution.v*2; %scale it to 6mgdw/gdw*h growth and scale to umol/gdw*h fluxes
logsPerHour = 1;
stepsPerHour = 100000;
simFluxes = simFluxesPerHour/logsPerHour;
simDurationHours = 200;
%necessary simFluxes = simFluxes/12; %per 5min. instead of 1h

met_pools = 0.011*ones(length(model.mets),1); %default: 11nmol/gdw, for umol/gdw*h fluxes
% read from file, spread over compartments, oriented on met throughput per
% compartment
met_throughput = abs(model.S)*abs(simFluxes);
N_species_lit_conc_table = readtable('N_relevant_species_no_cmp_literature_concentration.txt', 'ReadVariableNames', false, 'Delimiter', '\t');
N_species_id = string(N_species_lit_conc_table.Var1);
N_species_conc = N_species_lit_conc_table.Var2 * 0.011; %lit values are nmol/gFW
for i = 1:length(N_species_id)
   total_throughput = sum(met_throughput(mets_no_cmp == N_species_id(i)));
   current_pools = find(mets_no_cmp == N_species_id(i));
   for j = 1:length(current_pools)
       if met_throughput(current_pools(j)) > 0
           met_pools(current_pools(j)) = (met_throughput(current_pools(j))/total_throughput)*N_species_conc(i);
       end
   end
end
met_pools(met_pools<0.011) = 0.011;

%200x100000
history_length = simDurationHours * logsPerHour;
[history_ratio, hist_N14, hist_N15] = simMFA(model, simFluxes, history_length, atom_names, atom_map_rxns, ...
    atom_map_mapping, N15_feed_rxn_ids, N15_feed_atoms, met_pools, stepsPerHour/logsPerHour, model_bio1_rxn_index);


model2.c(model_bio1_rxn_index)=0;
model2.c(model_bio2_rxn_index)=1;

[minFlux2, maxFlux2] = fluxVariability(model2,100);
opt_bio2 = maxFlux2(model_bio2_rxn_index);
model2.lb(model_bio2_rxn_index)=opt_bio2;
variable_rxns = (minFlux2 < -0.01 & maxFlux2 > 0.01);
model2.c=ones(size(model2.rxns));
solution2 = optimizeCbModel(model2, 'min', 'zero');

last_met_throughput2 = abs(model.S)*abs(solution2.v);

% TODO N-limited means a lower growth - how much?
simFluxesPerHour2 = solution2.v*2; %scale it to 6mgdw/gdw*h growth and scale to umol/gdw*h fluxes
simFluxes2 = simFluxesPerHour2/logsPerHour;

met_pools2 = 0.011*ones(length(model.mets),1); %default: 11nmol/gdw, for umol/gdw*h fluxes
% read from file, spread over compartments, oriented on met throughput per
% compartment
met_throughput2 = abs(model.S)*abs(simFluxes2);
for i = 1:length(N_species_id)
   total_throughput = sum(met_throughput2(mets_no_cmp == N_species_id(i)));
   current_pools = find(mets_no_cmp == N_species_id(i));
   for j = 1:length(current_pools)
       if met_throughput2(current_pools(j)) > 0
           met_pools2(current_pools(j)) = (met_throughput2(current_pools(j))/total_throughput)*N_species_conc(i);
       end
   end
end
met_pools2(met_pools2<0.011) = 0.011;

[history_ratio2, hist2_N14, hist2_N15] = simMFA(model, simFluxes2, history_length, atom_names, atom_map_rxns, ...
    atom_map_mapping, N15_feed_rxn_ids, N15_feed_atoms, met_pools2, stepsPerHour/logsPerHour, model_bio2_rxn_index);


example_atoms=[311,298,641,284];
x_values = (0:history_length);
indices = x_values+1;
figure;
plot(x_values,history_ratio(example_atoms(1),indices), '-r', ...
    x_values,history_ratio2(example_atoms(1),indices), '-.r', ...
    x_values,history_ratio(example_atoms(2),indices), '-b', ...
    x_values,history_ratio2(example_atoms(2),indices), '-.b', ...
    x_values,history_ratio(example_atoms(3),indices),  '-g',...
    x_values,history_ratio2(example_atoms(3),indices), '-.g', ...
    x_values,history_ratio(example_atoms(4),indices),  '-c',...
    x_values,history_ratio2(example_atoms(4),indices), '-.c');

title("Simulated enrichment");
ylabel('Ratio of N15 of some N-atom-pools');
xlabel("simulated time [h]");
legend_strings=[(atom_names(example_atoms)+" (opt. bio)"), (atom_names(example_atoms)+" (N lim.)")];
legend(reshape(legend_strings',[8,1]));

atoms_with_N15_zero_for_bio2 = [308, 312, 639, 640];

ratioratio=ones(size(history_ratio));
ratioratio(history_ratio>0.01)=(history_ratio2(history_ratio>0.01)./history_ratio(history_ratio>0.01));
[ratioratio_max, max_i] = max(abs(1-ratioratio),[],2);

ratrioratio_maxdivergencevalues = zeros(size(ratioratio,1),1);
for i=1:size(ratioratio,1)
    ratrioratio_maxdivergencevalues(i) = ratioratio(i,max_i(i));
end
figure0 = figure;
histogram(ratrioratio_maxdivergencevalues(ratioratio_max>0 & ratioratio_max <1),0.8+(0:100)*0.01)

top_rxns_table = readtable('top_rxn.txt', 'ReadVariableNames', false);
top_rxns = string(top_rxns_table.Var1);
top_met_table = readtable('top_met.txt', 'ReadVariableNames', false);
top_met = extractAfter(string(top_met_table.Var1), atom_name_prefix_length);

colors = [1,1,0; 0.8, 0.8, 0.2; 0.5, 0.5, 0.5; 0.2, 0.2, 0.8; 0, 0, 1];
%max_for_pools = max([max(history_ratio,[],2), max(history_ratio2,[],2)],[],2);
max_for_pools = max(max(history_ratio,[],2),[],2);
[~,pools_sorted_indices] = sort(max_for_pools);

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

top_pools_sorted_indices = pools_sorted_indices;%(end-223:end);
top_mets = mets_no_cmp(atom_met_inx(top_pools_sorted_indices));
mets_sorted = mets_no_cmp(atom_met_inx(pools_sorted_indices));
top_aa_idx = ismember(top_mets, amino_acids);
sorted_aa_idx = find(ismember(mets_sorted, amino_acids));


indices = [11,21,51,101,201];
figure;
hold on
for index=1:5
    scatter(1:length(top_pools_sorted_indices), history_ratio(top_pools_sorted_indices,indices(index)), 150, 'x', 'MarkerEdgeColor',colors(index,:));
    %scatter(find(top_aa_idx), history_ratio(top_pools_sorted_indices(top_aa_idx),indices(index)), 200, 'd', 'MarkerEdgeColor',colors(index,:));
    scatter(1:length(top_pools_sorted_indices), history_ratio2(top_pools_sorted_indices,indices(index)), 150, 'o', 'MarkerEdgeColor',colors(index,:));
    %scatter(find(top_aa_idx), history_ratio2(top_pools_sorted_indices(top_aa_idx),indices(index)), 200, 'filled','o', 'MarkerEdgeColor',colors(index,:));
end
hold off
xlim([0 length(top_pools_sorted_indices)]); %225])
ylim([-0.01 1.01])
set(gca,'FontSize',20)
set(gca,'xtick',[])
lgd = legend({'opt. biomass after 10h'; 'N limited biomass after 10h'; 'opt. biomass after 20h'; 'N limited biomass after 20h'; ...
    'opt. biomass after 40h'; 'N limited biomass after 40h'; 'opt. biomass after 100h'; 'N limited biomass after 100h'; ...
    'opt. biomass after 200h'; 'N limited biomass after 200h' });
%lgd = legend({'opt. biomass after 10h'; 'opt. biomass after 20h';  ...
%    'opt. biomass after 40h'; 'opt. biomass after 100h';  ...
%    'opt. biomass after 200h';  });
lgd.NumColumns = 1;
lgd.Orientation = 'horizontal';
lgd.Box=false;
lgd.Location = 'northwest';
ylabel('Ratio of 15N for individual atom pool');
xlabel('individual N-atom pools');

figure;
hold on
for index=1:5
    scatter(52:length(sorted_aa_idx), history_ratio(pools_sorted_indices(sorted_aa_idx(52:end)),indices(index)), 150, 'x', 'MarkerEdgeColor',colors(index,:));
    scatter(52:length(sorted_aa_idx), history_ratio2(pools_sorted_indices(sorted_aa_idx(52:end)),indices(index)), 150, 'o', 'MarkerEdgeColor',colors(index,:));
end
hold off
xlim([51 length(sorted_aa_idx)+1])
ylim([-0.01 1.01])
ylabel('Ratio of 15N for individual atom pool');
xlabel('individual N-atom pools that reach > 1% 15N enrichment during 200h');
xticklabels(atom_names(pools_sorted_indices(sorted_aa_idx(52:end))));
xticks(52:length(sorted_aa_idx));

atom_names_no_cmp = extractBefore(atom_names,"[") + extractAfter(atom_names,"]");
history_ratio_aa = zeros(length(aa_atoms), history_length+1);
history_ratio2_aa = history_ratio_aa;
for aa_atom_i = 1:length(aa_atoms)
    aa_atom = aa_atoms(aa_atom_i);
    current_atoms_idx = find(ismember(atom_names_no_cmp, aa_atom));
    current_totals = met_pools(atom_met_inx(current_atoms_idx));
    current_partitioning = current_totals / (sum(current_totals));
    current_totals2 = met_pools2(atom_met_inx(current_atoms_idx));
    current_partitioning2 = current_totals2 / (sum(current_totals2));
    history_ratio_aa(aa_atom_i, : ) = history_ratio(current_atoms_idx,:)'*current_partitioning;
    history_ratio2_aa(aa_atom_i, : ) = history_ratio2(current_atoms_idx,:)'*current_partitioning2;
end

[~, sorted_aa_index] = sort(history_ratio_aa(:,201));

figure;
hold on
for index=1:5
    scatter(1:length(aa_atoms), history_ratio2_aa(:,indices(index)) - history_ratio_aa(:,indices(index)), ...
        150, 'x', 'MarkerEdgeColor',colors(index,:));
end
hold off
xlim([0 length(aa_atoms)+1])
ylim([-0.0 0.2])
lgd = legend({'after 10h'; 'after 20h'; 'after 40h'; 'after 100h'; 'after 200h' });

ylabel('Diff of Ratio of 15N for atom pool N-limited vs optimal biomass');
xlabel('Non-compartmentalized N-atom pools of AAs');
xticklabels(aa_atoms);
xticks(1:length(aa_atoms));

figure;
hold on
for index=1:5
    scatter(1:length(aa_atoms), history_ratio_aa(sorted_aa_index,indices(index)), 150, 'x', 'MarkerEdgeColor',colors(index,:));
    scatter(1:length(aa_atoms), history_ratio2_aa(sorted_aa_index,indices(index)), 150, 'o', 'MarkerEdgeColor',colors(index,:));
end
hold off
xlim([0 length(aa_atoms)+1])
ylim([-0.01 1.01])
ylabel('Ratio of 15N for atom pool');
xlabel('Non-compartmentalized N-atom pools of AAs');
xticklabels(aa_atoms(sorted_aa_index));
xticks(1:length(aa_atoms));

%aa_atoms = ["Asp:N#1"; "Glu:N#1"; "Ser:N#1"; "Cys:N#1"; "Thr:N#1"; "Gly:N#1";
%    "Met:N#1"; "Pro:N#1"; "Ala:N#1"; "Ile:N#1"; "Leu:N#1"; "Phe:N#1"; "Tyr:N#1"; "Val:N#1"
%    "Gln:N#1"; "Gln:N#2"; "Lys:N#1"; "Lys:N#2"; "Asn:N#1"; "Asn:N#2"; "Trp:N#1"; "Trp:N#2"; 
%    "His:N#1"; "His:N#2"; "His:N#3"; "Arg:N#1"; "Arg:N#2"; "Arg:N#3"; "Arg:N#4" ];
plotspecs = [ "-r"; "-y"; "-g"; "-.r"; "-.y"; "-.g";
            "-b"; "-.b"; "-c"; "-.c"; "-m"; "-.m"; "-k"; "-.k";
            "--r"; "--r"; "--g"; "--g"; "--y"; "--y"; "--b"; "--b";
            "--m"; "--m"; "--m"; "--k"; "--k"; "--k"; "--k"];

plotspecs = [ "-r"; "-r"; "-g"; "-.r"; "-.y"; "-.g";
            "-b"; "-.b"; "-c"; "-r"; "-r"; "-r"; "-r"; "-r";
            "-r"; "-r"; "--g"; "--g"; "--y"; "--y"; "--b"; "--b";
            "--m"; "--m"; "--m"; "--k"; "--k"; "--k"; "--k"];

plotindices = [ 1; 2; 10; 11; 12; 13; 14; 15; 16; 
            3; 4; 5; 6; 7; 8; 9; 17; 18; 19; 20;
            21; 22; 23; 24; 25; 26; 27; 28; 29 ];
        


x_values = (0:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values,history_ratio_aa(plotindices(index),:), plotspecs(plotindices(index)));
end
hold off
ylabel('Ratio of 15N for atom pools');
xlabel('Time after source switch [hours]');
legend(aa_atoms(plotindices));

x_values = (0:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values,history_ratio2_aa(plotindices(index),:), plotspecs(plotindices(index)));
end
hold off
ylabel('Ratio of 15N for atom pools');
xlabel('Time after source switch [hours]');
legend(aa_atoms(plotindices));


x_values = (0:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values,history_ratio_aa(plotindices(index),:)-history_ratio_aa(26,:), plotspecs(plotindices(index)));
end
hold off
ylabel('Difference of 15N Ratio for atom pools compared to Arg:N#1');
xlabel('Time after source switch [hours]');
legend(aa_atoms(plotindices));

x_values = (0:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values,history_ratio2_aa(plotindices(index),:)-history_ratio2_aa(26,:), plotspecs(plotindices(index)));
end
hold off
ylabel('Difference of 15N Ratio for atom pools compared to Arg:N#1');
xlabel('Time after source switch [hours]');
legend(aa_atoms(plotindices));

x_values = (0:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values,history_ratio2_aa(plotindices(index),:)-history_ratio_aa(plotindices(index),:), plotspecs(plotindices(index)));
end
hold off
ylabel({'Difference of 15N Ratio for atom pools'; 'between N-limited and optimal biomass'});
xlabel('Time after source switch [hours]');
lgd = legend(aa_atoms(plotindices));
lgd.NumColumns = 2;

x_values = (1:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values,history_ratio_aa(plotindices(index),x_values+1)-history_ratio_aa(plotindices(index),x_values), plotspecs(plotindices(index)));
end
hold off
ylabel('Increase of Ratio of 15N for atom pools over the last hour');
xlabel('Time after source switch [hours]');
legend(aa_atoms(plotindices));

x_values = (2:200);
figure;
hold on
for index=1:length(aa_atoms)
plot(x_values, ...
    (history_ratio_aa(plotindices(index),x_values+1)-history_ratio_aa(plotindices(index),x_values)) - ...
    (history_ratio_aa(plotindices(index),x_values)-history_ratio_aa(plotindices(index),x_values-1)), plotspecs(plotindices(index)));
end
hold off
ylabel('Increase of Increase of Ratio of 15N for atom pools over the last hour');
xlabel('Time after source switch [hours]');
legend(aa_atoms(plotindices));

atom_names(history_ratio(:,201)==0 & (bins == bins(nh3_import_id))')

photosystem_rxn_names = ["RBO_h"; "PGP_h"; "GCAO_p"; "CAT_p"; "GGAT_p"; "GlyDH_m"; "AMT_m";
                    "DHLDH1_m"; "GlyHMT_m"; "SGAT_p"; "GCEADH_p"; "MalDH1_p"; "GCEAK_h"];
photosystem_rxns = find(ismember(model.rxns, photosystem_rxn_names));
photosystem_mappings = find(ismember(S_N_rxns, photosystem_rxn_names));
photosystem_cleave_cond_mappings = [190; 191; 478; 694; 816; 817; 826; 1380]; %DHLDH1 478, GlyHMT 826, SGAT 1380 maybe not
% GlyDH_m !   GGAT_p !
unique(S_N_rxns(photosystem_cleave_cond_mappings));

