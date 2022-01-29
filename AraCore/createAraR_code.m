if 0
if 0
    %addpath(fullfile('~/apps/gurobi903','linux64', 'matlab'));
    addpath('~/apps/cobratoolbox/');
    initCobraToolbox(false);
    %changeCobraSolver('gurobi', 'all');
end
if 1
    model = readSBML('~/Documents/Promotion/nfm/data/AraCore2/ArabidopsisCoreModel.xml', 1000);
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
end
end
% goal: r1     = list("R"=c("Sout", "Sin"), "C"=c(-1, 1),     "E"="v1",   "T"=c("A", "A"))

atom_name_prefix_length = 2;

atom_N_id_table = readtable('~/Documents/Promotion/nfm/data/AraCore2/all_atoms.N.sorted.txt', 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_N_id_table.Var1),atom_name_prefix_length);

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

AtomTransitionRDT_table = readtable("~/Documents/Promotion/nfm/data/AraCore2/all_mapping.N.sorted.txt","Delimiter"," ", 'ReadVariableNames', false);
atom_map_rxns = string(AtomTransitionRDT_table.Var1);
atom_map_mapping = string(AtomTransitionRDT_table.Var2);

S_N = zeros(length(atom_cpds), length(atom_map_mapping));
S_N_rxns = string(zeros(size(S_N,2),1));
S_N_rxns_i = zeros(size(S_N,2),1);

for map_i = 1:length(atom_map_rxns) %284:284 %
    rxn = atom_map_rxns(map_i);
    sourceAtom = extractAfter(extractBefore(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    destAtom = extractAfter(extractAfter(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    source_i = find(atom_names == sourceAtom);
    dest_i = find(atom_names == destAtom);
    S_N(source_i,map_i) = -1;
    S_N(dest_i,map_i) = 1;
    rxn_i = find(model.rxns == rxn);
    if model.S(atom_met_inx(source_i), rxn_i) > 0
        S_N(:,current_single_transition) = -1*S_N(:,map_i);
    end
    S_N_rxns(map_i) = rxn;
    S_N_rxns_i(map_i) = rxn_i;
end

output = string(1:length(model.rxns))';
for rxn_i = 1:length(model.rxns)
    base_flux = solution.v(rxn_i);
    pattern = 'r%d     = list("R"=c(%s), "C"=c(%s),     "E"="v%d",   "T"=c(%s))';
    mets_i = find(model.S(:,rxn_i)~=0);
    mets = string(model.mets(mets_i));
    stoichio = model.S(mets_i, rxn_i);
    if base_flux < 0
        stoichio = -stoichio;
    end
    rxn = model.rxns(rxn_i);
    current_mappings = atom_map_mapping(atom_map_rxns == rxn);
    if ~isempty(current_mappings)
        atom_mapping = string(1:length(mets_i));
        name_mapping = containers.Map();
        currentID = 'A';
        for current_mapping_i = 1:length(current_mappings)
            current_mapping = current_mappings(current_mapping_i);
            sourceAtom = extractAfter(extractBefore(current_mapping,"="),atom_name_prefix_length);
            if ~ismissing(sourceAtom)
                destAtom = extractAfter(extractAfter(current_mapping,"="),atom_name_prefix_length);
                nextID = currentID;
                if ~isKey(name_mapping, cellstr(sourceAtom))
                    name_mapping(sourceAtom) = sprintf("%c",currentID);
                    nextID = currentID + 1;
                end
                if ~isKey(name_mapping, cellstr(destAtom))
                    name_mapping(destAtom) = sprintf("%c",currentID);
                    nextID = currentID + 1;
                end
                currentID = nextID;
            end
        end
        for met_i = 1:length(mets_i)
            met_atom_names = cellstr(atom_names(atom_met_inx == mets_i(met_i)));
            atom_mapping(met_i) = sprintf("%s",string(values(name_mapping,intersect(keys(name_mapping), met_atom_names))));
        end
    else
        atom_mapping = "";
    end
    output(rxn_i) = sprintf(pattern, rxn_i, printQuotedStrArray(mets), printNumArray(stoichio), rxn_i, printQuotedStrArray(atom_mapping));
end

out_text = sprintf('%s,\n', output);
fileID = fopen('rxnlist.R','w');
fprintf(fileID,'rxn <- list(%s\n)',out_text(1:end-2));
fclose(fileID);

function str = printQuotedStrArray(vec)
    str = sprintf(',"%s" ',vec);
    str = str(2:end-1);
end

function str = printNumArray(vec)
    str = sprintf(',%d ',vec);
    str = str(2:end-1);
end

