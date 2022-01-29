function [history_ratio, history_N14, history_N15] = simMFA(model, last_fluxes, ...
    history_length, atom_names, atom_map_rxns, atom_map_mapping, N15_feed_rxn_ids, ...
    N15_feed_atoms, met_pools, flux_part, export_rxn_ids)
% The plan: create the atom transition matrix for the N atoms
% rows: the N-atoms of all compounds
% columns: every atom transition, for one single N atom transferred

if exist('history_length', 'var') % Process arguments and set up problem
    if isempty(history_length)
        history_length = 1000;
    end
else
    history_length = 1000;
end

[S_N, S_N_rxns, S_N_rxns_i, atom_met_inx] = create_S_N(model, atom_names, ...
    atom_map_rxns, atom_map_mapping);

% simple workaround: very fast metabolites get a larger pool
met_throughput = abs(model.S)*abs(last_fluxes);
met_ratios = met_throughput ./ met_pools;
met_pools(met_ratios>(flux_part/2)) = met_throughput(met_ratios>(flux_part/2))/(flux_part/2);

N14_atoms = met_pools(atom_met_inx);
N15_atoms = zeros(size(N14_atoms));

history_N14 = zeros(size(N14_atoms,1),history_length+1);
history_N15 = zeros(size(N14_atoms,1),history_length+1);
history_N14(:,1) = N14_atoms;
history_N15(:,1) = N15_atoms;

N_import=[];
N_import_rxns=[];

minimal_flux = 1E-10;
%sanitize fluxes
last_fluxes(abs(last_fluxes)<minimal_flux)=0;
last_met_throughput = abs(model.S)*abs(last_fluxes);
last_atom_throughput = last_met_throughput(atom_met_inx)/flux_part;

%precalculating some values, which are reused every step
base_rxn_flux = last_fluxes/flux_part;
base_flux = last_fluxes(S_N_rxns_i)/flux_part;
flux_direction = sign(base_flux);
base_flux = abs(base_flux);
S_N_corrected = S_N' .* flux_direction;
S_N_source = S_N_corrected';
S_N_source(S_N_source>0)=0;
S_N_dest = S_N_corrected';
S_N_dest(S_N_dest<0)=0;
trans_source = zeros(size(base_flux));
trans_dest = zeros(size(base_flux));
for transition_i = 1 : size(S_N,2)
    source = find(S_N_source(:,transition_i)<0);
    dest = find(S_N_dest(:,transition_i)>0);
    if length(source) ~= 1 || length(dest) ~= 1
        source = length(atom_names)+1; % pseudo atom pool to dump
        dest = source;
    end
    trans_source(transition_i) = source;
    trans_dest(transition_i) = dest;
end

N15_atom_import = zeros(size(N14_atoms));
for feed_rxn_i = 1:length(N15_feed_rxn_ids)
    dest=N15_feed_atoms(feed_rxn_i);
    % relying on selection of feed rxns to have pos. flux
    N15_atom_import(dest) = N15_atom_import(dest) + base_rxn_flux(N15_feed_rxn_ids(feed_rxn_i));
    if base_rxn_flux(N15_feed_rxn_ids(feed_rxn_i)) < 0
        "negative feed! " + model.rxns(N15_feed_rxn_ids(feed_rxn_i))
    end
end

export_atom_flux = zeros(size(N14_atoms));
other_atom_flux = export_atom_flux;
for rxn_i = 1:length(export_rxn_ids)
    met_flux = model.S(:,export_rxn_ids(rxn_i))*base_rxn_flux(export_rxn_ids(rxn_i));
    atom_flux = met_flux(atom_met_inx);
    export_atom_flux(atom_flux < 0) = export_atom_flux(atom_flux < 0) + atom_flux(atom_flux < 0);
    other_atom_flux(atom_flux > 0) = other_atom_flux(atom_flux > 0) + atom_flux(atom_flux > 0);
end

% history_length steps with flux_part mini steps of 1/flux_part poolsize
for history_step = 1:history_length
    history_step
    for mini_step = 1:flux_part
        do_repeat = true;
        while do_repeat
            do_repeat = false;
            N15_ratios = [N15_atoms./history_N14(:,1); 0];

            % first calculate all outgoing atoms
            N15_atoms_out = S_N_source * (N15_ratios(trans_source) .* base_flux);
            N15_atom_export = export_atom_flux .* N15_ratios(1:end-1);
            N15_atoms_out = N15_atoms_out + N15_atom_export;
            N14_atoms_out = S_N_source * base_flux + export_atom_flux - N15_atoms_out;
            N15_atom_gain = S_N_dest * (N15_ratios(trans_source).*base_flux);
            N14_atom_gain = S_N_dest * base_flux - N15_atom_gain;

            N15_atom_gain = N15_atom_gain + N15_atom_import;

            %workaround: in AraCore, only ACP is involved, which gets no 15N
            N14_atom_gain = N14_atom_gain+other_atom_flux;

            N15_atom_change = N15_atom_gain + N15_atoms_out;
            N14_atom_change = N14_atom_gain + N14_atoms_out;

            if sum(last_atom_throughput < abs(N15_atom_gain) + abs(N14_atom_gain) + abs(N15_atoms_out) + abs(N14_atoms_out) -1e-5)
                find(last_atom_throughput < abs(N15_atom_gain) + abs(N14_atom_gain) + abs(N15_atoms_out) + abs(N14_atoms_out) -1e-5)
                do_repeat = true;
            elseif sum((N14_atoms + N14_atom_change) < 0) || sum((N15_atoms + N15_atom_change) < 0)
                do_repeat = true;
            elseif  sum(abs(N14_atoms+N14_atom_change+N15_atoms+ N15_atom_change-history_N14(:,1)) > 1E-5) > 0
                find(abs(N14_atoms+N15_atoms-history_N14(:,1)) > 1E-5)
                do_repeat = true;
            else
                N14_atoms = N14_atoms + N14_atom_change;
                N15_atoms = N15_atoms + N15_atom_change;
            end
        end
    end
    history_N14(:,history_step+1) = N14_atoms;
    history_N15(:,history_step+1) = N15_atoms;
end

history_ratio = zeros(size(history_N15));
history_ratio(history_N14 > 0) = history_N15(history_N14 > 0)  ./ (history_N14(history_N14 > 0) + history_N15(history_N14 > 0));


end