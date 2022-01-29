model = {};
model.S=[1,-1,0,0,0,0,0,0,-1,0,0,0; 
    0,1,-1,0,0,0,0,0,0,0,0,0;
    0,0,1,-1,0,0,0,0,0,0,0,0;
    0,0,0,1,-1,0,0,0,0,0,0,0;
    0,0,0,0,1,-1,0,0,0,0,0,0;
    0,0,0,0,0,1,-1,0,0,0,0,0;
    0,0,0,0,0,0,1,-1,0,0,0,1;
    0,0,0,0,0,0,0,0,1,-1,1,-1;
    0,0,0,0,0,0,0,0,0,1,-1,0];
model.rxns=["R_im"; "R_1";"R_2";"R_3"; "R_4"; "R_5";"R_6"; "R_ex"; "R_1b"; "R_cycle1"; "R_cycle2"; "R_1c" ];
model.mets=["M1"; "M2"; "M3"; "M4"; "M5"; "M6"; "M7"; "M2b"; "M2c" ];

atom_name_prefix_length = 2;

atom_names = ["M1:N#1"; "M2:N#1"; "M3:N#1"; "M4:N#1"; "M5:N#1"; "M6:N#1"; "M7:N#1" ; "M2b:N#1"; "M2c:N#1" ];

atom_map_rxns = ["R_1";"R_2";"R_3"; "R_4"; "R_5"; "R_6"; "R_1b"; "R_cycle1"; "R_cycle2"; "R_1c"];
atom_map_mapping = ["M_M1:N#1=M_M2:N#1"; "M_M2:N#1=M_M3:N#1"; "M_M3:N#1=M_M4:N#1"; 
    "M_M4:N#1=M_M5:N#1"; "M_M5:N#1=M_M6:N#1"; "M_M6:N#1=M_M7:N#1";
    "M_M1:N#1=M_M2b:N#1"; "M_M2b:N#1=M_M2c:N#1"; "M_M2c:N#1=M_M2b:N#1"; "M_M2b:N#1=M_M7:N#1"];

N15_feed_rxn_ids = 1;
N15_feed_atoms = 1;

%fluxes = [ 1; 1; 1; 1; 1; 1; 1; 1; 0; 0.5; 0.5; 0]; 
fluxes = [ 1.5; 1; 1; 1; 1; 1; 1; 1.5; 0.5; 0.5; 0.5; 0.5]; 

history_length = 2000;
history_points = [201, 401, 1001, 2001];
history_ratio = simMFA(model, fluxes, history_length, atom_names, ...
    atom_map_rxns, atom_map_mapping, N15_feed_rxn_ids, N15_feed_atoms);

x_values = (0:30);
indices = x_values*10+1;
figure;
plot(x_values,history_ratio(1,indices), '-r', ...
    x_values,history_ratio(2,indices), '-.r', ...
    x_values,history_ratio(3,indices), '-b', ...
    x_values,history_ratio(4,indices), '-.b', ...
    x_values,history_ratio(5,indices),  '-g',...
    x_values,history_ratio(6,indices), '-.g', ...
    x_values,history_ratio(7,indices),  '-c', ...
    x_values,history_ratio(8,indices),  '-.c', ...
    x_values,history_ratio(9,indices),  '-.y');

figure;
colors = [1,1,0; 0.8, 0.8, 0.2; 0.5, 0.5, 0.5; 0.2, 0.2, 0.8; 0, 0, 1];
linestyles = ["-xy"; "-xm"; "-xc"; "-xr"; "-xg"];
max_for_pools = max(history_ratio(:,1:110),[],2);
[~,pools_sorted_indices] = sort(max_for_pools);
indices = [10,20,40,70,110];
hold on
for index=1:5
    %scatter(1:9, history_ratio(pools_sorted_indices,indices(index)), 150, 'x', 'MarkerEdgeColor',colors(index,:));
    plot(1:9, history_ratio(pools_sorted_indices,indices(index)), linestyles(index));
end
hold off
xlim([0 10])
ylim([-0.01 1.01])
xticks(1:9);
xticklabels(model.mets(pools_sorted_indices));
