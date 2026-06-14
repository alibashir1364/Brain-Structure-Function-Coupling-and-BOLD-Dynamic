%% Load Glasser TSV file
filename = 'atlas-Glasser_dseg.tsv'; 
T = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');

%% Extract relevant columns

indices = T.index;                     
communities = T.community_economo;     % community_economo

%% Define the list of networks of interest
network_list = {'limbic', 'motor', 'sensory 1', 'sensory 2', 'association 1'};

%% Initialize cell array to hold indices for each network
network_indices = cell(length(network_list),1);

%% Loop through each network and collect indices
for i = 1:length(network_list)
    net_name = network_list{i};
    network_indices{i} = indices(strcmpi(communities, net_name));
end

%% Assign each network to a separate variable (optional)
limbic = network_indices{strcmpi(network_list,'limbic')};
motor = network_indices{strcmpi(network_list,'motor')};
sensory1 = network_indices{strcmpi(network_list,'sensory 1')};
sensory2 = network_indices{strcmpi(network_list,'sensory 2')};
association1 = network_indices{strcmpi(network_list,'association 1')};

%% Display sizes (optional)
disp('Number of regions per network:');
% disp(table(network_list', cellfun(@length, network_indices)', 'VariableNames', {'Network','NumRegions'}));

%% Now you have separate matrices: limbic, motor, sensory1, sensory2, association1
%% Load Glasser atlas TSV file
filename = 'atlas-Glasser_dseg.tsv';
T = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');

%% Extract columns
% فرض بر این است که ستون‌ها به این نام‌ها وجود دارند
glasser_idx = T.index;                    % ایندکس نواحی Glasser (1–360)
community_mesulam = T.community_mesulam;  % لیبل Mesulam

%% Define Mesulam communities of interest
mesulam_labels = {'unimodal', 'paralimbic', 'idiotypic', 'heteromodal'};

%% Initialize cell array
mesulam_indices = struct();

%% Collect indices for each Mesulam community
for i = 1:length(mesulam_labels)
    label = mesulam_labels{i};
    mesulam_indices.(label) = ...
        glasser_idx(strcmpi(community_mesulam, label));
end

%% Assign outputs to separate matrices
unimodal     = mesulam_indices.unimodal;
paralimbic  = mesulam_indices.paralimbic;
idiotypic   = mesulam_indices.idiotypic;
heteromodal = mesulam_indices.heteromodal;

%% Summary (optional)
fprintf('Mesulam community sizes:\n');
fprintf('Unimodal     : %d regions\n', length(unimodal));
fprintf('Paralimbic  : %d regions\n', length(paralimbic));
fprintf('Idiotypic   : %d regions\n', length(idiotypic));
fprintf('Heteromodal : %d regions\n', length(heteromodal));

%% Save results (optional)
save('Glasser_community_mesulam_indices.mat', ...
     'unimodal', 'paralimbic', 'idiotypic', 'heteromodal');
