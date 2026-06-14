clear; clc;

% ================================
% General parameters
% ================================
TR = 0.72;                 % Repetition time (seconds)
fs = 1/TR;                 % Sampling frequency
Nroi = 360;                % Number of brain regions (Glasser)
% ================================
% Load preprocessed BOLD time series
% BOLD: [ROI x Time x Subjects]
% ================================
load('ROI_ts_Glasser_360_HCP.mat');  
% Variable name assumed: BOLD
BOLD = roiTS_all;
[Nroi, T, Nsub] = size(BOLD);
Nsub = 12
% ================================
% Load structural connectivity matrix
% ================================
load('SC.mat');
% Variable assumed: connectome [ROI x ROI x Subjects]

% Compute group-average structural connectivity
SC = mean(SC_all,3);

% Normalize SC for numerical stability
SC = SC ./ max(SC(:));
% ================================
% Preprocess BOLD signals
% ================================
for s = 1:Nsub
    
    % Extract subject-specific BOLD data
    X = squeeze(BOLD(:,:,s));     % [ROI x Time]
    
    % Remove linear trends from each ROI
    X = detrend(X')';             % Detrend along time
    
    % Band-pass filter (resting-state range)
    X = bandpass(X',[0.01 0.1],fs)';
    
    % Z-score each ROI time series
    % (mean=0, std=1 across time)
    Xz(:,:,s) = zscore(X,0,2);
end
% ================================
% Define all unique functional edges
% ================================
edge_mask = triu(ones(Nroi),1);   % Upper triangular mask
edge_idx  = find(edge_mask);      % Linear indices of edges
Ne = length(edge_idx);            % Number of edges
% ================================
% Compute edge time series
% ================================
for s = 1:Nsub
    
    % Extract z-scored BOLD for subject s
    Z = Xz(:,:,s);                % [ROI x Time]
    
    % Preallocate edge time series matrix
    % Rows: edges, Columns: time
    EdgeTS = zeros(Ne, T);
    
    % Loop over edges
    for e = 1:Ne
        
        % Convert linear index to ROI pair (i,j)
        [i,j] = ind2sub([Nroi Nroi], edge_idx(e));
        
        % Instantaneous edge co-fluctuation
        % e_ij(t) = z_i(t) * z_j(t)
        EdgeTS(e,:) = Z(i,:) .* Z(j,:);
    end
    
    % Store edge time series for subject s
    Edge_TS(:,:,s) = EdgeTS;      % [Edges x Time x Subjects]
end
% ================================
% Estimate intrinsic timescale of each edge
% ================================
maxLag = 50;   % Number of lags (~36s)

for s = 1:Nsub
    
    for e = 1:Ne
        
        % Extract edge time series
        ets = Edge_TS(e,:,s);
        
        % Compute autocorrelation
        ac = autocorr(ets, 'NumLags', maxLag);
        
        % Integral timescale (sum of autocorrelation)
        edge_tau(e,s) = sum(ac);
    end
end
% ================================
% Edge-wise variance across time
% ================================
edge_var = squeeze(var(Edge_TS,[],2));   % [Edges x Subjects]
% ================================
% Mean absolute co-fluctuation amplitude
% ================================
edge_energy = squeeze(mean(abs(Edge_TS),2));
% ================================
% Structural connectivity mapped to edges
% ================================
SC_edge = SC(edge_idx);

% Compute node strengths
node_strength = sum(SC,2);

% Map node-level strength to edges
for e = 1:Ne
    [i,j] = ind2sub([Nroi Nroi], edge_idx(e));
    SC_edge_mean(e) = (node_strength(i) + node_strength(j)) / 2;
end
% ================================
% Average edge metrics across subjects
% ================================
tau_mean = mean(edge_tau,2);
var_mean = mean(edge_var,2);

% Partial correlation (optional control variables)
[r_tau,p_tau] = corr(tau_mean, SC_edge_mean','Type','Spearman');
[r_var,p_var] = corr(var_mean, SC_edge_mean','Type','Spearman');
