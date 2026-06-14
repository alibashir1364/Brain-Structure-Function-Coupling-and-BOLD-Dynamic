%% ============================================================
%  GSP-based decomposition of BOLD signals
%  Coupled (low-pass) vs Decoupled (high-pass)
%  Following:
%   - van De Ville et al., 2019 (Graph Signal Processing)
%   - Fulcher et al., 2019 (Intrinsic timescales)
% ============================================================

clear; clc;

%% =======================
%  Basic parameters
% =======================
TR = 0.72;                 % Repetition time (seconds)
fs = 1/TR;                 % Sampling frequency
Nroi = 360;                % Glasser atlas
low_freq_thr = 0.14;       % RLFP threshold (Fulcher)
graph_cut = 0.35;          % Fraction of low graph frequencies (van de Ville)

%% =======================
%  Load data
% =======================
% BOLD: [ROI x Time x Subjects]
load('ROI_ts_Glasser_360_HCP.mat');
BOLD = roiTS_all;
[Nroi, T, Nsub] = size(BOLD);

% Structural connectivity
load('HCPMMP1_acpc_connectome_data.mat');
SC = mean(SC_all,3);           % Group-average SC
SC = SC ./ max(SC(:));             % Normalize weights

% ROI volume (control variable)
load('Glasser360_ROI_volume.mat'); % ROI_volume [360 x 1]


%% =======================
%  Structural graph Laplacian
% =======================
% Degree matrix
D = diag(sum(SC,2));

% Combinatorial Laplacian
L = D - SC;

% Eigendecomposition (Graph Fourier basis)
[U, Lambda] = eig(L);
lambda = diag(Lambda);

% Sort eigenvalues/eigenvectors (ascending graph frequency)
[lambda, idx] = sort(lambda,'ascend');
U = U(:,idx);

%% =======================
%  Define graph spectral filters
% =======================
% Number of low-frequency graph modes
K = round(graph_cut * Nroi);

% Low-pass graph filter (coupled component)
H_low = U(:,1:K) * U(:,1:K)';

% High-pass graph filter (decoupled component)
H_high = U(:,K+1:end) * U(:,K+1:end)';

%% =======================
%  Preallocate results
% =======================
RLFP_low  = zeros(Nroi,Nsub);
RLFP_high = zeros(Nroi,Nsub);

tau_low   = zeros(Nroi,Nsub);
tau_high  = zeros(Nroi,Nsub);

SC_strength = sum(SC,2);   % Structural strength (weighted degree)

%% =======================
%  Main loop over subjects
% ========0===============
for s = 1:Nsub
    
    % -----------------------------------------
    % Extract and preprocess BOLD
    % -----------------------------------------
    X = squeeze(BOLD(:,:,s));      % [ROI x Time]
    X = detrend(X')';              % Remove linear trend
    
    % -----------------------------------------
    % GSP filtering (van de Ville, 2019)
    % -----------------------------------------
    % Coupled signal: smooth on SC
    X_low  = H_low  * X;
    
    % Decoupled signal: high graph variation
    X_high = H_high * X;
    
    % -----------------------------------------
    % Loop over ROIs
    % -----------------------------------------
    for r = 1:Nroi
        
        %% ===== Low-pass (coupled) signal =====
        x = X_low(r,:);
        
        % Power spectrum
        Xf = fft(x);
        P = abs(Xf).^2 / T;
        f = (0:T-1)*fs/T;
        
        % RLFP
        RLFP_low(r,s) = sum(P(f<low_freq_thr)) / sum(P);
        
        % Autocorrelation timescale
        ac = autocorr(x,'NumLags',50);
        tau_low(r,s) = sum(ac);
        
        
        %% ===== High-pass (decoupled) signal =====
        x = X_high(r,:);
        
        % Power spectrum
        Xf = fft(x);
        P = abs(Xf).^2 / T;
        
        % RLFP
        RLFP_high(r,s) = sum(P(f<low_freq_thr)) / sum(P);
        
        % Autocorrelation timescale
        ac = autocorr(x,'NumLags',50);
        tau_high(r,s) = sum(ac);
        
    end
end

%% =======================
%  Group averaging
% =======================
RLFP_low_mean  = mean(RLFP_low,2);
RLFP_high_mean = mean(RLFP_high,2);

tau_low_mean   = mean(tau_low,2);
tau_high_mean  = mean(tau_high,2);

%% =======================
%  Statistical analysis (Fulcher-style)
% =======================
% Partial Spearman correlations controlling for ROI volume

% Coupled component
[r_low,p_low] = partialcorr( ...
    SC_strength, RLFP_low_mean, roi_volume, ...
    'Type','Spearman');

% Decoupled component
[r_high,p_high] = partialcorr( ...
    SC_strength, RLFP_high_mean, roi_volume, ...
    'Type','Spearman');

%% =======================
%  Visualization
% =======================
figure;

subplot(1,2,1)
scatter(SC_strength, RLFP_low_mean,40,'filled')
xlabel('SC strength')
ylabel('RLFP (coupled)')
title(sprintf('Coupled: r = %.2f',r_low))
grid on

subplot(1,2,2)
scatter(SC_strength, RLFP_high_mean,40,'filled')
xlabel('SC strength')
ylabel('RLFP (decoupled)')
title(sprintf('Decoupled: r = %.2f',r_high))
grid on
%% ============================================================
%  Inter-individual analysis (subject-level correlations)
%  Following Fulcher et al.
% ============================================================

% Preallocate correlation results
r_low_ind  = zeros(Nsub,1);   % Coupled (low-pass)
p_low_ind  = zeros(Nsub,1);

r_high_ind = zeros(Nsub,1);   % Decoupled (high-pass)
p_high_ind = zeros(Nsub,1);

for s = 1:Nsub
    
    % -----------------------------------------
    % Coupled (low-pass) signal
    % -----------------------------------------
    % RLFP vector for subject s (360 x 1)
    rlfp_s_low = RLFP_low(:,s);
    
    % Partial Spearman correlation
    [r_low_ind(s), p_low_ind(s)] = partialcorr( ...
        SC_strength, rlfp_s_low, roi_volume, ...
        'Type','Spearman');
    
    
    % -----------------------------------------
    % Decoupled (high-pass) signal
    % -----------------------------------------
    rlfp_s_high = RLFP_high(:,s);
    
    [r_high_ind(s), p_high_ind(s)] = partialcorr( ...
        SC_strength, rlfp_s_high, roi_volume, ...
        'Type','Spearman');
end
% Summary statistics
mean_r_low  = mean(r_low_ind);
std_r_low   = std(r_low_ind);

mean_r_high = mean(r_high_ind);
std_r_high  = std(r_high_ind);

fprintf('Coupled (low-pass):  mean r = %.3f ± %.3f\n', ...
        mean_r_low, std_r_low);

fprintf('Decoupled (high-pass): mean r = %.3f ± %.3f\n', ...
        mean_r_high, std_r_high);
figure;

subplot(1,2,1)
histogram(r_low_ind,20)
xlabel('Partial Spearman r')
ylabel('Number of subjects')
title('Inter-individual correlation (Coupled)')
grid on

subplot(1,2,2)
histogram(r_high_ind,20)
xlabel('Partial Spearman r')
ylabel('Number of subjects')
title('Inter-individual correlation (Decoupled)')
grid on
%% ============================================================
%  Paired statistical test: Coupled vs Decoupled
% ============================================================

% Remove NaNs if any
valid_idx = ~isnan(r_low_ind) & ~isnan(r_high_ind);

r_low_valid  = r_low_ind(valid_idx);
r_high_valid = r_high_ind(valid_idx);

% -----------------------------------------
% Paired non-parametric test (Wilcoxon)
% -----------------------------------------
[p_paired,~,stats] = signrank(r_low_valid, r_high_valid);

fprintf('\nPaired test (Wilcoxon signed-rank):\n');
fprintf('p-value = %.4e\n', p_paired);

% -----------------------------------------
% Effect size (Cohen''s d for paired samples)
% -----------------------------------------
diff_r = r_low_valid - r_high_valid;
cohen_d = mean(diff_r) / std(diff_r);

fprintf('Effect size (Cohen''s d) = %.3f\n', cohen_d);
figure;
boxplot([r_low_valid r_high_valid], ...
        {'Coupled','Decoupled'});
ylabel('Partial Spearman r');
title(sprintf('Paired comparison (p = %.2e)', p_paired));
grid on;
%% ============================================================
%  Null model: Randomized graph spectrum
% ============================================================

Nnull = 100;   % Number of null realizations

r_low_null  = zeros(Nnull,1);
r_high_null = zeros(Nnull,1);

for n = 1:Nnull
    
    % -----------------------------------------
    % Randomize graph Fourier basis
    % -----------------------------------------
    % Random orthonormal matrix
    [Q,~] = qr(randn(Nroi));
    
    % Construct null Laplacian
    L_null = Q * diag(lambda) * Q';
    
    % Eigen-decomposition (for safety)
    [U_null,~] = eig(L_null);
    
    % Define null filters
    H_low_null  = U_null(:,1:K) * U_null(:,1:K)';
    H_high_null = U_null(:,K+1:end) * U_null(:,K+1:end)';
    
    % -----------------------------------------
    % Apply null filters (group-level)
    % -----------------------------------------
    RLFP_low_null  = zeros(Nroi,1);
    RLFP_high_null = zeros(Nroi,1);
    
    for s = 1:Nsub
        
        X = squeeze(BOLD(:,:,s));
        X = detrend(X')';
        
        X_low_null  = H_low_null  * X;
        X_high_null = H_high_null * X;
        
        for r = 1:Nroi
            
            % ----- Low-pass null -----
            x = X_low_null(r,:);
            Xf = fft(x);
            P  = abs(Xf).^2 / T;
            f  = (0:T-1)*fs/T;
            RLFP_low_null(r) = RLFP_low_null(r) + ...
                sum(P(f<low_freq_thr))/sum(P);
            
            % ----- High-pass null -----
            x = X_high_null(r,:);
            Xf = fft(x);
            P  = abs(Xf).^2 / T;
            RLFP_high_null(r) = RLFP_high_null(r) + ...
                sum(P(f<low_freq_thr))/sum(P);
        end
    end
    
    % Average across subjects
    RLFP_low_null  = RLFP_low_null  / Nsub;
    RLFP_high_null = RLFP_high_null / Nsub;
    
    % Correlation with SC strength
    [r_low_null(n), p_high_ind(n)] = partialcorr( ...
        SC_strength, RLFP_low_null, roi_volume, ...
        'Type','Spearman');
    
    [r_high_null(n), p_high_ind(n)] = partialcorr( ...
        SC_strength, RLFP_high_null, roi_volume, ...
        'Type','Spearman');
end
figure;
histogram(r_low_null,20,'FaceAlpha',0.5); hold on;
xline(r_low,'r','LineWidth',2);
xlabel('Partial Spearman r');
ylabel('Count');
title('Null distribution (Coupled)');
legend('Null','Observed');
grid on;
%% ============================================================
%  Network-specific analysis for Hypothesis 6
%  (Decoupled timescales in association cortex)
% ============================================================
idx_sensory1=sensory1;
idx_sensory2=sensory2;
idx_motor=motor;
idx_limbic=limbic;
idx_association1=association1;

network_names = {'Sensory1','Sensory2','Motor','Limbic','Association'};
network_idx   = {idx_sensory1, idx_sensory2, idx_motor, ...
                 idx_limbic, idx_association1};

Nnet = numel(network_names);

% Preallocate
mean_RLFP_low_net  = zeros(Nnet,1);
mean_RLFP_high_net = zeros(Nnet,1);

for n = 1:Nnet
    idx = network_idx{n};
    
    mean_RLFP_low_net(n)  = mean(RLFP_low_mean(idx));
    mean_RLFP_high_net(n) = mean(RLFP_high_mean(idx));
end
%% ------------------------------------------------------------
%  Statistical comparison: Association vs others (Decoupled)
% ------------------------------------------------------------

assoc_RLFP = RLFP_high_mean(idx_association1);

p_values = zeros(Nnet-1,1);
labels   = {};

cnt = 1;
for n = 1:Nnet-1
    other_RLFP = RLFP_high_mean(network_idx{n});
    
    % Non-parametric test (distribution-free)
    p_values(cnt) = ranksum(assoc_RLFP, other_RLFP);
    labels{cnt}   = sprintf('Association vs %s', network_names{n});
    cnt = cnt + 1;
end

% Display results
fprintf('\nHypothesis 6: Decoupled RLFP (Association cortex)\n');
for i = 1:numel(p_values)
    fprintf('%s: p = %.3e\n', labels{i}, p_values(i));
end
%% ------------------------------------------------------------
%  Effect size (Cohen's d)
% ------------------------------------------------------------

cohen_d = zeros(Nnet-1,1);

cnt = 1;
for n = 1:Nnet-1
    x = assoc_RLFP;
    y = RLFP_high_mean(network_idx{n});
    
    pooled_std = sqrt((var(x)+var(y))/2);
    cohen_d(cnt) = (mean(x)-mean(y)) / pooled_std;
    cnt = cnt + 1;
end
%% ------------------------------------------------------------
%  Figure: Network-wise RLFP (Coupled vs Decoupled)
% ------------------------------------------------------------

figure;
bar([mean_RLFP_low_net mean_RLFP_high_net]);
set(gca,'XTickLabel',network_names,'XTickLabelRotation',30);
ylabel('Mean RLFP');
legend({'Coupled (Low-pass)','Decoupled (High-pass)'});
title('Network-specific intrinsic timescales');
grid on;
%% ------------------------------------------------------------
%  Figure: Decoupled RLFP distribution per network
% ------------------------------------------------------------

data = [];
group = [];

for n = 1:Nnet
    data  = [data; RLFP_high_mean(network_idx{n})];
    group = [group; n*ones(numel(network_idx{n}),1)];
end

figure;
boxplot(data, group, 'Labels', network_names);
ylabel('Decoupled RLFP');
title('Structure-decoupled intrinsic timescales');
grid on;
%% ------------------------------------------------------------
%  SC–RLFP correlation within networks (Decoupled)
% ------------------------------------------------------------

fprintf('\nSC–Decoupled RLFP correlations within networks:\n');

for n = 1:Nnet
    idx = network_idx{n};
    
    [r,p] = partialcorr( ...
        SC_strength(idx), ...
        RLFP_high_mean(idx), ...
        roi_volume(idx), ...
        'Type','Spearman');
    
    fprintf('%s: r = %.2f\n', network_names{n}, r);
    fprintf('%s: p = %.2f\n', network_names{n}, p);
end


%% ============================================================
%  Hypothesis 6: Network-wise analysis (Glasser hierarchical networks)
%  Networks: heteromodal, idiotypic, paralimbic, unimodal
% =============================================================
idx_heteromodal=heteromodal;
idx_idiotypic=idiotypic;
idx_paralimbic=paralimbic;
idx_unimodal=unimodal;
network_names = {'Heteromodal','Idiotypic','Paralimbic','Unimodal'};
network_idx   = {idx_heteromodal, idx_idiotypic, idx_paralimbic, idx_unimodal};

nNet = numel(network_names);

% Preallocate results
mean_RLFP_high   = zeros(nNet,1);
r_SC_RLFP_high   = zeros(nNet,1);
p_SC_RLFP_high   = zeros(nNet,1);

%% ------------------------------------------------------------
%  Loop over networks
% ------------------------------------------------------------
for n = 1:nNet
    
    idx = network_idx{n};
    
    % Extract network-specific data
    SC_n   = SC_strength(idx);
    RLFP_n = RLFP_high_mean(idx);
    Vol_n  = roi_volume(idx);
    
    % -----------------------------
    % Mean timescale (RLFP_high)
    % -----------------------------
    mean_RLFP_high(n) = mean(RLFP_n);
    
    % -----------------------------
    % Partial correlation
    % SC vs RLFP_high | volume
    % -----------------------------
    [r,p] = partialcorr(SC_n, RLFP_n, Vol_n, 'Type','Spearman');
    
    r_SC_RLFP_high(n) = r;
    p_SC_RLFP_high(n) = p;
end

%% ============================================================
%  Figure 1: Mean RLFP_high across networks
% ============================================================
figure;
bar(mean_RLFP_high);
set(gca,'XTickLabel',network_names,'FontSize',12)
ylabel('Mean RLFP_{high} (decoupled)');
title('Structure-decoupled timescales across cortical hierarchies');
grid on;

%% ============================================================
%  Figure 2: SC–RLFP_high partial correlations
% ============================================================
figure;
bar(r_SC_RLFP_high);
set(gca,'XTickLabel',network_names,'FontSize',12)
ylabel('Partial Spearman r');
title('SC–timescale coupling (decoupled signal)');
yline(0,'--');
grid on;

%% ============================================================
%  Statistical comparison: Fisher Z between networks
% ============================================================
Z = atanh(r_SC_RLFP_high);   % Fisher Z-transform

fprintf('\nFisher Z comparison (decoupled SC–timescale coupling):\n');
for i = 1:nNet
    for j = i+1:nNet
        
        z_diff = (Z(i)-Z(j)) / sqrt(1/(numel(network_idx{i})-3) + ...
                                    1/(numel(network_idx{j})-3));
        p_diff = 2*(1-normcdf(abs(z_diff)));
        
        fprintf('%s vs %s: Z=%.2f, p=%.4f\n', ...
            network_names{i}, network_names{j}, z_diff, p_diff);
    end
end

%% ============================================================
%  Scatter plots for illustration (reviewer-facing)
% ============================================================
figure;
for n = 1:nNet
    subplot(2,2,n)
    idx = network_idx{n};
    
    scatter(SC_strength(idx), RLFP_high_mean(idx), 30, 'filled')
    lsline
    xlabel('SC strength')
    ylabel('RLFP_{high}')
    title(network_names{n})
    grid on
end
sgtitle('Decoupled dynamics: SC vs timescale (network-wise)');

%% ============================================================
%  Hypothesis 6: Heteromodal vs Non-heteromodal cortex
%  Structure-decoupled timescales (RLFP_high)
% ============================================================
idx_heteromodal=heteromodal;        % indices of heteromodal regions
idx_nonheteromodal = setdiff(1:360, idx_heteromodal);
% ------------------------------------------------------------
% Define networks
% ------------------------------------------------------------
network_names = {'Heteromodal','Non-heteromodal'};
idx_network   = {idx_heteromodal, setdiff(1:360, idx_heteromodal)};

nNet = 2;

% Preallocate
mean_RLFP_high = zeros(nNet,1);
r_SC_RLFP_high = zeros(nNet,1);
p_SC_RLFP_high = zeros(nNet,1);

%% ------------------------------------------------------------
% Network-wise analysis
% ------------------------------------------------------------
for n = 1:nNet
    
    idx = idx_network{n};
    
    % Extract data
    SC_n   = SC_strength(idx);
    RLFP_n = RLFP_high_mean(idx);
    Vol_n  = roi_volume(idx);
    
    % Mean decoupled timescale
    mean_RLFP_high(n) = mean(RLFP_n);
    
    % Partial Spearman correlation (control for volume)
    [r,p] = partialcorr(SC_n, RLFP_n, Vol_n, 'Type','Spearman');
    
    r_SC_RLFP_high(n) = r;
    p_SC_RLFP_high(n) = p;
end

%% ============================================================
%  Figure 1: Mean RLFP_high
% ============================================================
figure;
bar(mean_RLFP_high)
set(gca,'XTickLabel',network_names,'FontSize',12)
ylabel('Mean RLFP_{high}')
title('Decoupled timescales: heteromodal vs non-heteromodal')
grid on;

%% ============================================================
%  Figure 2: Partial correlation SC–RLFP_high
% ============================================================
figure;
bar(r_SC_RLFP_high)
set(gca,'XTickLabel',network_names,'FontSize',12)
ylabel('Partial Spearman r')
title('SC constraint on decoupled dynamics')
yline(0,'--')
grid on;

%% ============================================================
%  Statistical comparison: Fisher Z-test
% ============================================================
Z = atanh(r_SC_RLFP_high);

n1 = numel(idx_network{1});
n2 = numel(idx_network{2});

z_diff = (Z(1)-Z(2)) / sqrt(1/(n1-3) + 1/(n2-3));
p_diff = 2*(1-normcdf(abs(z_diff)));

fprintf('\nFisher Z-test (SC–RLFP_high):\n');
fprintf('Heteromodal vs Non-heteromodal: Z=%.2f, p=%.4f\n', z_diff, p_diff);

%% ============================================================
%  Scatter plots (for illustration)
% ============================================================
figure;

subplot(1,2,1)
scatter(SC_strength(idx_network{1}), RLFP_high_mean(idx_network{1}), 35, 'filled')
lsline
xlabel('SC strength')
ylabel('RLFP_{high}')
title('Heteromodal')
grid on

subplot(1,2,2)
scatter(SC_strength(idx_network{2}), RLFP_high_mean(idx_network{2}), 35, 'filled')
lsline
xlabel('SC strength')
ylabel('RLFP_{high}')
title('Non-heteromodal')
grid on

sgtitle('Structure–decoupled dynamics across cortical hierarchy');
