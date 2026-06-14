function GSCI = compute_GSCI(BOLD, SC, eigVec, eigVal, params)
% ---------------------------------------------------------
% Graph-Spectral Causality Index (GSCI)
% %
% INPUTS:
% BOLD   : [N x T x S] BOLD signals (N regions, T time, S subjects)
% SC     : [N x N] structural connectivity matrix
% eigVec : [N x N] Laplacian eigenvectors of SC
% eigVal : [N x 1] Laplacian eigenvalues
% params :
%   .low_idx   : indices of low graph frequencies
%   .high_idx  : indices of high graph frequencies
%   .n_null    : number of null realizations
%   .regionVol : [N x 1] region volumes (covariate)
%
% OUTPUT:
% GSCI structure with fields:
%   .r_emp_low, .r_emp_high
%   .r_null_low, .r_null_high
%   .GSCI_low, .GSCI_high
% ---------------------------------------------------------

[N, T, S] = size(BOLD);

% Preallocate
r_emp_low  = zeros(S,1);
r_emp_high = zeros(S,1);

r_null_low  = zeros(S, params.n_null);
r_null_high = zeros(S, params.n_null);

% SC strength
SC_strength = sum(SC,2);

% Loop over subjects
for s = 1:S

    X = squeeze(BOLD(:,:,s));   % [N x T]

    % ---------- Graph Fourier Transform ----------
    X_hat = eigVec' * X;        % [N x T]

    % Low- and High-frequency reconstruction
    X_low  = eigVec(:,params.low_idx)  * X_hat(params.low_idx,:);
    X_high = eigVec(:,params.high_idx) * X_hat(params.high_idx,:);

    % ---------- RLFP computation ----------
    RLFP_low  = mean(abs(X_low),  2);
    RLFP_high = mean(abs(X_high), 2);

    % ---------- Empirical partial correlations ----------
    r_emp_low(s)  = partialcorr(SC_strength, RLFP_low,  params.regionVol, 'Type','Spearman');
    r_emp_high(s) = partialcorr(SC_strength, RLFP_high, params.regionVol, 'Type','Spearman');

    % ---------- Null model ----------
    for k = 1:params.n_null

        % Randomize spectrum (permute eigenvectors)
        perm_idx = randperm(N);
        U_null   = eigVec(:,perm_idx);

        Xh_null  = U_null' * X;

        Xl_null  = U_null(:,params.low_idx)  * Xh_null(params.low_idx,:);
        Xh_null2 = U_null(:,params.high_idx) * Xh_null(params.high_idx,:);

        RLFP_l_null = mean(abs(Xl_null), 2);
        RLFP_h_null = mean(abs(Xh_null2),2);

        r_null_low(s,k)  = partialcorr(SC_strength, RLFP_l_null, params.regionVol, 'Type','Spearman');
        r_null_high(s,k) = partialcorr(SC_strength, RLFP_h_null, params.regionVol, 'Type','Spearman');
    end
end

% ---------- Compute GSCI ----------
GSCI_low  = (r_emp_low  - mean(r_null_low, 2)) ./ abs(r_emp_low);
GSCI_high = (r_emp_high - mean(r_null_high,2)) ./ abs(r_emp_high);

% ---------- Output ----------
GSCI.r_emp_low   = r_emp_low;
GSCI.r_emp_high  = r_emp_high;
GSCI.r_null_low  = r_null_low;
GSCI.r_null_high = r_null_high;
GSCI.GSCI_low    = GSCI_low;
GSCI.GSCI_high   = GSCI_high;

end
