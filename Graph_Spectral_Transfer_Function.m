%% ============================================================
% Graph Spectral Transfer Function (Subject-wise Implementation)
% ============================================================

[Nnodes, T, Nsub] = size(BOLD);
Nbands = 30;

%% ------------------------
% 1. Structural graph Laplacian
% ------------------------
% W = SC;
% D = diag(sum(W,2));
% L = D - W;
% 
% [U, Lambda] = eig(L);
% lambda = diag(Lambda);
% 
% % sort frequencies
% [lambda, idx] = sort(lambda,'ascend');
% U = U(:,idx);

%% ------------------------
% 2. Structural strength
% ------------------------
% SC_strength = sum(W,2);

%% ------------------------
% 3. Define spectral bands
% ------------------------
edges = round(linspace(1, Nnodes+1, Nbands+1));
lambda_mid = zeros(Nbands,1);

%% ------------------------
% 4. Allocate outputs
% ------------------------
H_lambda_subject = zeros(Nbands, Nsub);   % subject-wise transfer
H_lambda_group   = zeros(Nbands,1);        % group mean

%% ------------------------
% 5. Loop over spectral bands
% ------------------------
for m = 1:Nbands
    
    band_idx = edges(m):(edges(m+1)-1);
    lambda_mid(m) = mean(lambda(band_idx));
    
    fprintf('Processing spectral band %d / %d\n', m, Nbands)
    
    for s = 1:Nsub
        
        % ------------------------
        % Extract subject signal
        % ------------------------
        X = BOLD_proc(:,:,s);   % [360 x T]
        
        % ------------------------
        % Graph Fourier transform
        % ------------------------
        X_hat = U' * X;
        
        % ------------------------
        % Band-limited reconstruction
        % ------------------------
        X_band_hat = zeros(size(X_hat));
        X_band_hat(band_idx,:) = X_hat(band_idx,:);
        
        X_band = U * X_band_hat;   % [360 x T]
        
        % ------------------------
        % Compute RLFP
        % ------------------------
        RLFP_band = compute_RLFP(X_band, TR);
        
        % ------------------------
        % Partial Spearman correlation
        % ------------------------
        H_lambda_subject(m,s) = partialcorr( ...
            SC_strength, ...
            RLFP_band, ...
            roi_volume, ...
            'Type','Spearman');
    end
    
    % ------------------------
    % Group-level average
    % ------------------------
    H_lambda_group(m) = mean(H_lambda_subject(m,:), 'omitnan');
end
%%
figure('Color','w','Position',[200 200 720 420])
plot(lambda_mid, H_lambda_group,'-o','LineWidth',2,'MarkerSize',6)
xlabel('Graph frequency (Laplacian eigenvalue)','FontSize',11)
ylabel('Partial Spearman correlation','FontSize',11)
title('Graph Spectral Transfer Function (Group-level)','FontSize',12)
grid on
box off
%%
mean_H = mean(H_lambda_subject, 2, 'omitnan');
std_H  = std(H_lambda_subject,  [], 2, 'omitnan');

figure('Color','w','Position',[200 200 720 420])
hold on

% Shaded area (mean ± std)
fill([lambda_mid; flipud(lambda_mid)], ...
     [mean_H + std_H; flipud(mean_H - std_H)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.4);

% Mean curve
plot(lambda_mid, mean_H, '-k', 'LineWidth',2);

xlabel('Graph frequency (Laplacian eigenvalue)','FontSize',11)
ylabel('Partial Spearman correlation','FontSize',11)
title('Inter-individual variability of spectral transfer','FontSize',12)
box off
grid on
