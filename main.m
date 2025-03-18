clear all
close all
clc

%% Parameters setup
% System model: y_m = Ax_m + w_m, m = 1,...,M
M = 5;              % Number of snapshot vectors (columns)
N = 300;            % Signal length (rows)
sparsity_level = 20/N;  % Fraction of non-zero rows for 'random' sparsity
num_blocks = 3;     % Number of non-zero blocks for 'block' sparsity
sparsity_type = 'Hybrid';  % Options: 'Random', 'Block', or 'Hybrid'
block_size = 5;     % Size of each non-zero block
isolated_ratio = 0.25;  % Fraction of isolated non-zero rows for 'hybrid' sparsity
pilot_length = 30;  % Measurement dimension
threshold = 0.1;    % Threshold for support recovery

% SNR settings
snr_db = 0:4:16;    % SNR values in dB
SNR_length = length(snr_db);

% Number of Monte Carlo iterations
num_iterations = 300;

%% Initialize result matrices
% NMSE metrics
mse_adaptive_SBL = zeros(SNR_length, num_iterations);
mse_SBL = zeros(SNR_length, num_iterations);
mse_SBL_Dol = zeros(SNR_length, num_iterations);
mse_PC_SBL = zeros(SNR_length, num_iterations);
D_joint_MMSE = zeros(SNR_length, num_iterations);

% Support recovery metrics (F1 score)
adaptive_SBL_SRR = zeros(SNR_length, num_iterations);
SBL_SRR = zeros(SNR_length, num_iterations);
SBL_dol_SRR = zeros(SNR_length, num_iterations);
PC_SRR = zeros(SNR_length, num_iterations);

fprintf('Starting simulation with %d iterations...\n', num_iterations);
tic

%% Main simulation loop
for it = 1:num_iterations
    if mod(it, 10) == 0
        fprintf('Iteration %d of %d\n', it, num_iterations);
    end
    
    % Generate sparse signal
    X = generate_sparse_signal(N, M, sparsity_type, sparsity_level, num_blocks, block_size, isolated_ratio);
    nn_locations = find(vecnorm(X, 2, 2) > 0);
    Beta = ones(N, 1);
    
    % Generate sensing matrix
    L = pilot_length;
    Phi_B = double((randn(L*2, N) < 0));  % Binary Phi
    A = (reshape(lteSymbolModulate(Phi_B(:), 'QPSK'), L, N) / sqrt(L));
    
    % Calculate effective channel
    sig = A(:, nn_locations) * X(nn_locations, :);
    sig_Power = var(sig(:));
    
    % Loop over SNR values
    for e = 1:SNR_length
        % Add noise
        sig2e = sig_Power / (10 ^ (snr_db(e) / 10));
        Noise = sqrt(0.5 * sig2e) * (randn(L, M) + 1j * randn(L, M));
        Y = A * X + Noise;
        
        %% Oracle MMSE calculation
        R_diag = kron(diag(Beta(nn_locations)), eye(M));
        PHI = kron(A(:, nn_locations), eye(M));
        X_vec = X(nn_locations, :).';
        X_vec = X_vec(:);
        Noise_col = Noise.';
        B_vec = PHI * X_vec + Noise_col(:);
        
        X_joint_MMSE = R_diag * PHI' * ((PHI * R_diag * PHI' + sig2e * eye(M*L)) \ B_vec);
        D_joint_MMSE(e, it) = norm(X_vec - X_joint_MMSE)^2 / norm(X_vec)^2;
        
        %% M-SBL
        [SBL_xhat, ~, mse_SBLC, SBL_SRR(e, it), warm_init] = sbl_mmv(A, Y, N, X, nn_locations, sig2e, L);
        mse_SBL(e, it) = norm(SBL_xhat - X, 'fro')^2 / norm(X, 'fro')^2;
        SBL_SRR(e, it) = F1_score(SBL_xhat, nn_locations, threshold);
        
        %% pattern_coupled_MSBL 
        [X_pc, ~, mse_PC] = pattern_coupled_MSBL(A, Y, N, X, nn_locations, sig2e, L, warm_init);
        mse_PC_SBL(e, it) = norm(X_pc - X, 'fro')^2 / norm(X, 'fro')^2;
        PC_SRR(e, it) = F1_score(X_pc, nn_locations, threshold);
        
        %% SBL-DoL
        [mu, ~, ~, ~] = sparse_learning_dol_tv_even_odd_update_projection(A, Y, sig2e, 2000, 4, 0, 0, 1, 1);
        mse_SBL_Dol(e, it) = norm(mu - X, 'fro')^2 / norm(X, 'fro')^2;
        SBL_dol_SRR(e, it) = F1_score(mu, nn_locations, threshold);
        
        %% Proposed solution (Adaptive TV-SBL)
        [SBL_xhat_log, err] = Adaptive_TV_SBL(A, Y, N, sig2e, L, X, warm_init);
        mse_adaptive_SBL(e, it) = norm(SBL_xhat_log - X, 'fro')^2 / norm(X, 'fro')^2;
        adaptive_SBL_SRR(e, it) = F1_score(SBL_xhat_log, nn_locations, threshold);
    end
end

%% Process and display results
simulation_time = toc;
fprintf('Simulation completed in %.2f seconds\n', simulation_time);

% Calculate performance metrics with trimmed mean (removing extreme outliers)
percent_to_keep = 0.95;
num_keep = floor(percent_to_keep * num_iterations);

% Calculate NMSE in dB
NMSE_PC = 10 * log10(mean(mink(squeeze(mse_PC_SBL), num_keep, 2), 2))'
NMSE_proposed = 10 * log10(mean(mink(squeeze(mse_adaptive_SBL), num_keep, 2), 2))'
NMSE_MSBL = 10 * log10(mean(mink(squeeze(mse_SBL), num_keep, 2), 2))'
NMSE_DoL = 10 * log10(mean(mink(squeeze(mse_SBL_Dol), num_keep, 2), 2))'
NMSE_Oracle_MMSE = 10 * log10(mean(mink(squeeze(D_joint_MMSE), num_keep, 2), 2))'

% Calculate Support Recovery Rate (F1 Score)
SR_Proposed = mean(maxk(squeeze(adaptive_SBL_SRR), num_keep, 2), 2)'
SR_DOL = mean(maxk(squeeze(SBL_dol_SRR), num_keep, 2), 2)'
SR_SBL = mean(maxk(squeeze(SBL_SRR), num_keep, 2), 2)'
SR_PC = mean(maxk(squeeze(PC_SRR), num_keep, 2), 2)'

% Save results to MAT file --Comment if you want to save the results
%save_filename = sprintf('sparse_recovery_results_%s_sparsity.mat', lower(sparsity_type));
%save(save_filename, 'snr_db', 'NMSE_PC', 'NMSE_DOL', 'NMSE_MSBL', 'NMSE_aditia', ...
 %   'NMSE_Oracle_MMSE', 'SR_Proposed', 'SR_DOL', 'SR_SBL', 'SR_PC', 'sparsity_type', ...
  %  'sparsity_level', 'num_blocks', 'block_size', 'isolated_ratio', 'M', 'N', 'pilot_length');

%% Plotting results
% Create a function for consistent plot styling
plotResults(snr_db, sparsity_type, ...
    NMSE_proposed, NMSE_MSBL, NMSE_DoL, NMSE_PC, NMSE_Oracle_MMSE, ...
   SR_Proposed, SR_SBL, SR_DOL, SR_PC);

