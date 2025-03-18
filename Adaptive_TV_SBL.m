function [SBL_xhat, err] = Adaptive_TV_SBL(A, y, N, sig2e, Tau_p, x, sbl_start)

    % Initialization parameters
    k_max = 500;  % EM maximum number of iterations
    t_max = 20;   % ADMM maximum number of iterations
    M = size(y, 2);
    rho_0 = 0.01; % Regularizer for ADMM dual variable
    epsilon = 1e-3; % Stopping criterion for EM and ADMM loops

    % Initialize variables
    Gamma = zeros(N, N, k_max + 1);
    Gamma(:, :, 1) = diag(sbl_start);.1* eye(N);
    SBL_SRR = zeros(1);
    sbl_max_iter = 1;
    l_k = max(1, round(sqrt(M) - 1)); % Number of neighbors in Omega
    lambda = zeros(N, N); % ADMM dual variable
    c = zeros(N, N); % Splitting Variable
    gamma = real(diag(Gamma(:, :, 1)));

    % Iterative Process
    for k = 1:k_max

        % Compute Posterior Mean and Covariance (EM loop)
        F1 = Gamma(:, :, k) * A' / (A * Gamma(:, :, k) * A' + sig2e * eye(Tau_p));
        mu_x = F1 * y;
        Sigma_x = Gamma(:, :, k) - F1 * A * Gamma(:, :, k);

        % Compute S matrix
        for m=1:M
        St(:,:,m)=mu_x(:,m)*mu_x(:,m)';
    end
    S=sum(St,3)/M+Sigma_x;
        if k < sbl_max_iter
            gamma = real(diag(S));
        else
            % ADMM loop
            rho = max(1e-4, 0.5 * rho_0 * (1 + cos(k * pi / k_max))); % Cosine decay rule

            % Compute Beta coefficients
            beta = zeros(N, N);
            for i = 1:N
                for j = 1:i-1
                    beta(i, j) = exp(-(log(gamma(i)) - log(gamma(j)))^2);
                    beta(j, i) = beta(i, j);
                end
            end

            % ADMM Inner Loop
            for t = 1:t_max
                bdiag = gamma;

                for i = 1:N
                    Omega_bar = max(1, i - l_k):i-1;
                    Omega_tilde = i+1:min(i + l_k, N);

                    % Normalize beta weights
                   beta(i, [Omega_bar, Omega_tilde]) = beta(i, [Omega_bar, Omega_tilde]) / sum(beta(i, [Omega_bar, Omega_tilde]));

                    % Compute gamma updates
                    a_bar = sum(c(i, Omega_bar)' +  log(bdiag(Omega_bar)) + lambda(i, Omega_bar)' / rho);
                    a_tilde = sum(c(Omega_tilde, i) -log(bdiag(Omega_tilde)) + lambda(Omega_tilde, i) / rho);

                    beta_square = length(beta(i, [Omega_bar, Omega_tilde]));
                    gamma(i) = real(S(i, i) / (1 + rho * beta_square - rho * (a_bar - a_tilde)));

                    % Stability pruning
                    gamma(i) = max(gamma(i), 1e-9);
                    bdiag(i) = gamma(i);  %commnet that if you want pararlel update for gamma vector
                end

                % Update variables
                for i = 1:N
                    for j = 1:i-1
                        gmbar = real(log(gamma(i)) - log(gamma(j)) - lambda(i, j) / rho);
                        c(i, j) = real(sign(gmbar) * max(0, abs(gmbar) - beta(i, j) / (M * rho)));
                        c(j, i) = -c(i, j);
                        lambda(i, j) = real(lambda(i, j) + rho * (c(i, j) -  (log(gamma(i)) - log(gamma(j)))));
                        lambda(j, i) = -lambda(i, j);
                    end
                end

                % Convergence check
                GM(:,t)=gamma;
                if t > 2 && (norm(gamma - GM(:, t-1))^2 / norm(GM(:, t))^2) < 0.01
                    break;
                end
            end
        end

        % Update Gamma
        Gamma(:, :, k + 1) = diag(gamma);
        xhat_sbl(:, :, k) = mu_x;
        residual_norm(k) = norm(y - A * mu_x, 'fro')^2; % Check for best solution
        Err(k) = norm(mu_x - x, 'fro')^2 / norm(x, 'fro')^2;
        Ey=(A*Gamma(:,:,k)*A'+sig2e*eye(Tau_p));
        obj(k)=real(log(det(Ey))+trace(y*y'*inv(Ey)));
        % Stopping criterion
        if k > 20
            if norm(xhat_sbl(:, :, k) - xhat_sbl(:, :, k - 1), 'fro') < epsilon || isnan(norm(gamma))
                break;
            end
        end
    end

    % Select best solution based on residual error
  %  [~, I] = max(residual_norm);
    SBL_xhat = xhat_sbl(:, :, end);
    err = min(Err);
end
