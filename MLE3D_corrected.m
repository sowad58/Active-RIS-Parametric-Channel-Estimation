function [var_amp_g_est, var_phas_g_est, g_est, Azidx, Elidx, Didx] = MLE3D_corrected( ...
        y, L, B, Dh, a_range, lambda, M_V, M_H, d_V, d_H, ...
        varphi_range, theta_range, Pp)

% MLE3D   Maximum-likelihood estimate of the near-field RIS channel g
%   Implements the three closed-form expressions:


% ---------------- Pull basic grid sizes ----------------- %
% For NF, contains 2 anlges+distance.
thetaSRes  = size(a_range,3);
varphiSRes = size(a_range,2);
distSRes   = size(a_range,4);

sigma1_sq  = 4*10^-12;                                   % thermal-Noise part
sigma2_sq  = 4*10^-12;                                   % hardware-noise part

% ---------------- Pre-compute F⁻¹ once (no angle dependence) -------------- %
Sigma = sigma1_sq * (B * Dh) * (B * Dh)' + sigma2_sq * eye(L);
Finv  = inv(Sigma);
oneL  = ones(L,1);
d0    = real(oneL' * Finv * oneL);      % guaranteed real
t1    = oneL' * Finv * y;               % complex scalar

% Allocate cost matrix for grid search
psi_cost = zeros(varphiSRes, thetaSRes, distSRes);

for Didx = 1:distSRes
    for Elidx_ = 1:thetaSRes
        A   = a_range(:,:,Elidx_,Didx);   % M × varphiSRes
        S   = B * Dh * A;                  % L × varphiSRes  (s = B D_h a)
        T   = Finv * S;                    % L × varphiSRes  (t = F^{-1} s)

        % Numerator: | y^H F^{-1} B D_h a(ψ) |^2
        num = abs(y' * T).^2;              % 1 × varphiSRes

        % Denominator: a(ψ)^H D_h^H B^H F^{-1} B D_h a(ψ) = s^H t
        den = real(sum(conj(S) .* T, 1));  % 1 × varphiSRes (force real for stability)

        % Final cost (vectorized over azimuth/ψ)
        psi_cost(:,Elidx_,Didx) = (num ./ max(den, eps)).';  % varphiSRes × 1
    end
end




% ---------------- Pick best direction / distance ------------------------ %
[~,maxind]          = max(psi_cost,[],'all');
[Azidx,Elidx,Didx]  = ind2sub(size(psi_cost),maxind);
a                   = a_range(:,Azidx,Elidx,Didx);   % best steering vector
BDHa                = B * Dh * a;                    % L × 1

% ---------------- Scalars reused in β̂ and ω̂ ---------------------------- %
% --- per-image amplitude estimate β̂ ---
s = B * Dh * a;         % L×1 : s = B D_h a(ψ)
t = Finv * s;           % L×1 : t = F^{-1} s

% Numerator: | y^H F^{-1} B D_h a |^2
num_beta = abs(y' * t)^2;

% Denominator: (a^H D_h^H B^H F^{-1} B D_h a)^2 = (s^H t)^2
den_inner = real(s' * t);              % force real (numerical safety)
den_beta  = Pp * max(den_inner, 0)^2;  % nonnegative guard

% Amplitude estimate
var_amp_g_est = num_beta / max(den_beta, eps);


% ---------------- ω̂  (phase)  Eq. 19 ------------------------------------ %% single ψ
var_phas_g_est = -angle( y' * (Finv * (B*Dh*a)) );
    
% ---------------- Assemble ĝ ------------------------------------------- %
g_est = sqrt(var_amp_g_est) * exp(1j*var_phas_g_est) * a;


%% For calculation of the direct channel, d between UE-BS. (OPTIONAL)
% ---------------- d̂ = sqrt(alpha) * exp(j*vartheta) -------------------- %
% Eq. (phase):   v̂ = arg(1ᵀ(F⁻¹ y − √Pp F⁻¹ B D_h ĝ))
% Eq. (amplitude): sqrt(α̂) = |1ᵀ(F⁻¹ y − √Pp F⁻¹ B D_h ĝ)| / ( √Pp * (1ᵀ F⁻¹ 1) )

% r               = oneL' * (Finv*y - sqrt(Pp) * (Finv * (B*Dh*g_est)));
% vartheta_hat    = angle(r);                                  % ∈ (−π, π]
% sqrt_alpha_hat  = abs(r) / (sqrt(Pp) * d0);                  % d0 = oneL' * Finv * oneL (computed above)
% alpha_hat       = sqrt_alpha_hat^2;                          % optional, if you also want α̂
% d_est           = sqrt_alpha_hat * exp(1j * vartheta_hat);   % final scalar estimate


end
