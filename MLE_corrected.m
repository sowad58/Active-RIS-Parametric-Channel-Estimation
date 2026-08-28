function [var_amp_g_est, var_phas_g_est, g_est, Azidx, Elidx] = ...
         MLE_corrected(y, L, B, Dh, a_FarAppx_range, Pp)
% Implements closed-form MLE for (Az, El) grid using the equations in the images.

% Grid sizes
% For FF, we  only need the angles, not the distance
ElRes = size(a_FarAppx_range,3);
AzRes = size(a_FarAppx_range,2);

% Noise/interference covariance
sigma1_sq = 4*10^-12;
sigma2_sq = 4*10^-12;
Sigma = sigma1_sq * (B*Dh) * (B*Dh)' + sigma2_sq * eye(L);   % L×L
Finv  = inv(Sigma);

% Cost over the (Az,El) grid
psi_cost = zeros(AzRes, ElRes);

for Elidx = 1:ElRes
    A   = a_FarAppx_range(:,:,Elidx);     % M×AzRes  (each column is a(ψ))
    S   = B * Dh * A;                       % L×AzRes  (s = B D_h a)
    T   =    Finv * S;                    % L × varphiSRes  (t = F^{-1} s)

        % Numerator: | y^H F^{-1} B D_h a(ψ) |^2
        num = abs(y' * T).^2;              % 1 × varphiSRes

        % Denominator: a(ψ)^H D_h^H B^H F^{-1} B D_h a(ψ) = s^H t
        % den = real(sum(conj(S) .* T, 1));  % 1 × varphiSRes (force real for stability)
       den = sum(conj(S).*T, 1);

        % Final cost (vectorized over azimuth/ψ)
        psi_cost(:,Elidx) = (num ./ max(den, eps)).';  % varphiSRes × 1
end

% Argmax over the grid
[~,maxind]    = max(psi_cost(:));
[Azidx,Elidx] = ind2sub(size(psi_cost), maxind);

% Recompute amplitude/phase at the best direction
a = a_FarAppx_range(:,Azidx,Elidx);   % M×1
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
     % ∈(-π,π]

% ---------------- Assemble ĝ ------------------------------------------- %
g_est = sqrt(var_amp_g_est) * exp(1j*var_phas_g_est) * a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r               = oneL' * (Finv*y - sqrt(Pp) * (Finv * (B*Dh*g_est)));
% vartheta_hat    = angle(r);                                  % ∈ (−π, π]
% sqrt_alpha_hat  = abs(r) / (sqrt(Pp) * d0);                  % d0 = oneL' * Finv * oneL (computed above)̂
% d_est           = sqrt_alpha_hat * exp(1j * vartheta_hat);   % final scalar estimate

end
