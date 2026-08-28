clc; clear; close all;

%% ----------------------- Environment -----------------------------------
freq   = 28e9;                              % Carrier
lambda = physconst('LightSpeed')/freq;      % Wavelength
NFConf     = true;                          % near-field case
porpose    = true;                          % run proposed estimator
Far_approx = true;                          % run far-field approx (Set to true to see comparison)
LSConf     = true;                          % run LS estimator (Set to true to see comparison)
Rate   = 15;                                % distance between UE-ARIS
Racian = false;                             % LOS only unless true

% UPA Setup
M_H = 32; M_V = 32; M = M_H*M_V;
d_H = 1/2; d_V = 1/2;                       % in wavelengths
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D     = sqrt(Hsize^2 + Vsize^2);            % array diagonal

d_fraun = 2*(D^2)/lambda;                   % Fraunhofer distance
d_NF    = d_fraun/10;                       % upper NF threshold
d_bjo   = 2*D;                              % Bjornson lower NF threshold

disp(['The near-field upper threshold is ' num2str(d_NF) ' (m)']);
disp(['Bjornson distance is ' num2str(d_bjo) ' (m)']);

if d_NF < d_bjo
    disp('The code will not work since Bjornson distance is greater than upper near field distance.');
    return;
end

% Element coordinates
i_idx = @(m) mod(m-1,M_H);          
j_idx = @(m) floor((m-1)/M_H);      
U = zeros(3,M);
for m = 1:M
    ym = (-(M_H-1)/2 + i_idx(m))*d_H*lambda;
    zm = (-(M_V-1)/2 + j_idx(m))*d_V*lambda;
    U(:,m) = [0; ym; zm];
end

%% ----------------------- Estimation params ------------------------------
varphiSRes = 128;
thetaSRes  = 128;
distSRes   = 64;
Plim = 20;                     
sigma1_sq = 4e-12;             
sigma2_sq = 4e-12;             
P_total   = 200;               
alpha  = 0.25;                 
P_ris  = alpha * P_total;      
P_user = (1-alpha)*P_total;    
Pd = P_user/11;                
Pp = 10*Pd;                    
SNR1 = Pd/sigma1_sq;           
SNR2 = P_ris/sigma2_sq;        

% True Channel h (BS <-> RIS)
varphi_BS = -pi/6;
theta_BS  = 0;
h  = UPA_Evaluation(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H,Rate);
h  = h(:);
Dh = diag(h);

nbrOfAngleRealizations = 50; % Lowered for faster local testing (adjust as needed)
nbrOfNoiseRealizations = 20;

%% ----------------------- Phase 1: Pre-compute h_hat_bank ----------------
% Define Sensing Array (Semi-Passive Elements as per Hu et al.)
N_s_H = 2; N_s_V = 2;
firstSensorCol = floor((M_H-N_s_H)/2) + 1;
firstSensorRow = floor((M_V-N_s_V)/2) + 1;
sensorCols = firstSensorCol:(firstSensorCol+N_s_H-1);
sensorRows = firstSensorRow:(firstSensorRow+N_s_V-1);
[sensorRowGrid,sensorColGrid] = ndgrid(sensorRows,sensorCols);
sensorIndexMatrix = (sensorRowGrid-1)*M_H + sensorColGrid;

% Calibration parameters
L_h = 16; 
P_h_cal = 1;
hCalSNR_dB = 20;
hCalPilots = ones(1,L_h);
h_hat_bank = complex(zeros(M, nbrOfNoiseRealizations));

fprintf('--- Running Phase 1: BS-ARIS 2D-ESPRIT Estimation ---\n');
for n2_cal = 1:nbrOfNoiseRealizations
    [Y_h, ~] = generateBSARISCalibrationData(h, sensorIndexMatrix, hCalPilots, P_h_cal, hCalSNR_dB);
    [h_hat_temp, ~, ~, ~, ~, ~] = estimateBSARISChannel2DESPRIT(...
        Y_h, hCalPilots, P_h_cal, sensorIndexMatrix, lambda, M_V, M_H, d_V, d_H, Rate);
    h_hat_bank(:, n2_cal) = h_hat_temp;
end
fprintf('Phase 1 Complete.\n\n');

phase1_NMSE_linear = mean(vecnorm(h_hat_bank - h).^2) / norm(h)^2;
phase1_NMSE_dB = 10 * log10(phase1_NMSE_linear);
fprintf('Average Phase 1 BS-RIS NMSE: %.2f dB\n', phase1_NMSE_dB);

%% ----------------------- Storage Arrays ---------------------------------
capacity          = zeros(1, nbrOfAngleRealizations);
rate_proposed     = NaN(Plim-1, nbrOfAngleRealizations, nbrOfNoiseRealizations);
Far_rate_proposed = NaN(Plim-1, nbrOfAngleRealizations, nbrOfNoiseRealizations);
rate_LS           = NaN(Plim-1, nbrOfAngleRealizations, nbrOfNoiseRealizations);

g_NMSE_proposed     = NaN(Plim-1, nbrOfAngleRealizations, nbrOfNoiseRealizations);
Far_g_NMSE_proposed = NaN(Plim-1, nbrOfAngleRealizations, nbrOfNoiseRealizations);
g_NMSE_LS           = NaN(Plim-1, nbrOfAngleRealizations, nbrOfNoiseRealizations);

%% ----------------------- Codebooks --------------------------------------
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,0);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);

try
    load("WideTwobeam32.mat");                     
    beamresponses = [beamresponses, firsttarget, secondtarget];
catch
    disp('Warning: WideTwobeam32.mat not found. Using standard codebook.');
end
Farbeamresponses = beamresponses;

%% ----------------------- Realization loops ------------------------------
fprintf('--- Running Phase 2: Online User Tracking ---\n');
for n1 = 1:nbrOfAngleRealizations
    fprintf('Processing Angle Realization %d / %d...\n', n1, nbrOfAngleRealizations);
    
    if NFConf
        d_t = unifrnd(d_bjo, d_NF, 1);
    else
        d_t = unifrnd(d_fraun, 5*d_fraun, 1);
    end
    azimuth   = unifrnd(-pi/3, pi/3, 1);
    elevation = unifrnd(-pi/3, pi/3, 1);
    
    varphi_range = linspace(azimuth  - pi/12, azimuth  + pi/12, varphiSRes);
    theta_range  = linspace(elevation- pi/24, elevation+ pi/24, thetaSRes);
    dist_range   = zeros(1,distSRes);
    
    if NFConf
        mind = max([d_bjo, d_t - d_bjo/8]);
        maxd = min([d_fraun, d_t + d_fraun/8]);
    else
        mind = max([d_NF, d_t - d_fraun/4]);
        maxd = min([5*d_fraun, d_t + d_fraun/4]);
    end
    dist_range(:) = linspace(mind, maxd, distSRes);
    
    if porpose || LSConf
        a_range = zeros(M, varphiSRes, thetaSRes, distSRes);
        for l = 1:distSRes
            d_red = repelem(dist_range(l), 1, varphiSRes);
            parfor it = 1:thetaSRes
                a_range(:,:,it,l) = nearFieldChan( ...
                    d_red, varphi_range, repelem(theta_range(it),1,varphiSRes), U, lambda);
            end
        end
    end
    
    if Far_approx
        a_FarAppx_range = zeros(M, varphiSRes, thetaSRes);
        parfor it = 1:thetaSRes
            a_FarAppx_range(:,:,it) = UPA_Evaluate( ...
                lambda, M_V, M_H, varphi_range, repelem(theta_range(it),1,varphiSRes), d_V, d_H);
        end
    end
    
    % True RIS-user channel g
    if ~Racian
        g = nearFieldChan_withAMplitude(d_t, azimuth, elevation, U, lambda);
    else
        g = sqrt(K/(K+1))*nearFieldChan_withAmplitude(d_t,azimuth,elevation,U,lambda) + ...
            sqrtm(R)* sqrt(1/(K+1)/2) * (randn(M,1) + 1i*randn(M,1));
    end
    
    % Capacity (upper bound uses TRUE channels)
    abs_h1_sq = abs(g).^2;
    abs_h2_sq = abs(h).^2;
    numerator   = abs_h2_sq .* abs_h1_sq;
    denominator = abs_h1_sq*SNR1 + abs_h2_sq*SNR2 + 1;
    result = sum(numerator ./ denominator) * SNR1 * SNR2;
    capacity(n1) = log2(1 + result);
    
    %% --------------------- Noise realizations (PARFOR) ------------------
    parfor n2 = 1:nbrOfNoiseRealizations
        rate_proposed_loc       = NaN(Plim-1,1);
        Far_rate_proposed_loc   = NaN(Plim-1,1);
        rate_LS_loc             = NaN(Plim-1,1);
        g_NMSE_proposed_loc     = NaN(Plim-1,1);
        Far_g_NMSE_proposed_loc = NaN(Plim-1,1);
        g_NMSE_LS_loc           = NaN(Plim-1,1);
        
        V_ML   = sqrt(sigma1_sq/2) * (randn(M,Plim) + 1i*randn(M,Plim));
        noiseL = sqrt(sigma2_sq/2) * (randn(Plim,1) + 1i*randn(Plim,1));
        
        h_hat = h_hat_bank(:, n2);
        Dh_hat = diag(h_hat);
        Dh_hat_angles = diag(h_hat ./ max(abs(h_hat), eps));
        
        phases_all_NF_loc = Dh_hat_angles * beamresponses;
        phases_all_FF_loc = Dh_hat_angles * Farbeamresponses;
        
        % ====== Proposed (near-field) ===================================
        if porpose
            utilize = false(size(phases_all_NF_loc,2),1); utilize(end-1:end) = true;
            amplitude = sqrt(P_ris/M);
            B_all = zeros(Plim, M);
            B_all(1:2,:) = (amplitude * phases_all_NF_loc(:,utilize)).';
            
            for L = 2:Plim
                idx = L-1;
                B_L = B_all(1:L,:);
                B_Dh_V = diag( B_L * (Dh * V_ML(:,1:L)) );   
                y = sqrt(Pp)*(B_L*Dh*g) + B_Dh_V + noiseL(1:L);
                
                [~,~,g_est,~,~] = MLE3D_corrected(y, L, B_L, Dh_hat, a_range, ...
                    lambda, M_V, M_H, d_V, d_H, varphi_range, theta_range, Pp);
                
                c_align = (g_est' * g) / (g_est' * g_est);
                g_est_aligned = g_est * c_align;
                g_NMSE_proposed_loc(idx) = norm(g_est_aligned - g)^2 / max(norm(g)^2, eps);
                
                alpha2 = abs(g_est) .* abs(h_hat);
                beta  = abs(h_hat).^2;
                gamma = (abs(g_est).^2*SNR1 + 1)/SNR2;
                C     = sqrt(sum((abs(alpha2).^2 .* gamma)./(beta+gamma).^2))^(-1);
                p_k_est = (alpha2 ./ (beta + gamma)) * C;
                RISconfig = angle(Dh_hat * g_est);
                Phi = p_k_est .* exp(-1i*RISconfig);
                
                True_RIS_Power = sum( abs(Phi).^2 .* (abs(g).^2 * Pd + sigma1_sq) );
                Phi = Phi * sqrt(P_ris / True_RIS_Power);
                
                signal       = Phi.' * Dh * g;
                signal_power = abs(signal)^2;
                RIS_noise_proj = norm(Phi.'*Dh)^2;
                noise_power  = sigma2_sq + sigma1_sq*RIS_noise_proj;
                rate_proposed_loc(idx) = log2(1 + (Pd*signal_power)/max(noise_power,eps));
                
                if L < Plim
                    unusedIndices = find(~utilize);
                    phaseCands = phases_all_NF_loc(:,unusedIndices);
                    phi_ref   = p_k_est .* exp(-1i*RISconfig);
                    beam_gain = abs(phi_ref' * phaseCands);
                    [~,bestK] = max(beam_gain);
                    bestIdx   = unusedIndices(bestK);
                    utilize(bestIdx) = true;
                    Phi_new = p_k_est .* phases_all_NF_loc(:,bestIdx);
                    
                    True_Dict_Power = sum( abs(Phi_new).^2 .* (abs(g).^2 * Pd + sigma1_sq) );
                    Phi_new = Phi_new * sqrt(P_ris / True_Dict_Power);
                    B_all(L+1,:) = Phi_new.';
                end
            end
            g_NMSE_proposed(:, n1, n2) = g_NMSE_proposed_loc;
            rate_proposed(:, n1, n2)   = rate_proposed_loc;
        end
        
        % ====== Far-field approximation =================================
        if Far_approx
            utilize = false(size(phases_all_FF_loc,2),1); utilize(end-1:end) = true;
            amplitude = sqrt(P_ris/M);
            B_all = zeros(Plim, M);
            B_all(1:2,:) = (amplitude * phases_all_FF_loc(:,utilize)).';
            
            for L = 2:Plim
                idx = L-1;
                B_L = B_all(1:L,:);
                B_Dh_V = diag( B_L * (Dh * V_ML(:,1:L)) );   
                y = sqrt(Pp)*(B_L*Dh*g) + B_Dh_V + noiseL(1:L);
                
                [~,~,g_est,~,~] = MLE_corrected(y, L, B_L, Dh_hat, a_FarAppx_range, Pp);
                Far_g_NMSE_proposed_loc(idx) = norm(g_est - g)^2 / max(norm(g)^2,eps);
                
                alpha2 = abs(g_est) .* abs(h_hat);
                beta  = abs(h_hat).^2;
                gamma = (abs(g_est).^2*SNR1 + 1)/SNR2;
                C     = sqrt(sum((abs(alpha2).^2 .* gamma)./(beta+gamma).^2))^(-1);
                p_k_est = (alpha2 ./ (beta + gamma)) * C;
                RISconfig = angle(Dh_hat * g_est);
                Phi = p_k_est .* exp(-1i*RISconfig);
                
                True_RIS_Power = sum( abs(Phi).^2 .* (abs(g).^2 * Pd + sigma1_sq) );
                Phi = Phi * sqrt(P_ris / True_RIS_Power);
                
                signal       = Phi.' * Dh * g;
                signal_power = abs(signal)^2;
                RIS_noise_proj = norm(Phi.'*Dh)^2;
                noise_power  = sigma2_sq + sigma1_sq*RIS_noise_proj;
                Far_rate_proposed_loc(idx) = log2(1 + (Pd*signal_power)/max(noise_power,eps));
                
                if L < Plim
                    unusedIndices = find(~utilize);
                    phaseCands = phases_all_FF_loc(:,unusedIndices);
                    phi_ref   = p_k_est .* exp(-1i*RISconfig);
                    beam_gain = abs(phi_ref' * phaseCands);
                    [~,bestK] = max(beam_gain);
                    bestIdx   = unusedIndices(bestK);
                    utilize(bestIdx) = true;
                    Phi_new = p_k_est.* phases_all_FF_loc(:,bestIdx);
                    
                    True_Dict_Power = sum( abs(Phi_new).^2 .* (abs(g).^2 * Pd + sigma1_sq) );
                    Phi_new = Phi_new * sqrt(P_ris / True_Dict_Power);
                    B_all(L+1,:) = Phi_new.';
                end
            end
            Far_g_NMSE_proposed(:, n1, n2) = Far_g_NMSE_proposed_loc;
            Far_rate_proposed(:, n1, n2)   = Far_rate_proposed_loc;
        end
        
        % ====== LS Estimator ===========================================
        if LSConf
            DFT = fft(eye(M));
            randomOrdering = randperm(M);
            amp_train = sqrt(P_ris/M);
            
            for L = 2:Plim
                idx = L-1;
                B_L = amp_train * transpose(DFT(:, randomOrdering(1:L)));
                B_Dh_V = diag( B_L * (Dh * V_ML(:,1:L)) );   
                y = sqrt(Pp)*(B_L*Dh*g) + B_Dh_V + noiseL(1:L);
                
                A_mat = sqrt(Pp) * (B_L * Dh_hat);
                g_est_LS = pinv(A_mat) * y; 
                g_NMSE_LS_loc(idx) = norm(g_est_LS - g)^2 / max(norm(g)^2,eps);
                
                alpha2 = abs(g_est_LS).*abs(h_hat);
                beta  = abs(h_hat).^2;
                gamma = (abs(g_est_LS).^2*SNR1 + 1)/SNR2;
                C     = sqrt(sum((abs(alpha2).^2 .* gamma)./(beta+gamma).^2))^(-1);
                p_k_est = (alpha2 ./ (beta + gamma)) * C;
                
                RISconfig_raw = angle(Dh_hat*g_est_LS);
                Phi = p_k_est .* exp(-1i*RISconfig_raw);
                
                True_RIS_Power = sum( abs(Phi).^2 .* (abs(g).^2 * Pd + sigma1_sq) );
                Phi = Phi * sqrt(P_ris / True_RIS_Power);
                
                signal       = Phi.' * Dh * g;
                signal_power = abs(signal)^2;
                RIS_noise_proj = norm(Phi.'*Dh)^2;
                noise_power  = sigma2_sq + sigma1_sq*RIS_noise_proj;
                rate_LS_loc(idx) = log2(1 + (Pd*signal_power)/max(noise_power,eps));
            end
            g_NMSE_LS(:, n1, n2) = g_NMSE_LS_loc;
            rate_LS(:, n1, n2)   = rate_LS_loc;
        end
        
    end
end
disp('Simulation Complete!');

%% ========================================================================
% ======================== DATA AVERAGING & PLOTTING ======================
% ========================================================================

% Calculate averages over both angle and noise realizations
mean_capacity      = mean(capacity);
mean_rate_proposed = mean(rate_proposed, [2 3], 'omitnan');
mean_rate_far      = mean(Far_rate_proposed, [2 3], 'omitnan');
mean_rate_LS       = mean(rate_LS, [2 3], 'omitnan');

mean_NMSE_proposed = mean(g_NMSE_proposed, [2 3], 'omitnan');
mean_NMSE_far      = mean(Far_g_NMSE_proposed, [2 3], 'omitnan');
mean_NMSE_LS       = mean(g_NMSE_LS, [2 3], 'omitnan');

% Convert NMSE to dB
NMSE_proposed_dB = 10 * log10(mean_NMSE_proposed);
NMSE_far_dB      = 10 * log10(mean_NMSE_far);
NMSE_LS_dB       = 10 * log10(mean_NMSE_LS);

L_range = 2:Plim;

% --- Figure 1: Achievable Rate vs Number of Pilots ---
figure('Name','Achievable Rate','Color','w');
plot(L_range, mean_rate_proposed, '-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
if Far_approx
    plot(L_range, mean_rate_far, '-^', 'LineWidth', 2, 'MarkerSize', 6); 
end
if LSConf
    plot(L_range, mean_rate_LS, '-s', 'LineWidth', 2, 'MarkerSize', 6); 
end
yline(mean_capacity, '--k', 'Capacity Bound', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

grid on;
xlabel('Number of Pilots (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Achievable Rate (bps/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Active RIS Achievable Rate vs. Training Overhead', 'FontSize', 14);

% Dynamically build legend
leg_entries = {'Proposed Parametric MLE (Near-Field)'};
if Far_approx, leg_entries{end+1} = 'Far-Field Approximation'; end
if LSConf, leg_entries{end+1} = 'Unstructured LS Estimator'; end
legend(leg_entries, 'Location', 'best', 'FontSize', 11);

% --- Figure 2: NMSE vs Number of Pilots ---
figure('Name','Normalized Mean Square Error','Color','w');
plot(L_range, NMSE_proposed_dB, '-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
if Far_approx
    plot(L_range, NMSE_far_dB, '-^', 'LineWidth', 2, 'MarkerSize', 6); 
end
if LSConf
    plot(L_range, NMSE_LS_dB, '-s', 'LineWidth', 2, 'MarkerSize', 6); 
end

grid on;
xlabel('Number of Pilots (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Channel Estimation Error vs. Training Overhead', 'FontSize', 14);
legend(leg_entries, 'Location', 'best', 'FontSize', 11);



