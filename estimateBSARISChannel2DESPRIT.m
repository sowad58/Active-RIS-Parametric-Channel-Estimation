function [h_hat,varphiHat,thetaHat,muH,muV,bestResidual] = estimateBSARISChannel2DESPRIT(Y_h,pilotSymbols,pilotPower,sensorIndexMatrix,lambda,M_V,M_H,d_V,d_H,Rate)
    [N_s_V,N_s_H] = size(sensorIndexMatrix);
    [U_s,~,~] = svd(Y_h,'econ');
    dominantSubspace = U_s(:,1);
    U_grid = reshape(dominantSubspace,N_s_V,N_s_H);
    
    U_H1 = U_grid(:,1:end-1);
    U_H2 = U_grid(:,2:end);
    psiH = (U_H1(:)'*U_H2(:)) / max(real(U_H1(:)'*U_H1(:)),eps);
    muH = angle(psiH);
    
    U_V1 = U_grid(1:end-1,:);
    U_V2 = U_grid(2:end,:);
    psiV = (U_V1(:)'*U_V2(:)) / max(real(U_V1(:)'*U_V1(:)),eps);
    muV = angle(psiV);
    
    pilotEnergy = real(pilotSymbols*pilotSymbols');
    z_h = Y_h*pilotSymbols' / max(sqrt(pilotPower)*pilotEnergy,eps);
    
    signCandidates = [-1, -1; -1, 1; 1, -1; 1, 1];
    bestResidual = inf; bestShape = []; bestGain = 0; varphiHat = 0; thetaHat = 0;
    
    for candidateID = 1:size(signCandidates,1)
        signH = signCandidates(candidateID,1);
        signV = signCandidates(candidateID,2);
        
        sinThetaCandidate = min(max(real(signV*muV/(2*pi*d_V)),-1),1);
        thetaCandidate = asin(sinThetaCandidate);
        cosThetaCandidate = cos(thetaCandidate);
        if abs(cosThetaCandidate) < 1e-10, continue; end
        
        sinVarphiCandidate = min(max(real(signH*muH/(2*pi*d_H*cosThetaCandidate)),-1),1);
        varphiCandidate = asin(sinVarphiCandidate);
        
        h_shape_candidate = UPA_Evaluation(lambda,M_V,M_H,varphiCandidate,thetaCandidate,d_V,d_H,Rate);
        h_shape_candidate = h_shape_candidate(:);
        sensedShape = h_shape_candidate(sensorIndexMatrix);
        sensedShape = sensedShape(:);
        
        complexGainCandidate = (sensedShape'*z_h) / max(real(sensedShape'*sensedShape),eps);
        residualCandidate = norm(z_h-complexGainCandidate*sensedShape)^2;
        
        if residualCandidate < bestResidual
            bestResidual = residualCandidate; bestShape = h_shape_candidate; bestGain = complexGainCandidate;
            varphiHat = varphiCandidate; thetaHat = thetaCandidate;
        end
    end
    if isempty(bestShape), error('2-D ESPRIT failed.'); end
    h_hat = bestGain*bestShape;
end