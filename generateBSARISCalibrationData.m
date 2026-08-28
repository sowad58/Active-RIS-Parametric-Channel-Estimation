function [Y_h,noiseVariance] = generateBSARISCalibrationData(h_true,sensorIndexMatrix,pilotSymbols,pilotPower,hCalSNR_dB)
    h_s_true = h_true(sensorIndexMatrix);
    h_s_true = h_s_true(:);
    noiselessSignal = sqrt(pilotPower) * h_s_true * pilotSymbols;
    receivedSignalPower = mean(abs(noiselessSignal(:)).^2);
    
    if isinf(hCalSNR_dB)
        noiseVariance = 0;
        calibrationNoise = complex(zeros(size(noiselessSignal)));
    else
        noiseVariance = receivedSignalPower / max(10^(hCalSNR_dB/10),eps);
        calibrationNoise = sqrt(noiseVariance/2) * (randn(size(noiselessSignal)) + 1i*randn(size(noiselessSignal)));
    end
    Y_h = noiselessSignal + calibrationNoise;
end
