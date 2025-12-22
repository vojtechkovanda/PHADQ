% Chambolle-Pock dequantization with Phase aware prior
%
% Vojtěch Kovanda
% Brno University of Technology, 2025


addpath('phase_correction');
addpath('dataset');

method = 'consistent'; % 'inconsistent'

%% input signal
audiofile = 'dataset/EBU_SQAM/8.wav';
[x, param.fs] = audioread(audiofile);

% signal length
param.L = length(x);

% normalization
maxval = max(abs(x));
x = x/maxval;

%% generate observation y

% setting conversion parameters
param.delta = 6;

% quantization
xq = quant(x, param.delta);

%% parameters

winLen = 8192;
shiftLen = winLen/4;
FFTnum = 2*winLen;

% parameter setting for B-PHADQ
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;

paramsolver.tau = 1;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.rho = 1/3;  % relaxation parameter

paramsolver.lambda = [1, 0.1, 0.1, 0.001, 0.01, 0.001, 0.0001, 0.0001, 0.00001, 0.00001, 0.00001, 0.00001];  % threshold (regularization parameter)

paramsolver.I = 100;

%% iPC DGT

insig = xq;

a = param.a;
M = param.M;
w = param.w;

[win, ~] = generalizedCosWin(w, 'hann');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/w);
diff_win = numericalDiffWin(tight_win);
    
zeroPhaseFlag = true;
rotateFlag = true;

insig = zeroPaddingForDGT(xq, a, M);

[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(length(insig), w, a, M);

% DGTs (original, its adjoint, and one with the differentiated window)
G = @(x) FDGT(x, tight_win, sigIdx, M, rotIdx, zeroPhaseFlag);
G_adj = @(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*w;
G_diff = @(x) FDGT(x, diff_win, sigIdx, M, rotIdx, zeroPhaseFlag);

% function to calculate the instantaneous frequency of the input signal
omega = @(x) calcInstFreq(G(x), G_diff(x), M, w, rotateFlag);

% operator to correct phase rotation and its adjoint
R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

% time-directional difference
D = @(z) z(:,1:end-1) - z(:,2:end);
D_adj = @(z) [z(:,1), (z(:,2:end) - z(:,1:end-1)), -z(:,end)];

% iPC-DGT
hatG = @(x, omega) D(R(G(x), omega));
hatG_adj = @(u, omega) G_adj(R_adj(D_adj(u), omega));




%%
    omega_y = omega(insig);
    param.L = @(x) hatG(x, omega_y);
    param.L_adj = @(u) hatG_adj(u, omega_y);

    paramsolver.x0 = insig;
    paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));

        switch method
            case 'consistent'

        x_hat = CP(param, paramsolver, insig);
            
            case 'inconsistent'

        x_hat = CP_incons(param, paramsolver, insig);

        end

    outsig = x_hat(1:length(x));

   
SDR = 20*log10(norm(x,2)./norm(x-outsig, 2));
%SDRq = 20*log10(norm(x,2)./norm(x-insig, 2));

%fprintf('SDR of the quantized signal is %4.3f dB.\n', SDRq);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);
