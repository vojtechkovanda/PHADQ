% Chambolle-Pock sparsity based dequantization with Phase aware prior
%
% Vojtěch Kovanda
% Brno University of Technology, 2025

addpath('phase_correction');
addpath('dataset');
addpath('PEMO-Q');


%% input signal
audiofile = 'a42_accordion.wav';
[x, param.fs] = audioread(audiofile);

% signal length
param.L = length(x);

% normalization
maxval = max(abs(x));
x = x/maxval;

%% generate observation y

% setting conversion parameters
param.delta = 8;

% quantization
xq = quant(x, param.delta);

%% parameters

winLen = 2048;
shiftLen = winLen/4;
FFTnum = 2*winLen;

% parameter setting for PHAINs
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;

paramsolver.epsilon = 0.001;  % for stopping criterion

paramsolver.tau = 1;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.alpha = 1;  % relaxation parameter

paramsolver.lambda = [1, 1, 0.1, 0.01, 0.001, 0.001, 0.0001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001];  % threshold (regularization parameter)

paramsolver.I = 200;
paramsolver.J = 2;
paramsolver.K = 200;

%% iPC DGT

insig = xq;

a = param.a;
M = param.M;
w = param.w;

[win, ~] = generalizedCosWin(w, 'hanning');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/w);
diff_win = numericalDiffWin(tight_win);
    
zeroPhaseFlag = true;
rotateFlag = true;

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


% SPADQ settings
param.Ls = length(x);
param.lam = [0.1; 0.0038; 0.0023; 0.0012; 0.000094; 0.000032; 0.000013; 0.0000055];

    param.F = frametight(frame('dgtreal', {'hanning', param.w}, param.a, param.M));
    param.F = frameaccel(param.F, param.Ls);  % precomputation for a fixed signal length

    param.rho = 1;

%%

soft = @(z, lambda) sign(z).*max(abs(z) - lambda(param.delta), 0);

    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    param.prox = @(z) soft(z, lambda/sigma);

    %omega_y = omega(insig);
    omega_y = omega(xq);
    param.L = @(x) hatG(x, omega_y);
    param.L_adj = @(u) hatG_adj(u, omega_y);

    paramsolver.x0 = insig;
    paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
    x_old = insig;

    SDR_in_time_all = zeros(1, (paramsolver.J-1) * (paramsolver.I + paramsolver.K));


    for j = 1:paramsolver.J

        if j == 1

        [x_hat, SDR_in_time] = CP(param, paramsolver, xq, x);
        paramsolver.x0 = x_hat;

        SDR_in_time_all(1:paramsolver.I) = SDR_in_time;

        else

        [x_hat, SDR_in_time] = cp_alg(insig, param, paramsolver, x);

        SDR_in_time_all(paramsolver.I+1:end) = SDR_in_time;

        end       

    end



    outsig = x_hat;

    figure;
    plot(SDR_in_time_all);
   
[~, ~, ODG] = audioqual(x, outsig, param.fs);
SDR = 20*log10(norm(x,2)./norm(x-outsig, 2));

[~, ~, ODGq] = audioqual(x, insig, param.fs);
SDRq = 20*log10(norm(x,2)./norm(x-insig, 2));

fprintf('SDR of the quantized signal is %4.3f dB.\n', SDRq);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);
fprintf('ODG of the quantized signal is %4.3f.\n', ODGq);
fprintf('ODG of the reconstructed signal is %4.3f.\n', ODG);