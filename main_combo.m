% Chambolle-Pock sparsity based dequantization with Phase aware prior
%
% Vojtěch Kovanda
% Brno University of Technology, 2025

%% input signal
audiofile = 'clarinet_reference.wav';
[x, param.fs] = audioread(audiofile);

% signal length
param.L = length(x);

% normalization
maxval = max(abs(x));
x = x/maxval;

%% generate observation y

% setting conversion parameters
param.delta = 10;

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

paramsolver.tau = 1/2;  % step size
paramsolver.sigma = 1/2;  % step size
paramsolver.alpha = 1;  % relaxation parameter

paramsolver.lambda = 0.001;  % threshold (regularization parameter)

paramsolver.I = 100;
paramsolver.J = 1;

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



%%

soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);

    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    param.prox = @(z) soft(z, lambda/sigma);

    %omega_y = omega(insig);
    omega_y = omega(x);
    param.L = @(x) hatG(x, omega_y);
    param.L_adj = @(u) hatG_adj(u, omega_y);

    param.L1 = @(x) G(x);
    param.L1_adj = @(u) G_adj(u);

    paramsolver.x0 = insig;
    paramsolver.u1 = zeros(size(param.L1(zeros(length(insig), 1))));
    paramsolver.u2 = zeros(size(param.L(zeros(length(insig), 1))));
    x_old = insig;


        [x_hat] = CV(param, paramsolver, xq);
        

    outsig = x_hat;

[~, ~, ODG] = audioqual(x, outsig, param.fs);
SDR = 20*log10(norm(x,2)./norm(x-outsig, 2));