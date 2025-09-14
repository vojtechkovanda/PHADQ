% Chambolle-Pock dequantization with Phase aware prior
%
% Vojtěch Kovanda
% Brno University of Technology, 2025

function [SDR, ODG] = cp_main_oracle(audiofile, wordlength, method)

%wordlength = 6;

addpath('phase_correction');
addpath('dataset');
addpath('PEMO-Q');
addpath('Sounds');

%method = 'consistent'; % 'inconsistent'

%% input signal
% audiofile = 'a66_wind_ensemble_stravinsky.wav';
[x, param.fs] = audioread(audiofile);

%[x, param.fs] = audioread('a08_violin.wav');

% signal length
param.L = length(x);

% normalization
maxval = max(abs(x));
x = x/maxval;

%% generate observation y

% setting conversion parameters
param.delta = wordlength;

% quantization
xq = quant(x, param.delta);

%% parameters

winLen = 8192;
shiftLen = winLen/4;
FFTnum = 2*winLen;

% parameter setting for PHAINs
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;

paramsolver.epsilon = 0.001;  % for stopping criterion

paramsolver.tau = 1;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.alpha = 1/3;  % relaxation parameter

paramsolver.lambda = [2, 0.1, 0.1, 0.001, 0.01, 0.001, 0.0001, 0.0001, 0.00001, 0.00001, 0.00001, 0.00001];  % threshold (regularization parameter)

% if method == 'inconsistent'
%     paramsolver.lambda(8) = 0.00001;
% end
paramsolver.I = 200;
paramsolver.J = 1;

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

soft = @(z, lambda) sign(z).*max(abs(z) - lambda(param.delta), 0);

    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    param.prox = @(z) soft(z, lambda/sigma);

    %omega_y = omega(insig);
    omega_y = omega(x);
    param.L = @(x) hatG(x, omega_y);
    param.L_adj = @(u) hatG_adj(u, omega_y);

    paramsolver.x0 = insig;
    paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
    x_old = insig;

    SDR_in_time_all = zeros(1, paramsolver.J * paramsolver.I);


    for j = 1:paramsolver.J


        switch method
            case 'consistent'

        [x_hat, SDR_in_time, ODG_in_time] = CP(param, paramsolver, xq, x);
            
            case 'inconsistent'

        [x_hat, SDR_in_time, ODG_in_time] = CP_incons(param, paramsolver, xq, x);

        end

        if norm(x_old - x_hat) < paramsolver.epsilon
            break
        end

        omega_x_hat = omega(x_hat);
        param.L = @(x) hatG(x, omega_x_hat);
        param.L_adj = @(u) hatG_adj(u, omega_x_hat);
        paramsolver.x0 = x_hat;

        x_old = x_hat;

        SDR_in_time_all((j-1)*paramsolver.I+1:j*paramsolver.I) = SDR_in_time;

    end

    outsig = x_hat;

    figure;
    plot(SDR_in_time_all)
    hold on

    %figure;
    %plot(ODG_in_time);

    [ODG, idxODG] = max(ODG_in_time);
    [SDR, idxSDR] = max(SDR_in_time);

   
% [~, ~, ODG] = audioqual(x, outsig, param.fs);
% SDR = 20*log10(norm(x,2)./norm(x-outsig, 2));

[~, ~, ODGq] = audioqual(x, insig, param.fs);
SDRq = 20*log10(norm(x,2)./norm(x-insig, 2));

% fprintf('SDR of the quantized signal is %4.3f dB.\n', SDRq);
% fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);
% fprintf('ODG of the quantized signal is %4.3f.\n', ODGq);
% fprintf('ODG of the reconstructed signal is %4.3f.\n', ODG);