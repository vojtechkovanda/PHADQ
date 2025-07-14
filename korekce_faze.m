% Ukázka korekce fáze

addpath('phase_correction');

% generování vstupního signálu (například sinusovky, nebo wav souboru)
fs = 48000;
t = 0:1/fs:2;
x = 0.5*sin(2*pi*1186*t);


audiofile = 'clarinet_reference.wav';
[x, fs] = audioread(audiofile);

% příprava parametrů okna pro DGT
winlen = 2024;
wintype = 'hanning';
a = winlen/4;
M = winlen*2;


%% DGT + iPC (PHAIN)

insig = x;

[win, ~] = generalizedCosWin(winlen, 'hanning');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/winlen);
diff_win = numericalDiffWin(tight_win);

zeroPhaseFlag = true;
rotateFlag = true;

[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(length(insig), winlen, a, M);

G = @(x) FDGT(x, tight_win, sigIdx, M, rotIdx, zeroPhaseFlag);
G_adj = @(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*w;
G_diff = @(x) FDGT(x, diff_win, sigIdx, M, rotIdx, zeroPhaseFlag);

omega = @(x) calcInstFreq(G(x), G_diff(x), M, winlen, rotateFlag);

R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

% korekce fáze je použití operátoru R na spektrum dané operátorem G


%% Použití a vykreslení

% spektrogram signálu
X = G(x);

% výpočet okamžité frekvence
omegasin = omega(x);

% korekce fáze podle okamžité frekvence
X_c = R(X, omegasin);

% vykreslení
subplot(2,1,1);
phase_sgram(X);

subplot(2,1,2);
phase_sgram(X_c);

