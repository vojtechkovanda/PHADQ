addpath('audio_dequantization-main\Algorithms');
addpath('audio_dequantization-main\Tools');
addpath('audio_dequantization-main');

sounds = {'a08_violin.wav', 'a16_clarinet.wav', 'a18_bassoon.wav', 'a25_harp.wav', 'a35_glockenspiel.wav', 'a41_celesta.wav', ...
          'a42_accordion.wav', 'a58_guitar_sarasate.wav', 'a60_piano_schubert.wav', 'a66_wind_ensemble_stravinsky.wav'};

ODG = zeros(1, 8);
idxODG = zeros(1, 8);
ODGq = zeros(1, 8);
SDR = zeros(1, 8);
idxSDR = zeros(1, 8);
SDRq = zeros(1, 8);

SDRs = zeros(1, 8);
ODGs = zeros(1, 8);

ODGavg = zeros(1, 8);
idxODGavg = zeros(1, 8);
ODGqavg = zeros(1, 8);
SDRavg = zeros(1, 8);
idxSDRavg = zeros(1, 8);
SDRqavg = zeros(1, 8);

SDRsavg = zeros(1, 8);
ODGsavg = zeros(1, 8);

for j=2:10

for i=1:7

[SDR(i+1), idxSDR(i+1), SDRq(i+1), ODG(i+1), idxODG(i+1), ODGq(i+1)] = cp_main(sounds{j}, i+1);
[SDRs(i+1), ODGs(i+1)] = dequantization_main(sounds{j}, i+1);

end

ODGavg = ODGavg + ODG;
idxODGavg = idxODGavg + idxODG;
ODGqavg = ODGqavg + ODGq;
SDRavg = SDRavg + SDR;
idxSDRavg = idxSDRavg + idxSDR;
SDRqavg = SDRqavg + SDRq;

SDRsavg = SDRsavg + SDRs;
ODGsavg = ODGsavg + ODGs;


end

    ODGavg = ODGavg / 10;
    idxODGavg = idxODGavg / 10;
    ODGqavg = ODGqavg / 10;
    SDRavg = SDRavg / 10;
    idxSDRavg = idxSDRavg / 10;
    SDRqavg = SDRqavg / 10;

    SDRsavg = SDRsavg / 10;
    ODGsavg = ODGsavg / 10;

figure
plot(ODGavg);
hold on
plot(ODGqavg);
plot(ODGsavg);

figure
plot(SDRavg);
hold on
plot(SDRqavg);
plot(SDRsavg);

figure
plot(idxSDRavg);
hold on
plot(idxODGavg*5-10);