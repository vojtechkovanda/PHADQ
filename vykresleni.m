set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

figure
bar([ODGavg',ODGiavg', ODGsavg', ODGqavg']);
xlabel('word length (bps)');
ylabel('PEMO-Q ODG');
legend('B-PHADQ (consisitent)', 'B-PHADQ (inconsistent)', 'CP (sparsity based)', 'quantized');

figure
bar([SDRavg',SDRiavg', SDRsavg', SDRqavg']);
xlabel('word length (bps)');
ylabel('SDR (dB)');
legend('B-PHADQ (consisitent)', 'B-PHADQ (inconsistent)', 'CP (sparsity based)', 'quantized');

figure
bar([idxSDRavg',idxSDRiavg', idxODGavg', idxODGiavg']);
xlabel('word length (bps)');
ylabel('number of iteration');
legend('B-PHADQ (consisitent) SDR', 'B-PHADQ (inconsisitent) SDR', 'B-PHADQ (consisitent) ODG', 'B-PHADQ (inconsisitent) ODG');

