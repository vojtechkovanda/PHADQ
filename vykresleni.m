

figure
bar([ODGavg',ODGiavg', ODGsavg', ODGqavg']);
title('Average Values of ODG');
xlabel('bit depth');
ylabel('Average Values');
legend('PHADQ cons', 'PHADQ incons', 'CP sparsity based', 'quantized');

figure
bar([SDRavg',SDRiavg', SDRsavg', SDRqavg']);
title('Average Values of SDR');
xlabel('bit depth');
ylabel('Average Values');
legend('PHADQ cons', 'PHADQ incons', 'CP sparsity based', 'quantized');

figure
bar([idxSDRavg',idxSDRiavg', idxODGavg', idxODGiavg']);
title('Average number of best iteration ODG or SDR');
xlabel('bit depth');
ylabel('Average Values');
legend('PHADQ cons best SDR', 'PHADQ incons best SDR', 'PHADQ cons best ODG', 'PHADQ incons best ODG');

