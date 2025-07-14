% vykreslení fáze u spektrogramu
% phaseplot

function phase_sgram(X)

phase = angle(X);
amplitude = abs(X);
level = 10*log10(amplitude);

bound = level > 1000 | level < -1000;
level(bound) = 0;

maxval = max(max(level(2:end, :)));

difc = 30-maxval;

level = level + difc;

maxval = max(max(level(2:end, :)));

level = level./maxval;

value = level;

maxvalhue = max(max(abs(phase)));
hue = phase./(2*maxvalhue) + 0.5;

saturation = ones(size(level));

hsv = cat(3, hue, saturation, value);

image(hsv2rgb(hsv));

