inputFolder = 'dataset/EBU_SQAM';   % vstupní složka
outputFolder = 'dataset/EBU_SQAM';  % výstupní složka (může být i jiná)

for k = 9:10
    % načti soubor
    filename = fullfile(inputFolder, sprintf('%d.wav', k));
    if ~isfile(filename)
        warning('Soubor %s neexistuje, přeskočeno.', filename);
        continue;
    end
    
    [y, Fs] = audioread(filename);
    
    % vezmi jen 1. kanál (mono)
    if size(y,2) > 1
        y = y(:,1);
    end
    
    % počet vzorků pro 7 sekund
    nSamples = min(length(y), round(7 * Fs));
    
    % zkrať signál
    y_short = y(1:nSamples);
    
    % ulož výsledek jako "X u.wav"
    outname = fullfile(outputFolder, sprintf('%du.wav', k));
    audiowrite(outname, y_short, Fs);
end

disp('Hotovo!');