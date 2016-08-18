function [fACor] = AutocorEq(stP, ftrial, fAutoT, fpathN, fa, fL, fangle, fHc, fHt, fname)
% calculate the autocorrelation from any intial states

% Construct a recorder
fACorT = zeros(ftrial,fAutoT+1);

% for n = 1:ftrial
parfor n = 1:ftrial
    p = stP(n,:);
    [Pnew, ftP, HTd] = fasttwistEquilRand(fpathN, fAutoT, n, p, fa, fL, fangle, fHc, fHt);
    [fACT] = AutocorEnd(ftP);
    fACorT = [fACorT; fACT];
end

fACor = mean(fACorT);

fig = figure;
plot(0:fAutoT, fACor);
title(fname);
parsaveas(gcf, strcat(fpathN, fname, '.png'),'png');
close gcf;

save(strcat(fpathN, fname, '.txt'), 'fACor', '-ascii');

end