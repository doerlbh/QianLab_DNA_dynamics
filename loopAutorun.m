
% by Baihan Lin, August 2016

function [finP, stP, pTf, pEf, pDf] = loopAutorun(fnode, stP, ftrial,fpathN, fa, fL, fangle, fHc, fHt)
% simulate the looping and calculating autocorrelation for tau

% Construct a recorder

pTf = zeros(1,ftrial); % record times to form a loop
pEf = zeros(1,ftrial); % record energy to form a loop
pDf = zeros(1,ftrial); % record head-tail distances to form a loop
finP = zeros(ftrial,fnode); % record end states
fACorT = zeros(ftrial,fAutoT);

% for n = 1:ftrial
parfor n = 1:ftrial

% Generate a polymer

% initial direction = positive x
% default change = counterclockwise angle
% in this vector, first and last node are zero
% from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)

    p = stP(n,:);
    % To twist till looped
    [Pnew, ftP, HTd, fin] = fasttwistLoopRand(fpathN, n, p, fa, fL, fangle, fHc, fHt);
    pTf(n) = fin;
    pEf(n) = pE(Pnew, fHc, fHt);
    pDf(n) = HTd;
    finP(n,:) = Pnew;
    [fACT] = AutocorEnd(ftP);
    fACTt = [1:length(fACT);fACT];
    fACorT = [fACorT fACTt];
end

% plot histograms

fig1 = figure;
histogram(pTf, 'BinWidth', 50);
title(strcat('Time Histogram for N', num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)));
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pTf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pTf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Loop-T-Hist-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Loop-T-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
parsave(filename, pTf, '-ascii');

fig2 = figure;
histogram(pEf, 'BinWidth', 1);
% line([pEf(1) pEf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('Energy Histogram for N', num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)));
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pEf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pEf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Loop-E-Hist-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Loop-E-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
parsave(filename, pEf, '-ascii');

fig3 = figure;
histogram(pDf, 'BinWidth', 0.2);
% line([pDf(1),pDf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('HTdistance Histogram for N', num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pDf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pDf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Loop-D-Hist-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Loop-D-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
parsave(filename, pDf, '-ascii');

fnamela = strcat('Loop-Auto-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle));
fACor = unevenmean(fACorT);
fACor = fACor/fACor(1);

fig = figure;
plot(0:length(fACor)-1, fACor);
title(fnamela);
parsaveas(gcf, strcat(fpathN, fnamela, '.png'),'png');
close gcf;

parsave(strcat(fpathN, fnamela, '.txt'), fACor, '-ascii');

fnamelt = strcat('Loop-relaxTau-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)); 
relaxTau = Cor2tau(LACor, fpathN, fnamelt);

end


