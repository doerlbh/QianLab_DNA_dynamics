% by Baihan Lin, July 2016

function [finP, stP, pEf, pDf] = EquilSimulation1poly(fnode, ftrial,ftwist, fpathN, fa, fL, fangle, fHc, fHt)
% simulate the equilibrium distribution from one single state

% Generate a polymer

% initial direction = positive x
% default change = counterclockwise angle
% in this vector, first and last node are zero
% from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)

p = createRandPolymer(fnode); % randomly generate

% Construct a recorder

pEf = zeros(1,ftrial); % record energy in equilibrium
pDf = zeros(1,ftrial); % record head-tail distances in equilibrium
finP = zeros(ftrial,fnode); % record end states
stP = zeros(ftrial,fnode); % record inital states

% for n = 1:ftrial
parfor n = 1:ftrial
    stP(n,:) = p;
    [Pnew, ftP, HTd] = fasttwistEquilRand(fpathN, ftwist, n, p, fa, fL, fangle, fHc, fHt);
    pEf(n) = pE(Pnew, fHc, fHt);
    pDf(n) = HTd;
    finP(n,:) = Pnew;
end

% plot histograms

fig1 = figure;
histogram(pEf, 'BinWidth', 0.1);
% line([pEf(1) pEf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('Energy Histogram in Equilibrium for N', num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)));
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pEf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pEf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Equil-E-Hist-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Equil-E-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pEf', '-ascii');

fig2 = figure;
histogram(pDf, 'BinWidth', 10);
% line([pDf(1),pDf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('HTdistance Histogram in Equilibrium for N', num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pDf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pDf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Equil-D-Hist-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Equil-D-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pDf', '-ascii');

end