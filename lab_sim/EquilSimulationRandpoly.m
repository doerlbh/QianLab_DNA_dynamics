% developed by Baihan Lin, August 2016

function [finP, stP, pEf, pDf] = EquilSimulationRandpoly(fnode, ftrial,ftwist, fpathN, fa, fL, fangle, fHc, fHt)
% simulate the equilibrium distribution from a random state

% Construct a recorder

pEf = zeros(1,ftrial); % record energy in equilibrium
pDf = zeros(1,ftrial); % record head-tail distances in equilibrium
finP = zeros(ftrial,fnode); % record end states
stP = zeros(ftrial,fnode); % record inital states

% for n = 1:ftrial
parfor n = 1:ftrial
    % Generate a polymer
    
    % initial direction = positive x
    % default change = counterclockwise angle
    % in this vector, first and last node are zero
    % from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)
    
    p = createRandPolymer(fnode); % randomly generate
    stP(n,:) = p;
    [Pnew, ftP, HTd] = fasttwistEquilRand(fpathN, ftwist, n, p, fa, fL, fangle, fHc, fHt);
    pEf(n) = pE(Pnew, fHc, fHt);
    pDf(n) = HTd;
    finP(n,:) = Pnew;
end

% plot histograms

filename = strcat(fpathN, 'Equil-finP-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
parsave(filename, finP, '-ascii');

fig1 = figure;
histogram(pEf, 'BinWidth', 0.2);
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
parsave(filename, pEf, '-ascii');

fig2 = figure;
histogram(pDf, 'BinWidth', 0.2);
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
parsave(filename, pDf, '-ascii');

end
