% Histogram reader
%
% by Baihan Lin
% August 2016

clear all;
close all;

%% pEf

pEf = dlmread('/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/lab_sim/data/data_20160818/');

fig1 = figure;
histogram(pEf, 'BinWidth', 1);
title('Energy Histogram in rEquilibrium for N  -t  -a  -l -r  ');
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pEf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pEf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'rEquil-E-Hist-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

%% pDf

pDf = dlmread('/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/lab_sim/data/data_20160818/');

fig2 = figure;
histogram(pDf, 'BinWidth', 20);
title('HTdistance Histogram in rEquilibrium for N  -t  -a  -l -r  ');
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pDf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pDf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'rEquil-D-Hist-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

%% rAuto

AutoT = 2000;

Cor = dlmread('/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/lab_sim/data/data_20160818/');

fig = figure;
plot(0:AutoT, Cor);
title(fname);
parsaveas(gcf, strcat('   .png'),'png');
close gcf;

%% rTau
% run after rAuto

n = length(Cor);
x = (1:n).';
X = [ones(n,1) x];
y = log(Cor.');
b = X\y;
yCalc = X*b;

fig1 = figure;
scatter(x,y);
hold on
plot(x,yCalc);
xlabel('t')
ylabel('log(Autocorrelation)')
title(strcat('Autocorrelation to find relaxation tau(',fname,')'));
legend('Data','Slope & Intercept','Location','best');
grid on

tau = 1/b(2);

xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl = yc(1)*0.2+yc(2)*0.8;
text(xl,yl,strcat('tau=',num2str(tau)),'Color','red','FontSize',12);

parsaveas(gcf, strcat(fpathN, fname, '.png'),'png');
close gcf;

save(strcat(fpathN, fname, '.txt'), 'tau', '-ascii');
