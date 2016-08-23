% by Baihan Lin, August 2016

function tau = Cor2tau(Cor, fpathN, fname)
% find tau from autocorrelation

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

tau = -1/b(2);

xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl = yc(1)*0.2+yc(2)*0.8;
text(xl,yl,strcat('tau=',num2str(tau)),'Color','red','FontSize',12);

parsaveas(gcf, strcat(fpathN, fname, '.png'),'png');
close gcf;

save(strcat(fpathN, fname, '.txt'), 'tau', '-ascii');

% tau = zeros(1,n);
% for c = 1:n
%     tau(c) = -log(Cor(c))/c;
% end
% 
% fig1 = figure;
% plot(tau);
% title(fname);
% parsaveas(gcf, strcat(fpathN, fname, '.png'),'png');
% close gcf;
% 
% taur = real(tau);
% taui = imag(tau);
% 
% fig2 = figure;
% plot(taur);
% title(fname);
% parsaveas(gcf, strcat(fpathN, fname, '-r.png'),'png');
% close gcf;
% 
% fig3 = figure;
% plot(taui);
% title(fname);
% parsaveas(gcf, strcat(fpathN, fname, '-i.png'),'png');
% close gcf;
% 
% save(strcat(fpathN, fname, '-r.txt'), 'taur', '-ascii');
% save(strcat(fpathN, fname, '-i.txt'), 'taui', '-ascii');

end
