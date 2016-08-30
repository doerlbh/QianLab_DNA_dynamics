
Cor = [ 1 10 5 3 3 2 10 10 3 7];

n = length(Cor);
x = (1:n).';
X = [ones(n,1) x];
y = real(log(Cor)).';
b = X\y;
yCalc = X*b;

fig1 = figure;
scatter(x,y);
hold on
plot(x,yCalc);
xlabel('t')
ylabel('log(Autocorrelation)')
% title(strcat('Autocorrelation to find relaxation tau(',fname,')'));
legend('Data','Slope & Intercept','Location','best');
grid on


% n = 20;
% x = (4:23).';
% X = [ones(n,1) x];
% y = (1:20).';
% b = X\y;
% yCalc = X*b;
% 
% scatter(x,y);
% hold on
% plot(x,yCalc);
% xlabel('t')
% ylabel('log(Autocorrelation)')
% legend('Data','Slope & Intercept','Location','best');
% grid on
% 
% tau = 1/b(2);
% 
% xc = xlim;
% xl = xc(1)*0.2+xc(2)*0.8;
% yc = ylim;
% yl = yc(1)*0.2+yc(2)*0.8;
% text(xl,yl,strcat('tau=',num2str(tau)),'Color','red','FontSize',12);
% 
% 
