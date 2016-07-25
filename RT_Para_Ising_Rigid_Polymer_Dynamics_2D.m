% Realtime_Ising_Rigid_Polymer_Dynamics_2D_in_parallel
% To simulate the dynamics of a rigid polymer by Sling model in 2D
% by Baihan Lin, Qian Lab
% July 2016

clear all;
close all;

%% Initialization

rng(1234);                 % randomizer

trial = 1000;              % trials
%twist = 200;            % change of state
node = 400;               % nodes of rigid polymer

global pathN;
pathN = strcat('/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/data/',num2str(node),'/');

system(['mkdir ' pathN]);

angle = 0.2;           % in rad, angle changed in each twist
L = 1;                  % length of each segment of rigid polymer
a = 10;                 % threshold to form loop

%T = 300;                % temperature (K)
%k = 1.38064852e-23;     % Boltzmann constant (J/K)
%b = 1/(k*T);           % thermodynamic beta
%kT = 4.11e-21;          % at 25°C (298 K)
%b = 1/kT;               % thermodynamic beta

Hc = 1.0;   % in unit of kT, energy level of cis rigid configuration
Ht = 0.9;   % in unit of kT, energy level of trans rigid configuration
b = 1;      % redefined beta based on H

Pc = exp(-b*Hc)/(exp(-b*Hc)+exp(-b*Ht));       % Probablity of cis change
Pt = exp(-b*Ht)/(exp(-b*Hc)+exp(-b*Ht));       % Probablity of trans change

%% Generate a polymer

% initial direction = positive x
% default change = counterclockwise angle
% in this vector, first and last node are zero
% from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)

p = createRandPolymer(node); % randomly generate
%p = createPolymer(node,Pc); % naturally generate

%% Construct a recorder

pTf = zeros(1,trial+1); % record times to form a loop
pEf = zeros(1,trial+1); % record energy to form a loop
pDf = zeros(1,trial+1); % record head-tail distances to form a loop

%% Simulation of twisting

pTf(1) = 0;
pEf(1) = pE(p, Hc, Ht);
pDf(1) = HTdist(p, L, angle);

parfor n = 2:trial+1
    % To twist till looped
    [Pnew, HTd, fin] = twistLoopRand(pathN, n, p, Pc, Pt, a, L, angle, Hc, Ht);
    pTf(n) = fin;
    pEf(n) = pE(Pnew, Hc, Ht);
    pDf(n) = HTd;
end


%% plot histograms

fig1 = figure;
histogram(pTf(2:trial+1), 'BinWidth', 20);
title(strcat('Time Histogram for ', num2str(node),'node-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pTf(2:trial+1)))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pTf(2:trial+1)))),'Color','red','FontSize',12);
filename = strcat(pathN, 'T-Hist-N',num2str(length(p)),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle),'.png');
saveas(gcf, filename,'png');
%close gcf;

fig2 = figure;
histogram(pEf(2:trial+1), 'BinWidth', 0.05);
% line([pEf(1) pEf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('Energy Histogram for ', num2str(node),'node-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)));
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pEf(2:trial+1)))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pEf(2:trial+1)))),'Color','red','FontSize',12);
filename = strcat(pathN, 'E-Hist-N',num2str(length(p)),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle),'.png');
saveas(gcf, filename,'png');
%close gcf;

fig3 = figure;
histogram(pDf(2:trial+1), 'BinWidth', 0.05);
% line([pDf(1),pDf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('HTdistance Histogram for ', num2str(node),'node-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pDf(2:trial+1)))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pDf(2:trial+1)))),'Color','red','FontSize',12);
filename = strcat(pathN, 'D-Hist-N',num2str(length(p)),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle),'.png');
saveas(gcf, filename,'png');
%close gcf;


%% Local functions

function fp = createNatPolymer(fnode,fPc)
% To create a natural polymer with node nodes based on probability

if rand() > 0.5
    fp(2:fnode-1) = 1;      % clockwise
else
    fp(2:fnode-1) = -1;     % counterclockwise
end

for no = 3:fnode-1
    if rand() > fPc
        fp(no:fnode-1) = -fp(no:fnode-1);      % trans
    end
end

disp(strcat('createNatPolymer: ', num2str(fp)));

end

function fp = createRandPolymer(fnode)
% To create a random polymer with node nodes

fp = zeros(1,fnode);
parfor no = 2:fnode-1
    if rand() > 0.5
        fp(no) = 1;      % clockwise
    else
        fp(no) = -1;     % counterclockwise
    end
end

disp('------Simulation-Starts-------');
disp(strcat('createRandPolymer: ', num2str(fp)));
disp('--trials--');
end

function [fPnew, fHTd, ffin] = twistLoopRand(fpath, ft, fp, fPc, fPt, fa, fL, fangle, fHc, fHt)
% To twist randomly till formed a loop

disp(strcat('T-',num2str(ft),'-------------'));

fPnew = fp;
ffin = 1;

stair = length(fp)-1;

fig = figure;
fv = visV(buildV(fPnew, fL, fangle));
xt = fv(1,:);
yt = fv(2,:);
plot(xt(1:stair), yt(1:stair));
xlmin = min(xt);
xlmax = max(xt);
ylmin = min(yt);
ylmax = max(yt);
axis([ xlmin, xlmax, ylmin, ylmax]);
% quiver(0, 0, xt(stair), yt(stair),0,'r');
grid;
xlabel 'x';
ylabel 'y';
title(strcat('Simulation of ',num2str(length(fp)),' node rigid polymer dynamics'));
text(0,0,strcat('ffin',num2str(ffin)));

% disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));

while HTdist(fPnew, fL, fangle) > fa
    
    no =  randsample(2:length(fp)-1,1);
    
    Phypo = [fPnew(1:no-1), -fPnew(no), fPnew(no+1:end)];    % hypothetical change
    Pchg = pE(Phypo, fHc, fHt)/(pE(Phypo, fHc, fHt)+pE(fPnew, fHc, fHt));
    if rand() < Pchg
        fPnew = Phypo;      % change state
    end
    if HTdist(fPnew, fL, fangle) < fa
        break;
    end
    
    ffin = ffin+1;
    
    stair = length(fp)-1;
    fv = visV(buildV(fPnew, fL, fangle));
    xt = fv(1,:);
    yt = fv(2,:);
    
    plot(xt(1:stair), yt(1:stair));
    title(strcat('Simulation of T',num2str(ft),'-',num2str(length(fp)),' node rigid polymer dynamics'));
    grid;
    xlmin = min(min(xt),xlmin);
    xlmax = max(max(xt),xlmax);
    ylmin = min(min(yt),ylmin);
    ylmax = max(max(yt),ylmax);
    axis([ xlmin, xlmax, ylmin, ylmax]);
%     quiver(0, 0, xt(stair), yt(stair),0,'r');
    xc = xlim;
    xl = xc(1)*0.2+xc(2)*0.8;
    yc = ylim;
    yl1 = yc(1)*0.17+yc(2)*0.83;
    yl2 = yc(1)*0.23+yc(2)*0.77;
    text(xl,yl1,strcat('ffin=',num2str(ffin)),'Color','red','FontSize',12);
    text(xl,yl2,strcat('HTdist=',num2str(HTdist(fPnew, fL, fangle))),'Color','red','FontSize',12);
    drawnow;
    
    %     pause(0.6)
    
    disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));
    
end

fHTd = HTdist(fPnew, fL, fangle);

filename = strcat(fpath, 'N',num2str(length(fp)),'-T',num2str(ft),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

disp(strcat('final state: ', num2str(fPnew)));
disp(strcat('finish time: ', num2str(ffin)));
disp(strcat('finish dist: ', num2str(HTdist(fPnew, fL, fangle))));

end


function fE = pE(fp, fHc, fHt)
% To calculate energy of a certain polymer state

fE = 0;
fs = diag(fp(2:length(fp)-2).'*fp(3:length(fp)-1));

parfor no = 1:length(fp)-3
    if fs(no) == 1
        
        fE = fE + fHc;      % cis
    else
        fE = fE + fHt;     % trans
    end
end

end

function fD = HTdist(fp, fL, fangle)
% To calculate the head-tail distance of a polymer

fm = buildV(fp, fL, fangle);
fD = norm(sum(fm,2));

end

function fm = buildV(fp, fL, fangle)
% To build a vector set (matrix) based on given polymer states

fx = ones(1, length(fp)-1)*fL;
fy = zeros(1, length(fp)-1);
fm = [fx;fy];
rot1 = [cos(fangle) -sin(fangle); sin(fangle) cos(fangle)]; % clockwise
rot2 = [cos(-fangle) -sin(-fangle); sin(-fangle) cos(-fangle)]; % counterclockwise

for no = 2:length(fp)-1
    if fp(no) == 1
        fm(:,no:end) = rot1*fm(:,no:end);       % clockwise
    else
        fm(:,no:end) = rot2*fm(:,no:end);       % counterclockwise
    end
end

end

function fv = visV(fm)
% To build a vector set (matrix) based on given polymer states

fv = fm;
for no = 2:length(fm)
    fv(:,no) = fv(:,no-1)+fv(:,no);
end

end

