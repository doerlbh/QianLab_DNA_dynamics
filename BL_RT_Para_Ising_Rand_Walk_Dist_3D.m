% Realtime_Ising_Rigid_Polymer_Dynamics_3D_using_parallel_Rand_Walk_Dist
% To simulate the dynamics of a rigid polymer by Ising model in 3D
% by Baihan Lin, Qian Lab
% August 2016

clear all;
close all;

%% Initialization

rng(378);                 % randomizer

trial = 1000;              % trials
twist = 2000;            % change of set state changes
node = 500;               % nodes of rigid polymer

global pathN;
pathN = strcat('/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/data/3D-Eq-',num2str(node),'/');
system(['mkdir ' pathN]);

angle = 0.1;            % in rad, angle changed in each twist
L = 1;                  % length of each segment of rigid polymer
a = 20;                 % threshold to form loop

Hc = 1.0;   % in unit of kT, energy level of cis rigid configuration
Ht = 0.9;   % in unit of kT, energy level of trans rigid configuration
% b = 1;      % redefined beta based on H

% Pc = exp(-b*Hc)/(exp(-b*Hc)+exp(-b*Ht));       % Probablity of cis change
% Pc2 = exp(-b*Hc)/(2*exp(-b*Hc)+exp(-b*Ht));
% Pt = exp(-b*Ht)/(exp(-b*Hc)+exp(-b*Ht));       % Probablity of trans change
% Pt2 = exp(-b*Ht)/(2*exp(-b*Hc)+exp(-b*Ht));

%% Main functions

tic;
[pEf, pDf] = EquilSimulation1poly(node, trial, twist, pathN, a, L, angle, Hc, Ht);
tEq = toc;

tic;
[pEfr, pDfr] = EquilSimulationRandpoly(node, trial, twist, pathN, a, L, angle, Hc, Ht);
tEqr = toc;

disp(tEq);
disp(tEqr);

% [pTf, pEf, pDf] = loopSimulation(node, trial, pathN, a, L, angle, Hc, Ht);

%% Local functions

function [pEf, pDf] = EquilSimulation1poly(fnode, ftrial,ftwist, fpathN, fa, fL, fangle, fHc, fHt)
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

% for n = 1:ftrial
parfor n = 1:ftrial
    [Pnew, HTd] = twistEquilRand(fpathN, ftwist, n, p, fa, fL, fangle, fHc, fHt);
    pEf(n) = pE(Pnew, fHc, fHt);
    pDf(n) = HTd;
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
saveas(gcf, filename,'png');
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
saveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Equil-D-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pDf', '-ascii');

end

function [pEf, pDf] = EquilSimulationRandpoly(fnode, ftrial,ftwist, fpathN, fa, fL, fangle, fHc, fHt)
% simulate the equilibrium distribution from a random state

% Construct a recorder

pEf = zeros(1,ftrial); % record energy in equilibrium
pDf = zeros(1,ftrial); % record head-tail distances in equilibrium

% for n = 1:ftrial
parfor n = 1:ftrial
    % Generate a polymer
    
    % initial direction = positive x
    % default change = counterclockwise angle
    % in this vector, first and last node are zero
    % from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)
    
    p = createRandPolymer(fnode); % randomly generate
    
    [Pnew, HTd] = twistEquilRand(fpathN, ftwist, n, p, fa, fL, fangle, fHc, fHt);
    pEf(n) = pE(Pnew, fHc, fHt);
    pDf(n) = HTd;
end

% plot histograms

fig1 = figure;
histogram(pEf, 'BinWidth', 0.1);
% line([pEf(1) pEf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('Energy Histogram in rEquilibrium for N', num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)));
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pEf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pEf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Equil-E-Hist-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
saveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'rEquil-E-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pEf', '-ascii');

fig2 = figure;
histogram(pDf, 'BinWidth', 10);
% line([pDf(1),pDf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('HTdistance Histogram in rEquilibrium for N', num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pDf))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pDf))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Equil-D-Hist-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
saveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'rEquil-D-N',num2str(fnode),'-t',num2str(ftwist),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pDf', '-ascii');

end

function [pTf, pEf, pDf] = loopSimulation(fnode, ftrial,fpathN, fa, fL, fangle, fHc, fHt)
% simulate the looping events
% Generate a polymer

% initial direction = positive x
% default change = counterclockwise angle
% in this vector, first and last node are zero
% from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)

p = createRandPolymer(fnode); % randomly generate

% Construct a recorder

pTf = zeros(1,ftrial+1); % record times to form a loop
pEf = zeros(1,ftrial+1); % record energy to form a loop
pDf = zeros(1,ftrial+1); % record head-tail distances to form a loop

% Simulation of twisting

pTf(1) = 0;
pEf(1) = pE(p, fHc, fHt);
pDf(1) = HTdist(p, fL, fangle);

% for n = 2:ftrial+1
parfor n = 2:ftrial+1
    % To twist till looped
    [Pnew, HTd, fin] = twistLoopRand(fpathN, n, p, fa, fL, fangle, fHc, fHt);
    pTf(n) = fin;
    pEf(n) = pE(Pnew, fHc, fHt);
    pDf(n) = HTd;
end

% plot histograms

fig1 = figure;
histogram(pTf(2:ftrial+1), 'BinWidth', 50);
title(strcat('Time Histogram for N', num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pTf(2:ftrial+1)))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pTf(2:ftrial+1)))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Loop-T-Hist-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
saveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Loop-T-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pTf', '-ascii');

fig2 = figure;
histogram(pEf(2:ftrial+1), 'BinWidth', 0.2);
% line([pEf(1) pEf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('Energy Histogram for N', num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)));
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pEf(2:ftrial+1)))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pEf(2:ftrial+1)))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Loop-E-Hist-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
saveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Loop-E-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pEf', '-ascii');

fig3 = figure;
histogram(pDf(2:ftrial+1), 'BinWidth', 0.2);
% line([pDf(1),pDf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('HTdistance Histogram for N', num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)))
xc = xlim;
xl = xc(1)*0.2+xc(2)*0.8;
yc = ylim;
yl1 = yc(1)*0.17+yc(2)*0.83;
yl2 = yc(1)*0.23+yc(2)*0.77;
text(xl,yl1,strcat('mean=',num2str(mean(pDf(2:ftrial+1)))),'Color','red','FontSize',12);
text(xl,yl2,strcat('var=',num2str(var(pDf(2:ftrial+1)))),'Color','red','FontSize',12);
filename = strcat(fpathN, 'Loop-D-Hist-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
saveas(gcf, filename,'png');
close gcf;

filename = strcat(fpathN, 'Loop-D-N',num2str(fnode),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.txt');
save(filename, 'pDf', '-ascii');

end

function fp = createRandPolymer(fnode)
% To create a random polymer with node nodes

fp = zeros(1,fnode);
fp(2) = 0;

% for no = 3:fnode-1
parfor no = 3:fnode-1
    r = rand();
    if r < 1/3
        fp(no) = 0;      % not flip, stay trans
    else
        if 1/3 < r < 2/3
            fp(no) = 1;     % flip to cis CCW
        else
            fp(no) = -1;     % flip to cis CW
        end
    end
end

disp('------Simulation-Starts-------');
disp(strcat('createRandPolymer: ', num2str(fp)));
disp('--trials--');
end

function [fPnew, fHTd, ffin] = twistLoopRand(fpath, ft, fp, fa, fL, fangle, fHc, fHt)
% To twist randomly till formed a loop

disp(strcat('T-',num2str(ft),'-------------'));

fPnew = fp;
ffin = 1;

stair = length(fp)-1;

fig = figure;
fv = visV(buildV(fPnew, fL, fangle));
xt = fv(1,:);
yt = fv(2,:);
zt = fv(3,:);
plot3(xt(1:stair),yt(1:stair),zt(1:stair));
xlmin = -length(fp)*fL;
xlmax = length(fp)*fL;
ylmin = -length(fp)*fL;
ylmax = length(fp)*fL;
zlmin = -length(fp)*fL;
zlmax = length(fp)*fL;
axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
grid;
xlabel 'x';
ylabel 'y';
zlabel 'z';
title(strcat('3D_Simulation_of_Node_',num2str(length(fp)),'_Trial_',num2str(ft),'_rigid_polymer_dynamics'));axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);

disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));

while HTdist(fPnew, fL, fangle) > fa
    
    no =  randsample(2:length(fp)-1,1);
    
    Pno = randFlip(fPnew(no));
    Phypo1 = [fPnew(1:no-1), Pno(1), fPnew(no+1:end)];    % hypothetical change
    Phypo2 = [fPnew(1:no-1), Pno(2), fPnew(no+1:end)];    % hypothetical change
    
    Pch1 = pE(Phypo1, fHc, fHt)/(pE(Phypo1, fHc, fHt)+pE(Phypo2, fHc, fHt)+pE(fPnew, fHc, fHt));
    Pch2 = pE(Phypo2, fHc, fHt)/(pE(Phypo1, fHc, fHt)+pE(Phypo2, fHc, fHt)+pE(fPnew, fHc, fHt));
    rtt = rand();
    if rtt < Pch1
        fPnew = Phypo1;      % change state
    else
        if rtt < Pch1+Pch2
            fPnew = Phypo2;
        end
    end
    if HTdist(fPnew, fL, fangle) < fa
        break;
    end
    
    ffin = ffin+1;
    
    stair = length(fp)-1;
    
    fv = visV(buildV(fPnew, fL, fangle));
    xt = fv(1,:);
    yt = fv(2,:);
    zt = fv(3,:);
    
    plot3(xt(1:stair), yt(1:stair), zt(1:stair));
    grid;
    xlmin = -length(fp)*fL;
    xlmax = length(fp)*fL;
    ylmin = -length(fp)*fL;
    ylmax = length(fp)*fL;
    zlmin = -length(fp)*fL;
    zlmax = length(fp)*fL;
    axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
    xlabel 'x';
    ylabel 'y';
    zlabel 'z';
    title(strcat('3D_Simulation_of_Node_',num2str(length(fp)),'_Trial_',num2str(ft),'_rigid_polymer_dynamics'));axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
    xc = xlim;
    xl = xc(1)*0.2+xc(2)*0.8;
    yc = ylim;
    yl1 = yc(1)*0.14+yc(2)*0.86;
    yl2 = yc(1)*0.26+yc(2)*0.74;
    zc = zlim;
    zl = zc(1)*0.2+zc(2)*0.8;
    
    text(xl,yl1,zl, strcat('ffin=',num2str(ffin)),'Color','red','FontSize',12);
    text(xl,yl2,zl, strcat('HTdist=',num2str(HTdist(fPnew, fL, fangle))),'Color','red','FontSize',12);
    drawnow;
    
    %       pause(0.6)
    
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

function [fPnew, fHTd] = twistEquilRand(fpath, ftwist, ft, fp, fa, fL, fangle, fHc, fHt)
% To twist randomly till formed a loop

disp(strcat('T-',num2str(ft),'-------------'));

fPnew = fp;

stair = length(fp)-1;

fig = figure;
fv = visV(buildV(fPnew, fL, fangle));
xt = fv(1,:);
yt = fv(2,:);
zt = fv(3,:);
plot3(xt(1:stair),yt(1:stair),zt(1:stair));
xlmin = -length(fp)*fL;
xlmax = length(fp)*fL;
ylmin = -length(fp)*fL;
ylmax = length(fp)*fL;
zlmin = -length(fp)*fL;
zlmax = length(fp)*fL;
axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
grid;
xlabel 'x';
ylabel 'y';
zlabel 'z';
title(strcat('3D_Simulation_of_Node_',num2str(length(fp)),'_Trial_',num2str(ft),'_rigid_polymer_dynamics'));axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);

disp(strcat('Debug0: ', num2str(fPnew)));

for t = 1:ftwist
    
    no =  randsample(2:length(fp)-1,1);
    
    Pno = randFlip(fPnew(no));
    Phypo1 = [fPnew(1:no-1), Pno(1), fPnew(no+1:end)];    % hypothetical change
    Phypo2 = [fPnew(1:no-1), Pno(2), fPnew(no+1:end)];    % hypothetical change
    
    Pch1 = pE(Phypo1, fHc, fHt)/(pE(Phypo1, fHc, fHt)+pE(Phypo2, fHc, fHt)+pE(fPnew, fHc, fHt));
    Pch2 = pE(Phypo2, fHc, fHt)/(pE(Phypo1, fHc, fHt)+pE(Phypo2, fHc, fHt)+pE(fPnew, fHc, fHt));
    rtt = rand();
    if rtt < Pch1
        fPnew = Phypo1;      % change state
    else
        if rtt < Pch1+Pch2
            fPnew = Phypo2;
        end
    end
    if HTdist(fPnew, fL, fangle) < fa
        break;
    end
    
    stair = length(fp)-1;
    
    fv = visV(buildV(fPnew, fL, fangle));
    xt = fv(1,:);
    yt = fv(2,:);
    zt = fv(3,:);
    
    plot3(xt(1:stair), yt(1:stair), zt(1:stair));
    grid;
    xlmin = -length(fp)*fL;
    xlmax = length(fp)*fL;
    ylmin = -length(fp)*fL;
    ylmax = length(fp)*fL;
    zlmin = -length(fp)*fL;
    zlmax = length(fp)*fL;
    axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
    xlabel 'x';
    ylabel 'y';
    zlabel 'z';
    title(strcat('3D_Simulation_of_Node_',num2str(length(fp)),'_Trial_',num2str(ft),'_rigid_polymer_dynamics'));axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
    xc = xlim;
    xl = xc(1)*0.2+xc(2)*0.8;
    yc = ylim;
    yl1 = yc(1)*0.14+yc(2)*0.86;
    yl2 = yc(1)*0.26+yc(2)*0.74;
    zc = zlim;
    zl = zc(1)*0.2+zc(2)*0.8;
    
    text(xl,yl1,zl, strcat('t=',num2str(t)),'Color','red','FontSize',12);
    text(xl,yl2,zl, strcat('HTdist=',num2str(HTdist(fPnew, fL, fangle))),'Color','red','FontSize',12);
    drawnow;
    
    %       pause(0.6)
    
    disp(strcat('Debug ',num2str(t),': ', num2str(fPnew)));
    
end

fHTd = HTdist(fPnew, fL, fangle);

filename = strcat(fpath, 'N',num2str(length(fp)),'-T',num2str(ft),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
parsaveas(gcf, filename,'png');
close gcf;

disp(strcat('final state: ', num2str(fPnew)));
disp(strcat('finish dist: ', num2str(HTdist(fPnew, fL, fangle))));

end

function fPno = randFlip(fn)
% To randomly flip

if fn == 0
    fPno = [1, -1];
else
    fPno = [0, -fn];
end

end

function fE = pE(fp, fHc, fHt)
% To calculate energy of a certain polymer state

fE = 0;

% for no = 3:length(fp)-1
parfor no = 3:length(fp)-1
    if abs(fp(no)) == 1
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

ftrans = ones(1, length(fp)-2);
for count = 2:length(fp)-2
    ftrans(count) = -ftrans(count-1);
end
fm2D = build2DV([0, ftrans, 0], fL, fangle);

fx = fm2D(1,:);
fy = fm2D(2,:);
fz = zeros(1, length(fp)-1);
% length(fx)
% length(fy)
% length(fz)

fm = [fx;fy;fz];

for no = 2:length(fp)-1
    if fp(no) == 0
    else
        ax = fm(:,no-1)/norm(fm(:,no-1));
        u = ax(1);
        v = ax(2);
        w = ax(3);
        
        if fp(no) == 1
            ang = 2*pi/3; % CCW
        else
            ang = 4*pi/3; % CW
        end
        rot = [u^2+(1-u^2)*cos(ang), u*v*(1-cos(ang))-w*sin(ang), u*w*(1-cos(ang))+v*sin(ang);
            u*v*(1-cos(ang))+w*sin(ang), v^2+(1-v^2)*cos(ang), v*w*(1-cos(ang))-u*sin(ang);
            u*w*(1-cos(ang))-v*sin(ang), w*v*(1-cos(ang))+u*sin(ang), w^2+(1-w^2)*cos(ang)];
        fm(:,no:end) = rot*fm(:,no:end);
    end
end

end

function fm = build2DV(fp, fL, fangle)
% To build a 2D vector set (matrix) based on given polymer states

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

