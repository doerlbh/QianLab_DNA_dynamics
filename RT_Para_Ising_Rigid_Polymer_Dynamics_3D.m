% Realtime_Ising_Rigid_Polymer_Dynamics_3D_using_parallel
% To simulate the dynamics of a rigid polymer by Sling model in 3D
% by Baihan Lin, Qian Lab
% July 2016

clear all;
close all;

%% Initialization

rng(378);                 % randomizer

trial = 500;              % trials
%twist = 200;            % change of state
node = 100;               % nodes of rigid polymer

global pathN;
pathN = strcat('/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/data/3D-',num2str(node),'/');
system(['mkdir ' pathN]);

angle = 0.5;           % in rad, angle changed in each twist
L = 1;                  % length of each segment of rigid polymer
a = 10;                 % threshold to form loop

Hc = 1.0;   % in unit of kT, energy level of cis rigid configuration
Ht = 0.9;   % in unit of kT, energy level of trans rigid configuration
b = 1;      % redefined beta based on H

Pc = exp(-b*Hc)/(exp(-b*Hc)+exp(-b*Ht));       % Probablity of cis change
Pc2 = exp(-b*Hc)/(2*exp(-b*Hc)+exp(-b*Ht));
Pt = exp(-b*Ht)/(exp(-b*Hc)+exp(-b*Ht));       % Probablity of trans change
Pt2 = exp(-b*Ht)/(2*exp(-b*Hc)+exp(-b*Ht));

%% Generate a polymer

% initial direction = positive x
% default change = counterclockwise angle
% in this vector, first and last node are zero
% from node 2 to node n-1, they are either 1 (CW) or -1 (CCW)

p = createRandPolymer(node); % randomly generate

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
histogram(pTf(2:trial+1), 'BinWidth', 50);
title(strcat('Time Histogram for N', num2str(node),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)))
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
histogram(pEf(2:trial+1), 'BinWidth', 0.1);
% line([pEf(1) pEf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('Energy Histogram for N', num2str(node),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)));
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
histogram(pDf(2:trial+1), 'BinWidth', 0.1);
% line([pDf(1),pDf(1)],get(axes,'YLim'),'Color',[1 0 0],'LineWidth',3);
title(strcat('HTdistance Histogram for N', num2str(node),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)))
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

function fp = createRandPolymer(fnode)
% To create a random polymer with node nodes

fp = zeros(1,fnode);
fp(2) = 0;

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
zt = fv(3,:);
plot3(xt(1:stair),yt(1:stair),zt(1:stair));
xlmin = min(xt);
xlmax = max(xt);
ylmin = min(yt);
ylmax = max(yt);
zlmin = min(zt);
zlmax = max(zt);
axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
grid;
xlabel 'x';
ylabel 'y';
zlabel 'z';
title(strcat('3D Simulation of ',num2str(length(fp)),' node rigid polymer dynamics'));

disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));

while HTdist(fPnew, fL, fangle) > fa
    
    no =  randsample(2:length(fp)-1,1);
    
    Pno = randFlip(fPnew(no));
    Phypo = [fPnew(1:no-1), Pno, fPnew(no+1:end)];    % hypothetical change
    
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
    zt = fv(3,:);
    
    plot3(xt(1:stair), yt(1:stair), zt(1:stair));
    grid;
    xlmin = min(min(xt),xlmin);
    xlmax = max(max(xt),xlmax);
    ylmin = min(min(yt),ylmin);
    ylmax = max(max(yt),ylmax);
    zlmin = min(min(zt),zlmin);
    zlmax = max(max(zt),zlmax);
    axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
    xlabel 'x';
    ylabel 'y';
    zlabel 'z';
    title(strcat('3D Simulation of ',num2str(length(fp)),' node rigid polymer dynamics'));axis([ xlmin, xlmax, ylmin, ylmax, zlmin, zlmax]);
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

function fPno = randFlip(fn)
% To randomly flip

r = rand();
if fn == 0
    if r < 0.5
        fPno = 1;
    else
        fPno = -1;
    end
else
    if r < 0.5
        fPno = 0;
    else
        fPno = -fn;
    end
end

end

function fE = pE(fp, fHc, fHt)
% To calculate energy of a certain polymer state

fE = 0;

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
        ax = fm(:,no-1);
        u = ax(1);
        v = ax(2);
        w = ax(3);
        
        if fp(no) == 1
            ang = pi/3; % CCW
        else
            ang = 2*pi/3; % CW
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

