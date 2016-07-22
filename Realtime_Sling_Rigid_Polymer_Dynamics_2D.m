% Realtime_Sling_Rigid_Polymer_Dynamics_2D
% To simulate the dynamics of a rigid polymer by Sling model in 2D
% by Baihan Lin, Qian Lab
% July 2016

clear all;
close all;

global pathN;
pathN = '/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/';

%% Initialization

rng(2);                 % randomizer

trial = 10;              % trials
%twist = 200;            % change of state
node = 15;               % nodes of rigid polymer

angle = 0.05;           % in rad, angle changed in each twist
L = 1;                  % length of each segment of rigid polymer
a = 10;                 % threshold to form loop

%T = 300;                % temperature (K)
%k = 1.38064852e-23;     % Boltzmann constant (J/K)
%b = 1/(k*T);           % thermodynamic beta
%kT = 4.11e-21;          % at 25�C (298 K)
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

pF = zeros(1,trial); % record times to form a loop

%% Simulation of twisting

for n = 1:trial
    % To twist till looped
    [Pnew, fin] = twistLoopSeries(n, p, Pc, Pt, a, L, angle, Hc, Ht);
    pF(n) = fin;
end

fig = figure;
hist(pF);
title(strcat('Hist for ', num2str(length(fp)),'node,'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle)))
pathName ='/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/Hist-';
filename = strcat(pathName, num2str(length(fp)),'node,'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
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
for no = 2:fnode-1
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

function [fPnew, ffin] = twistLoopSeries(ft, fp, fPc, fPt, fa, fL, fangle, fHc, fHt)
% To twist node by node till formed a loop

disp(strcat('T-',num2str(ft),'-------------'));

fPnew = fp;
ffin = 1;

fig = figure;
fv = visV(buildV(fPnew, fL, fangle));
xt = fv(1,:);
yt = fv(2,:);
plot(xt(1:ffin), yt(1:ffin));
xlmin = min(xt);
xlmax = max(xt);
ylmin = min(yt);
ylmax = max(yt);
axis([ xlmin, xlmax, ylmin, ylmax]);
grid;
xlabel 'x';
ylabel 'y';
title(strcat('Simulation of ',num2str(length(fp)),' node rigid polymer dynamics'));
text(0,0,strcat('ffin',num2str(ffin)));

disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));

while HTdist(fPnew, fL, fangle) > fa
    
    for no = 2:length(fp)-1
        Phypo = [fPnew(1:no-1), -fPnew(no:end)];    % hypothetical change
        Pchg = pE(Phypo, fHc, fHt)/(pE(Phypo, fHc, fHt)+pE(fPnew, fHc, fHt));
        if rand() < Pchg
            fPnew(no:end) = -fPnew(no:end);      % change state
        end
        if HTdist(fPnew, fL, fangle) < fa
            break;
        end
        
        ffin = ffin+1;
        stair = ffin;
        if ffin > length(fp)-1
            stair = length(fp)-1;
        end
        
        fv = visV(buildV(fPnew, fL, fangle));
        xt = fv(1,:);
        yt = fv(2,:);
        
        plot(xt(1:stair), yt(1:stair));
        title(strcat('Simulation of ',num2str(length(fp)),' node rigid polymer dynamics'));
        grid;
        xlmin = min(min(xt),xlmin);
        xlmax = max(max(xt),xlmax);
        ylmin = min(min(yt),ylmin);
        ylmax = max(max(yt),ylmax);
        axis([ xlmin, xlmax, ylmin, ylmax]);
        xc = xlim;
        xl = xc(1)*0.2+xc(2)*0.8;
        yc = ylim;
        yl = yc(1)*0.2+yc(2)*0.8;
        text(xl,yl,strcat('ffin=',num2str(ffin)),'Color','red','FontSize',12);
        drawnow;
        %pause(0.002)
        
        disp(strcat('Debug ',num2str(ffin),': ', num2str(fPnew)));
    end
    
end

path ='/Users/DoerLBH/Dropbox/git/QianLab_DNA_dynamics/';
filename = strcat(path, num2str(length(fp)),'node-T',num2str(ft),'-a',num2str(fa),'-l',num2str(fL),'-r',num2str(fangle),'.png');
saveas(gcf, filename,'png');
close gcf;

disp(strcat('final state: ', num2str(fPnew)));
disp(strcat('finish time: ', num2str(ffin)));
disp(strcat('finish dist: ', num2str(HTdist(fPnew, fL, fangle))));

end


% function [fPnew, fE] = twistPoly(fp, fPc, fPt, fa, fL, fangle)
% % To twist at specific twist number
%
% fE = pE(fPnew, fHc, fHt);
%
% end

function fE = pE(fp, fHc, fHt)
% To calculate energy of a certain polymer state

fE = 0;
fs = diag(fp(2:length(fp)-2).'*fp(3:length(fp)-1));

for no = 1:length(fp)-3
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
rot1 = [cos(fangle) -sin(fangle); sin(fangle) cos(fangle)];
rot2 = [cos(-fangle) -sin(-fangle); sin(-fangle) cos(-fangle)];

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

