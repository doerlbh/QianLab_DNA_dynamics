% Sling_Rigid_Polymer_Dynamics_2D
% To simulate the dynamics of a rigid polymer by Sling model in 2D
% by Baihan Lin, Qian Lab
% July 2016

clear all;
close all;

%% Initialization

rng(1);                 % randomizer

trial = 1;              % trials
%twist = 200;            % change of state
global node = 3;               % nodes of rigid polymer

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
    [Pnew, fin] = twistLoopSeries(p, node, Pc, Pt, a, L, angle)
    
end

%% Local functions

function fp = createNatPolymer(node,fPc)
% To create a natural polymer with node nodes based on probability

if rand() > 0.5
    fp(2:node-1) = 1;      % clockwise
else
    fp(2:node-1) = -1;     % counterclockwise
end

for no = 3:node-1
    if rand() > fPc
        fp(no:node-1) = -fp(no:node-1);      % trans
    end
end

end

function fp = createRandPolymer(node)
% To create a random polymer with node nodes

for no = 2:node-1
    if rand() > 0.5
        fp(no) = 1;      % clockwise
    else
        fp(no) = -1;     % counterclockwise
    end
end

end

function [fPnew, ffin] = twistLoopSeries(fp, node, fPc, fPt, fa, fL, fangle)
% To twist node by node till formed a loop

fPnew = fp;
ffin = 0;
while HTdist(fPnew, fL, fangle) < a
    for no = 2:node-1
        Phypo = [fPnew(1:no-1), -fPnew(no:end)];    % hypothetical change
        if rand() < pE(Phypo)/(pE(Phypo)+pE(fPnew))
            fPnew(no:end) = -fPnew(no:end);      % change state
        end
    end
    ffin = ffin+1;
end

end

function [fPnew, fE] = twistPoly(fp, node, fPc, fPt, fa, fL, fangle)
% To twist at specific twist number

fE = pE(fPnew);

end

function fE = pE(fp)
% To calculate energy of a certain polymer state

fE = pEfull(fp, Hc, Ht);

end

function fE = pEfull(fp, Hc, Ht)
% To calculate energy of a certain polymer state

fE = 0;
fs = diag(fp(2:node-2).'*fp(3:node-1));
for no = 1:node-2
    if fs(no) == 1
        fE = fE + Hc;      % cis
    else
        fE = fE + Ht;     % trans
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

for no = 2:node-1
    if fp(no) == 1
        fm(:,no:end) = rot1*fm(:,no:end);       % clockwise
    else
        fm(:,no:end) = rot2*fm(:,no:end);       % counterclockwise
    end
end

end
