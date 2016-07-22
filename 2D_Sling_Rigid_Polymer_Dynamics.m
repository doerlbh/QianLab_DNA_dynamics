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
node = 8;               % nodes of rigid polymer

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
% from node 2 to node n-1, they are either 1 (cis) or -1 (trans)

p = createRandPolymer(node); % randomly generate
%p = createPolymer(node,Pc); % naturally generate

%% Construct a recorder

pF = zeros(1,trial); % record times to form a loop

%% Simulation of twisting

for n = 1:trial
    
    % To twist till looped
    [Pnew, fin] = twistLoopseries(p, node, Pc, Pt, a, L, angle)
    
end

%% Local functions

function fp = createPolymer(fnode,fPc)
% To create a natural polymer with fnode nodes based on probability

for no = 2:fnode-1
    if rand() < fPc
        fp(no) = 1;      % cis
    else
        fp(no) = -1;     % trans
    end
end

end

function fp = createRandPolymer(fnode)
% To create a random polymer with fnode nodes

for no = 2:fnode-1
    if rand() > 0.5
        fp(no) = 1;      % cis
    else
        fp(no) = -1;     % trans
    end
end

end

function [fPnew, ffin] = twistLoopseries(fp, fnode, fPc, fPt, fa, fL, fangle)
% To twist node by node till formed a loop

fPnew = fp;
ffin = 0;
while headtaildist(fPnew) < a {
    for no = 2:fnode-1
        Phypo = [fPnew(1:no-1), -fPnew(no:end)];
        if rand() < pE(Phypo)/(pE(Phypo)+pE(fPnew))
            fPnew(no:end) = -fPnew(no:end);      % change state
        end
    end
    ffin = ffin+1;
end
end

end

function [fPnew, fE] = twistPoly(fp, fnode, fPc, fPt, fa, fL, fangle)
% To twist at specific twist number
fE = pE(fPnew)

end

function fE = pE(fp)
% To calculate energy of a certain polymer state

fE = pE(fp, Hc, Ht)

end

function fE = pE(fp, Hc, Ht)
% To calculate energy of a certain polymer state

fE = 0
for no = 2:fnode-1
    if fp(no) = 1
        fE = fE + Hc;      % cis
    else
        fE = fE + Ht;     % trans
    end
end

end




