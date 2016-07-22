% Sling_Rigid_Polymer_Dynamics_2D
% To simulate the dynamics of a rigid polymer by Sling model in 2D
% by Baihan Lin, Qian Lab
% July 2016

%% Initialization

trial = 1;              % training
twist = 200;            % change of state
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

p = createRandPolymer(node); % randomly generate
%p = createPolymer(node,Pc,Pt); % naturally generate

%% Construct a recorder

pF = []




%% Simulation of twisting

for n = 1:trial
    
    
    % To twist till looped
    [Pnew, fin] = twistLoop(p, node, Pc, Pt, a, L, angle)
    
end


%% Local functions

function fp = createPolymer(fnode,fPc,fPt)
% To create a natural polymer with fnode nodes based on probability

end

function fp = createRandPolymer(fnode)
% To create a random polymer with fnode nodes

end

function [fPnew, ffin] = twistLoop(fp, fnode, fPc, fPt, fa, fL, fangle)
% To twist till formed a loop


end

function [fPnew, fE] = twistPoly(fp, fnode, fPc, fPt, fa, fL, fangle)
% To twist at specific twist number

    fE = polyenergy(fPnew)
    
end

