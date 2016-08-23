% Realtime_Ising_Rigid_Protein_Dynamics_3D_using_parallel_Rand_Walk_Dist_Auto
% To simulate the dynamics of a rigid protein by Ising model in 3D plus
% autocorrelation
%
% by Baihan Lin, Baker Lab
% August 2016

clear all;
close all;

%% Initialization

rng(110);                 % randomizer

trial = 1000;             % trials
twist = 2000;             % change of set state changes
% node = 500;               % nodes of rigid polymer
AutoT = 1000;             % autocorrelation run time

global pathN;
pathN = '/Users/sunnylinL/Dropbox/Sim1/data/Loop2_20160822/';
% pathN = '/gscratch/stf/sunnylin/other/sim/data/';

system(['mkdir ' pathN]);

angle = 0.05;            % in rad, angle changed in each twist
L = 1;                  % length of each segment of rigid polymer
a = 20;                 % threshold to form loop

Hc = 1.0;   % in unit of kT, energy level of cis rigid configuration
Ht = 0.9;   % in unit of kT, energy level of trans rigid configuration

% parfor it = 1:5
% it = 1;
for it = 1:5
    
    node = it*200;
    
    %% Main functions for Equilibrium
    
%     [finPr, stPr, pEfr, pDfr] = EquilSimulationRandpoly(node, trial, twist, pathN, a, L, angle, Hc, Ht);
    
    [NfinPr, NstPr, NpEfr, NpDfr] = EquilNatpoly(node, trial, twist, pathN, a, L, angle, Hc, Ht);
        
    [finP,stP,pTf,pEf,pDf] = loopAutorun(node,NfinPr,trial,pathN,a,L,angle,Hc,Ht);
    
end
