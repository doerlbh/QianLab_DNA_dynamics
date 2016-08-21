% Realtime_Ising_Rigid_Protein_Dynamics_3D_using_parallel_Rand_Walk_Dist_Auto
% To simulate the dynamics of a rigid protein by Ising model in 3D plus
% autocorrelation
%
% by Baihan Lin, Baker Lab
% August 2016

clear all;
close all;

%% Initialization

rng(1111);                 % randomizer

trial = 2000;             % trials
twist = 1000;             % change of set state changes
% node = 500;               % nodes of rigid polymer
AutoT = 2000;             % autocorrelation run time

global pathN;
% pathN = '/Users/sunnylinL/Dropbox/Sim1/data/data_20160817/';
% pathN = '/gscratch/stf/sunnylin/other/sim/data/';
pathN = '/gscratch/stf/sunnylin/other/sim/data/data_20160821/';

system(['mkdir ' pathN]);

angle = 0.02;            % in rad, angle changed in each twist
L = 1;                  % length of each segment of rigid polymer
a = 20;                 % threshold to form loop

Hc = 1.0;   % in unit of kT, energy level of cis rigid configuration
Ht = 0.9;   % in unit of kT, energy level of trans rigid configuration

% parfor it = 1:5
it = 1;
% for it = 1:5
    
    node = it*100;
    
    %% Main functions for Equilibrium
    
    [finPr, stPr, pEfr, pDfr] = EquilSimulationRandpoly(node, trial, twist, pathN, a, L, angle, Hc, Ht);
    
    %% Autocorrelation for Equilibrium
    
    fnameer = strcat('Equil-Auto-N',num2str(node),'-A',num2str(AutoT),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle));
    [RACorr] = AutocorEq(finPr, trial, AutoT, pathN, a, L, angle, Hc, Ht, fnameer);
    
    fnameer = strcat('Equil-Tau-N',num2str(node),'-A',num2str(AutoT),'-a',num2str(a),'-l',num2str(L),'-r',num2str(angle)); 
    rtaur = Cor2tau(RACorr, pathN, fnameer);
    
    %% Main functions for Looping
    
%     [lfinPr, lstPr, lpTfr, lpEfr, lpDfr] = loopSimulationRandpoly(node, finPr, trial, pathN, a, L, angle, Hc, Ht);
    
% end
