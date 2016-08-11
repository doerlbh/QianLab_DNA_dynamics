% Rand_Walk_3D
% Modified by Baihan Lin, Qian Lab, July 2016
% Originally developed by Felix Ye, Qian Lab, ?? ??

%% This is based on Flory's book 
clear all
close all
clc
 
T       =@(theta, phi) [cos(theta),          sin(theta),            0; 
                        sin(theta)*cos(phi), -cos(theta)*cos(phi), sin(phi);
                        sin(theta)*sin(phi), -cos(theta)*sin(phi), -cos(phi)];
theta   = 0.1;
l       = 1;
L_p     = l/(1-cos(theta));  %the Kuhn length is b=2L_p=400l
M       = 10^5;
% for case 2
theta_k = 0.5; 
gamma   = cos(theta_k)/cos(theta);
 
% for case 3 %This is the probability 
sigma   = 0.9;
s(1)    = 1/(1+2*sigma);     %phi=0
s(2)    = sigma/(1+2*sigma); %phi=2pi/3
s(3)    = sigma/(1+2*sigma); %phi=-2pi/3
f       = @(x) 2/3*pi*heaviside(x-s(1))-4/3*pi*heaviside(x-s(1)-s(2)); % sampling function
 
rng default
 
%% Case 1: L_c\approx  L_p, homogenous
% k_select = 2:2:20; %L_c= 1:10 Kuhn length 
% 
% parfor k_n = 1:10
%     k = k_select(k_n);
%     N = k*round(L_p); %N=200*k
%     R_e = zeros(3,M);
    
   % auto_correlation=zeros(1,N);
    
%    for j=1:M
 
%         T_select=cell(1,N);
%         T_select{1}=T(theta, rand*2*pi);
%         
%         for i=1:N-1
%             T_select{i+1}=T_select{i}*T(theta, rand*2*pi);
%         end
%         
%         R = zeros(3,N);
%         R(:,1) = [1;0;0];
%         for i=1:N-1
%             R(:,i+1)=T_select{i}*R(:,1);
%         end
%         
%         R_e(:,j)=sum(R,2);
%       %  auto_correlation=auto_correlation+R(:,1)'*R(:,:);
%     end
%    % auto_correlation=auto_correlation./M;
%     parsave(sprintf('output%d.mat', k),  R_e);
%   
% 
% % 
%  end
 
%% Case 2: L_c\approx L_p, there is a kink in the middle
k_select = 2:2:20;
 
parfor  k_n = 1:10
     k=k_select(k_n);
    N   = k*round(L_p); %N=200
    R_e = zeros(3,M);
    
  %  auto_correlation = zeros(1,N);
    for j=1:M
        T_select = cell(1,N);      
        T_select{1} = T(theta, rand*2*pi);
        for i=1:N-1
            if i~=round(N/2)
                T_select{i+1} = T_select{i}*T(theta, rand*2*pi);
            else
                T_select{i+1} = T_select{i}*T(theta_k, rand*2*pi); %kink at the middle
            end
        end
        
        R = zeros(3,N);
        R(:,1) = [1;0;0];
        for i = 1:N-1
            R(:,i+1) = T_select{i}*R(:,1);
        end
        R_e(:,j) = sum(R,2);
       % auto_correlation = auto_correlation+R(:,1)'*R(:,:);
    end
    %auto_correlation = auto_correlation./M;
    parsave(sprintf('output_kink%d_gamma%d.mat', k,gamma), R_e);
  
 
% 
 end
 
%% Case 3: Reversible kink, phi=0, 2pi/3, -2pi/3
 
 
 
%       
% 
% parfor k = 1:7
%      N   = k*round(L_p); %N=200
%      R_e = zeros(3,M);
%     
%     for j = 1:M
%         T_select = cell(1,N);
%         T_select{1} = T(theta, f(rand));
%         for i = 1:N-1
%             T_select{i+1} = T_select{i}*T(theta, f(rand));
%         end
%         
%         R = zeros(3,N);       
%         R(:,1) = [1;0;0];
%         for i = 1:N-1
%             R(:,i+1) = T_select{i}*R(:,1);
%         end
%         R_e(:,j) = sum(R,2);
%         
%     end
%       parsave(sprintf('output_reversible_kink%d_sigma%d.mat', k,sigma), R_e);
%     
% end
