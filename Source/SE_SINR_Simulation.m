%This Matlab script can be used to reproduce Figures 5.4(b), 5.5, 5.6(b), 5.7, and 5.8 in the monograph:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.1 (Last edited: 2021-11-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Empty workspace and close figures
close all
clear all
clc

%% Define simulation setup

Nsc=12;
N_symb=14;

%Number of Monte-Carlo setups
nbrOfSetups = 100;

%Number of channel realizations per setup
nbrOfRealizations = N_symb*1e3;

%Number of APs 
L = 100;

%Number of antennas per AP
N = 4;

%Number of UEs in the network
K = 20;

%Length of coherence block (12*14 = 168)
tau_c = Nsc*N_symb;     

%Length of pilot sequences
tau_p = 10;

%Compute the prelog factor assuming only uplink data transmission
prelogFactor = (1-tau_p/tau_c);


%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Prepare to save simulation results
SE_P_MMSE_Total = zeros(K,nbrOfSetups); %P-MMSE(DCC)
SE_small_MMSE_Total = zeros(K,nbrOfSetups); %Small-Cell MMSE
SINR_P_MMSE = zeros(K,nbrOfRealizations,nbrOfSetups);
SINR_small_MMSE = zeros(K,nbrOfRealizations,nbrOfSetups);

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p,1,0);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    % Elimino ASD_varphi y ASD_theta
    

    %% Cell-Free Massive MIMO with DCC
    
    %Compute SE using combiners and results in Section 5 for centralized
    %and distributed uplink operations for DCC
    [SE_P_MMSE, SE_small_MMSE, SINR_P_MMSE_INST, SINR_small_MMSE_INST] ...
          = functionComputeSE_uplink(Hhat,H,D,D_small,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);    
    
    %Save SE values
    SE_P_MMSE_Total(:,n) = SE_P_MMSE;
    SE_small_MMSE_Total(:,n) = SE_small_MMSE;
    SINR_P_MMSE(:,:,n)=SINR_P_MMSE_INST;
    SINR_small_MMSE(:,:,n)=SINR_small_MMSE_INST;
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end



Results_file=['Results_SE_SINR_K' num2str(K) '_L' num2str(L) '_N' num2str(N) '.mat'];
% save(fullfile('Results',Results_file),'SE_P_MMSE_Total','SE_small_MMSE_Total','SINR_P_MMSE','SINR_small_MMSE','nbrOfRealizations','nbrOfSetups');

nbrOfTTIs = nbrOfRealizations/N_symb;

SINR_P_MMSE_bloc = zeros(K,nbrOfTTIs,nbrOfSetups);
SINR_small_bloc = zeros(K,nbrOfTTIs,nbrOfSetups);

for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    for bloc=1:nbrOfTTIs
        SINR_P_MMSE_bloc(:,bloc,n) = mean(real(SINR_P_MMSE(:,1+(bloc-1)*N_symb:bloc*N_symb,n)),2);
        SINR_small_bloc(:,bloc,n) = mean(real(SINR_small_MMSE(:,1+(bloc-1)*N_symb:bloc*N_symb,n)),2);
    end

end

save(fullfile('Results',Results_file),'SE_P_MMSE_Total','SE_small_MMSE_Total','SINR_P_MMSE_bloc','SINR_small_bloc','nbrOfRealizations','nbrOfTTIs','nbrOfSetups');



% Plot Spectral Efficiency results
figure;
hold on; box on;
set(gca,'fontsize',16); % Tamaño y fuente de los ejes

plot(sort(SE_P_MMSE_Total(:)),linspace(0,1,K*nbrOfSetups),'b','LineWidth',2);
plot(sort(SE_small_MMSE_Total(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'CellFree mMIMO','Small-Cell'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);
