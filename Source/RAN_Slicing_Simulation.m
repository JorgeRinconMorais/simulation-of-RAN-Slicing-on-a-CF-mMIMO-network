%%%%%%%%%%%%%%% RAN Slicing %%%%%%%%%%%%%%%

close all
clear all
% clc



deployment='Cell_Free';
% deployment='Small_Cell';

% category = 'eMBB';
category = 'uRLLC';



%% Define simulation setup

%Number of APs 
L = 100;

%Number of antennas per AP
N = 1;

%Number of UEs in the network
K = 20;

%Length of coherence block
tau_c = 168; %12*14

%Length of pilot sequences
tau_p = 10;

%Compute the prelog factor assuming only uplink data transmission
prelogFactor = (1-tau_p/tau_c);



% MCS:
Qm = zeros(28,1);
Qm(1:5)=2;      % QPSK
Qm(6:11)=4;      % 16QAM
Qm(12:20)=6;      % 64QAM
Qm(21:28)=8;      % 256QAM

Kcod= [120 193 308 449 602 378 434 490 553 616 658 466 517 567 616 666 719 772 822 873 682.5 711 754 797 841 885 916.5 948]';
Rcod = Kcod./1024;
gamma = 2.^(Qm.*Rcod) - 1;

TBS = zeros(28,1);
switch category

    case 'eMBB'
    TBS = floor(prelogFactor.*tau_c.*Qm.*Rcod);

    case 'uRLLC'   
    TBS_aux = floor(prelogFactor.*tau_c.*Qm.*Rcod);
    TBS(1:4) = [8 12 18 29];
    TBS(5:end) = TBS_aux(1:end-4);

end
% [[0:27]' 10*log10(gamma) Qm Rcod TBS]



Results_SINR_file=['Results_SE_SINR_K' num2str(K) '_L' num2str(L) '_N' num2str(N) '.mat'];
load(fullfile('Results',Results_SINR_file));

[[0:27]' 10*log10(gamma) Qm Rcod TBS]


%%%%%%%%%%%%%%% Preparation %%%%%%%%%%%%%%%

% -	RT Printing Machine.
    
PayLoad = 2000;     % 250 bytes -> Message size from 40 to 250 bytes [bits].
Tp = 2*10^-3;       % Transmission period [s].
D = 1*10^-3;        % 1 ms -> Latencia [s].
TTI = D;            % TTI = 1 ms (TTI <= D)[s].
u = 0;
BLER=1e-4;

%%%%%%%%%%%%%%% End Preparation %%%%%%%%%%%%%%%



SINR = zeros(K,nbrOfTTIs);

N_PRB_Prov_mean = zeros(nbrOfSetups,1);
N_PRB_Prov_5p = zeros(nbrOfSetups,1);

col_mean = zeros(nbrOfSetups,1);
col_5p = zeros(nbrOfSetups,1);

Pc_simul_mean = zeros(nbrOfSetups,1);
Pc_simul_5p = zeros(nbrOfSetups,1);

for n = 1:nbrOfSetups
    
    % Display simulation progress
    % disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    switch deployment

        case 'Cell_Free'
        SINR = SINR_P_MMSE_bloc(:,:,n);

        case 'Small_Cell'   
        SINR = SINR_small_bloc(:,:,n);

    end
        
    %%%%%%%%%%%%%%% Provitioning - Mean %%%%%%%%%%%%%%%
    
    SINR_mean = mean(SINR,2);
    N_PRB_mean = zeros(K,1);

    for k=1:K

        I_MCS = find(SINR_mean(k)>=gamma, 1, 'last' );
        if isempty(I_MCS)
            I_MCS=1;
        end

        I_MCS = I_MCS-1;    % Better Provisioning
        if I_MCS==0
            I_MCS = 1;
        end

        N_PRB_mean(k) = ceil(PayLoad./TBS(I_MCS));

    end
    
    N_PRB_Prov_mean(n) = max(N_PRB_mean);
    
    
    
    %%%%%%%%%%%%%%% Provitioning - 5 Percentile %%%%%%%%%%%%%%%
    
    SINR_5p = zeros(K,1);
    N_PRB_5p = zeros(K,1);

    for k=1:K

        [Fx,x]=ecdf(SINR(k,:));
        for il=1:length(Fx)-1
    	    if Fx(il)==Fx(il+1)
    	    Fx(il)=Fx(il)-0.01/il; 
     	    end
            if x(il)==x(il+1)
    	    x(il)=x(il)-0.01/il; 
            end
        end
    
        SINR_5p(k) = interp1(Fx,x,5/100);
        
        I_MCS = find(SINR_5p(k)>=gamma, 1, 'last' );
        if isempty(I_MCS)
            I_MCS=1;
        end

        I_MCS = I_MCS-1;    % Better Provisioning
        if I_MCS==0
            I_MCS = 1;
        end

        N_PRB_5p(k) = ceil(PayLoad./TBS(I_MCS));
    end
    
    N_PRB_Prov_5p(n) = max(N_PRB_5p);


    
    %%%%%%%%%%%%%%% System Level Simulation %%%%%%%%%%%%%%%
    
    SINR_Simul = repmat(SINR,1,1,length(gamma));

    gamma_Simul  = zeros(K,nbrOfTTIs,length(gamma));
    for g=1:length(gamma)
        gamma_Simul(:,:,g)=gamma(g);
    end

    [~,p]=min(abs(gamma_Simul-SINR_Simul),[],3);
    lgc = gamma(p)>SINR;
    I_MCS_aux = p-lgc;
    I_MCS = I_MCS_aux + (I_MCS_aux==0);
    N_PRB_inst = ceil(PayLoad./TBS(I_MCS));
    
    
    
    for k=1:K
        col_mean(n) = col_mean(n) + sum(N_PRB_inst(k,:)>N_PRB_Prov_mean(n));  % Mean

        col_5p(n) = col_5p(n) + sum(N_PRB_inst(k,:)>N_PRB_Prov_5p(n));  % 5 Percentile
    end
        
    % Pc for Mean
    Pc_simul_mean(n) = col_mean(n)./(K*nbrOfTTIs);
    
    % Pc for 5 Percentile
    Pc_simul_5p(n) = col_5p(n)./(K*nbrOfTTIs);

end



Results_file=['Results_RAN_Slicing_' category '_' deployment '_K' num2str(K) '_L' num2str(L) '_N' num2str(N) '.mat'];
save(fullfile('Results',Results_file),'N_PRB_Prov_mean','N_PRB_Prov_5p','col_mean','col_5p','Pc_simul_mean','Pc_simul_5p','nbrOfTTIs','nbrOfSetups');



[sum(Pc_simul_mean<1e-2) sum(Pc_simul_5p<1e-2) sum(Pc_simul_mean<1e-3) sum(Pc_simul_5p<1e-3)]
[N_PRB_Prov_mean N_PRB_Prov_5p col_mean col_5p Pc_simul_mean<1e-2 Pc_simul_5p<1e-2 Pc_simul_mean<1e-3 Pc_simul_5p<1e-3]

format long g;
[[0:27]' 10*log10(gamma) Qm Rcod TBS]