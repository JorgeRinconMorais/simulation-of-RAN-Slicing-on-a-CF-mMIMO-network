function [SE_P_MMSE, SE_small_MMSE, SINR_P_MMSE_INST, SINR_small_MMSE_INST] = functionComputeSE_uplink(Hhat, H, D, D_small, B, C, tau_c, tau_p, nbrOfRealizations, N, K, L, p, R, pilotIndex)

% La función functionComputeSE_uplink calcula la eficiencia espectral (SE) del enlace ascendente 
% utilizando esquemas de combinación P-MMSE y MMSE de small cells.

%Compute uplink SE for P-MMSE and small-cell MMSE combining schemes using the capacity
%bound in Theorem 5.1 for the centralized schemes and the capacity bound
%in Theorem 5.4 for the distributed schemes. Compute the genie-aided uplink
%SE from Corollary 5.9 for the centralized operation and from Corollary 5.10 
%for the distributed operation. 
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%D                 = DCC matrix for cell-free setup with dimension L x K 
%D_small           = DCC matrix for small-cell setup with dimension L x K
%B                 = Matrix with dimension N x N x L x K
%C                 = Matrix with dimension N x N x L x K
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs 
%L                 = Number of APs
%p                 = Uplink transmit power per UE
%R                 = Matrix with dimension N x N x L x K
%pilotIndex        = Vector containing the pilot assigned to each UE
%
%OUTPUT:
%SE_P_MMSE         = SEs achieved with P-MMSE combining in (5.16)
%SE_small_MMSE     = SEs achieved with L-MMSE combining for small-cell setup

%Compute the prelog factor assuming only uplink data transmission
prelogFactor = (1-tau_p/tau_c);
% Para ajustar la SE en función del tiempo dedicado a la transmisión de datos 
% en comparación con el tiempo total del bloque de coherencia.

%Prepare to store simulation results
SE_P_MMSE = zeros(K,1);
SE_small_MMSE = zeros(K,1);
SINR_P_MMSE_INST = zeros(K,nbrOfRealizations);
SINR_small_MMSE_INST = zeros(K,nbrOfRealizations);
% Matrices para almacenar resultados

%% Go through all channel realizations
for n = 1:nbrOfRealizations  % Bucle de 1 a nº de realizaciones

    %Consider the centralized schemes
    
    %Go through all UEs
    for k = 1:K  % Bucle sobre todos los usuarios
        
        % Identificación de los APs que sirven al UE k en configuraciones cell-free y small-cell:
        %Determine the set of serving APs for UE k
        servingAPs = find(D(:,k)==1); %cell-free setup
        servingAP_small = find(D_small(:,k)==1); %small-cell setup
        
        %Compute the number of APs that serve UE k
        La = length(servingAPs);
        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        % Inicialización de matrices para almacenar realizaciones de canal y matrices de correlación de errores de estimación:
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        % Extracción de realizaciones de canal y matrices de correlación de errores para los APs que sirven al UE k:
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
        %Compute P-MMSE combining according to (5.16)
        v = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p+C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        SINR_P_MMSE_INST(k,n) = numerator/denominator;

        %Compute the SE by computing the instantaneous SE for one
        %channel realization according to (5.4)
        SE_P_MMSE(k) = SE_P_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        %Extract the required channel estimates and estimation correlation
        %matrices for the SE computation of UE k in small-cell setup
        Hhatallj_active_small = reshape(Hhat((servingAP_small-1)*N+1:servingAP_small*N,n,:),[N K]);
        C_tot_blk_small = reshape(sum(C(:,:,servingAP_small,:),4),[N,N]);
        %Compute L-MMSE combining for small-cell setup according to (5.29)
        v_small = p*((p*(Hhatallj_active_small*Hhatallj_active_small')+p*C_tot_blk_small+eye(N))\Hhatallj_active_small(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)
        %for L-MMSE combining in small-cell setup
        numerator = p*abs(v_small'*Hhatallj_active_small(:,k))^2;
        denominator = p*norm(v_small'*Hhatallj_active_small)^2 + v_small'*(p*C_tot_blk_small+eye(N))*v_small - numerator;
        
        SINR_small_MMSE_INST(k,n) = numerator/denominator;

        %Compute the SE by computing the instantaneous SE for one
        %channel realization according to (5.4) in small-cell setup
        SE_small_MMSE(k) = SE_small_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
    end
    
end

end
