imulation for ch estimation in LTE 
% clc 
clear all; 
close all; 
%%system parameter 
Trans_Antenna_Num = 1; 
Recei_Antenna_Num = 1; 
CellNum = 3; 
Nfft = 2048; 
NusedSc = 12*111; 
Ts = 1/(15000*2048); 
Ncp_normal = round(4.6875e-6/Ts); 
Ncp_extend = round(16.67e-6/Ts); 
N_RB_min = 6; 
N_RB_max = 111; 
N_RBsc = 12; 
Nsym_normal = 7; 
Nsym_extend = 6; 
Nslot_perFrame = 20; 
%%User parameter 
Current_cell = 1; 
CP_idx = 0;    % 0 for normal CP, and 1 for extend CP 
Num_UsedSlot = 2; 
Num_UsedRB = 10; 
UsedRB_idx = 1:10;%sort(randint(1,Num_UsedRB,111)+1); 
if CP_idx == 0 
    Ncp = Ncp_normal; Nsym_slot = Nsym_normal; 
elseif CP_idx == 1 
    Ncp = Ncp_extend; Nsym_slot = Nsym_extend; 
end 
 
PowerList = [10 0 0];   %the power from different Bs 
CH_profile{1} = [1 1/6 1/6];   %channel delay for Bs_1 
CH_profile{2} = [1 1/6 1/6];  %channel delay for Bs_2 
CH_profile{3} = [1 1/6 1/7];  %channel delay for Bs_3 
 
%%%main function %%%%%% 
for cellIdx = 1:CellNum 
    Tr_F_data = zeros(Nfft,Nsym_slot*Num_UsedSlot); 
    UsedF_data = zeros(N_RBsc*N_RB_max,Nsym_slot*Num_UsedSlot); 
    for RB_idx = UsedRB_idx 
        source{cellIdx,RB_idx} = qammod(randint(N_RBsc,Nsym_slot*Num_UsedSlot,4),4)/sqrt(2); 
        UsedF_data((RB_idx-1)*N_RBsc+1:RB_idx*N_RBsc,:) = source{cellIdx,RB_idx};       
    end 
    Tr_F_data(1024-12*111/2+1:1024+12*111/2,:) = UsedF_data;          %mapping the frequency domain data 
    Tr_T_data{cellIdx} = ifft(ifftshift(Tr_F_data))*sqrt(2048);               
    Tr_T_data{cellIdx} = [Tr_T_data{cellIdx}(end-Ncp+1:end,:);Tr_T_data{cellIdx}];         %add cp 
    Tr_T_data{cellIdx} = reshape(Tr_T_data{cellIdx} ,1,[]); 
    Tr_T_data{cellIdx} = Tr_T_data{cellIdx}.*10^(PowerList(cellIdx)/20);    %trans power 
    Re_T_data{cellIdx} = filter(CH_profile{cellIdx},1,Tr_T_data{cellIdx});    %pass ch; 
end 
%%receiver 
Re_T_signal = Re_T_data{1}+Re_T_data{2}+Re_T_data{3};   %combin the signal from different Bs 
     
%%%ideal ch estimation%%% 
Rec_Signal = Re_T_data{Current_cell}; 
Rec_Signal_tmp = reshape(Rec_Signal,Nfft+Ncp,Nsym_slot*Num_UsedSlot); 
Rec_Signal_tmp(1:Ncp,:)=[]; 
Rec_F_Signal = fftshift(fft(Rec_Signal_tmp))/sqrt(Nfft); 
Rec_Used_Signal = Rec_F_Signal(1024-12*111/2+1:1024+12*111/2,:); 
for RB_idx = UsedRB_idx 
    Ideal_H{RB_idx} = Rec_Used_Signal((RB_idx-1)*N_RBsc+1:RB_idx*N_RBsc,:)./source{Current_cell,RB_idx};     
    DemodSignal =  Rec_Used_Signal((RB_idx-1)*N_RBsc+1:RB_idx*N_RBsc,:)./Ideal_H{RB_idx}; 
%     figure,plot(DemodSignal,'ro');grid on;title('Ideal H; domodulation signal'); 
%     pause; 
    %%mse 
end 
%%%%%%%%%%%%%%%% 
% Noise = sqrt(12*10/2048)*sqrt(0.5)*(randn(size(Re_T_signal))+sqrt(-1)*randn(size(Re_T_signal))); 
% Re_T_signal = Re_T_signal+Noise; 
SIR =  10*log10(mean(abs(Re_T_data{1}).^2)/(mean(abs(Re_T_data{2}+Re_T_data{3}).^2))), 
%%Estch with the combined signal 
Rec_Signal = Re_T_signal; 
Rec_Signal_tmp = reshape(Rec_Signal,Nfft+Ncp,Nsym_slot*Num_UsedSlot); 
Rec_Signal_tmp(1:Ncp,:)=[]; 
Rec_F_Signal = fftshift(fft(Rec_Signal_tmp))/sqrt(Nfft); 
Rec_Used_Signal = Rec_F_Signal(1024-12*111/2+1:1024+12*111/2,:); 
 
for RB_idx = UsedRB_idx 
    %ch estimation 
    Pilot_f_idx = [06,12,03,09,06,12,03,09]; 
    Pilot_t_idx = [01,01,05,05,08,08,12,12]; 
    Rec_RB_Data = Rec_Used_Signal((RB_idx-1)*N_RBsc+1:RB_idx*N_RBsc,:); 
    Tra_RB_Data = source{Current_cell,RB_idx}; 
    EstH = zeros(size(Tra_RB_Data)); 
    for rs_idx = 1:length(Pilot_f_idx)       
        EstH(Pilot_f_idx(rs_idx),Pilot_t_idx(rs_idx)) = Rec_RB_Data(Pilot_f_idx(rs_idx),Pilot_t_idx(rs_idx))./Tra_RB_Data(Pilot_f_idx(rs_idx),Pilot_t_idx(rs_idx)); 
    end 
    %%interpolation 
    T_delta01 = (EstH([3,9],12)-EstH([3,9],5))/7; 
    T_delta02 = (EstH([6,12],8)-EstH([6,12],1))/7; 
    EstH([3 9],:) = kron(EstH([3,9],5),ones(1,14))+diag(T_delta01)*kron([-4:9],ones(2,1)); 
    EstH([6 12],:) = kron(EstH([6 12],1),ones(1,14))+diag(T_delta02)*kron([-4:9],ones(2,1)); 
    Est_H{RB_idx} = [EstH(3,:);EstH(3,:);interp1([3:3:12],EstH(3:3:12,:),[3:1:12])]; 
    DemodSignal =Rec_RB_Data./Est_H{RB_idx}; 
%     figure,plot(DemodSignal,'bo');grid on;title('Est H; domodulation signal'); 
%     pause; 
    %%mse 
    MeanSqrt = sqrt(mean(mean(abs(Ideal_H{RB_idx}).^2))); 
    Mse(RB_idx) = (sum(sum((abs((Ideal_H{RB_idx}- Est_H{RB_idx})./MeanSqrt).^2))))/(N_RBsc*Nsym_slot*Num_UsedSlot); 
end 
Mean_Mse = mean(Mse(UsedRB_idx)), 
 
% keyboard; 
