% clc;
% clear;
% 
% s=9;
% N_ofdm=10;
% 
% N_carrier= 100;
% nsamp=20;
% N_fft=nsamp*N_carrier;
% 
% M=4;
% 
% m=0.25;
% length_cp=round(m*N_fft);
% 
% ebn0_db = [0:1:10];
% ebn0=10.^(0.1*ebn0_db);
% esn0=log2(M)*ebn0;
% esn0_db=10*log10(esn0);
% snr_db=esn0_db-10*log10(nsamp);
% 
% ber_thertical = [];
% ber=[];
%      
%     rng(s);
%     serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
% %      serial_bit=zeros(1,N_carrier*N_ofdm*log2(M));
%     serial_symbol = qammod(serial_bit',M,'InputType','bit');
%     TxConstellation = serial_symbol;
% 
%     parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
%     
%     %zero_padding
%     offset_1 = round((N_fft-N_carrier)/2); 
%     offset_2= N_fft-N_carrier-offset_1;
%     
%     parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm); parallel_symbol;zeros(offset_2,N_ofdm)];                    
%     
%     OFDM_symbol_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));
% 
% 
%     cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_ofdm);
%     OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];    
%     first_ofdm_symbol= OFDM_symbol_with_cp(:,1);
%       
% figure;
% [psd,f] = periodogram(first_ofdm_symbol, rectwin(length(first_ofdm_symbol)), N_fft*2, ...
%                       10, 'centered');
% plot(f,10*log10(psd),'linewidth',1.5);
% xlabel('Normalized frequency');
% ylabel('PSD (dBW/Hz)');

clc;
clear;

s=9;
N_ofdm=100;

N_carrier= 10;
nsamp=5;
N_fft=nsamp*N_carrier;

m=0;
length_cp=round(m*N_fft);

M=4;

ebn0_db = [0];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);
     
 rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
    
    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
     
%     for i=2:N_carrier
%         parallel_symbol(i,:)=0;
%     end
    
    
    parallel_symbol_1=parallel_symbol(1:round(N_carrier/2),:);
    parallel_symbol_2=parallel_symbol((round(N_carrier/2)+1):N_carrier,:);
    
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
    offset=N_fft-N_carrier;
    
   
    parallel_symbol_zeropadding = [parallel_symbol_1;zeros(offset,N_ofdm);parallel_symbol_2];   
     
    OFDM_symbol_with_zeros=ifft(parallel_symbol_zeropadding);
    
    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_ofdm);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];  
    OFDM_serial=reshape(OFDM_symbol_with_cp,[1,(N_fft+length_cp)*N_ofdm]);
    
    received=awgn(OFDM_serial,ebn0_db,'measured');
    
      
figure;
% 
% f=((1:length(OFDM_serial))-length(OFDM_serial)/2)/length(OFDM_serial); % PSD at the transmitter end
% plot(f,20*log10(abs(fftshift(fft(OFDM_serial))))-max(20*log10(abs(fftshift(fft(OFDM_serial))))))
f=((1:length(received))-length(received)/2)/length(received); % PSD at the receiver end
plot(f,20*log10(abs(fftshift(fft(received))))-max(20*log10(abs(fftshift(fft(received))))))
ylim([-60 0]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)');






