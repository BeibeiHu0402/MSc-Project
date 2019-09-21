clc;
clear;

s=9;
N_ofdm=10000;

N_carrier= 10;
    

%% no oversampling
nsamp=1;
N_fft=nsamp*N_carrier;

M=4;

m=0;
length_cp=round(m*N_fft);

papr=[];
papr_db=[];
  
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
    
    serial_symbol = qammod(serial_bit',M,'InputType','bit');

    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
    
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
    
    parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm); parallel_symbol;zeros(offset_2,N_ofdm)];                    
    OFDM_symbol_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));

    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_ofdm);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];
    length_ofdm=length(OFDM_symbol_with_cp);     
                                           
    OFDM_serial=reshape(OFDM_symbol_with_cp,[1,(N_fft+length_cp)*N_ofdm]);
      



%% oversampled by 5
nsamp_5=5;
N_fft=nsamp_5*N_carrier;

M=4;

m=0;
length_cp=round(m*N_fft);



papr=[];
papr_db=[];
  
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
    
    serial_symbol = qammod(serial_bit',M,'InputType','bit');

    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
    
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
    
    parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm); parallel_symbol;zeros(offset_2,N_ofdm)];                    
    OFDM_symbol_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));

    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_ofdm);
    OFDM_symbol_with_cp_5=[cp;OFDM_symbol_with_zeros];
    length_ofdm=length(OFDM_symbol_with_cp_5);     
                                           
    OFDM_serial=reshape(OFDM_symbol_with_cp_5,[1,(N_fft+length_cp)*N_ofdm]);
      
%PAPR

PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW'); % no oversampling
[~,~,paprOFDM] = PAPR2(OFDM_symbol_with_cp);
% for i=1:N_ofdm
%     symbol=OFDM_symbol_with_cp(:,i);
%     papr(i)=max(abs(symbol).^2)/mean(abs(symbol).^2);
%     papr_db(i)=10*log10(papr(i));
% end
PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW'); % oversampled by 5
[~,~,paprOFDM_5] = PAPR2(OFDM_symbol_with_cp_5);
% for i=1:N_ofdm
%     symbol=OFDM_symbol_with_cp_5(:,i);
%     papr(i)=max(abs(symbol).^2)/mean(abs(symbol).^2);
%     papr_db(i)=10*log10(papr(i));
% end
[cdf1, PAPR1] = ecdf(paprOFDM);
[cdf2, PAPR2] = ecdf(paprOFDM_5);
semilogy(PAPR1,1-cdf1,'linewidth',1.5,'color','b');
hold on
semilogy(PAPR2,1-cdf2,'linewidth',1.5,'color','r');
xlabel('PAPR0(dB)');
ylabel('PAPR>PAPR0');









