clc;
clear;

s=9;
N_ofdm=100;

N_carrier= 100;
nsamp=20;
N_fft=nsamp*N_carrier;

M=4;

m=0.25;
length_cp=round(m*N_fft);

ebn0_db = [0:1:10];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

ber_thertical = [];
ber=[];
     
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
    
    serial_symbol = qammod(serial_bit',M,'InputType','bit');
    TxConstellation = serial_symbol;

    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
    
    %zero_padding
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
    
    parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm); parallel_symbol;zeros(offset_2,N_ofdm)];                    
    
    OFDM_symbol_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));

%     first_ofdm_symbol= OFDM_symbol_with_zeros(:,1);

    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_ofdm);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];    
    first_ofdm_symbol= OFDM_symbol_with_cp(:,1);
      
figure;
[psd,f] = periodogram(first_ofdm_symbol, rectwin(length(first_ofdm_symbol)), N_fft*2, ...
                      10, 'centered');
plot(f,10*log10(psd),'linewidth',1.5);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)');

% PSD with/without CP?



