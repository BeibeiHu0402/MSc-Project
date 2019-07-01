clc;
clear;

s=9;
N_ofdm=10;

N_carrier= 256;
nsamp=100;
N_fft=nsamp*N_carrier;

M=4;

m=0.2;
length_cp=round(m*N_fft);

ebn0_db = [0:1:10];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

ber_thertical = [];
ber=[];

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
      
%PAPR
PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprOFDM] = PAPR2(OFDM_symbol_with_cp);

for i=1:N_ofdm
    symbol=OFDM_symbol_with_cp(:,i);
    papr(i)=max(abs(symbol).^2)/mean(abs(symbol).^2);
    papr_db(i)=10*log10(papr(i));
end

%oversampling better shows PAPR




