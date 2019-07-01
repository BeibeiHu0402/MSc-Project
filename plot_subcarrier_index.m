clc;
clear;

s=9;
N_ofdm=1;

N_carrier= 4;
nsamp=2;
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

n_zoom=50;

    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
    
    serial_symbol = qammod(serial_bit',M,'InputType','bit');

    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
     
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;    

    parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm); parallel_symbol;zeros(offset_2,N_ofdm)];                                      
 
    first_ofdm_symbol_before_ifft=parallel_symbol_zeropadding(:,1);% the first symbol before IFFT
    diagonal=diag(first_ofdm_symbol_before_ifft); % plot overlapping subcarrier spectrum.
    diagonal_ifft=ifft(diagonal);
    diagonal_fft = fft(diagonal_ifft,N_carrier*n_zoom);%zoom out nzoom resolution
    
figure;
plot(-N_fft/2:N_fft/(N_carrier*n_zoom):(N_fft-N_fft/(N_carrier*n_zoom))/2,abs(diagonal_fft),'linewidth',1.5);
grid on;
xlabel('Subcarrier Index');
ylabel('Amplitude');


