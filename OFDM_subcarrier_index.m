clc;
clear;

s=9;
N_ofdm=1000;

N_carrier= 10;
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
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
    
    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);
    parallel_symbol_1=parallel_symbol(1:round(N_carrier/2),:);
    parallel_symbol_2=parallel_symbol((round(N_carrier/2)+1):N_carrier,:);
    
%     offset_1 = round((N_fft-N_carrier)/2); 
%     offset_2= N_fft-N_carrier-offset_1;
    offset=N_fft-N_carrier;
    
   
    parallel_symbol_zeropadding = [parallel_symbol_1;zeros(offset,N_ofdm);parallel_symbol_2];                                
 
    first_ofdm_symbol_before_ifft=parallel_symbol_zeropadding(:,1);% the first symbol before IFFT
    diagonal=diag(first_ofdm_symbol_before_ifft); % plot overlapping subcarrier spectrum.
    diagonal_ifft=ifft(ifftshift(diagonal));
    diagonal_fft = fft(diagonal_ifft,N_carrier*n_zoom);%zoom out nzoom resolution
    
figure;
plot(-N_fft/2:N_fft/(N_carrier*n_zoom):(N_fft-N_fft/(N_carrier*n_zoom))/2,abs(diagonal_fft),'linewidth',1.5);
grid on;
xlabel('Subcarrier Index');
ylabel('Amplitude');


