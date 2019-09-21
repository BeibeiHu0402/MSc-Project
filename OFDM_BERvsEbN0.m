clc;
clear;

s=9;
N_symbol=100000;
N_carrier= 10;
nsamp=2;
N_fft=nsamp*N_carrier;
m=0.25;
length_cp=round(m*N_fft);

M=4;
ebn0_db = [0:1:10];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

ber_thertical_4 = [];
ber_ofdm_4=[];

for i=1:length(snr_db)
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_symbol*log2(M)));% generate bits
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);% QAM modulation
    parallel_symbol=reshape(serial_symbol,[N_carrier,N_symbol]); % serial-to-parallel   
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
    parallel_symbol_zeropadding = [zeros(offset_1,N_symbol); parallel_symbol;zeros(offset_2,N_symbol)];%oversampling                                 
    OFDM_symbol_with_zeros=ifft(parallel_symbol_zeropadding);% ifft
    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_symbol);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];% insert cyclic prefix                                              
    OFDM_serial=reshape(OFDM_symbol_with_cp,[1,(N_fft+length_cp)*N_symbol]);% parallel-tp-serial    
    received=awgn(OFDM_serial,snr_db(i),'measured');% AWGN channel   
    received_parallel=reshape(received,[N_fft+length_cp,N_symbol]); % serial-to-parallel   
    received_parallel(1:length_cp,:)=[]; %remove cp(now parallel)         
    received_fft_with_zeros=fft(received_parallel); % fft 
    received_fft=received_fft_with_zeros(offset_1+1:N_carrier+offset_1,1:N_symbol);% remove zeros
    received_demod=qamdemod(received_fft,M,'OutputType','bit');% QAM demodulation
    received_serial=reshape(received_demod,[1,N_carrier*N_symbol*log2(M)]); % parallel-to-serial
    [number,ratio] = biterr(serial_bit,received_serial);% compute BER   
    
    ber_ofdm_4(i)=ratio;
    ber_thertical_4(i) = berawgn(ebn0_db(i),'qam',M);
end




M=8;
ebn0_db = [0:1:10];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

% ber_thertical_4 = [];
ber_thertical_8 = [];
% ber_ofdm_4=[];
ber_ofdm_8=[];


for i=1:length(snr_db)
    
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_symbol*log2(M)));
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
    
    parallel_symbol=reshape(serial_symbol,[N_carrier,N_symbol]);
    
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;

    parallel_symbol_zeropadding = [zeros(offset_1,N_symbol); parallel_symbol;zeros(offset_2,N_symbol)];                                
 
    OFDM_symbol_with_zeros=ifft(parallel_symbol_zeropadding);

    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_symbol);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];
                                              
    OFDM_serial=reshape(OFDM_symbol_with_cp,[1,(N_fft+length_cp)*N_symbol]);
    
    received=awgn(OFDM_serial,snr_db(i),'measured');
    
    received_parallel=reshape(received,[N_fft+length_cp,N_symbol]);
    
    %remove cp(now parallel)
    received_parallel(1:length_cp,:)=[];
          
    %fft operation on parallel matrix
    received_fft_with_zeros=fft(received_parallel);  
    received_fft=received_fft_with_zeros(offset_1+1:N_carrier+offset_1,1:N_symbol);

    received_demod=qamdemod(received_fft,M,'OutputType','bit');
    received_serial=reshape(received_demod,[1,N_carrier*N_symbol*log2(M)]);
    
    [number,ratio] = biterr(serial_bit,received_serial);   
    disp(ratio); 
    ber_ofdm_8(i)=ratio;
    ber_thertical_8(i) = berawgn(ebn0_db(i),'qam',M);
end

semilogy(ebn0_db,ber_ofdm_4,'k-o','linewidth',1.5,'color','b');
hold on;
semilogy(ebn0_db,ber_thertical_4,'k-o','linewidth',1.5,'color','r');
hold on;
semilogy(ebn0_db,ber_ofdm_8,'k-o','linewidth',1.5,'color','k');
hold on;
semilogy(ebn0_db,ber_thertical_8,'k-o','linewidth',1.5,'color','m');



