clc;
clear;

s=9;
N_ofdm=100000;

N_carrier= 4;
nsamp=2;
N_fft=nsamp*N_carrier;

M=4;

m=0.25;
length_cp=round(m*N_fft);

ebn0_db = [10];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

ber_thertical = [];
ber=[];
    
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_ofdm*log2(M)));
    
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',true);
    TxConstellation = serial_symbol;

    parallel_symbol=reshape(serial_symbol,[N_carrier,N_ofdm]);

    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
     
    parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm); parallel_symbol;zeros(offset_2,N_ofdm)];                    
    OFDM_symbol_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));

    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_ofdm);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];
                                              
    OFDM_serial=reshape(OFDM_symbol_with_cp,[1,(N_fft+length_cp)*N_ofdm]);
    
    received=awgn(OFDM_serial,snr_db,'measured');
    
    received_parallel=reshape(received,[N_fft+length_cp,N_ofdm]);
     
    %remove cp(now parallel)
    received_parallel(1:length_cp,:)=[];
          
    %fft operation on parallel matrix
    received_fft_with_zeros=fftshift(fft(received_parallel));  
    received_fft=received_fft_with_zeros(offset_1+1:N_carrier+offset_1,1:N_ofdm);
    RxConstellation = reshape(received_fft,[1,N_carrier*N_ofdm]);
    received_demod=qamdemod(received_fft,M,'OutputType','bit');
    received_serial=reshape(received_demod,[1,N_carrier*N_ofdm*log2(M)]);
    [number,ratio] = biterr(serial_bit,received_serial);   
%     disp(ratio); 
    ber=ratio;

%consteallation digaram
%  scatterplot(TxConstellation);
%  hold on;
%  line([-2,2],[-1,-1],'linestyle',':','color','r','linewidth',1.5);
%  line([-2,2],[1,1],'linestyle',':','color','r','linewidth',1.5);
%  line([-1,-1],[-2,2],'linestyle',':','color','r','linewidth',1.5);
%  line([1,1],[-2,2],'linestyle',':','color','r','linewidth',1.5);
%  hold off;
 
 scatterplot(RxConstellation);
 hold on;
 line([-15,15],[-1,-1],'linestyle',':','color','r','linewidth',1.5);
 line([-15,15],[1,1],'linestyle',':','color','r','linewidth',1.5);
 line([-1,-1],[-15,15],'linestyle',':','color','r','linewidth',1.5);
 line([1,1],[-15,15],'linestyle',':','color','r','linewidth',1.5);
 hold off;




