clc;
clear;

%% parameters
s=99;
N_fbmc=10000;
N_carrier=10;
nsamp=1;
N_fft=nsamp*N_carrier;

M=4; %4-QAM
K=4; %overlapping factor

ebn0_db = [0];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

ber_thertical = [];
ber=[];

FBMC_serial=zeros(1,K*N_fft+(N_fbmc-1)*N_fft/2);
parallel_symbol=zeros(N_carrier,N_fbmc);

PAPR_parallel=[];
%% overlapping factor and coefficients
switch K
    case 2
        Hk_1=sqrt(2)/2;
    case 3
        Hk_1=0.911438;
        Hk_2=0.411438;
    case 4
        Hk_1=0.971960;
        Hk_2=sqrt(2)/2;
        Hk_3=0.235147;
    otherwise
        return
end

%% PSD measurement

%     rng(s);
%     serial_bit=round(randi([0 1],1,N_carrier*N_fbmc*log2(M)));
%     serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
%    
%     parallel_symbol=reshape(serial_symbol,[N_carrier,N_fbmc]);
%     
% %     for i=1:N_carrier/2
% %         parallel_symbol(i,:)=0;
% %     end
% %     
% %     for i=N_carrier/2+2:N_carrier
% %         parallel_symbol(i,:)=0;
% %     end
%     
%     parallel_symbol_1=parallel_symbol(1:round(N_carrier/2),:);
%     parallel_symbol_2=parallel_symbol((round(N_carrier/2)+1):N_carrier,:);
%     
%     offset_1 = round((N_fft-N_carrier)/2); 
%     offset_2= N_fft-N_carrier-offset_1;
%     offset=N_fft-N_carrier;
%     
%      parallel_symbol_zeropadding = [parallel_symbol_1;zeros(offset,N_fbmc);parallel_symbol_2];   
% %     parallel_symbol_zeropadding = [zeros(offset_1,N_ofdm);parallel_symbol;zeros(offset_2,N_ofdm)]; 
%      
%     FBMC_symbols_with_zeros=ifft(parallel_symbol_zeropadding);
% %     FBMC_symbols_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));
   
 
    %transmitter end
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier/2*N_fbmc*log2(M)));
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
   
    parallel_symbol_before_oqam=reshape(serial_symbol,[N_carrier/2,N_fbmc]);
    
%     for i=2:N_carrier/2
%         parallel_symbol_before_oqam(i,:)=0;
%     end
    
    
    x=real(parallel_symbol_before_oqam);
    y=1i*imag(parallel_symbol_before_oqam);
    
    for n=1:N_fbmc
    if rem(n,2)==1     % Odd symbols
        parallel_symbol(1:2:N_carrier,n) = x(:,n);
        parallel_symbol(2:2:N_carrier,n) = y(:,n);
    else               % Even symbols
        parallel_symbol(1:2:N_carrier,n)= y(:,n);
        parallel_symbol(2:2:N_carrier,n) =x(:,n);
   
    end
    end
    
    parallel_symbol_1=parallel_symbol(1:round(N_carrier/2),:);
    parallel_symbol_2=parallel_symbol((round(N_carrier/2)+1):N_carrier,:);    
    offset=N_fft-N_carrier;
    
    parallel_symbol_zeropadding = [parallel_symbol_1;zeros(offset,N_fbmc);parallel_symbol_2];   
%     parallel_symbol_zeropadding=[-1+0j;0-j;-1+0j;0-j;1+0j;0-j;-1+0j;0+j];

    FBMC_symbols_with_zeros=ifft(parallel_symbol_zeropadding);


%        
     for i=0:K*N_fft-1
        switch K
            case 2
                h(i+1)=1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)));
            case 3
                h(i+1)=1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)))+2*(Hk_2*cos(2*pi*2*i/(K*N_fft)));
            case 4
                h(i+1)=1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)))+2*(Hk_2*cos(2*pi*2*i/(K*N_fft)))-2*(Hk_3*cos(2*pi*3*i/(K*N_fft)));
            otherwise
         return
        end
     end
       
     for n=1:N_fbmc
        a=FBMC_symbols_with_zeros(:,n);
        b=[];
        for i=1:K
            b=[b,a'];
        end
        c=b.*h;
        PAPR_parallel=[PAPR_parallel,c.'];
     end
   
     received=awgn(PAPR_parallel,0,'measured');
     
     
% figure;
% f=((1:length(PSD_serial))-length(PSD_serial)/2)/length(PSD_serial);
% plot(f,20*log10(abs(fftshift(fft(PSD_serial))))-max(20*log10(abs(fftshift(fft(PSD_serial))))));



% f=((1:length(received))-length(received)/2)/length(received);
% plot(f,20*log10(abs(fftshift(fft(received))))-max(20*log10(abs(fftshift(fft(received))))));


PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprFBMC] = PAPR2(PAPR_parallel);


% for i=1:N_ofdm
%     symbol=OFDM_symbol_with_cp(:,i);
%     papr(i)=max(abs(symbol).^2)/mean(abs(symbol).^2);
%     papr_db(i)=10*log10(papr(i));
% end

[cdf1, PAPR1] = ecdf(paprFBMC);
semilogy(PAPR1,1-cdf1,'linewidth',1.5,'color','b');
xlabel('PAPR0(dB)');

ylabel('PAPR>PAPR0');





