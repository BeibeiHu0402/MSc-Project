clc;
clear;

%% parameters
s=99;
N_symbol=10000;
N_carrier=10;
nsamp=1;
N_fft=nsamp*N_carrier;

M=4; %4-QAM
K=4; 
m=0;
length_cp=round(m*N_fft);

FBMC_serial=zeros(1,K*N_fft+(N_symbol-1)*N_fft/2);
parallel_symbol=zeros(N_carrier,N_symbol);

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

%% OQAM-FBMC
   
   %transmitter end
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier/2*N_symbol*log2(M)));
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
   
    parallel_symbol_before_oqam=reshape(serial_symbol,[N_carrier/2,N_symbol]);
     
    x=real(parallel_symbol_before_oqam);
    y=1i*imag(parallel_symbol_before_oqam);
    
    for n=1:N_symbol
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
    
    parallel_symbol_zeropadding = [parallel_symbol_1;zeros(offset,N_symbol);parallel_symbol_2];   

    FBMC_symbols_with_zeros=ifft(parallel_symbol_zeropadding);
      
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
       
     for n=1:N_symbol
        a=FBMC_symbols_with_zeros(:,n);
        b=[];
        for i=1:K
            b=[b,a'];
        end
        c=b.*h;
        PAPR_parallel=[PAPR_parallel,c.'];
     end
   
     received=awgn(PAPR_parallel,0,'measured');
    
%% OFDM
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier*N_symbol*log2(M)));
    
    serial_symbol = qammod(serial_bit',M,'InputType','bit');

    parallel_symbol=reshape(serial_symbol,[N_carrier,N_symbol]);
    
    offset_1 = round((N_fft-N_carrier)/2); 
    offset_2= N_fft-N_carrier-offset_1;
    
    parallel_symbol_zeropadding = [zeros(offset_1,N_symbol); parallel_symbol;zeros(offset_2,N_symbol)];                    
    OFDM_symbol_with_zeros=ifft(ifftshift(parallel_symbol_zeropadding));

    cp=OFDM_symbol_with_zeros(N_fft-length_cp+1:N_fft,1:N_symbol);
    OFDM_symbol_with_cp=[cp;OFDM_symbol_with_zeros];
    length_ofdm=length(OFDM_symbol_with_cp);     
                                           
    OFDM_serial=reshape(OFDM_symbol_with_cp,[1,(N_fft+length_cp)*N_symbol]);
    
     
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprFBMC] = PAPR(PAPR_parallel);

PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprOFDM] = PAPR(OFDM_symbol_with_cp);


% for i=1:N_ofdm
%     symbol=OFDM_symbol_with_cp(:,i);
%     papr(i)=max(abs(symbol).^2)/mean(abs(symbol).^2);
%     papr_db(i)=10*log10(papr(i));
% end

[cdf1, PAPR1] = ecdf(paprFBMC);
[cdf2, PAPR2] = ecdf(paprOFDM);

semilogy(PAPR1,1-cdf1,'linewidth',1.5,'color','b');
hold on;
semilogy(PAPR2,1-cdf2,'linewidth',1.5,'color','r');

xlabel('PAPR0(dB)');
ylabel('PAPR>PAPR0');





