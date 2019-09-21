clc;
clear;

%% parameters
s=99;
N_symbol=10000;
N_carrier=10;
nsamp=2;
N_fft=nsamp*N_carrier;
K=4; %overlapping factor

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

for i=0:K*N_fft-1
   switch K
      case 2
          h(i+1)=(1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft))))/(1+2*(Hk_1));
      case 3
          h(i+1)=(1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)))+2*(Hk_2*cos(2*pi*2*i/(K*N_fft))))/(1+2*(Hk_1+Hk_2));
      case 4
          h(i+1)=(1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)))+2*(Hk_2*cos(2*pi*2*i/(K*N_fft)))-2*(Hk_3*cos(2*pi*3*i/(K*N_fft))))/(1+2*(Hk_1+Hk_2+Hk_3));
      otherwise
          return
   end
end

% %% noise 
% M=4; %QAM
% ebn0_db = [5];
% ebn0=10.^(0.1*ebn0_db);
% esn0=log2(M)*ebn0;
% esn0_db=10*log10(esn0);
% snr_db=esn0_db-10*log10(nsamp);
% 
% ber_thertical_4 = [];
% ber_fbmc_4=[];
% 
% FBMC_serial=zeros(1,K*N_fft+(N_symbol-1)*N_fft/2);
% parallel_symbol=zeros(N_carrier,N_symbol);
% 
%     %transmitter end
%     rng(s);
%     serial_bit=round(randi([0 1],1,N_carrier/2*N_symbol*log2(M)));
%     serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
%     TxConstellation=serial_symbol;
%     parallel_symbol_before_oqam=reshape(serial_symbol,[N_carrier/2,N_symbol]);
%     
%     x=real(parallel_symbol_before_oqam);
%     y=1i*imag(parallel_symbol_before_oqam);
%     
%     for n=1:N_symbol
%     if rem(n,2)==1     % Odd symbols
%         parallel_symbol(1:2:N_carrier,n) = x(:,n);
%         parallel_symbol(2:2:N_carrier,n) = y(:,n);
%     else               % Even symbols
%         parallel_symbol(1:2:N_carrier,n)= y(:,n);
%         parallel_symbol(2:2:N_carrier,n) =x(:,n);
%    
%     end
%     end
%     
%     parallel_symbol_1=parallel_symbol(1:round(N_carrier/2),:);
%     parallel_symbol_2=parallel_symbol((round(N_carrier/2)+1):N_carrier,:);    
%     offset=N_fft-N_carrier;
%     
%     parallel_symbol_zeropadding = [parallel_symbol_1;zeros(offset,N_symbol);parallel_symbol_2];   
% %     parallel_symbol_zeropadding=[-1+0j;0-j;-1+0j;0-j;1+0j;0-j;-1+0j;0+j];
% 
%     FBMC_symbols_with_zeros=ifft(parallel_symbol_zeropadding);
%            
%     for n=1:N_symbol
%         one_FBMC_symbol=FBMC_symbols_with_zeros(:,n);
%         duplication=[];
%         for i=1:K
%             duplication=[duplication,one_FBMC_symbol.'];
%         end
%         filtered_output=duplication.*h;
%         FBMC_serial(1+(n-1)*N_fft/2:(n-1)*N_fft/2+K*N_fft)=FBMC_serial(1+(n-1)*N_fft/2:(n-1)*N_fft/2+K*N_fft)+filtered_output;
%     end
%  
% 
%  for i=1:length(snr_db)
%      received_bit=[];
%      RxConstellation=[];
%     received=awgn(FBMC_serial,snr_db(i),'measured');
% %      received=FBMC_serial;
%      
%      for n=1:N_symbol
%         R=received(1+(n-1)*N_fft/2:(n-1)*N_fft/2+K*N_fft).*h; % We apply the filter
% 
%         U=zeros(1,N_fft);
%          for k=1:4
%             U=U+R(1,1+(k-1)*N_fft:k*N_fft);                  % Summation
%          end
% 
%         U=U.';
%         SEST=fft(U);
%         SEST=SEST/0.6863;                           % Normalization
%         
%         SEST_1=SEST(1:round(N_carrier/2),:);
%         SEST_2=SEST(round(N_carrier/2)+offset+1:N_fft,:);
%         SEST=[SEST_1;SEST_2];
%         
%    
%         if rem(n,2)==1
%             receiver_real=real(SEST(1:2:N_carrier));
%             receiver_imaginary=imag(SEST(2:2:N_carrier));
%             receiver_demod = complex(receiver_real, receiver_imaginary);
%         else
%             receiver_imaginary=imag(SEST(1:2:N_carrier));
%             receiver_real=real(SEST(2:2:N_carrier));
%             receiver_demod = complex(receiver_real, receiver_imaginary);
%         end
%         RxConstellation=[RxConstellation,receiver_demod.'];
%         bit_demod=qamdemod(receiver_demod,M,'OutputType','bit')';% OQAM demodulation
%         received_bit=[received_bit,bit_demod];
%      end
%      
%        [number,ratio] = biterr(serial_bit,received_bit);   
%      disp(ratio); 
%     ber_fbmc_4(i)=ratio;
%     ber_thertical_4(i) = berawgn(ebn0_db(i),'qam',M);
%            
%  end

 %% noise 
M=8; %QAM
ebn0_db = [5];
ebn0=10.^(0.1*ebn0_db);
esn0=log2(M)*ebn0;
esn0_db=10*log10(esn0);
snr_db=esn0_db-10*log10(nsamp);

ber_thertical_8 = [];
ber_fbmc_8=[];

FBMC_serial=zeros(1,K*N_fft+(N_symbol-1)*N_fft/2);
parallel_symbol=zeros(N_carrier,N_symbol);

    %transmitter end
    rng(s);
    serial_bit=round(randi([0 1],1,N_carrier/2*N_symbol*log2(M)));
    serial_symbol = qammod(serial_bit',M,'InputType','bit','PlotConstellation',false);
    TxConstellation=serial_symbol;
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
%     parallel_symbol_zeropadding=[-1+0j;0-j;-1+0j;0-j;1+0j;0-j;-1+0j;0+j];

    FBMC_symbols_with_zeros=ifft(parallel_symbol_zeropadding);
           
    for n=1:N_symbol
        one_FBMC_symbol=FBMC_symbols_with_zeros(:,n);
        duplication=[];
        for i=1:K
            duplication=[duplication,one_FBMC_symbol.'];
        end
        filtered_output=duplication.*h;
        FBMC_serial(1+(n-1)*N_fft/2:(n-1)*N_fft/2+K*N_fft)=FBMC_serial(1+(n-1)*N_fft/2:(n-1)*N_fft/2+K*N_fft)+filtered_output;
    end
 

 for i=1:length(snr_db)
     received_bit=[];
     RxConstellation=[];
    received=awgn(FBMC_serial,snr_db(i),'measured');
%      received=FBMC_serial;
     
     for n=1:N_symbol
        R=received(1+(n-1)*N_fft/2:(n-1)*N_fft/2+K*N_fft).*h; % We apply the filter

        U=zeros(1,N_fft);
         for k=1:4
            U=U+R(1,1+(k-1)*N_fft:k*N_fft);                  % Summation
         end

        U=U.';
        SEST=fft(U);
        SEST=SEST/0.6863;                           % Normalization
        
        SEST_1=SEST(1:round(N_carrier/2),:);
        SEST_2=SEST(round(N_carrier/2)+offset+1:N_fft,:);
        SEST=[SEST_1;SEST_2];
        
   
        if rem(n,2)==1
            receiver_real=real(SEST(1:2:N_carrier));
            receiver_imaginary=imag(SEST(2:2:N_carrier));
            receiver_demod = complex(receiver_real, receiver_imaginary);
        else
            receiver_imaginary=imag(SEST(1:2:N_carrier));
            receiver_real=real(SEST(2:2:N_carrier));
            receiver_demod = complex(receiver_real, receiver_imaginary);
        end
        RxConstellation=[RxConstellation,receiver_demod.'];
        bit_demod=qamdemod(receiver_demod,M,'OutputType','bit')';% OQAM demodulation
        received_bit=[received_bit,bit_demod];
     end
     
       [number,ratio] = biterr(serial_bit,received_bit);   
     disp(ratio); 
    ber_fbmc_8(i)=ratio;
    ber_thertical_8(i) = berawgn(ebn0_db(i),'qam',M);
           
end

% scatterplot(TxConstellation);

scatterplot(RxConstellation);
 hold on;
 line([-15,15],[-1,-1],'linestyle',':','color','r','linewidth',1.5);
 line([-15,15],[1,1],'linestyle',':','color','r','linewidth',1.5);
 line([-1,-1],[-15,15],'linestyle',':','color','r','linewidth',1.5);
 line([1,1],[-15,15],'linestyle',':','color','r','linewidth',1.5);
 line([-3,-3],[-15,15],'linestyle',':','color','r','linewidth',1.5);
 line([3,3],[-15,15],'linestyle',':','color','r','linewidth',1.5);
 hold off;



