clc;
clear;

H1=0.971960;
H2=sqrt(2)/2;
H3=0.235147;
HkOneSided = [H1 H2 H3];
Hk = [fliplr(HkOneSided) 1 HkOneSided];

N_fft=10;

for K=2:4
for i=0:K*N_fft-1
   switch K
      case 2
          Hk_1=sqrt(2)/2;
          a(i+1)=(1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft))))/(1+2*(Hk_1));
      case 3
          Hk_1=0.911438;
          Hk_2=0.411438;
          b(i+1)=(1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)))+2*(Hk_2*cos(2*pi*2*i/(K*N_fft))))/(1+2*(Hk_1+Hk_2));
      case 4
          Hk_1=0.971960;
          Hk_2=sqrt(2)/2;
          Hk_3=0.235147;
          c(i+1)=
          (1-2*(Hk_1*cos(2*pi*1*i/(K*N_fft)))+
          2*(Hk_2*cos(2*pi*2*i/(K*N_fft)))-
          2*(Hk_3*cos(2*pi*3*i/(K*N_fft))))/(1+2*(Hk_1+Hk_2+Hk_3));
      otherwise
          return
   end 
end
end

fvtool(c,1,b,1,a,1);
legend('K=4','K=3','K=2');




