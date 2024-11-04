function [stima,Fourier_coeff,LR] = japan_coefficientFFT1(P,x,h,K,L,T)      
      
% Calcola i coefficienti di Fourier per l=1,...,L usando la FFT:

n=length(P);
Fourier_coeff=zeros(1,length(x)); LR=zeros(1,length(x));
ret=diff(P); P=P(1:n-1,:);
      
 for q=1:length(x) 
     R = sqrt((1/h)*K((P-x(q))/h)).*ret ;
     fft_v=fft(R); % FFT
     idx=L+1:-1:2; ff=fft_v(idx);
     fft_def=[conj(ff); fft_v(1:L+1)]; % Fourier coeff. of log-returns

     coeff=sum(fft_def.*fft_def(2*L+1:-1:1));

     Fourier_coeff(q)=coeff.*(1/(2*L+1)); % 0-th Fourier coeff. of variance
     
     LR(q)= (T/(n*h))*sum(K((P-x(q))/h)); % denominator
 end

stima=Fourier_coeff./LR;
     
 end

        