function [stima,Fourier_coeff,LR] = japan_coefficientFFT3(p,x,h,K,L,T)      
      
% Calcola i coefficienti di Fourier per l=1,...,L usando la FFT:

n=length(p);

P=p; X=x;
for q=1:length(x)-1
    P=[P,p];
end
for k=1:n-2
    X=[X;x];
end
ret=diff(P); P=P(1:n-1,:);
R=sqrt((1/h)*K((P-X)/h)).*ret;
fft_v=fft(R); % FFT
idx=L+1:-1:2; ff=fft_v(idx,:);
fft_def=[conj(ff); fft_v(1:L+1,:)]; % Fourier coeff. of log-returns

% idxk=-L:1:L; nshift=L+1;
% idxx=idxk+nshift;
% Capp=fft_def(idxx,:);
% coeff=sum(Capp.*fft_def(nshift-L:nshift+L,:));

coeff=sum(fft_def.*fft_def(2*L+1:-1:1,:));

Fourier_coeff=coeff.*(1/(2*L+1)); % 0-th Fourier coeff. of variance
     
LR= (T/(n*h))*sum(K((P-X)/h)); % denominator

stima=Fourier_coeff./LR;
     
 end

        