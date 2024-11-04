function [spot] = FE_spot_vol_FFT(P,T,N,M)

% Computes the spot variance by the Fourier estimator with Dirichlet kernel and FFT

% Input variables:
% P vector of the observed log-prices
% N maximum Fourier frequency for price returns
% M maximum Fourier frequency for spot variance

% Output variables:
% spot vector of spot variance at the time grid tau

r=diff(P); fft_v=fft(r);
idx=M+N+1:-1:2; ff=fft_v(idx);
fft_def=[conj(ff) fft_v(1:M+N+1)];
fft_def=fft_def./T; % Fourier coeff. of log-returns

idxk=-N:1:N; nshift=M+N+1;
for kk=-M:M
    idxx=idxk+nshift+kk;
    Capp=fft_def(idxx);
    coeff(M+kk+1)=Capp*fft_def(nshift-N:nshift+N)';
end
C=coeff.*(T/(2*N+1)); % Fourier coeff. of variance

k=(-M:1:M);
f=C.*(1-abs(k)/M);
Fsum=(2*M+1)*ifft(f);
Fsum=exp(-1i*2*pi*M*(k+M)/(2*M+1)).*Fsum;
spot=real(Fsum);