function spot = FE_spot_vol(P,t,tau,T,N,M)

% Computes the spot variance by the Fourier estimator with Dirichlet kernel

% Input variables:
% P vector of the observed log-prices
% t vector of the observation times
% tau vector of the times where the volatility is estimated
% N maximum Fourier frequency for price returns
% M maximum Fourier frequency for spot variance

% Output variables:
% spot vector of spot variance at the time grid tau
% c_s Fourier coefficients of the spot variance

n=max(size(P)); nv=max(size(tau)); const=2*pi/T;
c_pp=zeros(N+M,1); c_p=zeros(2*N+2*M+1,1);
r=diff(P); spot=zeros(1,nv);
c_0=sum(r); c_s=zeros(1,2*M+1);

for k=1:N+M
    c_pp(k)=sum(exp(-1i*const*k*t(1:end-1)).*r);
end
for j=1:N+M
    c_p(j)=conj(c_pp(N+M+1-j))/T;
end
    c_p(N+M+1)=c_0/T;
for j=1:N+M
    c_p(N+M+1+j)=c_pp(j)/T;
end

% Fourier coefficients of the spot variance in [0,T]
fact=T/(2*N+1);
nshift=N+M+1;
for k=-M:M
    c_s(k+M+1)=0.0;
    for l=-N:N
        c_s(k+M+1)=c_s(k+M+1)+fact*(c_p(l+nshift)*c_p(k-l+nshift));
    end
end
for it=1:nv
    spot(it)=0.0;
    for k=-M:M
        spot(it)=spot(it)+(1-abs(k)/M)*c_s(k+M+1)*exp(1i*tau(it)*const*k);
    end
end
spot=real(spot);