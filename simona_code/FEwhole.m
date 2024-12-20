function ker=FEwhole(P,h,x,K,t,T,n,m,Ndati)

% Fourier based kernel estimator (all-in-all):

ker=zeros(1,length(x)); 

% N=floor(Ndati/2); M=100; vol = FE_spot_vol(P,t,t(1:m:end-1),T,N,M); % spot volatility on each day using all intraday data
N=floor(Ndati/2); M=100; vol = FE_spot_vol1(P,t,t(1:m:end-1),T,N,M); % FFT-based spot volatility on each day using all intraday data

for q=1:length(x)
    num=K((x(q)-P(1:m:end-1))/h).*vol';
    den=K((x(q)-P(1:m:end-1))/h);
    
    %ker(q)=sum(num)*n/(T*sum(den));
    ker(q)=sum(num)/sum(den);
end
