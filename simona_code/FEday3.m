function ker=FEday3(P,h,x,K,t,T,n,dayindex)

% Fourier based kernel estimator (day by day):

ker=zeros(1,length(x)); vol=zeros(n,1);
for k=1:n
%     N=floor(m/2); M=3; vol(k) = FE_spot_vol(P((k-1)*m+1:min(k*m+1,length(P))),t((k-1)*m+1:min(k*m+1,length(P))),t((k-1)*m+1),T/n,N,M); % spot volatility on each day using intraday data
    p=P(dayindex(k):dayindex(k+1)); m=length(p); % intraday data
    N=floor(m/2); M=3; vol(k) = FE_spot_vol(p,t(dayindex(k):dayindex(k+1)),t(dayindex(k)),T/n,N,M); % spot volatility on each day using intraday data
end

P_daily=P(dayindex); % daily data
for q=1:length(x)
    num=K((x(q)-P_daily(1:end-1))/h).*vol;
    den=K((x(q)-P_daily(1:end-1))/h);
    
    %ker(q)=sum(num)*n/(T*sum(den));
    ker(q)=sum(num)/sum(den);
end
