function ker=FEday(P,h,x,K,t,T,n,m)

% Fourier based kernel estimator (day by day):

ker=zeros(1,length(x)); vol=zeros(n,1);
for k=1:n
    %vol(k) = sum(diff(P((k-1)*m+1:k*m+1)).^2); % integrated volatility on each day using intraday data
    N=floor(m/2); M=10; vol(k) = FE_spot_vol(P((k-1)*m+1:k*m+1),t((k-1)*m+1:k*m+1),t((k-1)*m+1),T/n,N,M); % spot volatility on each day using intraday data
end

for q=1:length(x)
    num=K((x(q)-P(1:m:end-1))/h).*vol;
    den=K((x(q)-P(1:m:end-1))/h);
    
    %ker(q)=sum(num)*n/(T*sum(den));
    ker(q)=sum(num)/sum(den);
end
