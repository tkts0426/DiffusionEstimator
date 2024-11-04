function ker=RV3(P,h,x,K,T,n,dayindex)

% Realized Volatility based kernel estimator:

ker=zeros(1,length(x)); vol=zeros(n,1);
for k=1:n
%     vol(k) =sum(diff(P((k-1)*m+1:min(k*m+1,length(P)))).^2); % integrated volatility on each day using intraday data
    vol(k) =sum(diff(P(dayindex(k):dayindex(k+1))).^2); % integrated volatility on each day using intraday data
end

P_daily=P(dayindex); % daily data
for q=1:length(x)
    num=K((x(q)-P_daily(1:end-1))/h).*vol;
    den=K((x(q)-P_daily(1:end-1))/h);
    
    ker(q)=sum(num)*n/(T*sum(den));
end
