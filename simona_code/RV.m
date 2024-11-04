function ker=RV(P,h,x,K,T,n,m)

% Realized Volatility based kernel estimator:

ker=zeros(1,length(x)); vol=zeros(n,1);
for k=1:n
    vol(k) =sum(diff(P((k-1)*m+1:k*m+1)).^2); % integrated volatility on each day using intraday data
end

for q=1:length(x)
    num=K((x(q)-P(1:m:end-1))/h).*vol;
    den=K((x(q)-P(1:m:end-1))/h);
    
    ker(q)=sum(num)*n/(T*sum(den));
end
