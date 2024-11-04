function ker=FZ(P,h,x,K,T,n)

% Florens-Zmirou kernel estimator:

ker=zeros(1,length(x));
vol =diff(P).^2/(T/n) ;% spot volatility at each node

for q=1:length(x)
    num=K((x(q)-P(1:n))/h).*vol;
    den=K((x(q)-P(1:n))/h);
    
    %ker=sum(num)*n/(T*sum(den));
    ker(q)=sum(num)/sum(den);
end