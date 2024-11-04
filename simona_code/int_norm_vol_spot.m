function INT=int_norm_vol_spot(X,x,sigma2,T,N,K,h)

dt=T/N; INT=zeros(1,length(x));

for q=1:length(x)
    INT(q)=sum(sigma2(1:end-1).*K((X(1:end-1)-x(q))/h)/h)*dt;
end