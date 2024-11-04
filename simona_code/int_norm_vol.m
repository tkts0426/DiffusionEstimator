function INT=int_norm_vol(X,x,eta,gamma,T,N,K,h)

dt=T/N; INT=zeros(1,length(x));

for q=1:length(x)
    INT(q)=sum((eta*X(1:end-1).^gamma).^2.*K((X(1:end-1)-x(q))/h)/h)*dt;
end