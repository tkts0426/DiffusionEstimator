function DEN=denorm_vol(X,x,T,N,K,h)

dt=T/N; DEN=zeros(1,length(x));

for q=1:length(x)
    DEN(q)=sum(K((X(1:end-1)-x(q))/h)/h)*dt;
end