function ker=Kernel_BPReference(P,vol,h,x,K,ND)

%n=length(P);
n=length(vol);

num=K((x-P)/h).*vol;
den=K((x-P)/h);
    
%ker=sum(num)*n/(T*sum(den));
ker=sum(num)/sum(den);