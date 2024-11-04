function [stima,Fourier_coeff,LR] = japan_coefficient2(P,x,t,h,K,L,T)      
      
% Calcola i coefficienti di Fourier per l=1,...,L:

n=length(P)-1;
somma1=zeros(1,L); somma2=zeros(1,L); stima=zeros(L,length(x)); s1=0; s2=0;
results=zeros(1,length(x)); Fourier_coeff=zeros(L,length(x)); LR=zeros(1,length(x));
ret=(P(2:n)-P(1:n-1));
      
 for q=1:length(x)   
     l=0;
     s1 = sum(exp(-2*pi*1i*l*(t(1:n-1))/T).*sqrt((1/h)*K((P(1:n-1)-x(q))/h)).*ret) ;
     s2 = sum(exp(2*pi*1i*l*(t(1:n-1))/T).*sqrt((1/h)*K((P(1:n-1)-x(q))/h)).*ret);
     results(q)=s1*s2;
     
     for l=1:L
            somma1(l) = sum(exp(-2*pi*1i*l*(t(1:n-1))/T).*sqrt((1/h)*K((P(1:n-1)-x(q))/h)).*ret) ;        
            somma2(l) = sum(exp(2*pi*1i*l*(t(1:n-1))/T).*sqrt((1/h)*K((P(1:n-1)-x(q))/h)).*ret);
            results(q)=results(q)+2*somma1(l).*somma2(l);
            Fourier_coeff(l,q)=results(q)*(1/(2*l+1));
            
            LR(q)= (T/(n*h))*sum(K((P(1:n-1)-x(q))/h));
            stima(l,q)=Fourier_coeff(l,q)./LR(q);
     end    
 end

        