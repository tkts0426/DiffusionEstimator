function [stima,Fourier_coeff,LR] = japan_coefficient3(p,x,t,h,K,L,T)      
      
% Calcola i coefficienti di Fourier per l=1,...,L:

n=length(p);
somma1=zeros(L,length(x)); somma2=zeros(L,length(x)); stima=zeros(L,length(x)); 
Fourier_coeff=zeros(L,length(x));

P=p; tt=t(1:n-1); X=x;
for q=1:length(x)-1
    P=[P,p]; tt=[tt,t(1:n-1)];
end
for k=1:n-2
    X=[X;x];
end
ret=diff(P); P=P(1:n-1,:);
      
     s1 = sum(sqrt((1/h)*K((P-X)/h)).*ret) ;
     s2 = sum(sqrt((1/h)*K((P-X)/h)).*ret);
     results=s1.*s2; % l=0
     
     LR= (T/(n*h))*sum(K((P-X)/h));
     
     for l=1:L
            somma1(l,:) = sum(exp(-2*pi*1i*l*tt/T).*sqrt((1/h).*K((P-X)/h)).*ret) ;        
%             somma2(l,:) = sum(exp(2*pi*1i*l*tt/T).*sqrt((1/h).*K((P-X)/h)).*ret);
            somma2(l,:) =conj(somma1(l,:));
            results=results+2*somma1(l,:).*somma2(l,:);
            Fourier_coeff(l,:)=results*(1/(2*l+1));
            
            stima(l,:)=Fourier_coeff(l,:)./LR;
     end    
 end

        