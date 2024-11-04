function [a] = japan_coefficient(P,x,t,h,k,H)      
      
%     Input: 
%     n=numero di dati
%     nf = numero di frequenze
%     u = serie storica di dati simulati
%     Output:
% h = bandwidth
%     a(1:nf+1), b(1:nf+1) coefficienti di Fourier di dp
%     N.B. a(k+1) e' il k-esimo coefficiente
n=length(P)-1;
T=2;
L=n/10;
M=n/20;
%ND=50;
%m=n/2;
delta=1/25;
somma1=zeros(1,L); somma2=zeros(1,L);
      
%            somma1(k)=1;
%            somma2(k)=1;


            
     for l=1:L
	   
            somma1(l) = sum(exp(-2*pi*i*l*(t(2:n))/T).*sqrt((1/h)*H(P(2:n)-x)/h).*(P(2:n)-P(1:n-1))) ;        
            somma2(l) = sum(exp(-2*pi*i*(k-l)*(t(2:n))/T).*sqrt((1/h)*H(P(2:n)-x)/h).*(P(2:n)-P(1:n-1)));
     end
     % den=H((x-P)/h);      
     results=somma1.*somma2;
     
     a=sum(results)*(1/(2*L+1));
     %a=sum((sin(delta*[1:M])./(delta*[1:M])).^2.*results(1:M))*(1/(2*L+1));
     
     %*ND)/(T*sum(den));
     
     
