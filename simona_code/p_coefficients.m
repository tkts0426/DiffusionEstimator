function [a,b] = p_coefficients(n,nf,u,t)      
      
%     Input: 
%     n=numero di dati
%     nf = numero di frequenze
%     u = serie storica di dati
%     Output:
%     a(1:nf+1), b(1:nf+1) coefficienti di Fourier di dp
%     N.B. a(k+1) e' il k-esimo coefficiente
      

%tmp=zeros(n,1);
%for i=1:n
%    tmp(i) = 2*pi*(t(i)-t(1))/(t(n)-t(1)); % normalizza gli istanti delle osservazioni su [0,2 pi]
%end
            
      for k=1:nf+1
	   
            somma1 = sum(u(1:n-1).*(cos((k-1)*t(2:n))-cos((k-1)*t(1:n-1))));            
            somma2 = sum(u(1:n-1).*(sin((k-1)*t(2:n))-sin((k-1)*t(1:n-1))));
                      
         a(k) = (1/pi)*(u(n)-u(1)-somma1);
         b(k) = -(1/pi)*somma2;

%         a(k) = (1/pi)*(u(n)-u(1)+somma1)
%         b(k) = +(1/pi)*somma2
                  
     end
            
     a(1) = a(1)/2.0;
     b(1) = 0.0;
