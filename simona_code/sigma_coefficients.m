function [av,bv] = sigma_coefficients(ndata,xx,tt,N,Nv,n0)
      
% Author: Roberto Reno', 2002

%     Computes the Nv coefficients of volatility (av,bv)
%     given N coefficients of price (a,b)
%     the convention is that the coefficients range from 1 to N+1
%     according to a_k = a(k+1), k =0,N
%     then volatility is given by fourier-fejer formula
%     n0 = see Malliavin-Mancino
%     DEVE ESSERE Nv <= N-n0


% INPUT: 
%     ndata (INTEGER): the number of data for the stock
%     xx(1:ndata): vettore colonna dei prezzi per lo stock 
%     tt(1:ndata): vettore colonna dei tempi per lo stock
%     N: the cut frequency for the stock
% OUTPUT:
%     intvol: integrated variance versus # of Fourier coefficients
%     coeffA,coeffB: i coefficienti (coeffA(1) = 2 * pi * a0(sigma^2), coeffB(1) = 0.)


%n0 = 1;
      
[a,b] = p_coefficients(ndata,N,xx,tt); % a(1:N+1),b(1:N+1) coefficienti di Fourier (N.B. a(k+1) e' il k-esimo coefficiente)

av(1)=0.5*sum(a(n0+1:N+1).^2+b(n0+1:N+1).^2)*pi/(N+1-n0); % a_0(sigma^2)
%av(1)=0.5*(2*a(1)^2+sum(a(n0+1:N+1).^2+b(n0+1:N+1).^2))*pi/(N+1/2); % a_0(sigma^2)
bv(1)=0;
for k=1:Nv
    % Usa N-k termini: E' IL MIGLIORE!
    av(k+1)=sum(a(n0+1:N-k+1).*a(n0+1+k:N+1)+b(n0+1:N-k+1).*b(n0+1+k:N+1))*pi/(N+1-n0-k); % a_k(sigma^2)
    bv(k+1)=sum(a(n0+1:N-k+1).*b(n0+1+k:N+1)-b(n0+1:N-k+1).*a(n0+1+k:N+1))*pi/(N+1-n0-k); % b_k(sigma^2)
    % Usa N-Nv termini: E' SBAGLIATO, MA FUNZIONA BENE!!! NO!!!!!!
%      av(k+1)=sum(a(n0+1:N-Nv+1).*a(n0+1+Nv:N+1)+b(n0+1:N-Nv+1).*b(n0+1+Nv:N+1))*pi/(N-Nv+1-n0); % a_k(sigma^2)
%      bv(k+1)=sum(a(n0+1:N-Nv+1).*b(n0+1+Nv:N+1)-b(n0+1:N-Nv+1).*a(n0+1+Nv:N+1))*pi/(N-Nv+1-n0); % b_k(sigma^2)
     % Usa N-Nv termini:
%      av(k+1)=sum(a(n0+1:N-Nv+1).*a(n0+1+k:N-Nv+1+k)+b(n0+1:N-Nv+1).*b(n0+1+k:N-Nv+1+k))*pi/(N-Nv+1-n0); % a_k(sigma^2)
%      bv(k+1)=sum(a(n0+1:N-Nv+1).*b(n0+1+k:N-Nv+1+k)-b(n0+1:N-Nv+1).*a(n0+1+k:N-Nv+1+k))*pi/(N-Nv+1-n0); % b_k(sigma^2)
     
    %av(k+1)=sum(a(n0+1:N-k+1).*a(n0+1+k:N+1))*2*pi/(N+1-n0); % a_k(sigma^2)
    %bv(k+1)=sum(a(n0+1:N-k+1).*b(n0+1+k:N+1))*2*pi/(N+1-n0); % b_k(sigma^2)
end

