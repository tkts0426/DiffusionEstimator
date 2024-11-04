function [Vol] = FourierExpansion(M,t,x,K,indici,n,n0)

%     Roberto Reno', 2002
%     Returns the expansion Fourier-Fejer with coefficients a,b in the
%     points t as Vol. Coefficients are for [0,2pi].
%     M is the number of coefficients, the convention is
%     a_k = a(k+1) , k =0,M
%     Vol is a M+1,n matrix containing al the partial sums of the Fourier
%     expansion up to M.
 
  N=min(M+1,(K-1)/2); % numero coefficienti del prezzo (ciò serve per non avere effetto di aliasing)
  %N=(K-1)/2; % numero coefficienti del prezzo (ciò serve per non avere effetto di aliasing)
%   M=N-n0; % maximum M allowed
%   Nv=N-n0; % numero coefficienti di Fourier
  M=N-n0; % maximum M allowed
  Nv=M; % numero coefficienti di Fourier
  %N=M+1
  [a,b] = sigma_coefficients(K,x,t,N,Nv,n0); % calcola Nv coefficienti di Fourier della volatilità mediante N coefficienti del prezzo

  for i=1:n
      Vol(1,i) = a(1);
      for m=1:M
          %tt=t(indici(i)); Vol(i)=Vol(i)+sum(a(2:M+1).*cos([1:M]*tt)+b(2:M+1).*sin([1:M]*tt)); % Fourier
          tt=t(i); Vol(m+1,i)=Vol(m,i)+(1-m/M)*(a(m+1)*cos(m*tt)+b(m+1)*sin(m*tt)); % Fourier-Fejer
      end
  end
      
