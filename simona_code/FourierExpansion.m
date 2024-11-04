function [Vol] = FourierExpansion(M,t,x,K,indici,n,n0,delta)

%     Roberto Reno', 2002
%     Returns the expansion Fourier-Fejer with coefficients a,b in the
%     points t as Vol. Coefficients are for [0,2pi].
%     M is the number of coefficients, the convention is
%     a_k = a(k+1) , k =0,M
 

  N=min(M+1,(K-1)/2); % numero coefficienti del prezzo (ciò serve per non avere effetto di aliasing)
  %N=(K-1)/2; % numero coefficienti del prezzo (ciò serve per non avere effetto di aliasing)
  M=N-n0; % maximum M allowed
%    Nv=N-n0; % numero coefficienti di Fourier
  %M=100; % maximum M allowed
  Nv=M; % numero coefficienti di Fourier
  %N=M+1
  [a,b] = sigma_coefficients(K,x,t,N,Nv,n0); % calcola Nv coefficienti di Fourier della volatilità mediante N coefficienti del prezzo

  for i=1:n
      Vol(i) = a(1);
      %tt=t(indici(i)); Vol(i)=Vol(i)+sum(a(2:M+1).*cos([1:M]*tt)+b(2:M+1).*sin([1:M]*tt)); % Fourier
      %tt=t(indici(i)); Vol(i)=Vol(i)+sum((1-[1:M]/M).*(a(2:M+1).*cos([1:M]*tt)+b(2:M+1).*sin([1:M]*tt))); % Fourier-Fejer
      tt=t(i); Vol(i)=Vol(i)+sum((sin(delta*[1:M])./(delta*[1:M])).^2.*(a(2:M+1).*cos([1:M]*tt)+b(2:M+1).*sin([1:M]*tt))); % Variant of Fourier-Fejer
  end
      
