h=0.0064/5; % 0.0064/6;

for k=1:n_traj           
            % Calcolo la funzione esatta a numeratore:
            INTEGRALE(k,:)=int_norm_vol(p2(:,k),x,eta,gamma,T,N,K,h);
            DENOM(k,:)=denorm_vol(p2(:,k),x,T,N,K,h);
            INTEGRALE_spot(k,:)=int_norm_vol_spot(p2(:,k),x,sigma2(:,k),T,N,K,h);
end
 
INT_esatto=mean(INTEGRALE);
DEN_esatto=mean(DENOM);

 % Confronto del quoziente:
   figure
   plot(x,mean(INTEGRALE_spot)./DEN_esatto,'-')
   hold on
   plot(x,INT_esatto./DEN_esatto,'r')
   %plot(x,extinf,'red')
   %plot(x,extsup,'red')
   plot(x,(eta*x.^gamma).^2,':')
   xlabel('r')
   ylabel('sigma^2 (r)')
   title('QUOZIENTE reale : esatto in rosso e con spot vol - (curva reale :)')

