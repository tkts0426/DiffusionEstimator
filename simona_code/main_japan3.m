clear
clc
disp('********************************************************************************')
disp('      Estimate of the spot volatility on N estimation points over ND days       ')
disp('                 by means of Realized Volatility type methods                   ')
disp('                            and Fourier estimator                               ')
disp('********************************************************************************')
disp(' ')

randn('seed',1)

% % FUNZIONA MALE (TROPPO POCHI DATI): ----------------------------------------------------
% n_traj=10; % numero di traiettorie
% fprintf('Computes the Fourier expansion on 1 day, for a total of n_traj = %i trajectories - intraday data \n',n_traj)
% ND=1; n_oss_intraday=3600; % 1 minute data
% N=n_oss_intraday*ND; % number of trading times
% T=ND; % Total period of time (expressed in days)
% D=1; % no sparse sampling
% %h=0.0064/5; 
% x=[0:0.01:0.1]; 
% % --------------------------------------------------------------------------

% % % VA UN PO' MEGLIO (MA ANCORA TROPPO POCHI DATI): ----------------------------------------------------
% n_traj=10; % numero di traiettorie
% fprintf('Computes the Fourier expansion on 1 day, for a total of n_traj = %i trajectories - intraday data \n',n_traj)
% ND=1; n_oss_intraday=7200; % 30 seconds frequency
% N=n_oss_intraday*ND; % number of trading times
% T=ND; % Total period of time (expressed in days)
% D=1; % no sparse sampling
% %h=0.0064/5; 
% x=[0:0.01:0.1]; 
% % --------------------------------------------------------------------------

% FUNZIONA MOLTO BENE: ----------------------------------------------------
n_traj=1; % numero di traiettorie
ND=250; n_oss_intraday=100; % (ND=500; D=n_oss_intraday/4;)
fprintf('Computes the Fourier expansion all in all (on %i days) - intraday data \n', ND)
N=n_oss_intraday*ND; % number of trading times (complessivo)
T=ND; % Total period of time (expressed in days)
D=1; % D=n_oss_intraday/2; % serve per eventuale sparse sampling (altrimenti D = 1)
%h=0.064/5; 
x=[0.04:0.01:0.5]; 
% --------------------------------------------------------------------------

% % FUNZIONA MOLTO BENE: stesso numero di dati del caso precedente, ma pensato su un giorno e non 500-------
% n_traj=10; % numero di traiettorie
% fprintf('Computes the Fourier expansion on 1 day, for a total of n_traj = %i trajectories - intraday data \n',n_traj)
% ND=1; n_oss_intraday=50000;
% N=n_oss_intraday*ND; % number of trading times
% T=ND; % Total period of time (expressed in days)
% D=25; % serve per eventuale sparse sampling (altrimenti D = 1)
% %h=0.064/5; 
% x=[0.03:0.01:0.1]; % BISOGNA RESTRINGERE L'INTERVALLO, PERCHE' IL RANGE DEI PREZZI E' LIMITATO!
% % --------------------------------------------------------------------------

disp('...')
   % Stima del coefficiente di diffusione usando tutti i prezzi intraday:-------------
   %h1=1.06; S= std(p1); h2=h1*S*(N+1)^(-1/5)
   L=(N/D)/2; % maximum Fourier frequency
   K=inline('exp(-x.^2/2)/sqrt(2*pi)');
    
        t=zeros(N+1,n_traj); p1=zeros(N+1,n_traj); sigma2=zeros(N+1,n_traj);
        vol=zeros(L,length(x)); Variance=zeros(L,length(x));
        INTEGRALE=zeros(n_traj,length(x)); INT_stimato=zeros(L,length(x));
        DENOM=zeros(n_traj,length(x)); DEN_stimato=zeros(1,length(x));
        INTEGRALE_spot=zeros(n_traj,length(x));
        ker_FR=zeros(L,length(x)); Fourier_int=zeros(L,length(x));

% Modello Chan: ---------------------------------------------
    alpha=0.079; beta=0.093; gamma=1.474; eta=0.794; r0=0.065;
    t=[]; p1=[]; sigma=[]; dt=T/N;
    %[t,p1,sigma2]=Chan(T,N,r0,alpha,beta,gamma,eta); % simulazione su [0,T]
% -----------------------------------------------------------
% % Modello Heston: ---------------------------------------------
%     dt=T/N;
%     alpha=2.0; beta=0.01; ni=0.1; sigma0=beta; mu=0.0; rho=0.0; % B-R model
%     S0=100; % stock parameters
%     [t,p1,sigma]=Heston(T,N,S0,sigma0,alpha,beta,ni,mu,rho); % simulazione N-F o B-R su [0,T]
%     sigma2=[sigma0;sigma].^2;
% % -----------------------------------------------------------
obs_prices=[];
for k=1:n_traj
    [t(:,k),p1(:,k),sigma2(:,k)]=Chan(T,N,r0,alpha,beta,gamma,eta);
%   sigma_eps=2*std(diff(p1(:,k)));
    sigma_eps=0; % microstructure effect
    p2(:,k)=p1(:,k)+(sigma_eps)^2*randn(size(p1(:,k)));%observed price
    obs_prices=[obs_prices;p2(:,k)];
end
figure; subplot(2,1,1), hist(obs_prices,200); title('Observed prices'), %axis([0,0.5,0,6e4])

h1=3; S= std(p2); h2=h1*S*N^(-1/5); % bandwidth parameter (kernel estimator)
% h1=3; S= std(p2); h2=h1*S*(N/D)^(-1/5); % bandwidth parameter (kernel estimator)

 for k=1:n_traj           
     h=h2(k); % bandwidth parameter
     
            % Calcolo la funzione esatta a numeratore:
            INTEGRALE(k,:)=int_norm_vol(p2(:,k),x,eta,gamma,T,N,K,h);
            DENOM(k,:)=denorm_vol(p2(:,k),x,T,N,K,h);
            INTEGRALE_spot(k,:)=int_norm_vol_spot(p2(:,k),x,sigma2(:,k),T,N,K,h);
            
            % Stima Fourier:
             K1=length(p2(:,k)); %sec1F=t(:,k)*86400;
            newx=p2(:,k); %[newx]=RemoveLinear(K1,p2(:,k),sec1F); % remove linear trend
            tau=t(:,k); %[tau]=normalize(sec1F,K1,2*pi); % normalizzo gli istanti delle osservazioni su [0,2 pi]
                       
%             [ker_FR,Fourier_int,LR]=japan_coefficient3(newx(1:D:end),x,tau(1:D:end),h,K,L,T);
            [ker_FR(L,:),Fourier_int(L,:),LR]=japan_coefficientFFT1(newx(1:D:end),x,h,K,L,T);
            vol=vol+ker_FR; INT_stimato=INT_stimato+Fourier_int;
            DEN_stimato=DEN_stimato+LR;
            Variance=Variance+ker_FR.^2;
            
%  
 end
         
%--------------------------------------------------------------

% Intervallo confidenza naive:  
vol=vol/n_traj; INT_stimato=INT_stimato/n_traj; DEN_stimato=DEN_stimato/n_traj;
Variance=(Variance-n_traj*vol.^2)/(n_traj-1); % varianza campionaria 
   
figure
   plot(x,vol(L,:),'-')
   hold on
   %plot(x,extinf,'red')
   %plot(x,extsup,'red')
   plot(x,(eta*x.^gamma).^2,':')
   plot(x,vol(L,:)+1.96*sqrt(Variance(L,:)/n_traj),'r')
   plot(x,vol(L,:)-1.96*sqrt(Variance(L,:)/n_traj),'r')
    xlabel('r')
 ylabel('sigma^2 (r)')
  title('Fourier estimator. Volatility function: real (:) and estimated (-)')
  
% % Altro intervallo confidenza [Jiang and Knight] su ultima traiettoria:  
%  %localtime=mean(LR); Variance=var(ker_FR); Vhat=vol.^2./((ND*h/T)*localtime); % SBAGLIATO!
%  Vhat=ker_FR(L,:).^2./((N/D)*h2(k)*LR);
% 
% figure
%    plot(x,ker_FR(L,:),'-')
%    hold on
%    plot(x,(eta*x.^gamma).^2,':')
%    plot(x,ker_FR(L,:)+1.96*sqrt(Vhat),'r')
%    plot(x,ker_FR(L,:)-1.96*sqrt(Vhat),'r')
%     xlabel('r')
%  ylabel('sigma^2 (r)')
%   title('ULTIMA TRAIETTORIA. Funzione volatilità esatta : e stimata -')

 % Confronto dei numeratori:
 if (n_traj == 1)
     INT_esatto=INTEGRALE;
 else
     INT_esatto=mean(INTEGRALE);
 end
  figure
 plot(x,INT_esatto,':')
 hold on
 plot(x,INT_stimato(L,:),'-')
 title('NUMERATOR exact (:) e estimated (-)')
 
  % Confronto dei denominatori:
  if (n_traj == 1)
     DEN_esatto=DENOM;
 else
     DEN_esatto=mean(DENOM);
 end
  figure
 plot(x,DEN_esatto,':')
 hold on
 plot(x,DEN_stimato,'-')
  title('DENOMINATORE esatto : e stimato -')

  % Confronto del quoziente:
  figure
 plot(x,INT_esatto./DEN_esatto,'--')
 hold on
 plot(x,INT_stimato(L,:)./DEN_stimato,'-')
 plot(x,(eta*x.^gamma).^2,':')
 title('QUOZIENTE esatto -- e stimato -')
 
%  % Confronto del quoziente:
%  figure
%  plot(x,mean(INTEGRALE_spot)./DEN_esatto,'-')
%    hold on
%    plot(x,INT_esatto./DEN_esatto,'r') % coincide con quella sopra
%    %plot(x,extinf,'red')
%    %plot(x,extsup,'red')
%    plot(x,(eta*x.^gamma).^2,':')
%    xlabel('r')
%    ylabel('sigma^2 (r)')
%    title('QUOZIENTE reale : esatto in rosso e con spot vol - (curva reale :)')

% -------------------------------------------------------------------------
% ALTRI STIMATORI: --------------------------------------------------------
disp('OTHER KERNEL ESTIMATORS:')
disp('Florens-Zmirou estimator - with daily data...')
volFZ=zeros(1,length(x));
for k=1:n_traj           
    h=h2(k); % bandwidth parameter
    P=p2(1:n_oss_intraday:end,k); % daily data
    FlZm=FZ(P,h,x,K,T,ND);
    volFZ=volFZ+FlZm;
end
volFZ=volFZ/n_traj;

disp('Realized Volatility based estimator - with intraday data...')
volRV=zeros(1,length(x));
for k=1:n_traj           
    h=h2(k); % bandwidth parameter
    P=p2(:,k); % high-frequency data
    RealVol=RV(P,h,x,K,T,ND,n_oss_intraday);
    volRV=volRV+RealVol;
end
volRV=volRV/n_traj;

disp('Fourier based estimator day-by-day - with intraday data...')
volF=zeros(1,length(x));
for k=1:n_traj           
    h=h2(k); % bandwidth parameter
    P=p2(:,k); % high-frequency data
    Fourier=FEday(P,h,x,K,t(:,k),T,ND,n_oss_intraday);
    volF=volF+Fourier;
end
volF=volF/n_traj;

disp('Fourier based estimator all-in-all - with intraday data...')
volF2=zeros(1,length(x));
for k=1:n_traj           
    h=h2(k); % bandwidth parameter
    P=p2(:,k); % high-frequency data
    Fourier2=FEwhole(P,h,x,K,t(:,k),T,ND,n_oss_intraday,N);
%     Fourier2=FEwholeFFT(P,h,x,K,t(:,k),T,ND,n_oss_intraday,N); % non del tutto corretto!
    volF2=volF2+Fourier2;
end
volF2=volF2/n_traj;

figure
   plot(x,volFZ,'-')
   hold on
   plot(x,volRV,'--')
   plot(x,volF,'-o')
   plot(x,volF2,'-*')
   plot(x,(eta*x.^gamma).^2,':')
%    plot(x,vol(L,:)+1.96*sqrt(Variance(L,:)/n_traj),'r')
%    plot(x,vol(L,:)-1.96*sqrt(Variance(L,:)/n_traj),'r')
    xlabel('r')
 ylabel('sigma^2 (r)')
  title('Florens-Zmirou (-) Realized Vol. (--) Fourier d. (-o) Fourier w. (-*)')

figure
estimates=[vol(L,:);volFZ;volRV;volF;volF2];
plot(x,estimates)
hold on
plot(x,(eta*x.^gamma).^2,'black')
   plot(x,vol(L,:)+1.96*sqrt(Variance(L,:)/n_traj),':')
   plot(x,vol(L,:)-1.96*sqrt(Variance(L,:)/n_traj),':')
xlabel('r')
ylabel('sigma^2 (r)')
