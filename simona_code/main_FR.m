clear
clc
disp('********************************************************************************')
disp('      Estimate of the spot volatility on N estimation points over ND days       ')
disp('                 by means of Realized Volatility type methods                   ')
disp('                            and Fourier estimator                               ')
disp('********************************************************************************')
disp(' ')

% ND=100; 
% N1g=(24*60*60)*(6/24); % number of seconds in one day, 6 hours trading
% N=N1g*ND; % number of trading times
% T=N/(24*60*60); % Total period of time (expressed in days)
% n_oss_intraday=N1g;
% D=n_oss_intraday;
% %N=D*ND; % number of trading times
%  % Total period of time (expressed in days)

ND=500; n_oss_intraday=100;
N=n_oss_intraday*ND; % number of trading times
T=ND; % Total period of time (expressed in days)
D=n_oss_intraday;

m=200;% Realized Volatility type estimators
D1=10; L=5; M=D1; %Ogawa


% Modello Chan: ---------------------------------------------
    alpha=0.079; beta=0.093; gamma=1.474; eta=0.794; r0=0.03;
    t=[]; p1=[]; sigma=[]; dt=T/N;
    [t,p1,sigma2]=Chan(T,N,r0,alpha,beta,gamma,eta); % simulazione su [0,T]
% -----------------------------------------------------------
% % Modello Heston: ---------------------------------------------
%     dt=T/N;
%     alpha=2.0; beta=0.01; ni=0.1; sigma0=beta; mu=0.0; rho=0.0; % B-R model
%     S0=100; % stock parameters
%     [t,p1,sigma]=Heston(T,N,S0,sigma0,alpha,beta,ni,mu,rho); % simulazione N-F o B-R su [0,T]
%     sigma2=[sigma0;sigma].^2;
% % -----------------------------------------------------------

% Stima della volatilità spot:
%RV_BP=RealVar_BP(p1,N,m,dt);% Realized Variance estimator using the m+1 returns around t_k (m deve essere pari!) (at nodes k+1=t_k, for k=1,2,...,N)
sigma_eps=2*std(diff(p1));
p2=p1+(sigma_eps)^2*randn(size(p1));%observed price 
   % FOURIER: -------------------------------------------
disp('Computes the Fourier expansion all in all (on ND days)')
Nmax=1000000; % Fourier estimator (case without noise)
%delta= 1/25;
%p1=p1'; % vettore colonna
%         sec1F=t*86400; % tempo in secondi
%        % t_k=sec1F(index); % instants of estimation (in seconds)
%         K1=length(p1); % numero di osservazioni 
%         %[sec,p,K]=choose(sec1F,p1(:,i),K1,frequenza); % sceglie solo certe osservazioni in base alla frequenza fissata
%                 n0=1;
%                 % Notazione trigonometrica: % ------------------------------
%                 [newx]=RemoveLinear(K1,p1,sec1F); % remove linear trend
%                 [tau]=normalize(sec1F,K1,2*pi); % normalizzo gli istanti delle osservazioni su [0,2 pi]
%                 if (Nmax > N)
%                     vol1=FourierExpansion(Nmax,tau,newx,K1,p1,N,n0,dt)*(2*pi)/T; % Fourier estimator (riscalo le volatilità sull'intervallo [0,2pi] 
%                                                                                      % (perchè il processo GARCH non è stato creato su [0,2pi]))
%                 else
%                     vol1=FourierExpansionM(Nmax,tau,newx,K1,p1,N,n0,dt)*(2*pi)/T; % Fourier estimator (riscalo le volatilità sull'intervallo [0,2pi] 
%                                                                                      % (perchè il processo GARCH non è stato creato su [0,2pi]))
%                 end
% 
%    IntSqErr=(T/ND)*sum((vol1(1:end)'-sigma2(1:end-1)).^2.)/(N-1);
   %-----------------------------------------------------------------------
  
   x=[0:0.01:0.2]; %h=1.06;S= 0.0626; h=h*S*(N+1)^(-1/5);
h=0.01; 
K=inline('exp(-x.^2/2)/sqrt(2*pi)')

n_traj=2; % numero di traiettorie
for k=1:n_traj
    [t,p1,sigma2]=Chan(T,N,r0,alpha,beta,gamma,eta);
     %RV_FZ1=Realvar_FZ1(p1,dt);
    sec1F=t*86400; % tempo in secondi
       % t_k=sec1F(index); % instants of estimation (in seconds)
        K1=length(p1); % numero di osservazioni 
        %[sec,p,K]=choose(sec1F,p1(:,i),K1,frequenza); % sceglie solo certe osservazioni in base alla frequenza fissata
                n0=1;
                % Notazione trigonometrica: % ------------------------------
                [newx]=RemoveLinear(K1,p1,sec1F); % remove linear trend
                [tau]=normalize(sec1F,K1,2*pi); % normalizzo gli istanti delle osservazioni su [0,2 pi]
                if (Nmax > N)
                    vol=FourierExpansion(Nmax,tau,newx,K1,p1,N,n0,dt)*(2*pi)/T; % Fourier estimator (riscalo le volatilità sull'intervallo [0,2pi] 
                                                                                     % (perchè il processo GARCH non è stato creato su [0,2pi]))
                else
                    vol=FourierExpansionM(Nmax,tau,newx,K1,p1,N,n0,dt)*(2*pi)/T; % Fourier estimator (riscalo le volatilità sull'intervallo [0,2pi] 
                                                                                     % (perchè il processo GARCH non è stato creato su [0,2pi]))
                end
    %RV_FR=Realvar_SO1(p1,N+1,D,L,M,T);
%     D1=1; L=80; M=60;
   % RV_SO3=RealVar_SOF(p1(1:end),N,D1,L,M,dt);
   
   % RV_BP=RealVar_BP(p1,N,m,dt);
   %RV_SO1=Realvar_Ogawa(p1,N,L,n,M,DELTA,DELTAsec,index);
for q=1:length(x)
    
    
   %ker_(k,q)=Kernel_BPReference(p1(1:D:end-1-L*D), RV_SO1',T,h,x(q),K);
   ker_FR(k,q)=Kernel_BPReference(p1(1:D:end-1), vol(1:D:end)',h,x(q),K);
   %ker_FZ1(k,q)=Kernel_BPReference(p1(1:D:end-1),RV_FZ1(1:D:end),h,x(q),K);
    
    
    LT(k,q)=(T/(N*h))*sum(K((x(q)-p1)/h));
end
end
% vol=mean(ker_BP);
vol=mean(ker_FR);
%vol3=mean(ker_FZ1);


localtime=mean(LT); Variance=var(ker_FR); Vhat=vol.^2./((N*h/T)*localtime); 
W=localtime;
%Variance3=var(ker_FZ1);Vhat3=vol3.^2./((N*h/T)*localtime);
for q=1:length(x)
    % Primo modo (di romuald):
   % extinf4(q)=(sqrt(N*h)*W.*vol3(q))/(1.96+sqrt(N*h*W));
    %extsup4(q)=(sqrt(N*h)*W.*vol3(q))/(-1.96+sqrt(N*h*W));
    
    % Secondo modo, Monte Carlo su varie traiettorie:
    %extinf3(q)=vol3(q)-1.96*sqrt(Variance(q))/sqrt(n_traj);
    %extsup3(q)=vol3(q)+1.96*sqrt(Variance(q))/sqrt(n_traj);
    
    % Terzo modo (Teorema 1 pag 624 [Jiang,Knight]):
    extinf(q)=vol(q)-1.96*sqrt(Vhat(q));
    extsup(q)=vol(q)+1.96*sqrt(Vhat(q));
%     extinf1(q)=vol3(q)-1.96*sqrt(Vhat3(q));
%     extsup1(q)=vol3(q)+1.96*sqrt(Vhat3(q));
    
end
% Mean=mean(vol)
% sd=std(vol)
% Meanl=mean(extinf)
% Meanu=mean(extsup)
% sdl=std(extinf)
% sdu=std(extsup)

%--------------------------------------------------------------
 figure(2)
 plot(x,vol,'-')
 hold on
 plot(x,extinf,'red')
 plot(x,extsup,'red')
  plot(x,(eta*x.^gamma).^2,':')
 xlabel('r')
 ylabel('?^2 (r)')
  title('Fourier Extimator ')
%     X=(vol'-sigma([index,N+1])).^2; % CORRETTO???
%     [IntSqErr,DevStd]=MeanVar(X,n,ND,T) % CORRETTO???
                  
%     figure(15)
%     subplot(2,1,1)
%     %plot([1,index],vol,'-')
%     plot(t(1:end-1)*ND/T,vol1,'magenta')
%     subplot(2,1,2)
%     plot(t(1:end-1)*ND/T,sigma2(1:end-1),'red')
%     ylabel('Spot Volatility')
%     xlabel('Day')
    %----------------------------------------------------------------------
%     figure(3)
%  plot(x,vol3,'-')
%  hold on
%  plot(x,extinf1,'red')
%  plot(x,extsup1,'red')
%   plot(x,(eta*x.^gamma).^2,':')
%  xlabel('r')
%  ylabel('?^2 (r)')
%   title('Florens_Zmirou ')
    
    
