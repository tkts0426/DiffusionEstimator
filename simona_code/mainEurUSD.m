clear
clc
disp('********************************************************************************')
disp('         Estimate of the diffusion coefficient on N estimation points           ')
disp('                 by means of Realized Volatility type methods                   ')
disp('                            and Fourier estimator                               ')
disp('********************************************************************************')
disp(' ')

% Questi comandi non caricano bene i tempi di contrattazione.
% DATI='5m'; [data,txt,raw] = xlsread('EURUSD_Candlestick_5_m_BID_29.02.2016-03.06.2016.xlsx',1,''); n_oss_intraday=288; x=[0.08:0.001:0.15]; % 3 months
% % DATI='15m'; [data,txt,raw] = xlsread('EURUSD_Candlestick_15_m_BID_04.01.2016-31.12.2016.xlsx',1,''); n_oss_intraday=240; x=[0.04:0.001:0.14]; % 1 year
% % Di quest'ultimo file la colonna dei volumi è in formato cella e non va bene:
% % DATI='1m'; [data,txt,raw] = xlsread('EURUSD_Candlestick_1_m_BID_04.01.2016-31.12.2016.xlsx',1,''); n_oss_intraday=1440; x=[0.04:0.004:0.14]; % 1 year

% Usa questi:
DATI='5m'; [~, ~, raw, dates] = xlsread('EURUSD_Candlestick_5_m_BID_29.02.2016-03.06.2016.xlsx','EURUSD_Candlestick_5_m_BID_29.0','A1:F27624','',@convertSpreadsheetExcelDates); n_oss_intraday=288; x=[0.08:0.001:0.15]; % 3 months
% DATI='15m'; [~, ~, raw, dates] = xlsread('EURUSD_Candlestick_15_m_BID_04.01.2016-31.12.2016.xlsx','EURUSD_Candlestick_15_m_BID_04.','A1:F34840','',@convertSpreadsheetExcelDates); n_oss_intraday=96; x=[0.04:0.001:0.14]; % 1 year

disp('Creating data structure...')
[txt,data] = conversion(raw);

switch DATI
    case '5m'
        lag=txt(289)-txt(1); % one-day duration
    case '15m'
        lag=txt(97)-txt(1); % one-day duration
end        

Volume=data(:,5); N=length(Volume); ind=[];
for i=1:N
    if (Volume(i) > 0)
        ind=[ind,i];
    end
end
P=log(data(ind,1)*1e-5); % open price
% P=log(data(ind,4)*1e-5); % close price
Time=txt(ind); Volume=Volume(ind); 
N=length(P)-1;

dayindex = indici(Time,lag); % Calcola gli indici di inizio giorno 
% OSSERVAZIONE: non noicidono sempre con la mezzanotte, a volte sono le 23:00
% o le 22:00.

% Analisi dei dati: -------------------------------------------------------
r=diff(P); % log-returns
figure, subplot(2,2,1), plot(r), ylabel('log-return'),  axis([1,length(r),-0.01,0.01]), %title('5 sec. return')
subplot(2,2,2), autocorr(r,[])

fprintf('Total number of quotes: %i \n', length(P))
fprintf('EUR/USD exchange rate: Mean: %f    Sdt. Dev.: %f    Min: %i    Max: %i  \n', mean(exp(P)),std(exp(P)),min(exp(P)),max(exp(P)))
fprintf('log-return (%%): Mean: %i    Sdt. Dev.: %i    Min: %i    Max: %i  \n', mean(r)*100,std(r)*100,min(r)*100,max(r)*100)
[ACF,lags,bounds] = autocorr(r,[])
% -------------------------------------------------------------------------

T=1; L=floor(N/2); 
K=inline('exp(-x.^2/2)/sqrt(2*pi)');
h1=3; S= std(P); h=h1*S*N^(-1/5); % bandwidth parameter (kernel estimator)
[vol,Fourier_int,LR]=japan_coefficientFFT1(P,x,h,K,L,T);

switch DATI
    case '5m'
        figure; 
        subplot(2,1,1), hist(P,200); title('EUR/USD 5-minute 29.02.2016-03.06.2016'), axis([0.07,0.16,0,600])
        subplot(2,1,2), plot(x,vol),axis([0.07,0.16,1e-3,5e-3])
        xlabel('r'), ylabel('sigma^2 (r)')
        title('Volatility function estimate: Fourier (-) and Fourier 2 (--)')
    case '15m'
        figure; 
        subplot(2,1,1), hist(P,200); title('EUR/USD 1-minute 04.01.2016-31.12.2016'), axis([0.02,0.16,0,500])
        subplot(2,1,2), plot(x,vol), axis([0.02,0.16,0.006,0.012])
        xlabel('r'), ylabel('sigma^2 (r)')
        title('Volatility function estimate: Fourier (-) and Fourier 2 (--)')
    case '1m'
        figure; 
        subplot(2,1,1), hist(P,200); title('EUR/USD 1-minute 04.01.2016-31.12.2016'), axis([0.03,0.15,0,400])
        subplot(2,1,2), plot(x,vol), axis([0.03,0.15,0.003,0.03])
        xlabel('r'), ylabel('sigma^2 (r)')
        title('Volatility function estimate: Fourier (-) and Fourier 2 (--)')
end

% -------------------------------------------------------------------------
% ALTRI STIMATORI: --------------------------------------------------------

% disp('OTHER KERNEL ESTIMATORS:')
% disp('Florens-Zmirou estimator - with daily data...')
% P_daily=P(1:n_oss_intraday:end); % daily data
% ND=length(P_daily)-1;
% volFZ=FZ(P_daily,h,x,K,T,ND);
% 
% disp('Realized Volatility based estimator - with intraday data...')
% volRV=RV2(P,h,x,K,T,ND,n_oss_intraday);
% 
% disp('Fourier based estimator day-by-day - with intraday data...')
% t=[0:1/(length(Time)-1):1]';
% volF=FEday2(P,h,x,K,t,T,ND,n_oss_intraday);
% 
% disp('Fourier based estimator all-in-all - with intraday data...')
% volF2=FEwhole2(P,h,x,K,t,T,ND,n_oss_intraday,N);

disp('OTHER KERNEL ESTIMATORS:')
disp('Florens-Zmirou estimator - with daily data...')
P_daily=P(dayindex); % daily data
ND=length(P_daily)-1;
volFZ=FZ(P_daily,h,x,K,T,ND);

disp('Realized Volatility based estimator - with intraday data...')
volRV=RV3(P,h,x,K,T,ND,dayindex);

disp('Fourier based estimator day-by-day - with intraday data...')
t=[0:1/(length(Time)-1):1]';
volF=FEday3(P,h,x,K,t,T,ND,dayindex);

disp('Fourier based estimator all-in-all - with intraday data...')
volF2=FEwhole3(P,h,x,K,t,T,ND,dayindex,N);

hold on, plot(x,volF2,'--'),

figure
   plot(x,volFZ,'-')
   hold on
   plot(x,volRV,'--')
   plot(x,volF,'-o')
   plot(x,volF2,'-*')
%    plot(x,vol(L,:)+1.96*sqrt(Variance(L,:)/n_traj),'r')
%    plot(x,vol(L,:)-1.96*sqrt(Variance(L,:)/n_traj),'r')
    xlabel('r')
 ylabel('sigma^2 (r)')
  title('Florens-Zmirou (-) Realized Vol. (--) Fourier d. (-o) Fourier w. (-*)')

figure
estimates=[vol;volFZ;volRV;volF;volF2];
plot(x,estimates)
xlabel('r')
ylabel('sigma^2 (r)')

figure
estimates=[vol;volF2];
plot(x,estimates)
xlabel('r')
ylabel('sigma^2 (r)')
