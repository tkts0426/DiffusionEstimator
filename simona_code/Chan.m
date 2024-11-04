function [t,r,sigma2]=Chan(T,N,r0,alpha,beta,gamma,eta)

%clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulazione traiettoria moto browniano geometrico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ND=1;
% N1g=(24*60*60)*(6/24); % number of seconds in one day, 6 hours trading
% N=N1g*ND; % number of trading times
% T=N/(24*60*60); % Total period of time (expressed in days)

% N=5000; % number of trading times
% dt=1/252;
% T=N*dt;

% alpha=0.079; beta=0.093; gamma=1.474; eta=0.794; r0=0.065;
    t=[]; p1=[]; sigma=[];


dt=T/N; r=zeros(N+1,1); t=zeros(N+1,1); r(1)=r0;

    for i=1:N
        t(i+1)=t(i)+dt;
        r(i+1)=r(i)+beta*(alpha-r(i))*dt+eta*r(i)^gamma*sqrt(dt)*randn;
    end
    sigma2=(eta*r.^gamma).^2;

% %  subplot(2,1,1)
% % plot(t,sigma2)
% %plot(sigma2)
%  xlabel('t')
%  ylabel('sigma2')
%  title('Chan model')

%garchplot(Innovations, sigma, r)
