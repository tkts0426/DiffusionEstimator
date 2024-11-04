function [tau]=normalize(sec,K,b)

% Normalizza gli istanti delle osservazioni su [0,b]

tau=zeros(K,1);
for j=1:K
  tau(j)=b*(sec(j)-sec(1))/(sec(K)-sec(1));
end
