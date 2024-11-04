%% Import the data, extracting spreadsheet dates in Excel serial date format
% [~, ~, raw, dates] = xlsread('EURUSD_Candlestick_15_m_BID_04.01.2016-31.12.2016.xlsx','EURUSD_Candlestick_15_m_BID_04.','A2:F34845','',@convertSpreadsheetExcelDates);
[~, ~, raw, dates] = xlsread('EURUSD_Candlestick_15_m_BID_04.01.2016-31.12.2016.xlsx','EURUSD_Candlestick_15_m_BID_04.','A1:F34844','',@convertSpreadsheetExcelDates);
[N,P]=size(raw);

% creazione matrici  
dati=NaN(N,P-1);
dati(1:N,1:P-1)=cell2mat(raw(1:N,2:P)); % prezzi e volumi


% conversione date
dat=raw(1:end,1);
%dat([3 6])='/';
dat1=dat(1:end,1);

 tmp=cell2mat(dat(1,1))
 tmp0=datetime(tmp(1:16),'Format','dd.MM.yyyy HH:mm');



for i=1:N
    tmp=cell2mat(dat(i,1));
    tmp=tmp(1:16);
    %disp(i)
    %datacol{i,1}=datetime(tmp,'Format','dd.MM.yyyy HH:mm');
    tmp1=datetime(tmp,'Format','dd.MM.yyyy HH:mm');
    tmp0=[tmp0; tmp1];
end

date=tmp0(2:end,1); % giorno e orario di contrattazione

