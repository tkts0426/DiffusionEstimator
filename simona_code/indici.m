function a = indici(time,lag)

% Calcola gli indici di inizio giorno

a=[1]; rif=time(1);
N=length(time);

for i=1:N-1
%     if ((time(i+1) - rif == lag) || (time(i+1) - rif == 2*lag) || (time(i+1) - rif == 3*lag))
    if (time(i+1) - rif >= lag)
        a=[a;i+1]; rif=time(i+1);
    end
end
a=[a;N];