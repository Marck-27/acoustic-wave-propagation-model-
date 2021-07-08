% función que aproxima la frecuencia dominate (o máxima) de una señal
% temporal.
% Esta frecuencia esta definida como el punto en el dominio de frecuencia
% angular, en donde el espectro de amplitud alcanza su máximo.

function fd=frecDominant(dt,signal)
    % Entrada: 'dt' y 'signal'
    % Salida: fd = frecuencia(angular) dominante
    
    N=length(signal);
    NQ=round(N/2+1);
    
    F=1/dt; %frecuencia lineal
    df=F/N; %espaciemiento dominio de frecuencia lineal
    f=(0:N-1)*df;% Dominio discreto de frecuencia lineal
%    w=2*pi*f; % Dominio discreto de frecuencia angular
    
    F_signal=dt*fft(signal); %Tranformada de Fourier (numérica) de signal(t)
%     plot(w,real(F_signal))% Grafica 
%     plot(w(1:NQ),real(F_signal(1:NQ)))% Grafica hasta NQ
%     title('Tranformada de Fourier')
%     xlabel('frecuencias w')
    
    Espect=abs(F_signal); % Espectro de amplitud o frecuencias
%     plot(w,Espect); % Grafica
%     plot(w(1:NQ),Espect(1:NQ))% Grafica hasta NQ
%     title('Espectro de frecuencias')
%     xlabel('frecuencias w')
    
    fmax=max(Espect(1:NQ));
    for i=1:NQ
        if Espect(i)==fmax
            k=i;
        end
    end
       
%    fd=w(k); % frecuencia domianate angular 'w' (rad/seg)
    fd=f(k); % frecuencia domianate lineal 'f' (Hz)

end
