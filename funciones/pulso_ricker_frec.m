
% pulso de ricker en el dominio del tiempo y en términos de la frecuencia

% IN:
% nt = numero de puntos en el tiempo
% fk = frecuencia (en Hz.) del pulso de ricker
% dt = espaciamiento temporal
% ts = fase o tiempo de llegada
% amp = amplitud de la forma de onda del pulso

function s = pulso_ricker_frec(nt,fk,dt,ts,amp)

    t=(0:(nt-1))*dt; % t(i+1)=t0+i*dt, i=0,1,...,(N-1)
    tp=1/fk;% Periodo mínimo = Altas frecuencias
    %ts=tp;%tp;% Tiempo de llegada
    %amp=1; % amplitud
    
    
    argp=pi*((t-ts)/tp);
    s=-(argp.^2 - 0.5).*exp(-(argp.^2)); % PULSO DE RICKER
    
    s=(1/max(s))*s; % normalizamos
    s=(amp)*s;
end