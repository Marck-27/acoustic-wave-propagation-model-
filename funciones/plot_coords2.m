function [] = plot_coords2( x,z,xp,zp,xst,zst,nbc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    nx0=nbc+1; % posicion real de x0
    nxf=length(x)-nbc; % posicion real de xf
    

    figure;
    zh1=z(zst(1)); % haltura z con respecto a los geófonos
    hold on
    view(0,-90)
    plot(x,zh1,'.k')% Eje x con fronteras CPML
    plot(x(1:nx0-1),zh1,'og')% Frontera Izquierda CPML
    plot(x(nxf+1:end),zh1,'og')% Frontera Derecha CPML
    plot(x(xp),z(zp),'*r','LineWidth',10)% Coordenadas de las fuentes
    plot(x(xst),z(zst),'ob')% Coordenadas de los geófonos
    title('Coordenadas de CPML, Fuentes y Receptores')
    xlabel('Eje X (m)')
    ylabel('Eje Z (m)')
    axis([x(1),x(end) z(1) z(end)])
    grid on


end

