% CASO SH con esquema de velocidad-esfuerzos y fronteras CPML:

% PARÁMETROS DE ENTRADA:
% vel(:,:) -> matriz 2D del modelo de velocidades sin considerar fronteras absorbentes
% nbc  -> # de nodos para fronteras absorbentes
% Nt   -> # de iteraciones en el tiempo
% dt   -> espaciamiento en el dominio del tiempo
% s(:) -> vector de valores para la onduleta de la fuente o explosión 
% x0   -> coordenada  inicial 'x' para el dominio de interés sin absorbencia
% dx   -> espaciamiento en el dominio del espacio (dz=dx)
% nx   -> # de nodos del eje 'x' (sin contar nodos absorbentes)
% z0   -> coordenada  inicial 'z' para el dominio de interés sin absorbencia
% nz   -> # de nodos del eje 'z' (sin contar nodos absorbentes)
% xp   -> nodo o posición en el eje x para la coordenada x de la fuente (dentro de la región sin absorbencia)
% zp   -> nodo o posición en el eje z para la coordenada z de la fuente (dentro de la región sin absorbencia)
% xst(:)  -> vector de nodos de las coordenadas x para los receptores (que guardan las trazas)
% zst(:)  -> vector de nodos de las coordenadas z para los receptores (que guardan las trazas)
% isFS -> quitar o poner superficie libre (false=No, true=Sí)
% animar -> muestra la animación de la propagación de la onda (false=No, true=Sí)
% sismograma -> muestra el sismograma final de la propagación de la onda (false=No, true=Sí)
% Rc <--- coeficiente teórico de reflexión para las CPML


% PARÁMETROS DE SALIDA:
% seis(:,:) -> matriz 2D de las trazasa registradas a travez del tiempo (sismograma)
%             (registro de una sóla línea de receptores paralela al eje x en la parte superior del domio de interés)

function [seis_Vy]=PS_Acustic_trazas(vel,nbc,Nt,dt,s,x0,dx,nx,z0,nz,xp,zp,xst,zst,isFS,animar,Rc)

%% Creamos video
crea_video=0; % 0<--- No   1<--- Si

if crea_video==1
    my_video=VideoWriter('video.avi');% Nombre del video
    my_video.FrameRate=10;
    my_video.Quality=100; % calidad del video
    open(my_video); % inicia grabació del video (HAY QUE CERRARLO MÁS ADELANTE)
end
%% 
%clear; close all; clc

%% Configuración de colores para gáficos
seiscolor=load('colorbar.txt');
if animar
    % PARAMETROS PARA FIJAR ESCALA DE COLORES:
    k_color=1e-5;% calibrador color 
    minVel=min(min(k_color*vel));
    maxVel=max(max(k_color*vel));
    if (minVel ~=   maxVel)% fijamos límites de color            
        color_min=minVel;
        color_max=maxVel;
    else
        color_min=minVel;
        color_max=1.1*minVel;                         
    end 
end
%%

%addpath((genpath(pwd))) 

%isFS=true; % ponemos superficie libre (false=No, true=Sí)
%animar=true;% se muestra la animación  (false=No, true=Sí)
%sismograma=true;% se muestra el sismograma  (false=No, true=Sí)

%% Numero de puntos en el dominio del tiempo (Número de iteraciones)
%Nt=2001;
%dt=0.0019;

%% ********* Carga de la onduleta para la fuente o explosión **************
%load('pulso')
% % s(:) -> Vector de tamaño 4001 X 1 que contiene la onduleta de la fuente

%% Datos de la malla
%nbc=20;% # nodos absorbentes CPML en las 4 fronteras

% ****************************** EJE X ******************************
%x0=0;% coordenada inicial real sin CPML en x
%dx=10;% Espaciemiento en x
%nx=151;%401;% Numero total de nodos en x sin absorbencia
Nx=nx+2*nbc;% Numero total de nodos en x, incluyendo nodos con absorbencia
nxCPML=nbc;% # nodos absorbentes con CPML en x
x=(x0-nxCPML*dx):dx:(x0+(Nx-nxCPML-1)*dx); % malla 

nx0=nxCPML+1; % posicion real de x0
nxf=Nx-nxCPML; % posicion real de xf
% % Zona de absorbencia izquierda: x(1):...:x(nx0-1)
% % Zona sin absorbencia: x(nx0):...:x(nxf)
% % Zona de absorbencia derecha: x(nxf+1):...:x(Nx)
xf=x(nxf);% coordenada final real sin CPML en x

% ****************************** EJE Z ******************************
%z0=0;% coordenada inicial real sin CPML en z
dz=dx;% Espaciemiento en z
%nz=151;% Numero total de nodos en z sin absorbencia
Nz=nz+2*nbc;% Numero total de nodos en z, incluyendo nodos con absorbencia
nzCPML=nxCPML;% # nodos absorbentes con CPML en z
z=(z0-nzCPML*dz):dz:(z0+(Nz-nzCPML-1)*dz); % malla


nz0=nzCPML+1; % posicion real de z0
nzf=Nz-nzCPML; % posicion real de zf
% % Zona de absorbencia superior: z(1):...:z(nz0-1)
% % Zona sin absorbencia: z(nz0):...:z(nzf)
% % Zona de absorbencia inferior: z(nzf+1):...:z(Nz)
zf=z(nzf);% coordenada final real sin CPML en z

%% Construcción dominio de número de onda (kx,kz)
% CONSTRUCCION DE SIMETRIA PARA LA FFT EN 2D:
% Dado que se esta aplicando la Trans. Fourier, entonces 
% F{Dx_v}=ikx*F{v} , F{Dz_v}=ikz*F{v}, etc. así que 
% construiremos  'kx' y 'kz' de tal forma que sean simétricos con respecto
% a la posición de las frecuencias de Nyquist (NQx,NQz)


% Discretización del dominio de numero de onda 'kx':
fx=1/dx; % fecuencia lineal en x
%fxMAX=fx/2; % frecuencia lineal maxima
dfx=fx/Nx;
if mod(Nx,2)==0
    NQx=round(Nx/2+1); % posicion frecuencia de Nyquist en x    
    fx1D=(0:NQx-1)*dfx; % frecuencia lineal discreta hasta NQx    
else
    NQx=round(Nx/2); % posicion frecuencia de Nyquist en x    
    fx1D=(0:NQx)*dfx; % frecuencia lineal discreta hasta NQx
end
kx=2*pi*fx1D; % Dominio discreto de numero de onda hasta NQx 
kx1D=[kx -kx(NQx-1:-1:2)]; % Dominio de numero de onda kx 'simétrico' en 1D

% % simetria de kx en 2D
% kx2D=zeros(Nz,Nx);
% for iz=1:Nz
%     kx2D(iz,:)=kx1D; 
% end


% Discretización del dominio de numero de onda 'kz':
fz=1/dz; % fecuencia lineal en z
%fzMAX=fz/2; % frecuencia lineal maxima
dfz=fz/Nz;
if mod(Nz,2)==0
    NQz=round(Nz/2+1); % posicion frecuencia de Nyquist en z    
    fz1D=(0:NQz-1)*dfz; % frecuencia lineal discreta hasta NQz    
else
    NQz=round(Nz/2); % posicion frecuencia de Nyquist en z    
    fz1D=(0:NQz)*dfz; % frecuencia lineal discreta hasta NQz
end
kz=2*pi*fz1D; % Dominio discreto de numero de onda hasta NQz 
kz1D=[kz -kz(NQz-1:-1:2)]; % Dominio de numero de onda kz 'simétrico' en 1D

% % simetria de kz en 2D
% kz2D=zeros(Nz,Nx);
% for ix=1:Nx
%     kz2D(:,ix)=(kz1D)'; 
% end

%% ******************* Carga de modelo de velocidades *******************
%load('mod_vel_marmousi2_151_401_10m')
% vel : Matriz (151 X 401) del modelo de velocidad VERDADERO (sin considerar nodos con absorbencia)

%vel=max(max(vel))*ones(size(vel));%Prueba con 'v' homogénea 

% AJUSTAMOS LAS DIMENSIONES DEL MODELO DE VELOCIDADES PARA PONER FRONTERAS ABOSRBENTES
%vel=[repmat(vel(:,1),1,nbc), vel, repmat(vel(:,end),1,nbc)];
%vel=[repmat(vel(1,:),nbc,1); vel; repmat(vel(end,:),nbc,1)];
vel=redim_bound(vel,nbc);

%% Coordenadas de la fuente y receptores

% NODOS DE LAS COORDENADAS DE LA FUENTES (DEBEN ESTAR DENTRO DE LA REGIÓN SIN FRONTERAS ABSORBENTES)
%xp=round(Nx/2);
%zp=nz0+1;%round(Nz/2);

% NODOS DE LAS COORDENADAS DE LAS ESTACIONES (PARALELAS AL EJE X )
%hst=2; % <--- Cada 'hst' nodos se ubicará un receptor a lo largo del eje 'x'
%xst=nx0:hst:nxf;% # del nodo de la coordenada 'x' de los receptores sin nodos absorbentes
%    %xst=round((nx0+nxf)*(1/4)); % Registro de una sólo receptor
%zst=(nz0+1)*ones(1,length(xst));     % # del nodo de la coordenada 'z=0 constante' de los receptores 
N_tr=length(xst);% <--- Numero de trazas a registrar


% NODOS DE LAS COORDENADAS DE LAS ESTACIONES (PARALELAS AL EJE Z )
% hst=2; % <--- Cada 'hst' nodos se ubicará un receptor a lo largo del eje 'z'
% zst=nz0:hst:nzf;% # del nodo de la coordenada 'z' de los receptores sin nodos absorbentes
%     %zst=round((nz0+nzf)*(1/4)); % Registro de una sólo receptor
% xst=(nx0+1)*ones(1,length(zst));     % # del nodo de la coordenada 'x=0 constante' de los receptores 
% N_tr=length(zst);% <--- Numero de trazas a registrar




%% Discretización para el dominio del tiempo
V_S=max(max(vel));

% Condición de estabilidad:
%st=(1/sqrt(2))*(dx/V_S); % para onda acústica
%st=6/(7*sqrt(2))*(dx/V_S); % estabilidad para 4o orden
st=2/(pi*sqrt(2))*(dx/V_S); % estabilidad para en el método pseudo-espectral

if dt >= st
    fprintf('\n Error, dt = %f y debe ser menor que %f \n',dt,st);
    fprintf('        o reducir valores de velocidad \n');
    stop;
end
if dt <= 0
    fprintf('\n Error, dt = %f y debe ser positivo\n',dt);
    stop;
end

%%

% AJUSTAMOS LA ONDULETA 's' DE LA FUENTE A CUALQUIER NUMERO DE PUNTOS
nn0=length(s);
if nn0<Nt
    pulso=zeros(Nt,1);
    pulso(1:nn0)=s;
else
    pulso=s;
end
%plot(t,pulso)

fd=frecDominant(dt,pulso);

%% PARÁMETROS PARA FRONTERAS C-PML:

% parámetros 'a' y 'b' fijos en las 4 fronteras:
[bx_L,ax_L,bx_R,ax_R,bz_T,az_T,bz_B,az_B]=CPML_Acustic_params_ab(Rc,fd,dt,V_S,x,dx,x0,xf,nxCPML,nx0,nxf,z,dz,z0,zf,nzCPML,nz0,nzf);

% inicialiazmos las funciones 'PSI' con ceros:
[psi_V_xL,psi_V_xR,psi_Syx_xL,psi_Syx_xR,psi_V_zT,psi_V_zB,psi_Syz_zT,psi_Syz_zB]=CPML_Acustic_psi_zero(Nz,Nx,nzCPML,nxCPML);

%% ************************************************************************
% Inician iteraciones en el tiempo:
%**************************************************************************
Vy=zeros(size(vel));   % Velocidades iniciales 
%Uy=zeros(size(vel));   % Desplazamientos iniciales
Syx=zeros(size(vel)); % Esfuerzos iniciales
Syz=zeros(size(vel));
seis_Vy=zeros(Nt,N_tr); % Aqui almacenaremos las trazas de cada receptor

% Argumentos de la ifft:
argX=(sqrt(-1)*kx1D);
argZ=(sqrt(-1)*kz1D);

% parámetros para reducir fenómeno de Nyquist
% (Ver artículo de T.Ozdenvar & McMechan 1996)
% (Carcione 1999, Staggered mesh for anisotropic and viscoelastic wave eq)
Ex= sqrt(-1)*( (dx/2)*kx1D );
Ez= sqrt(-1)*( (dz/2)*kz1D );

for it=1:Nt % # nodos en tiempo = # de iteraciones  
    
    %********************* INTRODUCIMOS FUENTE *****************************       
    %Syz(zp,xp)=Syz(zp,xp)-pulso(it);  
    Vy(zp,xp)=Vy(zp,xp)+pulso(it);
   
    % Calculamos los campos de onda espaciales del tiempo 'it':
    [Vy,Syz,Syx,...
    psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB]...
    ...
    =PS_Acustic_Field_2D_cpml...
    ...
    (dt,vel,argX,argZ,Ex,Ez,nbc,Vy,Syz,Syx,...
    ax_L,bx_L,ax_R,bx_R,az_T,bz_T,az_B,bz_B,...
    psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB,isFS);
    
    
   %********** Registro en línea de receptores paralelos al eje x *********
   % En la iteración 'it' del tiempo actual guardamos el valor de U en las
   % coordenadas de cada receptor   
    for ss=1:N_tr %recorremos el # de receptores (= # de trazas)
        % En cada iteración 'it' guardamos el valor de V en la coordenada
        % de la linea de receptores (x( xst(ss) ),z(zst)), con z(zst) fija 
        seis_Vy(it,ss)=Vy(zst(ss),xst(ss));
    end
        
    if animar        
        if mod(it,10)==0            
            figure(20);
            
%             % NO SE MUESTRAN LAS FRONTERAS ABSOR.
%             surf(Vy(nz0:nzf,nx0:nxf)+k_color*vel(nz0:nzf,nx0:nxf) );
%             text(nx,nz,0,['it= ',num2str( it )]);
            
            %SE MUESTRAN LAS FRONTERAS ABSORBENTES
            surf(Vy + k_color*vel ); 
            text(nx+2*nbc,nz+2*nbc,0,['it= ',num2str( it )]);
            
            shading interp;            

            %CONFIGURACIÓN DE POSICIÓN Y EJES DE GRAFICACIÓN
            axis equal;
            view(0,-90)            

            % CONFIGURACIÓN DE COLORES
            %colormap jet;
            colormap(seiscolor);               
            caxis([color_min,color_max]);
            
            %colorbar;
            xlabel('Eje x')
            ylabel('Eje z')          
            
            if crea_video==1
                % Vamos grabando el video: 
                ventana_actual=getframe(gcf);% gcf=get current frame
                writeVideo(my_video,ventana_actual);
            end

            pause(1e-6);
        end
    end   
    
end % end nodos en el tiempo

if crea_video==1
    % Cerramos el video
    close(my_video);
end

end