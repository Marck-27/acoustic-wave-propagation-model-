
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
% xp(:)-> vector de nodos o posiciones en el eje x para la coordenada x de la fuente (dentro de la región sin absorbencia)
% zp(:)-> vector de nodos o posiciones en el eje z para la coordenada z de la fuente (dentro de la región sin absorbencia)
% ns   -> # de fuentes o explosiones
% xst(:)  -> vector de nodos de las coordenadas x para los receptores (que guardan las trazas)
% zst(:)  -> vector de nodos de las coordenadas z para los receptores (que guardan las trazas)
% ng   -> # de geófonos o estaiones 
% isFS -> quitar o poner superficie libre (false=No, true=Sí)
% animar -> muestra la animación de la propagación de la onda (false=No, true=Sí)
% sismograma -> muestra el sismograma final de la propagación de la onda (false=No, true=Sí)


function [s,nt,dt,nbc,x0,dx,nx,z0,nz,xp,zp,ns,xst,zst,ng,isFS,animar,sismograma,Rc]=indata_Acustic_Marmousi


addpath((genpath(pwd)));

isFS=true; % ponemos superficie libre (false=No, true=Sí)
animar=true;% se muestra la animación  (false=No, true=Sí)
sismograma=false;% se muestra el sismograma  (false=No, true=Sí)

%% Numero de puntos en el dominio del tiempo (Número de iteraciones)

% No. de muestras en el tiempo
nt=5000; % for Large scale model

k=1;% Calibrador de espaciamientos 'dt', 'dx' y frecuencias

dt=k*0.001;

%% Datos de la malla
nbc=30;% # nodos absorbentes CPML en las 4 fronteras
Rc=1e-5; % coeficiente teórico de reflexión


% ****************************** EJE X ******************************
% % Cargamos el modelo de marmousi
load('./modelos_vel/vel_marmousi2_model/vel_model.mat');% vel(:,:)
[nz,nx]=size(vel);

x0=0;% coordenada inicial real sin CPML en x
dx=k*10;% Espaciemiento en x
%nx=nx;% Numero total de nodos en x sin absorbencia
Nxx=nx+2*nbc;% Numero total de nodos en x, incluyendo nodos con absorbencia
nxCPML=nbc;% # nodos absorbentes con CPML en x
x=(x0-nxCPML*dx):dx:(x0+(Nxx-nxCPML-1)*dx); % malla con nodos CPML 

nx0=nxCPML+1; % posicion real de x0
nxf=Nxx-nxCPML; % posicion real de xf
% % Zona de absorbencia izquierda: x(1):...:x(nx0-1)
% % Zona sin absorbencia: x(nx0):...:x(nxf)
% % Zona de absorbencia derecha: x(nxf+1):...:x(Nx)
%xf=x(nxf);% coordenada final real sin CPML en x

% ****************************** EJE Z ******************************
z0=0;% coordenada inicial real sin CPML en z
dz=dx;% Espaciemiento en z
%nz=nz;% Numero total de nodos en z sin absorbencia
Nzz=nz+2*nbc;% Numero total de nodos en z, incluyendo nodos con absorbencia
nzCPML=nxCPML;% # nodos absorbentes con CPML en z
z=(z0-nzCPML*dz):dz:(z0+(Nzz-nzCPML-1)*dz); % malla con nodos CPML

nz0=nzCPML+1; % posicion real de z0
%nzf=Nzz-nzCPML; % posicion real de zf
% % Zona de absorbencia superior: z(1):...:z(nz0-1)
% % Zona sin absorbencia: z(nz0):...:z(nzf)
% % Zona de absorbencia inferior: z(nzf+1):...:z(Nz)
%zf=z(nzf);% coordenada final real sin CPML en z

% zz=z(nz0:nzf);% <--- malla del eje z, sin CPML
% save('./Real_Data/Coord_aprox_real_data/mesh_z.mat','zz');
%% Nodos de Coordenadas de las fuentes
% NODOS DE LAS COORDENADAS DE LAS FUENTES (DEBEN ESTAR DENTRO DE LA REGIÓN SIN FRONTERAS ABSORBENTES)
% Vector de nodos de las coordenadas de la fuente o explosión
dis_abs=20;%50;% Numéro de nodos separados de las fornteras absorbentes
h_exp=1;% <--- Cada 'h_exp' nodos se ubicará una fuente a lo largo del eje 'x'
xp=(nx0+dis_abs):h_exp:(nxf-dis_abs);

%zp=(nz0+1)*ones(size(xp)); %<--- fuente en superficie
zp=(nz0+100)*ones(size(xp)); %<--- fuente en profundidad

ns=length(xp); % número de fuentes o explosiones

%% Nodos de Coordenadas de los geófonos
% NODOS DE LAS COORDENADAS DE LAS ESTACIONES (PARALELAS AL EJE X )
hst=1; % <--- Cada 'hst' nodos se ubicará un receptor a lo largo del eje 'x'
xst=nx0:hst:nxf;% # del nodo de la coordenada 'x' de los receptores sin nodos absorbentes
%xst=round((nx0+nxf)*(1/4)); % Registro de una sólo receptor

zst=(nz0+1)*ones(1,length(xst));     % # del nodo de la coordenada 'z=0 constante' de los receptores 
%N_tr=length(xst);% <--- Numero de trazas a registrar
ng=length(xst); % número de estaciones o geófonos

%% Graficamos todas las coordenadas
%plot_coords2( x,z,xp,zp,xst,zst,nbc);
  
%% ********* Carga de la onduleta para la fuente o explosión **************

% Frecuencia máxima en Hz soportada por las Diferencias Finitas:
fk_max=(1/k)*15;% Marmousi

ts=1/fk_max;amp=1;
s=pulso_ricker_frec(nt,fk_max,dt,ts,amp);

% ******************* grafica de la fuente *********************
figure;plot(s);title([num2str(fk_max*k),'Hz Source wavelet']);%hold on; plot(s1,'o')

end