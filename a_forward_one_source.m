
%%
% Autor: Marcos Bernal Romero (marck_rt@hotmail.com)
% Departamento de Matemáticas, Facultad de Ciencias, UNAM.

%% Modelado directo de la ecuación de onda para medios acusticos:
% Dtt_P - v²(Dxx_P + Dzz_P) = f 

% v(:,:) <-- modelo de velocidades

% usando la formulación velocidad-esfuerzos:
% Dt_V   = v²(Dx_Sxx + Dz_Szz) + f
% Dt_Sxx = Dx_V
% Dt_Szz = Dz_V

% Se implementan fronteras absorbentes tipo C-PML 

%%
clear; close all; clc;

addpath((genpath(pwd)));

%% PARAMETROS MALLA, FUENTES, ESTACIONES:

% Parámetros del modelo canadiense
dir_model='vel_marmousi2_model';
[s,nt,dt,nbc,x0,dx,nx,z0,nz,sx,sz,ns,gx,gz,ng,isFS,animar,sismograma,Rc]=indata_Acustic_Marmousi;
fprintf('\n====== Trabajando con el modelo canadiense ======\n')

is=round((10/20)*ns); % <-- número de fuente

%% Cargamos modelo de velocidad (SIN NODOS ABSORVENTES)

% modelo de marmousi
load(['./modelos_vel/',num2str(dir_model),'/vel_model.mat']);% vel(:,:)
%vel=vel(1:nz,1:nx);nt=4000;

% modelo homogéneo
%vel=3000*ones(nz,nx);nt=1500;

% modelo de capas
vel = linear_model(nz,nx,1,3000,3200,3600,3800,4200);nt=2000;

%% CONSTRUIMOS TRAZAS SINTETICAS
disp('***** Construyendo trazas observadas: *****');

% CONTROL DE ERRORES NUMÉRICOS:
Vmin=min(min(vel));
Vmax=max(max(vel));
fmax=frecDominant(dt,s);
control_error(Vmin,Vmax,fmax,dt,dx);       

disp(['Modelando, fuente ',num2str(is),' de ',num2str(ns)]);
tic;
[seis_real]=FD_Acustic_trazas(vel,nbc,nt,dt,s,x0,dx,nx,z0,nz,sx(is),sz(is),gx,gz,isFS,animar,Rc);  
%[seis_real]=PS_Acustic_trazas(0.7*vel,nbc,nt,dt,s,x0,dx,nx,z0,nz,sx(is),sz(is),gx,gz,isFS,animar,Rc);
toc;


%% Gráficas:

figure;
imagesc(seis_real);
colormap(gray);caxis([-0.001 0.001]);
title('sismograma')
ylabel('tiempo')
xlabel('Trazas')


figure;
subplot(1,4,1)
tr1=round( (2/10)*size(seis_real,2) );% <--- traza a graficar
hold on
plot(seis_real(1:nt,tr1),'LineWidth',2)
title(['traza No.',num2str(tr1)])
view(90,90)

subplot(1,4,2)
tr2=round( (4/10)*size(seis_real,2) );% <--- traza a graficar
hold on
plot(seis_real(1:nt,tr2),'LineWidth',2)
title(['traza No.',num2str(tr2)])
view(90,90)

subplot(1,4,3)
tr3=round( (7/10)*size(seis_real,2) );% <--- traza a graficar
hold on
plot(seis_real(1:nt,tr3),'LineWidth',2)
title(['traza No.',num2str(tr3)])
view(90,90)

subplot(1,4,4)
tr4=round( (8/10)*size(seis_real,2) );% <--- traza a graficar
hold on
plot(seis_real(1:nt,tr4),'LineWidth',2)
title(['traza No.',num2str(tr4)])
view(90,90)