
function [bx_L,ax_L,bx_R,ax_R,bz_T,az_T,bz_B,az_B]=CPML_Acustic_params_ab(Rc,fd,dt,V_S,x,dx,x0,xf,nxCPML,nx0,nxf,z,dz,z0,zf,nzCPML,nz0,nzf)
%% PARÁMETROS PARA FRONTERAS CPML
N=2;
kx=1;
kz=1;
%Rc=1e-10; % coeficiente teórico de reflexión
alpha_max=pi*fd;
Vp=V_S; % si el medio no es homogéneo se toma la velocidad máxima

% Numero de nodos en total para el eje 'x' y eje 'z'
Nx=length(x);
Nz=length(z);

% ***************** Construcción de CPMLs en eje x ************************
Lx=dx*nxCPML; % longitud de cada CPML horizontal
dx0= -(N+1)*Vp*log(Rc)/(2*Lx);% ~173.4 con nbc=20

% función de amortiguamiento 'dx(x)' izquierda:
dxx_L(1:nxCPML)=dx0*((x(1:nx0-1)-x0)/Lx).^N; % restamos x0 cuando x0~=0

% función de amortiguamiento 'dx(x)' derecha:
dxx_R(1:nxCPML)=dx0*((x(nxf+1:Nx)-xf)/Lx).^N;%dxx_L(nxCPML:-1:1);

% recta 'alpha_x' con pendiente + (en frontera izquierda)
mx_positiva=(alpha_max/(x(nx0-1)-x(1))) ;
alpha_xL(1:nxCPML)= mx_positiva * (x(1:nx0-1)-x(1));

% recta 'alpha_x' con pendiente - (en frontera derecha)
mx_negativa=(alpha_max/(x(nxf+1)-x(Nx)));
alpha_xR(1:nxCPML)= mx_negativa * (x(nxf+1:Nx)-x(Nx));


% ***************** Construcción de CPMLs en eje z ************************
Lz=dz*nzCPML; % longitud de cada CPML vertical
dz0= -(N+1)*Vp*log(Rc)/(2*Lz);

% función de amortiguamiento 'dz(z)' superior (Top):
dzz_T(1:nzCPML)=dz0*((z(1:nz0-1)-z0)/Lz).^N; % restamos z0 cuando z0~=0

% función de amortiguamiento 'dz(z)' inferior (Bottom):
dzz_B(1:nzCPML)=dz0*((z(nzf+1:Nz)-zf)/Lz).^N;%dzz_T(nzCPML:-1:1);

% recta 'alpha_z' con pendiente + (en frontera superior)
mz_positiva=(alpha_max/(z(nz0-1)-z(1))) ;
alpha_zT(1:nzCPML)= mz_positiva * (z(1:nz0-1)-z(1));

% recta 'alpha_z' con pendiente - (en frontera inferior)
mz_negativa=(alpha_max/(z(nzf+1)-z(Nz)));
alpha_zB(1:nzCPML)= mz_negativa * (z(nzf+1:Nz)-z(Nz));


% *****************
% Las funciones 'psi'para CPML estan dadas por:
%   psi_V_x(i)=bx(i).*psi_V_x(i)+ax(i)*Dx_V;
%   psi_V_z(i)=bz(i).*psi_V_z(i)+az(i)*Dz_V;
%   psi_Syx(i)=bx(i).*psi_Syx(i)+ax(i)*Dx_Syx; 
%   psi_Syz(i)=bz(i).*psi_Syz(i)+az(i)*Dz_Syz;    
% donde los coeficientes 'ax', 'bx', 'az' y 'bz' se calculan como sigue:

% coeficientes 'ax' y 'bx' en frontera izquierda
bx_L(1:nxCPML)=exp(-( dxx_L/kx + alpha_xL )*dt);
ax_L(1:nxCPML)=( dxx_L./(kx*(dxx_L + alpha_xL*kx) )).*(bx_L-1);

% coeficientes 'ax' y 'bx' en frontera derecha
bx_R(1:nxCPML)=exp(-( dxx_R/kx + alpha_xR )*dt);%bx_L(nxCPML:-1:1);
ax_R(1:nxCPML)=( dxx_R./(kx*(dxx_R + alpha_xR*kx) )).*(bx_R-1);%ax_L(nxCPML:-1:1);

% coeficientes 'az' y 'bz' en frontera superor
bz_T(1:nzCPML)=exp(-( dzz_T/kz + alpha_zT )*dt);
az_T(1:nzCPML)=( dzz_T./(kz*(dzz_T + alpha_zT*kz) )).*(bz_T-1);

% coeficientes 'az' y 'bz' en frontera inferior
bz_B(1:nxCPML)= exp(-( dzz_B/kz + alpha_zB )*dt);%bz_T(nzCPML:-1:1);
az_B(1:nxCPML)= ( dzz_B./(kz*(dzz_B + alpha_zB*kz) )).*(bz_B-1);%az_T(nzCPML:-1:1);


end