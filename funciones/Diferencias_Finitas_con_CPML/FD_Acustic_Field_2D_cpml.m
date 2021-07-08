


% Aproximación espacial 2D de los campos: Vy,Syz,Syx  
% y sus repectivas derivadas: Dx_Vy,Dz_Vy,Dx_Syx,Dz_Syz
% usando diferencias finitas de 4o orden con el esquema Stagered grid
% Velocidad-Esfuerzo, usando fronteras absorbentes tipo C-PML, cuyos
% parámetros son:
% ax_L,bx_L,ax_R,bx_R,az_T,bz_T,az_B,bz_B,
% psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB


function ...
    [Vy,Syz,Syx,...
    psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB]...
    ...
    =FD_Acustic_Field_2D_cpml...
    ...
    (dt,dx,dz,vel,nbc,Vy,Syz,Syx,...
    ax_L,bx_L,ax_R,bx_R,az_T,bz_T,az_B,bz_B,...
    psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB,isFS)
    
    

    %**********************************************************************
    [Nz,Nx]=size(vel);
    
    nzCPML=nbc;
    nxCPML=nbc;
    
    nx0=nxCPML+1; % posicion real de x0
    nxf=Nx-nxCPML; % posicion real de xf
    
    nz0=nzCPML+1; % posicion real de z0
    nzf=Nz-nzCPML; % posicion real de zf
    
    % Derivadas a usar:
    Dx_Vy=zeros(size(vel));
    Dz_Vy=zeros(size(vel));
    Dx_Syx=zeros(size(vel));
    Dz_Syz=zeros(size(vel));
    
    
    % Coeficientes para Diferencias finitas de 4o Orden:
    A=-1/24;
    B=9/8;
   
    %********************** Aproximación de esfuerzos Syx *****************
    
    % Aproximación dela derivada 'Dx_V' en toda la malla (en dirección x):
    i=2:Nx-2;   
    Dx_Vy(:,i)=( A*( Vy(:,i+2)-Vy(:,i-1) ) + B*( Vy(:,i+1)-Vy(:,i) ) )/dx;
    
    % Aplicamos CPML en la frontera izquierda (L):
    i=1:nx0-1;  
    psi_V_xL=bx_L.*psi_V_xL + ax_L.*Dx_Vy(:,i);
    Dx_Vy(:,i)=Dx_Vy(:,i) + psi_V_xL;
    
    % Aplicamos CPML en la frontera derecha (R):
    i=nxf+1:Nx;
    psi_V_xR=bx_R.*psi_V_xR + ax_R.*Dx_Vy(:,i);
    Dx_Vy(:,i)=Dx_Vy(:,i) + psi_V_xR;
    
    % Aproximamos esfuerzos Syx:
    Syx = Syx + dt*Dx_Vy; 


    %********************** Aproximación de esfuerzos Syz *****************  
    
    % Aproximación dela derivada 'Dz_V' en toda la malla (en dirección z):
    l=2:Nz-2;
    Dz_Vy(l,:)=( A*( Vy(l+2,:)-Vy(l-1,:) ) + B*( Vy(l+1,:)-Vy(l,:) ) )/dz;
    
    % Aplicamos CPML en la frontera superior (T):
    l=1:nz0-1;
    psi_V_zT=(bz_T').*psi_V_zT + (az_T').*Dz_Vy(l,:);
    Dz_Vy(l,:)=Dz_Vy(l,:) + psi_V_zT;
    
    % Aplicamos CPML en la frontera inferior (B):
    l=nzf+1:Nz;
    psi_V_zB=(bz_B').*psi_V_zB + (az_B').*Dz_Vy(l,:);
    Dz_Vy(l,:)=Dz_Vy(l,:) + psi_V_zB;
    
    % Aproximamos esfuerzos Syz:
    Syz = Syz + dt*Dz_Vy;
  
    
    % ya que aproximamos Syx y Syz, podemos aproximar las velocidades 'V'
    %************** Aproximación de velocidades ***************************
    
    % Aproximación dela derivada 'Dx_Syx' en toda la malla (en dirección x):
    i=3:Nx-1;   
    Dx_Syx(:,i) = ( A*(Syx(:,i+1)-Syx(:,i-2)) + B*( Syx(:,i)- Syx(:,i-1) ) )/dx;
    
    % Aplicamos CPML en la frontera izquierda (L):
    i=1:nx0-1;
    psi_Syx_xL=bx_L.*psi_Syx_xL + ax_L.*Dx_Syx(:,i);   
    Dx_Syx(:,i) = Dx_Syx(:,i) + psi_Syx_xL;
    
    % Aplicamos CPML en la frontera derecha (R):
    i=nxf+1:Nx;   
    psi_Syx_xR=bx_R.*psi_Syx_xR + ax_R.*Dx_Syx(:,i);
    Dx_Syx(:,i) = Dx_Syx(:,i) + psi_Syx_xR;
    
    
    % Aproximación dela derivada 'Dz_Syz' en toda la malla (en dirección z):
    l=3:Nz-1;
    Dz_Syz(l,:) = ( A*(Syz(l+1,:)-Syz(l-2,:)) + B*( Syz(l,:)- Syz(l-1,:) ) )/dz;
    
    % Aplicamos CPML en la frontera superior (T):
    l=1:nz0-1;
    psi_Syz_zT=(bz_T').*psi_Syz_zT + (az_T').*Dz_Syz(l,:);    
    Dz_Syz(l,:) = Dz_Syz(l,:) +  psi_Syz_zT;
    
    % Aplicamos CPML en la frontera inferior (B):
    l=nzf+1:Nz;
    psi_Syz_zB=(bz_B').*psi_Syz_zB + (az_B').*Dz_Syz(l,:);
    Dz_Syz(l,:) = Dz_Syz(l,:) +  psi_Syz_zB;

    
    % Aproximamos Vy:
    Vy = Vy + dt*(vel.^2).*( Dx_Syx + Dz_Syz ); 
    
    
    %************** Condición de superficie libre para Vy ****************
    if isFS        
        % superficie libre en: 
        zFS=nzCPML+1;
        
        %ponemos ceros en los nodos vacios (arriba de la superficie libre)
        Vy(1:zFS-1,:) = 0.0;
        
        % ponemos ceros sobre la superficiel libre y hacemos una simetria:
        Vy(zFS , :)=0.0;        
        Vy(zFS-1:-1:zFS-2 , :) = Vy(zFS+1:zFS+2 , :);% simetria par             
    end
    
    
    %Uy = Uy + dt*Vy; % desplazamiento 'Uy'
    
end