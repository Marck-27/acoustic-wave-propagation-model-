


% Aproximación espacial 2D de los campos: Vy,Syz,Syx  
% usando método Pseudo-Spectral (PS) con el esquema Stagered grid
% Velocidad-Esfuerzo, usando fronteras absorbentes tipo C-PML, cuyos
% parámetros son:
% ax_L,bx_L,ax_R,bx_R,az_T,bz_T,az_B,bz_B,
% psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB
% 'argX' y 'argZ' son argumentos de la ifft en el método Pseudo-Spectral
% 'Ex' y 'Ez' son parámetros de una exponencial que sirven para reducir el ruido de Nyquist 


function ...
    [Vy,Syz,Syx,...
    psi_V_xL,psi_V_xR,psi_V_zT,psi_V_zB,psi_Syx_xL,psi_Syx_xR,psi_Syz_zT,psi_Syz_zB]...
    ...
    =PS_Acustic_Field_2D_cpml...
    ...
    (dt,vel,argX,argZ,Ex,Ez,nbc,Vy,Syz,Syx,...
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
    

   
    %********************** Aproximación de esfuerzos Syx *****************  
    
    
    % Aproximacion de dreivada Dx_Vy en toda la malla:
    KXX=exp(-Ex);
    Dx_Vy=real( ifft2( KXX .* argX .*fft2(Vy) ) );   
    
    % Actualizamos las derivadas con CPML en la frontera izquierda (L) 
     i=1:nx0-1; 
        psi_V_xL=bx_L.*psi_V_xL+ax_L.*Dx_Vy(:,i);
        Dx_Vy(:,i)=Dx_Vy(:,i) + psi_V_xL;
        
     % Actualizamos las derivadas con CPML en la frontera derecha (R)   
     i=nxf+1:Nx;
        psi_V_xR=bx_R.*psi_V_xR + ax_R.*Dx_Vy(:,i);
        Dx_Vy(:,i)=Dx_Vy(:,i) + psi_V_xR;  
        
      % Aproximamos Syx:
      Syx = Syx + dt*Dx_Vy;
      
    %********************** Aproximación de esfuerzos Syz *****************  
    
    % Aproximacion de dreivada Dz_Vy en toda la malla:
    KZZ=exp(-Ez);
    Dz_Vy=real( ifft2( KZZ' .* argZ' .* fft2(Vy) ) );
    
    % Actualizamos las derivadas con CPML en la frontera superior (T)
    l=1:nz0-1;   
        psi_V_zT = bz_T'.*psi_V_zT + az_T'.*Dz_Vy(l,:);
        Dz_Vy(l,:)=Dz_Vy(l,:) + psi_V_zT;
        
    % Actualizamos las derivadas con CPML en la frontera inferior (B)
    l=nzf+1:Nz;       
        psi_V_zB=bz_B'.*psi_V_zB + az_B'.*Dz_Vy(l,:);
        Dz_Vy(l,:)=Dz_Vy(l,:) + psi_V_zB;   
        
    % Aproximamos Syz:
    Syz = Syz + dt*Dz_Vy;
    
    % ya que aproximamos Syx y Syz, podemos aproximar las velocidades 'V'
    %************** Aproximación de velocidades ***************************
    
    % Aproximacion de dreivada Dx_Syx en toda la malla:
    KXX=exp(Ex);
    Dx_Syx=real( ifft2( KXX .* argX .*fft2(Syx) ) );
    
        % Actualizamos las derivadas con CPML en la frontera izquierda (L) 
     i=1:nx0-1; 
        psi_Syx_xL=bx_L.*psi_Syx_xL+ax_L.*Dx_Syx(:,i);
        Dx_Syx(:,i)=Dx_Syx(:,i) + psi_Syx_xL;
        
     % Actualizamos las derivadas con CPML en la frontera derecha (R)   
     i=nxf+1:Nx;
        psi_Syx_xR=bx_R.*psi_Syx_xR + ax_R.*Dx_Syx(:,i);
        Dx_Syx(:,i)=Dx_Syx(:,i) + psi_Syx_xR;
       
    
    
    % Aproximacion de dreivada Dz_Syz en toda la malla:
    KZZ=exp(Ez);
    Dz_Syz=real( ifft2( KZZ' .* argZ' .*fft2(Syz) ) );
    
    % Actualizamos las derivadas con CPML en la frontera superior (T)
    l=1:nz0-1;   
        psi_Syz_zT = bz_T'.*psi_Syz_zT + az_T'.*Dz_Syz(l,:);
        Dz_Syz(l,:)=Dz_Syz(l,:) + psi_Syz_zT;
        
    % Actualizamos las derivadas con CPML en la frontera inferior (B)
    l=nzf+1:Nz;       
        psi_Syz_zB=bz_B'.*psi_Syz_zB + az_B'.*Dz_Syz(l,:);
        Dz_Syz(l,:)=Dz_Syz(l,:) + psi_Syz_zB;
    
    
    % Aproximamos Syz:
    Vy = Vy+ dt*(vel.^2).*(Dx_Syx + Dz_Syz);
    
    
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
    
   % Uy = Uy + dt*Vy; % desplazamiento 'U' 
    
end