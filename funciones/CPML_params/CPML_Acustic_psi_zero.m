function [psi_V_xL,psi_V_xR,psi_Syx_xL,psi_Syx_xR,psi_V_zT,psi_V_zB,psi_Syz_zT,psi_Syz_zB]=...
    CPML_Acustic_psi_zero(Nz,Nx,nzCPML,nxCPML)

    % valores iniciales par las 'psi' (que dependen de x):
    psi_V_xL=zeros(Nz,nxCPML);% psi(v(x)) izquierda 
    psi_V_xR=zeros(Nz,nxCPML);% psi(v(x)) derecha

    psi_Syx_xL=zeros(Nz,nxCPML); % psi(Syx(x)) izquierda
    psi_Syx_xR=zeros(Nz,nxCPML); % psi(Syx(x)) derecha

    % valores iniciales par las 'psi' (que dependen de z):
    psi_V_zT=zeros(nzCPML,Nx);% psi(v(z)) superior
    psi_V_zB=zeros(nzCPML,Nx);% psi(v(z)) inferior

    psi_Syz_zT=zeros(nzCPML,Nx); % psi(Syz(z)) superior
    psi_Syz_zB=zeros(nzCPML,Nx); % psi(Syz(z)) inferior

end