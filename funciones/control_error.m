function control_error(Vmin,Vmax,fmax,dt,dx,m)   
    
    if nargin > 5 % Estabilidad para filtro de onda S
        
        n=m;
        h=1;
        
    else % Estabilidad para onda acústica
                
        %n=12;% numero de muestras por longitud de onda para FD de 2o orden    
        %h=1; % para FD de 2o orden
        
        n=8; % numero de muestras por longitud de onda para FD de 4o orden        
        h=7/6; % para FD de 4o orden
        
    end
    

    c_stab=dx/(h*sqrt(2)*Vmax); % Cota para error de Estabilidad
    
    c_disp=Vmin/(n*fmax); % Cota para error de Dispersión
    
      
    if dt > c_stab
        
        disp('ERROR DE ESTABILIDAD')
        fprintf('dt= %f y debe ser menor o igual que %f \n',dt,c_stab)
        fprintf('Reducir Vmax ó incremente dx \n\n')
        stop;
    
    else % Hay Estabilidad
        
        if dx > c_disp % Hay error de dispersion
            disp('Hay estabilidad con error de dispersión')
            fprintf('dx= %f y debe ser menor o igual que %f \n',dx,c_disp)
            fprintf('Incrementar dt ó incremente periodo de la fuente \n\n')
            
        else % No hay error de dispersion
            disp('Hay estabilidad sin dispersión!')
        end
        
    end

end
