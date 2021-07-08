
% lilear model:

function [lin_vel] = linear_model(nz,nx,blur,v1,v2,v3,v4,v5)

    lin_vel=zeros(nz,nx);
    
    Int1=round(1/5 * nz); 
    Int2=round(2/5 * nz);
    Int3=round(3/5 * nz);
    Int4=round(4/5 * nz);
   
    
    
    lin_vel(1:Int1 , :)=v1;
    lin_vel(Int1:Int2 , :)=v2;
    lin_vel(Int2:Int3 , :)=v3;
    lin_vel(Int3:Int4 , :)=v4;
    lin_vel(Int4:end , :)=v5;
    
    lin_vel=filter_2Dfield(lin_vel,blur);

end