
% función que ajusta el modelo de velocidades 'vel' 
% para agregar nodos con condiciones de frontera absorbentes

function v=redim_bound(vel,nbc)
v=[repmat(vel(:,1),1,nbc), vel, repmat(vel(:,end),1,nbc)]; % copia columnas exteriores
v=[repmat(v(1,:),nbc,1); v; repmat(v(end,:),nbc,1)]; % copia renglones exteriores
end