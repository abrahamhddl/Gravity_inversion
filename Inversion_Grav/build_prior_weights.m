function W = build_prior_weights(mesh)

%% Perfil central

y0 = mean(mesh.yc);

%% Distancia normalizada

r = abs(mesh.yc - y0);

r = r/max(r);

%% Peso coseno

wY = 0.2 + 0.8*cos(pi*r/2).^2;

%% Expandir a todas las celdas

weights = zeros(mesh.nModel,1);

nxy = mesh.NX*mesh.NZ;

for iy = 1:mesh.NY

    ind = (iy-1)*nxy + (1:nxy);

    weights(ind) = wY(iy);

end

%% Matriz diagonal

W = spdiags(weights,0,mesh.nModel,mesh.nModel);

end