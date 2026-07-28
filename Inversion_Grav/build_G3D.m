function G = build_G3D(mesh,...
                       x_obs,...
                       y_obs,...
                       z_obs)
%==============================================================
% BUILD_G3D
%
% Construye la matriz de sensibilidad
%==============================================================

cacheFile = sprintf('G_%dx%dx%d_%d.mat',...
    mesh.NX,...
    mesh.NY,...
    mesh.NZ,...
    length(x_obs));

if exist(cacheFile,'file')

    fprintf('\n');
    fprintf('=====================================\n');
    fprintf('Cargando G desde disco...\n');

    S = load(cacheFile);

    G = S.G;
    
    disp(cacheFile)
    delete(cacheFile)
    return

end

%% Número de datos

nData = length(x_obs);

%% Número de parámetros

nModel = mesh.nModel;

%% Reservar memoria

G = zeros(nData,nModel,'single');

%% Geometría de prismas

x1 = mesh.xc - mesh.dx/2;
x2 = mesh.xc + mesh.dx/2;

y1 = mesh.yc - mesh.dy/2;
y2 = mesh.yc + mesh.dy/2;

z1 = mesh.zc - mesh.dz/2;
z2 = mesh.zc + mesh.dz/2;

%% Construcción

parfor c = 1:nModel

    [iz,ix,iy] = ind2sub([mesh.NZ mesh.NX mesh.NY],c);

    g = prism_gz_vector(...
        x_obs,...
        y_obs,...
        z_obs,...
        x1(ix),x2(ix),...
        y1(iy),y2(iy),...
        z1(iz),z2(iz),...
        1);

    G(:,c) = single(g);

    if mod(c,500)==0
        fprintf('Prisma %d de %d\n',c,nModel);
    end

end

fprintf('\n');
fprintf('Guardando G...\n');

save(cacheFile,'G','-v7.3');
end