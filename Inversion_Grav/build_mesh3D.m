function mesh = build_mesh3D(NX,NY,NZ,...
                             xmin,xmax,...
                             zmin,zmax,...
                             dy_profile)
%==============================================================
%
% BUILD_MESH3D
%
%==============================================================

%% Número de celdas

mesh.NX = NX;
mesh.NY = NY;
mesh.NZ = NZ;

%% Dominio

mesh.xmin = xmin;
mesh.xmax = xmax;

mesh.ymin = 0;
mesh.ymax = (NY-1)*dy_profile;

mesh.zmin = zmin;
mesh.zmax = zmax;

%% Tamaño de celda

mesh.dx = (xmax-xmin)/NX;
mesh.dy = dy_profile;
mesh.dz = (zmax-zmin)/NZ;

%% Centros

mesh.xc = xmin + (0.5:NX-0.5)*mesh.dx;

mesh.yc = mesh.ymin + (0.5:NY-0.5)*mesh.dy;

mesh.zc = zmin + (0.5:NZ-0.5)*mesh.dz;

%% Malla

[mesh.Z,...
 mesh.X,...
 mesh.Y] = ndgrid(mesh.zc,...
                  mesh.xc,...
                  mesh.yc);

%% Número de parámetros

mesh.nModel = NX*NY*NZ;

end