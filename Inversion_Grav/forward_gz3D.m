function gz = forward_gz3D( ...
                rho,...
                xmin,xmax,...
                ymin,ymax,...
                zmin,zmax,...
                x_obs,...
                y_obs,...
                z_obs)

%==========================================================
%
% FORWARD_GZ3D
%
% Forward gravimétrico 3D mediante prismas rectangulares
%
% INPUT
%
% rho      NZ x NX x NY
%
% dominio físico
%
% x_obs
% y_obs
% z_obs
%
% OUTPUT
%
% gz (mGal)
%
%==========================================================

[NZ,NX,NY] = size(rho);

%% Tamaños de celda

dx = (xmax-xmin)/NX;
dy = (ymax-ymin)/NY;
dz = (zmax-zmin)/NZ;

%% Centros

xc = xmin + (0.5:NX-0.5)*dx;
yc = ymin + (0.5:NY-0.5)*dy;
zc = zmin + (0.5:NZ-0.5)*dz;

%% Número de estaciones

Nobs = length(x_obs);

gz = zeros(Nobs,1);

%% ======================================================
% Suma de todos los prismas
%% ======================================================

for iobs = 1:Nobs

    g = 0;

    for iy = 1:NY

        y1 = yc(iy)-dy/2;
        y2 = yc(iy)+dy/2;

        for ix = 1:NX

            x1 = xc(ix)-dx/2;
            x2 = xc(ix)+dx/2;

            for iz = 1:NZ

                if rho(iz,ix,iy)==0
                    continue
                end

                z1 = zc(iz)-dz/2;
                z2 = zc(iz)+dz/2;

                g = g + prism_gz(...
                        x_obs(iobs),...
                        y_obs(iobs),...
                        z_obs(iobs),...
                        x1,x2,...
                        y1,y2,...
                        z1,z2,...
                        rho(iz,ix,iy));

            end

        end

    end

    gz(iobs)=g;

end

end