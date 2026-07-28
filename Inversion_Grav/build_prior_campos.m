function [rho_prior,alpha,priorInfo] = build_prior_campos(mesh,fileCampos)

Campos = readmatrix(fileCampos);

rho_ref = 2670;

Campos = Campos - rho_ref;

[NZc,NXc] = size(Campos);

%%---------------------------------------------
% Geometría física del modelo de Campos
%%---------------------------------------------

Lx_campos = 220e3;

dxC = Lx_campos/NXc;

itrench = 15;

x_trench = (itrench-0.5)*dxC;

xC = (0.5:NXc-0.5)*dxC - x_trench;

zC = mesh.zmin + ...
     (0.5:NZc-0.5)*(mesh.zmax-mesh.zmin)/NZc;

%%---------------------------------------------
% Interpolador
%%---------------------------------------------

Fc = griddedInterpolant({zC,xC},...
                        Campos,...
                        'nearest',...
                        'none');

perfil_interp = Fc(mesh.Z(:,:,1),...
                   mesh.X(:,:,1));

perfil_interp(isnan(perfil_interp)) = 0;

figure

imagesc(mesh.xc/1000,...
         mesh.zc/1000,...
         perfil_interp)

set(gca,'YDir','reverse')
axis image
colorbar

title('perfil_interp')


rho_prior = zeros(mesh.NZ,...
                  mesh.NX,...
                  mesh.NY);

for iy=1:mesh.NY

    rho_prior(:,:,iy)=perfil_interp;

end

%%---------------------------------------------
% Diagnóstico
%%---------------------------------------------

priorInfo.xmin = min(xC);
priorInfo.xmax = max(xC);

priorInfo.zmin = min(zC);
priorInfo.zmax = max(zC);

alpha = rho_prior(:);

end