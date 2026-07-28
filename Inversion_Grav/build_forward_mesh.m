function forward = build_forward_mesh(mesh,rho_prior,NXf,NYf,NZf)
%==============================================================
%
% BUILD_FORWARD_MESH
%
% Construye la malla fina e interpola el modelo previo.
% Posteriormente aplica el padding utilizado por forward_gz3D.
%
%==============================================================

%%------------------------------------------------------------
% Malla fina
%%------------------------------------------------------------

forward.NX = NXf;
forward.NY = NYf;
forward.NZ = NZf;

forward.dx = (mesh.xmax-mesh.xmin)/NXf;
forward.dy = (mesh.ymax-mesh.ymin)/NYf;
forward.dz = (mesh.zmax-mesh.zmin)/NZf;

forward.xf = mesh.xmin + (0.5:NXf-0.5)*forward.dx;
forward.yf = mesh.ymin + (0.5:NYf-0.5)*forward.dy;
forward.zf = mesh.zmin + (0.5:NZf-0.5)*forward.dz;

[forward.Zf,...
 forward.Xf,...
 forward.Yf] = ndgrid(forward.zf,...
                      forward.xf,...
                      forward.yf);

%%------------------------------------------------------------
% Interpolar modelo a la malla fina
%%------------------------------------------------------------

Fprior = griddedInterpolant({mesh.zc,...
                             mesh.xc,...
                             mesh.yc},...
                             rho_prior,...
                             'nearest',...
                             'nearest');

forward.rho_fine = Fprior(forward.Zf,...
                          forward.Xf,...
                          forward.Yf);

%%------------------------------------------------------------
% Padding
%%------------------------------------------------------------

forward.rho_pad = pad_xy_replicate(forward.rho_fine);

%%------------------------------------------------------------
% Dominio del padding
%%------------------------------------------------------------

Lx = mesh.xmax-mesh.xmin;
Ly = mesh.ymax-mesh.ymin;

forward.xmin_pad = mesh.xmin-Lx;
forward.xmax_pad = mesh.xmax+Lx;

forward.ymin_pad = mesh.ymin-Ly;
forward.ymax_pad = mesh.ymax+Ly;

end