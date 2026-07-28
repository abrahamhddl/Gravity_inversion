clc
clear

addpath('../functions')

%%========================================================
% Datos
%%========================================================
data = readmatrix('C:\Users\titis\Documents\Kyoto_University\Data_gsg\prof_Q\perfiles_paralelos.csv');

perfil = data(:,1);

x_obs = data(:,2)*1000;

dy_profile = 5000;

y_obs = perfil*dy_profile;

z_obs = zeros(size(x_obs));

d_obs = data(:,5);

d_obs = d_obs-mean(d_obs);

sigma_g = 1;
%%========================================================
% Dominio
%%========================================================
perfiles = unique(perfil);

NX = 50;
NZ = 30;
NY = numel(perfiles);

xmin = -200e3;
xmax = 200e3;

zmin = 0;
zmax = 25e3;

dy_profile = 5000;

mesh = build_mesh3D(...
            NX,NY,NZ,...
            xmin,xmax,...
            zmin,zmax,...
            dy_profile);


W = build_prior_weights(mesh);
%%========================================================
% Prior
%%========================================================

[rho_prior,alpha,priorInfo] = ...
    build_prior_campos(...
        mesh,...
        'C:\Users\titis\Documents\Kyoto_University\Data_gsg\array\Campos_array20x20.csv');


%%========================================================
% Forward mesh
%%========================================================

forward = build_forward_mesh(...
            mesh,...
            rho_prior,...
            68,...
            30,...
            71);

%%========================================================
% Sensibilidad
%%========================================================

disp('Construyendo G ...')

tic

G = build_G3D(...
        mesh,...
        x_obs,...
        y_obs,...
        z_obs);

toc

G = double(G);

%%========================================================
% Regularización
%%========================================================

clear L

L = generar_L3D(mesh.NX,...
                mesh.NY,...
                mesh.NZ);

L = sparse(L);

%% Crear estado nuevo

new.mesh = mesh;

new.rho_prior = rho_prior;

new.alpha = alpha;

new.G = G;

new.L = L;

compare_states(old,new)

%%========================================================
% MAP
%%========================================================

gamma = 1e-4;

tic
result = solve_MAP(...
            G,...
            L,...
            W,...
            alpha,...
            d_obs,...
            sigma_g,...
            gamma);
toc

rho_map = reshape(result.m_map,...
                  mesh.NZ,...
                  mesh.NX,...
                  mesh.NY);
%% ==========================================
% Diagnostico por perfil Y
%% ==========================================

perfiles = unique(perfil);

nXZ = mesh.NZ*mesh.NX;

r_prior = L*result.m_map - alpha;
r_data  = d_obs - result.g_map;

w = full(diag(W));

RMSprior = zeros(mesh.NY,1);
RMSpriorW = zeros(mesh.NY,1);
RMSEdata = zeros(mesh.NY,1);
wY_plot = zeros(mesh.NY,1);

for iy = 1:mesh.NY

    %% Modelo
    indM = (iy-1)*nXZ + (1:nXZ);

    RMSprior(iy) = ...
        sqrt(mean(r_prior(indM).^2));

    RMSpriorW(iy) = ...
        sqrt(mean(w(indM).*r_prior(indM).^2));

    wY_plot(iy) = mean(w(indM));

    %% Datos
    idxD = perfil == perfiles(iy);

    RMSEdata(iy) = ...
        sqrt(mean(r_data(idxD).^2));

end
%% ==========================================
% problema condicionado a sigma cuadrada
%% ==========================================

a0 = 2;
b0 = 1;

result_conj = solve_MAP_conjugate(...
                    G,...
                    L,...
                    W,...
                    alpha,...
                    d_obs,...
                    gamma,...
                    a0,...
                    b0);
%% ============================================================
% VISUALIZACION 3D DEL MODELO MAP
% Parte directamente de result_conj.m_map
%% ============================================================

%% 1. Reconstruir modelo 3D original
rho_map_conj = reshape(result_conj.m_map,...
                  mesh.NZ,...
                  mesh.NX,...
                  mesh.NY);
%% ============================================================
% 2. Crear interpolador 3D
%% ============================================================
Fmap_conj = griddedInterpolant(...
            {mesh.zc,...
             mesh.xc,...
             mesh.yc},...
            rho_map_conj,...
            'linear',...
            'nearest');
%% ============================================================
% 3. Malla fina SOLO PARA VISUALIZACION
%% ============================================================
NXplot = 120;
NYplot = 120;
NZplot = 60;

xq = linspace(mesh.xc(1),...
              mesh.xc(end),...
              NXplot);

yq = linspace(mesh.yc(1),...
              mesh.yc(end),...
              NYplot);

zq = linspace(mesh.zc(1),...
              mesh.zc(end),...
              NZplot);

[Zq,Xq,Yq] = ndgrid(zq,xq,yq);
%% Interpolar
rho_interp_conj = Fmap_conj(Zq,Xq,Yq);
%% ============================================================
% 4. Convertir a orden X-Y-Z para SLICE
%% ============================================================
V = permute(rho_interp_conj,[2 3 1]);

[Xs,Ys,Zs] = meshgrid(...
                    xq/1000,...
                    yq/1000,...
                    zq/1000);
%% ============================================================
% 5. POSICIONES DE LOS CORTES
%% ============================================================

% Cortes verticales X-Z
yslice = [50 200 350];

% Cortes verticales Y-Z
xslice = [];

% Cortes horizontales X-Y
zslice = [];
