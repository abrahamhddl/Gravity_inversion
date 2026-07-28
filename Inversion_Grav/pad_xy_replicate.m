function rho_pad = pad_xy_replicate(rho)
%==========================================================
%
% PAD_XY_REPLICATE
%
% Replica el modelo en X e Y utilizando el valor
% de las caras laterales.
%
% Entrada:
%   rho      -> NZ x NX x NY
%
% Salida:
%   rho_pad  -> NZ x (3*NX) x (3*NY)
%
%==========================================================

[NZ,NX,NY] = size(rho);

%% ========================================================
% Padding en X
%% ========================================================

left_pad  = repmat(rho(:,1,:),  [1 NX 1]);

right_pad = repmat(rho(:,end,:),[1 NX 1]);

rho_x = cat(2,left_pad,rho,right_pad);

%% ========================================================
% Padding en Y
%% ========================================================

front_pad = repmat(rho_x(:,:,1),  [1 1 NY]);

back_pad  = repmat(rho_x(:,:,end),[1 1 NY]);

rho_pad = cat(3,front_pad,rho_x,back_pad);

end