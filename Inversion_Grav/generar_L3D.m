function L = generar_L3D(nx,ny,nz)
%==============================================================
% GENERAR_L3D
%
% Operador GMRF 3D de 6 vecinos.
%
% Orden compatible con:
%
%   rho(NZ,NX,NY)
%   rho(:)
%   reshape(m,NZ,NX,NY)
%
% Para cada nodo:
%
%   (Nvecinos + 1)*m_i - sum(m_j)
%
% de modo que el prior sea:
%
%   ||L*m - alpha||^2
%
%==============================================================

N = nx*ny*nz;

% Máximo 7 vecinos no nulos + diagonal
L = spalloc(N,N,7*N);

for iy = 1:ny

    for ix = 1:nx

        for iz = 1:nz

            %--------------------------------------
            % Índice MATLAB correcto: Z-X-Y
            %--------------------------------------

            idx = sub2ind([nz nx ny],...
                          iz,ix,iy);

            n = 0;

            %======================================
            % Z
            %======================================

            if iz > 1

                id = sub2ind([nz nx ny],...
                             iz-1,ix,iy);

                L(idx,id) = -1;

                n = n+1;

            end

            if iz < nz

                id = sub2ind([nz nx ny],...
                             iz+1,ix,iy);

                L(idx,id) = -1;

                n = n+1;

            end

            %======================================
            % X
            %======================================

            if ix > 1

                id = sub2ind([nz nx ny],...
                             iz,ix-1,iy);

                L(idx,id) = -1;

                n = n+1;

            end

            if ix < nx

                id = sub2ind([nz nx ny],...
                             iz,ix+1,iy);

                L(idx,id) = -1;

                n = n+1;

            end

            %======================================
            % Y
            %======================================

            if iy > 1

                id = sub2ind([nz nx ny],...
                             iz,ix,iy-1);

                L(idx,id) = -1;

                n = n+1;

            end

            if iy < ny

                id = sub2ind([nz nx ny],...
                             iz,ix,iy+1);

                L(idx,id) = -1;

                n = n+1;

            end

            %======================================
            % Diagonal
            %
            % +1 por condición informativa alpha
            %======================================

            L(idx,idx) = n + 1;

        end

    end

end

end