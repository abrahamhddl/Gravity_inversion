function result = solve_MAP_conjugate(...
                    G,...
                    L,...
                    W,...
                    alpha,...
                    d_obs,...
                    gamma,...
                    a0,...
                    b0)

%==============================================================
% MAP para modelo Normal - Inversa Gamma
%
% d | m,sigma2  ~ N(Gm, sigma2 I)
%
% p(m | sigma2) propto
% exp[-gamma/(2 sigma2) *
%      (Lm-alpha)' W (Lm-alpha)]
%
% sigma2 ~ InvGamma(a0,b0)
%==============================================================

nData = length(d_obs);

%% ============================================================
% 1. Media posterior / MAP de m
% =============================================================

b = G'*d_obs ...
  + gamma*(L'*(W*alpha));

Afun = @(x) ...
      G'*(G*x) ...
    + gamma*(L'*(W*(L*x)));

tol   = 1e-6;
maxit = 500;

[m_map,flag,relres,iter] = ...
    pcg(Afun,b,tol,maxit);

%% ============================================================
% 2. Prediccion y residual de datos
% =============================================================

g_map = G*m_map;

r_data = d_obs-g_map;

SSE_data = r_data'*r_data;

rmse = sqrt(mean(r_data.^2));

%% ============================================================
% 3. Residual del prior
% =============================================================

r_prior = L*m_map-alpha;

SSE_prior = r_prior'*(W*r_prior);

%% ============================================================
% 4. Posterior de sigma^2
%
% sigma2 | d ~ InvGamma(an,bn)
% =============================================================

an = a0 + nData/2;

bn = b0 ...
   + 0.5*(SSE_data + gamma*SSE_prior);

%% Media y moda de sigma^2

if an > 1
    sigma2_mean = bn/(an-1);
else
    sigma2_mean = NaN;
end

sigma2_mode = bn/(an+1);

sigma_mean = sqrt(sigma2_mean);
sigma_mode = sqrt(sigma2_mode);

%% ============================================================
% 5. Guardar resultados
% =============================================================

result.m_map = m_map;
result.g_map = g_map;

result.residual = r_data;

result.rmse = rmse;

result.flag   = flag;
result.relres = relres;
result.iter   = iter;

result.SSE_data  = SSE_data;
result.SSE_prior = SSE_prior;

result.a_post = an;
result.b_post = bn;

result.sigma2_mean = sigma2_mean;
result.sigma2_mode = sigma2_mode;

result.sigma_mean = sigma_mean;
result.sigma_mode = sigma_mode;

end