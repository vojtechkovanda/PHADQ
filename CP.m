function x = CP(param, paramsolver, xq)

%% initialization

x = paramsolver.x0;
q = paramsolver.u0;
p_old = paramsolver.x0;


tau = paramsolver.tau;
sigma = paramsolver.sigma;
rho = paramsolver.rho;

clip = @(x) (sign(x).*min(abs(x), paramsolver.lambda(param.delta)));

%% iteration

for i = 1:paramsolver.I
    
    waitbar(i/paramsolver.I);

     q = clip(q + sigma.*param.L(x));
     u = p_old - tau.*param.L_adj(q);
     p = projection(u, xq, param.delta);
     x = p + rho.*(p-p_old);
     p_old = p;
    
end