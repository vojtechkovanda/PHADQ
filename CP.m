function [x, SDR] = CP(param, paramsolver, xq, in)

%% initialization

x = paramsolver.x0;
u = paramsolver.u0;
SDR = zeros(1, paramsolver.I);

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;


%% iteration

for i = 1:paramsolver.I
    
    waitbar(i/paramsolver.I);

    
    p = projection(x - tau*sigma*param.L_adj(u), xq, param.delta);
    v = u + param.L(2*p - x);
    q = v - param.prox(v);
    
    x = x + alpha*(p - x);
    u = u + alpha*(q - u);
    
    SDR(i) = 20*log10(norm(in,2)./norm(in-x, 2));

end