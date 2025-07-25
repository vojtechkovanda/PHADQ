function [x, SDR, ODG] = CP_incons(param, paramsolver, xq, in)

%% initialization

x = paramsolver.x0;
q = paramsolver.u0;
p_old = paramsolver.x0;

SDR = zeros(1, paramsolver.I);

ODG = zeros(1, floor(paramsolver.I/5));
ODG = [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, ODG];

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;

clip = @(x) (sign(x).*min(abs(x), paramsolver.lambda(param.delta)));


%% iteration

for i = 1:paramsolver.I
    
    waitbar(i/paramsolver.I);

     q = clip(q + sigma.*param.L(x));
     u = p_old - tau.*param.L_adj(q);
     p = 1/(1+tau)*(tau*projection(u, xq, param.delta) + u);
     x = p + alpha.*(p-p_old);
     p_old = p;
    
    SDR(i) = 20*log10(norm(in,2)./norm(in-x, 2));

        if mod(i, 5) == 0
        
       [~, ~, ODG(i/5+10)] = audioqual(in, x, param.fs);

         end

end