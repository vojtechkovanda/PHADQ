function [p, SDR] = cp_alg(y, param, paramsolver, in)
% Chambolle-Pock algorithm
%
% Vojtěch Kovanda
% Brno University of Technology, 2024

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x) (sign(x).*min(abs(x), 1));

zeta = param.lam(param.delta); % setting threshold for clipping
sig = 1/zeta;

%% initial values
i = 0;

x = paramsolver.x0;
p = x;
q = frana(param.F, x);
SDR = zeros(paramsolver.K, 1);

%% algorithm
while i < paramsolver.K

    i = i + 1;
     
waitbar(i/paramsolver.K);
    
     q = clip(q + sig.*frana(param.F, x));
     p_old = p;
     p1 = frsyn(param.F, q);
     p1 = p1(1:param.Ls);
     p = projection(p-zeta*p1, y, param.delta);
     x = p + param.rho*(p - p_old);

     SDR(i) = 20*log10(norm(in,2)./norm(in-p, 2));

end

