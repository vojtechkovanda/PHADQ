function [x] = CV(param, paramsolver, xq)

%% initialization

lam = 0.1; % setting threshold for clipping

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x) (sign(x).*min(abs(x), lam));

x = paramsolver.x0;

u3 = x;
u1 = paramsolver.u1;
u2 = paramsolver.u2;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;


%% iteration

for i = 1:paramsolver.I
    
    waitbar(i/paramsolver.I);

    U1 = param.L1_adj(u1);
    U2 = param.L_adj(u2);
    U3 = u3;
    
    x_tild = (x - tau * (U1+U2+U3));
    x = alpha * x_tild + (1 - alpha) * x;

    bL = 2*x_tild-x;

    p1 = u1 + sigma * param.L1(bL);
    u1_tild = clip(p1);
    u1 = alpha * u1_tild + (1 - alpha) * u1;

    p2 = u2 + sigma * param.L(bL);
    u2_tild = param.prox(p2);
    u2 = alpha * u2_tild + (1 - alpha) * u2;

    p3 = u3 + sigma * bL;
    u3_tild = p3 - sigma * projection(p3/sigma, xq, param.delta);
    u3 = alpha * u3_tild + (1 - alpha) * u3;
        
end