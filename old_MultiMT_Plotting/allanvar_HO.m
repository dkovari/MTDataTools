function V = allanvar_HO(tau,alpha,kappa,kbT)

if nargin<4
    kbT = 4.11; %pN*nm
end

V = 2*kbT*alpha^2/kappa^3./tau.^2.*(kappa*tau/alpha + 2*exp(-kappa*tau/alpha) - 1/2*exp(-2*kappa*tau/alpha) - 3/2);