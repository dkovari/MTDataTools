function [kappa,alpha,k_ci,a_ci] = FitAllanVariance_HO(av,t,Navg,varargin)
% Fit the Allan variance of a signal to a stochastic harmonic oscillator
% 
% Input:
%   av: Vector contraining Allan variance at varous durations
%   t: duration (tau) for each calculated allan variance
%   Navg: the number of separate pairs used to calculate each element in av
%   kappa_0: (optional) initial starting point for kappa search
%   alpha_0: (optional) initial starting point for alhpa search
%   kbT: (default=4.11 pN*nm) thermal energy
%        if kbT is not specified, the AV is assumed to have nm^2 units
%        Resulting kappa will have units pN/nm
%                  alpha will have pN*s/nm
%
% Output:
%   kappa: estimated harmonic oscillator stiffness (pN/nm default units)
%   alpha: estimated drag coefficient (pN*s/nm default units)
%   k_ci, a_ci: +/- 95% confidence limit for the estimates

p = inputParser;

p.CaseSensitive = false;

addParameter(p,'KappaStart',0,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'AlphaStart',0,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'kbT',4.11,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'AlphaLower',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'AlphaUpper',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'KappaLower',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'KappaUpper',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>=0);

parse(p,varargin{:});

%Default kbT=4.11 pN*nm

kbT = p.Results.kbT;
alpha_0 = p.Results.AlphaStart;
kappa_0 = p.Results.KappaStart;

aLow = p.Results.AlphaLower;
aHigh =p.Results.AlphaUpper;
kLow = p.Results.KappaLower;
kHigh = p.Results.KappaUpper;

if aLow>=aHigh
    error('AlphaLow must be less than AlphaUpper');
end

if kLow>=kHigh
    error('KappaLow must be less than KappaUpper');
end


av = av/(2*kbT); %rescale to simplify equations


%rescale to be close to range 1-10;
ScaleFact = round(log10(nanmean(av)));
av = av*10^-ScaleFact;
alpha_0 = alpha_0*10^ScaleFact;
kappa_0 = kappa_0*10^ScaleFact;
aLow = aLow*10^ScaleFact;
aHigh = aHigh*10^ScaleFact;
kLow = kLow*10^ScaleFact;
kHigh = kHigh*10^ScaleFact;


if alpha_0<= aLow
    if isinf(aHigh)
        alpha_0 = aLow+0.1;
    else
        alpha_0 = (aHigh+aLow)/2;
    end
end
if kappa_0<= kLow
    if isinf(kHigh)
        kappa_0 = kLow+0.1;
    else
        kappa_0 = (kHigh+kLow)/2;
    end
end

    function [LLF,G,H] = nlogFn(p)
        a= p(1);
        k =p(2);
        
        KT = k*t;
        KTa = KT/a;
        EnKTa = exp(-KTa);
        En2KTa = EnKTa.^2;
        
        invT2 = t.^(-2);
        
        V = a^2/k^3.*invT2.*(KTa + 2*EnKTa - 1/2*En2KTa - 3/2);
        avHO_V = av./V;
        LLF = sum(Navg.*( avHO_V + log(V)));
        
        if nargout>1
            dVa = (k^-3).*invT2.*( (-3*a+KT) -En2KTa.*(a+KT) +2*EnKTa.*(2*a+KT));
            dVk = 0.5*(k^-4).*invT2.*a.*( 9*a -4*KT + En2KTa.*(3*a+2*KT) -4*EnKTa.*(3*a+KT) );
            
            X = Navg.*(1-avHO_V)./V;
            G = [ sum(X.*dVa);...
                  sum(X.*dVk)];
              if nargout>2
                  K2T2 = KT.^2;
                  dVaa = (a^-2)*(k^-3)*invT2.*(-3*a^2 +2*EnKTa.*(2*a^2+2*a*KT+K2T2) -En2KTa.*(a^2+2*a*KT+2*K2T2) );
                  dVkk = (k^-5)*invT2.*(-6*a*(3*a-KT) +2*EnKTa.*(12*a^2+6*a*KT+K2T2) -2*En2KTa.*(3*a^2+3*a*KT+K2T2) );
                  dVka = (k^-4)/a*invT2.*( a*(9*a-2*KT) -2*EnKTa.*(6*a^2+4*a*KT+K2T2) +En2KTa.*(3*a^2+4*a*KT+2*K2T2) );
                  
                  Y = Navg.*(2*avHO_V-1).*(V.^-2);
                  H(2,2) = sum( (dVk.^2).*Y + X.*dVkk );
                  H(1,1) = sum( (dVa.^2).*Y + X.*dVaa );
                  H(2,1) = sum( dVa.*dVk.*Y + X.*dVka );
                  H(1,2) = H(2,1);
                  
              end
        end
    end

    function Hout = HesFcn(p,~)
        [~,~,Hout] = nlogFn(p);
    end

options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','notify','HessianFcn',@HesFcn,'Algorithm','trust-region-reflective');
[x,~,~,~,~,~,hessian] = fmincon(@nlogFn,[alpha_0,kappa_0],[],[],[],[],[aLow,kLow],[aHigh,kHigh],[],options);

kappa = x(2)*10^-ScaleFact;
alpha = x(1)*10^-ScaleFact;

a_ci = 1.96/sqrt((hessian(1,1)))*10^-ScaleFact;
k_ci = 1.96/sqrt((hessian(2,2)))*10^-ScaleFact;

end