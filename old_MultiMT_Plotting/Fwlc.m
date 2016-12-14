function F = Fwlc(Lo,P,x)
x(x>=Lo) = NaN;
x(x<=0) = NaN;
F = 4.11./P.*(1/4*(1-x./Lo).^(-2)-1/4+x./Lo);