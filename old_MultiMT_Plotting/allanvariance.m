function [av,M,Navg] = allanvariance(X,M)
% allanvariance(X,dT): Calculate the Allan Variance
% Calculates the Alan Variance of a time-ordered signal.
% The algorithm is based on the method presented in:
%   Lansdorp et al. Rev. Sci Inst. 83 (2012)
% The calculation is performed using:
%   av_m = 1/(2(N-2m)(m*tc)^2)*sum(cumX(k+2m)-2ComX(k+m)+CumX(k))^2
% and CumX is the cumulative sum of X
% 
% Input:
%  X: time-ordered signal
%     if X is a MxN matrix then it is interpreted as a set of N
%     time-ordered signals.
%  M: optionally specify the m values to calculate av over
%     default: M = (2.^(0:floor(log2(N)-1)))'; Octave Sampling
%              where N is the number of time-points in X
%
% Output:
%   av: the Allan Variance (sigma(m)^2)
%   M: the time-average width m used for calculating each simga(m)^2
%   Navg: the number of differences used to calculate each sigma(m)^2
%           Navg=(N-2*m+1)
%   Navg can be used as a weighting factor when calculating fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016, Daniel T. Kovari, Emory University
%% Change Log:
%   2016-12-08: Initial file creation
%               Algorithm is based on Lansdorp et al. Rev. Sci Inst. 83 (2012)
%               Except that Eq 18 in Lansdorp does not use the cumlative in
%               the sum, which it should.

%% 
if isrow(X)
    X = X';
end
N = size(X,1);

%M = [1,2:2:floor(N/2)]';
if nargin<2
    M = (2.^(0:floor(log2(N)-1)))';
end

X = [zeros(1,size(X,2));cumsum(X,1)];

av = NaN(numel(M),size(X,2));
for n=1:numel(M)
    m=M(n);
    av(n,:) = 1/2/m^2/(N-2*m+1)*...
        sum( ( X(2*m+1:N+1,:) -2*X(m+1:N-m+1,:) +X(1:N-2*m+1,:) ).^2,1);
end

if nargout>2
    Navg = N-2*M+1;
end