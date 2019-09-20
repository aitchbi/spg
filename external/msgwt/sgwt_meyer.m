% sgwt_meyer : Returns Meyer-like wavelet generating kernel
%
% function r=sgwt_meyer(x,varargin)
%
% defines function g(x) with 
% g(x) = sin(pi/2*V(x/w1-1)) for w1<x<2w1
% g(x) = cos(pi/2*V(x/2w1-1)) for 2w1<x<4w1
% where V is the Meyer wavelet auxiliary function
%
% Inputs :
% x : array of independent variable values
% t1 : transition region
%
% Outputs :
% r - result (same size as x)
%
% 2011, Nora Leonardi, Dimitri Van De Ville

function r=sgwt_meyer(x,varargin)
    control_params={'t1',1};
    argselectAssign(control_params);
    argselectCheck(control_params,varargin);
    argselectAssign(varargin);
    
    r=zeros(size(x));

    w1=2*t1/3;
    
    r1=find(x>=w1 & x<2*w1);
    r2=find(x>=2*w1 & x<2^2*w1);

    r(r1)=sin(pi/2*meyeraux(x(r1)/w1-1));
    r(r2)=cos(pi/2*meyeraux(x(r2)/(2*w1)-1));
end
