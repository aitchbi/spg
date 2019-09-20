% sgwt_meyer_end : Returns Meyer-like wavelet generating kernel
% keep filter constant at 1 at fine scale
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

function r=sgwt_meyer_end(x,varargin)
    control_params={'t1',1};
    argselectAssign(control_params);
    argselectCheck(control_params,varargin);
    argselectAssign(varargin);
    
    r=zeros(size(x));
    
    w1=2*t1/3;
    
    r1=find(abs(x)>=w1 & abs(x)<2*w1);
    r2=find(abs(x)>=2*w1);

    r(r1)=sin(pi/2*meyeraux(abs(x(r1))/w1-1));
    r(r2)=1;
end
