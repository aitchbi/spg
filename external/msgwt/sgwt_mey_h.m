% sgwt_mey_h : Returns Meyer-like scaling function
%
% function r=sgwt_mey_h(x,varargin)
%
% defines function h(x) with 
% h(x) = 1 for 0<x<w1
% h(x) = cos(pi/2*V(x/w1-1)) for w1<x<2w1
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

function r=sgwt_mey_h(x,varargin)
    control_params={'t1',1};
    argselectAssign(control_params);
    argselectCheck(control_params,varargin);
    argselectAssign(varargin);
    
    r=zeros(size(x));

    w1=2*t1/3;
    
    r1= x>=0 & x<w1;
    r2= x>=w1 & x<2*w1;

    r(r1)=1;
    r(r2)=cos(pi/2*meyeraux(x(r2)/w1-1));
end
