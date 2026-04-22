function [ytilde,Htilde2,Hctilde2] = whitenMeasModelConsiderCov(y,V,Htilde,Hctilde)
%{
Inputs:
    >y: unwhitened prefit residual
    >V: square root of measurement covariance matrix (upper triangular)
    >Htilde: linearized measurement matrix

Outputs:
    >ytilde: whitened measurement
    >Htilde2: whitened linearized measurement matrix
%}

ytilde = V\y;
Htilde2 = V\Htilde;
Hctilde2 = V\Hctilde;

end

