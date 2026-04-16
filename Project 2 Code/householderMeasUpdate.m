function [R,b] = householderMeasUpdate(Rbar,bbar,Htilde2,ytilde)
%{
Inputs:
    >Rbar: square root of information matrix
    >bbar: square root equivalent of the state update
    >ytilde: whitened measurement
    >Htilde2: whitened linearized measurement matrix

Outputs:
    >R: transformed Rbar
    >b: transformed bbar
%}

% Construct matrix
Atilde = [Rbar bbar(:);
          Htilde2 ytilde];

% QR decomposition
A = qr(Atilde);

% Extract components
Rsize = size(Rbar,1);
R = A(1:Rsize,1:Rsize);
b = A(1:Rsize,end);

end

