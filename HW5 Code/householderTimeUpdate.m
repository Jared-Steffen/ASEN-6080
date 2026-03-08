function [Rbaru,Rbarux,butilde,Rbar,bbar] = householderTimeUpdate(Ru,bbaru,Rtilde,Gamma,bhat)
%{
Inputs:
    >Ru: square root of process noise matrix Q
    >bbaru:
    >Rtilde: square root of information matrix
    >Gamma: SNC process noise mapping matrix
    >bhat: whitened measurement

Outputs:
    >Rbaru: transformed Ru neeed for SRIF smoother
    >Rbarux: needed for the smoother
    >butilde: needed for the smoother
    >Rbar: transformed Rtilde needed for measurement update
    >bbar: transformed bbar needed for measurement update
%}

% Construct matrix
Atilde = [Ru zeros(size(Ru,1),size(Rtilde,1)) bbaru(:);
          -Rtilde*Gamma Rtilde bhat(:)];

% QR decomposition
A = qr(Atilde);

% Extract components
Rusize = size(Ru,1);
Rtildesize = size(Rtilde,1);
Rbaru = A(1:Rusize,1:Rusize);
Rbarux = A(1:Rusize,Rusize+1:Rusize+Rtildesize);
butilde = A(1:Rusize,end);
Rbar = A(Rusize+1:end,Rusize+1:Rusize+Rtildesize);
bbar = A(Rusize+1:end,end);

end


