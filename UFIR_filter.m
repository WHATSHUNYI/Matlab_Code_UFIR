function [m, G] = UFIR_filter(m,A,H,B,u,G,Y)

Pm     = A * m + B * u;
%G     =  (H'* H + (A * G * A')^-1)^-1;

G = A * G *A'-A * G *A'*(eye(size(H'*H *A * G *A'))+ H'*H *A * G *A')^-1* H'* H *A * G *A';
Gain  =  G * H';

m     =  Pm +  Gain * (Y - H * Pm);

end