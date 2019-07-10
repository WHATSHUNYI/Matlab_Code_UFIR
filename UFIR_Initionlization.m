% Please cite the paper concerning "Bayesian state estimation on finite
% horizons: The case of linear state®Cspace model" when you use these code.
% A: state transaction matrix
% B: control input matrix
% U: Control Input
% Y: Measurment
% Q: Process noise covariance
% R: Measurment noise coariance
% H: Measurment matrix
% State-space mode
% x = A * x + B *u + w
% y = H * x + v
% Please note Q and R are only invovled in the computation of error covariance.
function [m,G,P] = UFIR_Initionlization(A,H,B,U,Y,Q,R)
    dim_x = size(A,1);
    dim_u = size(B,2);
    ExH = H * A; 
    ExB = B; 
    ExQ = Q;
    ExR = R;
    ExE = eye(dim_x);
    DiagH = H;

 for i = 2: dim_x
        ExH = [H*A^i; ExH];
        ExB = blkdiag(B,ExB);
        ExE = blkdiag(eye(dim_x),ExE);
        ExQ = blkdiag(Q,ExQ);
        ExR = blkdiag(R,ExR);
        DiagH = blkdiag(H,DiagH);
  end
    
 N= dim_x;  
 for j=1:N-1                                                          %Bæÿ’Û
    for i=1:j
        W=eye(dim_x);
        for h=i:j
            W= W * A;
        end
       ExB(:,j*dim_u+1:(j+1)*dim_u)=ExB(:,j*dim_u+1:(j+1)*dim_u)+[zeros((i-1)*dim_x,dim_u);W*B;zeros((N-i)*dim_x,dim_u)];   
    end  
 end  
 
N= dim_x;  
for j=1:N-1                                                          % Eæÿ’Û
    for i=1:j
        W=eye(dim_x);
        for h=i:j
            W= W * A;
        end
       ExE(:,j*dim_x+1:(j+1)*dim_x)=ExE(:,j*dim_x+1:(j+1)*dim_x)+[zeros((i-1)*dim_x,dim_x);W;zeros((N-i)*dim_x,dim_x)];   
    end  
end  
 
 
    ExEE = DiagH * ExE;                                                       % H matric in paper
    L = DiagH * ExB;   
    G = A ^(dim_x) * (ExH' * ExH) ^-1 * (A^(dim_x))';                    % Noise Power Gain
    Gain = A^(dim_x) * (ExH' * ExH) ^-1 * ExH';                            % Bias Correction 
    m    = Gain * Y + (ExB(1:dim_x,:) - Gain * L) * U;                       % State Estimate
    P    =  (Gain * ExEE  -ExE(1:dim_x,:) ) * ExQ * (Gain * ExEE  -ExE(1:dim_x,:))'+  Gain * ExR * Gain';                   % Noise Covaraince
end