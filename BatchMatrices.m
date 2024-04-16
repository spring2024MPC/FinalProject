function [Ab,Bb,Qb,Rb] = BatchMatrices(A,B,Q,R,P,N)
% function [Ab,Bb,Qb,Rb] = BatchMatrices(A,B,Q,R,P,N)
%
% This function computes the batch cost and dynamics matrices for
% unconstrained receeding horizon control with a quadratic objective and
% LTI dynamics.
%
% Inputs
%   A,B   = discrete time dynamics and input matrix for plant
%   Q,R   = running cost, x'Qx + u'Ru
%   P     = terminal cost x'Px
%   N     = horizon
%
% Outputs
%   Ab,Bb = batch dynamics matrices designed so that
%           xb = Ab x(t) + Bb ub,
%           where
%               xb      = [x_0'; x_1';...;x_N'],
%               ub      = [u_0'; u_1';...;u_(N-1)'],
%               x_0     = x(t), and
%               x_(i+1) = A x_i + B u_i.
%  Qb,Rb = batch cost matrices designed so that
%           xb'Qbxb + ub'Rub = x_N'Px_N + sum_(i=0)^(N-1) ( x_i'Qx_i + u_i'Ru_i )
 
n = size(A,1); % # states
m = size(B,2); % # inputs
 
%% Make the Q and R matrices
% kron and blkdiag are helpful here.
% Qb =  [Q 0 ...    0]
%       [0 Q ...    0]
%       [    ...     ]
%       [0   ...0 Q 0]
%       [0   ...0 0 P]
Qb = blkdiag(kron(eye(N), Q), P);
% Rb =  [R 0 ...  0]
%       [0 R ...  0]
%       [    ...   ]
%       [0   ...0 R]
Rb = kron(eye(N), R);
 
%% Make A and B matrices
Ab = zeros((N+1)*n,n);  % pre-allocate space for Ab
Ab(1:n,:) = eye(n);     % initialization
Bb = zeros((N+1)*n,N*m);% pre-allocate space for Bb
 
 
% Ab =  [ A^0 ], Bb =   [ 0         0   0 ...      0] <- block-row 0, rows 1 to n
%       [ A^1 ]         [ B         0   0 ...      0] <- block-row 1, rows n+1 to 2n
%       [ ... ]         [ AB        B   0 ...      0] <- block-row 2, rows 2n+1 to 3n
%       [ ... ]         [           ...             ]
%       [ A^N ]         [ A^(N-1)B  ...         AB B]\
Bb(n+1:2*n,1:m) = B;
Bprev = Bb(n+1:2*n, 1:m);
for iter = 1:N % loop over rows in Ab and Bb
    ri = iter*n + 1; % first row in Ab and Bb to be edited
    rf = (iter+1)*n; % last row in Ab to be edited
    % Write [[block-row]] iter of Ab in terms of previous block-row.
    Ab(ri:rf,:) = A^iter;
    % Write [[block-column]] iter of Bb in terms of previous block-column,
    % but go back to front.
    if iter ~= N
        %Bb((ri+n):(rf+n),:) = [A^(iter-1)*B,Bprev];
        Bb(ri:rf,1:(iter+1)*m) = [A^(iter-1)*B,Bprev];
        Bprev = Bb(ri:rf,1:(iter+1)*m);
    end
end

