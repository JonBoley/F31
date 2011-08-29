function x = lscov3(A,b,V)

[nobs,nvar] = size(A); % num observations, num predictor variables
nrhs = size(b,2);      % num right-hand sides
rankV = nobs;




% V given as a weight vector.
D = sqrt(full(V(:)));
b = D.*b;
C1=A(:,1).*D;
C2=A(:,2).*D;
A(:,1)=C1;
A(:,2)=C2;

[Q,R,perm] = qr(A,0);
z = Q'*b;

%    if issparse(R)
%        % Use R to remove dependent columns in A (c.f. mldivide).
%        diagR = abs(diag(R));
%        keepCols = (diagR > 20.*max(diagR).*eps(class(R)).*(nobs+nvar));
%    else
%        % Use the rank-revealing QR to remove dependent columns of A.
%        keepCols = (abs(diag(R)) > abs(R(1)).*max(nobs,nvar).*eps(class(R)));
%    end
%    rankA = sum(keepCols);
%    if rankA < nvar
%        warning('MATLAB:lscov:RankDefDesignMat', 'A is rank deficient to within machine precision.');
%        R = R(keepCols,keepCols);
%        z = z(keepCols,:);
%        perm = perm(keepCols);
%   end

    % Compute the LS coefficients, filling in zeros in elements corresponding
    % to rows of R that were thrown out.
    xx = R \ z;
    if issparse(xx)
        x = sparse(nvar,nrhs);
    else
        x = zeros(nvar,nrhs);
    end
    x(perm,1:nrhs) = xx;

