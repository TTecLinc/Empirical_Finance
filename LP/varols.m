function [u beta a0]=varols(z, p)
global T vars

%	Generate the Dependent Variable Vector

y = z((p+1):T,:);

%	Generate the lags
x = ones(T-p,1);

i = 1;
if i > p;
    display('i>p');
else
	while i<=p;
		x = [x z((p+1-i):(T-i),:)];
		i = i+1;
	end
end


%	OLS on the VAR

beta = (x'*x)\(x'*y);
% A\=inv(A)

%	Get Structural residuals from reduced form

u = y - x*beta;
vcv = u'*u;

a0 = zeros(vars,vars);
% Cholesky Decomposition
% eye: Identity matrix
%A = spdiags(B,d,A) replaces the diagonals specified by d with the columns of B. The output is sparse.
%A(logical(eye(size(A))))=B
%m_0=eye(vars);

a0 = chol(vcv)/diag(diag(chol(vcv)));
e = u*a0;
end