function [x, w] = gauss_quadrature(n, up, down, H, cH, weighting)
% Gives the nodes and weights for Gaussian approximation.
  i_ary = 0:(2*n);
  if strcmp(weighting,'Fraction')
     moment = cH./(i_ary-H+0.5).*up.^(i_ary-H+0.5);
  elseif strcmp(weighting,'Legendre')
     moment = 1./(i_ary + 1).*(up.^(i_ary+1) - down.^(i_ary+1));
  end

  h = zeros ( n + 1, n + 1 );
  for i = 1 : n + 1
    for j = 1 : n + 1
      h(i,j) = moment(i+j-1);
    end
  end

  r = chol ( h );
  alpha = zeros ( n, 1 );
  alpha(1) = r(1,2) / r(1,1);
  for i = 2 : n
    alpha(i) = r(i,i+1) / r(i,i) - r(i-1,i) / r(i-1,i-1);
  end

  beta = zeros ( n - 1, 1 );
  for i = 1 : n - 1
    beta(i) = r(i+1,i+1) / r(i,i);
  end
  jacobi = diag ( alpha, 0 ) + diag ( beta, -1 ) + diag ( beta, +1 );
  [ evec, eval ] = eig ( jacobi );

  x = diag ( eval );
  w = moment(1) * evec(1,:).^2;

  x = x(:);
  w = w(:);
  out = [x, w];

end

