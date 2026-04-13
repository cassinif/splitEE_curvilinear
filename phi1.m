function ret = phi1(z)
  idx = abs(z) < 1;
  ret = zeros(size(z));
  a = 1./factorial(17:-1:1);
  ret(idx) = polyval(a,z(idx));
  ret(~idx) = (exp(z(~idx))-1)./z(~idx);
end
