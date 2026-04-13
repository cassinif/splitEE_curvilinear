clear all
close all

% Example in Section 5.1

rng(2);

n = [25,50];  %nrho ntheta
%n = [40,80];

tstar = 1;

alpha1par = 0.899;
alpha2par = 0.2;
alpha3par = 0.2;
beta1par = -0.91;
beta2par = -alpha1par;
gammapar = 0.00387;
deltapar = 0.0075;
rhostar = 1;

hrho = rhostar/(n(1)-1/2);
rho1 = hrho/2;
rho = linspace(rho1,rhostar,n(1))';

theta = linspace(0,2*pi,n(2)+1)';
theta = theta(1:n(2));
htheta = 2*pi/n(2);

[RHO,THETA] = ndgrid(rho,theta);
X = RHO.*cos(THETA);
Y = RHO.*sin(THETA);

avec_rho = -2/hrho^2*ones(n(1),1);
bvec_rho = [2/hrho^2;1./hrho^2+1./((2*(2:n(1)-1)'-1)*hrho^2)];
cvec_rho = [2*(1:n(1)-2)'./((2*(1:n(1)-2)'+1)*hrho^2);2/hrho^2];
Arho = diag(avec_rho)+diag(bvec_rho,1)+diag(cvec_rho,-1);

Atheta = toeplitz([-2,1,zeros(1,n(2)-2)]/htheta^2);
Atheta(1,n(2)) = 1/htheta^2;
Atheta(n(2),1) = 1/htheta^2;

gu = @(u,v) alpha1par*u.*(1-alpha2par*v.^2)+v.*(1-alpha3par*u);
gv = @(u,v) beta1par*v.*(1+alpha1par*alpha2par/beta1par*u.*v)+u.*(beta2par+alpha3par*v);

U0 = rand(n)-1/2;
V0 = rand(n)-1/2;

deltavec = [1;sqrt(cumprod(cvec_rho./bvec_rho))];

S = (Arho.*(deltavec.'))./deltavec;

[Q,lambda] = eig(S,'vector');
Qd = deltavec.*Q;
Qi = (Q./deltavec).'; 

[Qth,lambdath] = eig(Atheta,'vector');
d1 = 1./(rho.^2);
Lmat = d1*lambdath';

mref = 50000;
mrange = [200:700:3000];

disp('Computing reference...')
m = mref;
disp(sprintf('m = %i',m))
tau = tstar/m;
U = U0;
V = V0;
t = 0;
PdirRu = Qd*(phi1((tau*gammapar)*lambda).*Qi);
PdirRv = Qd*(phi1((tau*deltapar)*lambda).*Qi);
pldu = phi1((tau*gammapar)*Lmat);
pldv = phi1((tau*deltapar)*Lmat);
for jj = 1:m
  gUn = gu(U,V);
  gVn = gv(U,V);

  Fn = gammapar*(Arho*U + d1.*(U*Atheta'))+gUn;
  U = U + tau*PdirRu*(pldu.*(Fn*Qth))*Qth';
  Fn = deltapar*(Arho*V + d1.*(V*Atheta'))+gVn;
  V = V + tau*PdirRv*(pldv.*(Fn*Qth))*Qth';
  t = t + tau;
end

Uref = U;
Vref = V;
nrmU = norm(Uref,'fro');
nrmV = norm(Vref,'fro');
disp('Reference computed!')

Lap = kron(speye(n(2)),sparse(Arho)) + kron(sparse(Atheta),spdiags(d1,0,n(1),n(1)));
Mu = gammapar*Lap;
Mv = deltapar*Lap;
Gfun = @(u,v) [gu(u,v);gv(u,v)];
Ffun = @(u,v) [Mu*u;Mv*v] + Gfun(u,v);
pn = prod(n);

counter = 0;
for m = mrange
  disp(sprintf('m = %i',m))
  counter = counter + 1;
  tau = tstar/m;

  disp('Forward Euler')
  tic
  uv = [U0(:);V0(:)];
  t = 0;
  for jj = 1:m
    uv = uv + tau*Ffun(uv(1:pn),uv(pn+1:2*pn));
    t = t + tau;
  end
  cpu_ee(counter) = toc;
  err_ee(counter) = norm([norm(Uref(:)-uv(1:pn),'fro')/nrmU,norm(Vref(:)-uv(pn+1:2*pn),'fro')/nrmV],2);

  disp('Exponential Euler')
  tic
  uv = [U0(:);V0(:)];
  t = 0;
  zz = zeros(2*pn,1);
  tol = 1e-6;
  for jj = 1:m
    Fn = Ffun(uv(1:pn),uv(pn+1:2*pn));
    uv = uv + kiops(tau,@(x) [Mu*x(1:pn);Mv*x(pn+1:2*pn)],[zz,Fn],tol,10,10,128,false);
    t = t + tau;
  end
  cpu_expe(counter) = toc;
  err_expe(counter) = norm([norm(Uref(:)-uv(1:pn),'fro')/nrmU,norm(Vref(:)-uv(pn+1:2*pn),'fro')/nrmV],2);

  disp('Split exponential Euler')
  tic
  U = U0;
  V = V0;
  t = 0;
  PdirRu = Qd*(phi1((tau*gammapar)*lambda).*Qi);
  PdirRv = Qd*(phi1((tau*deltapar)*lambda).*Qi);
  pldu = phi1((tau*gammapar)*Lmat);
  pldv = phi1((tau*deltapar)*Lmat);
  for jj = 1:m
    gUn = gu(U,V);
    gVn = gv(U,V);

    Fn = gammapar*(Arho*U + d1.*(U*Atheta'))+gUn;
    U = U + tau*PdirRu*(pldu.*(Fn*Qth))*Qth';
    Fn = deltapar*(Arho*V + d1.*(V*Atheta'))+gVn;
    V = V + tau*PdirRv*(pldv.*(Fn*Qth))*Qth';
    t = t + tau;
  end
  cpu_spexpe(counter) = toc;
  err_spexpe(counter) = norm([norm(Uref-U,'fro')/nrmU,norm(Vref-V,'fro')/nrmV],2);
end

figure
loglog(mrange,err_ee,'xb')
hold on
loglog(mrange,err_expe,'sm')
loglog(mrange,err_spexpe,'or')
loglog(mrange,err_spexpe(end)*(mrange/mrange(end)).^(-1),'--k')
legend('Forward Euler','Exponential Euler','Split Exponential Euler')
xlabel('m')
ylabel('error')

figure
loglog(cpu_ee,err_ee,'-xb')
hold on
loglog(cpu_expe,err_expe,'-sm')
loglog(cpu_spexpe,err_spexpe,'-or')
legend('Forward Euler','Exponential Euler','Split Exponential Euler')
xlabel('WC time')
ylabel('error')
