clear all
close all

% Example in Section 5.2

rng(2);

indic = true;

n = [40,80]; % nrho ntheta

tstar = 1600;

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

nsteps = 10000;
tau = tstar/nsteps;

U = U0;
V = V0;
t = 0;

tic

deltavec = [1;sqrt(cumprod(cvec_rho./bvec_rho))];

S = (Arho.*(deltavec.'))./deltavec;

[Q,lambda] = eig(S,'vector');
Qd = deltavec.*Q;
Qi = (Q./deltavec).'; 
PdirRu = Qd*(phi1((tau*gammapar)*lambda).*Qi);
PdirRv = Qd*(phi1((tau*deltapar)*lambda).*Qi);

[Qth,lambdath] = eig(Atheta,'vector');
d1 = 1./(rho.^2);
Lmat = d1*lambdath';
pldu = phi1((tau*gammapar)*Lmat);
pldv = phi1((tau*deltapar)*Lmat);

if indic == true
  counter = 1;
  wrho = hrho*[1/4+1/2,ones(1,n(1)-2),1/2];
  ar = pi*rhostar^2;
  meanval(counter) = htheta*(wrho*(rho.*sum(U,2)))/ar;
end

for jj = 1:nsteps
  gUn = gu(U,V);
  gVn = gv(U,V);

  Fn = gammapar*(Arho*U + d1.*(U*Atheta'))+gUn;
  U = U + tau*PdirRu*(pldu.*(Fn*Qth))*Qth';
  Fn = deltapar*(Arho*V + d1.*(V*Atheta'))+gVn;
  V = V + tau*PdirRv*(pldv.*(Fn*Qth))*Qth';

  t = t + tau;

  if indic == true
    counter = counter + 1;
    meanval(counter) = htheta*(wrho*(rho.*sum(U,2)))/ar;
  end

end

toc

figure;
pcolor(X,Y,U);
shading interp
axis equal
colorbar
title(sprintf('u at time t = %.2f',t))
xlabel('x')
ylabel('y')
drawnow

if indic == true
  figure
  plot(0:tau:tstar,meanval)
  title('<u>')
  xlabel('t')
  ylabel('Integral mean')
  drawnow
end
