clear all
close all

% Example in Section 5.3

rng(2);

indic = true;

n = [160,160]; % nrho ntheta

tstar = 2.5;

alpha1par = 500;
alpha2par = 0.14;
beta1par = 1.34;
deltapar = 50;

lambdapar = -1.95;

rhostar = 1;

hrho = rhostar/(n(1)+(1-lambdapar)/2);
rho1 = (1-lambdapar)*hrho/2;
rho = linspace(rho1,rhostar,n(1)+1)';
rho = rho(1:n(1));

theta = linspace(0,2*pi,n(2)+1)';
theta = theta(1:n(2));
htheta = 2*pi/n(2);

[RHO,THETA] = ndgrid(rho,theta);
X = RHO.*cos(THETA);
Y = RHO.*sin(THETA);

avec_rho = -2./(hrho^(2+lambdapar)*((1-lambdapar)/2+(1:n(1))'-1).^lambdapar);
bvec_rho = ((1:n(1)-1)'-lambdapar)./(hrho^(2+lambdapar)*((1-lambdapar)/2+(1:n(1)-1)'-1).^(1+lambdapar));
cvec_rho = (1:n(1)-1)'./(hrho^(2+lambdapar)*((1-lambdapar)/2+(1:n(1)-1)').^(1+lambdapar));
Arho = diag(avec_rho)+diag(bvec_rho,1)+diag(cvec_rho,-1);

Atheta = toeplitz([-2,1,zeros(1,n(2)-2)]/htheta^2);
Atheta(1,n(2)) = 1/htheta^2;
Atheta(n(2),1) = 1/htheta^2;

gu = @(u,v) alpha1par*(alpha2par-u+u.^2.*v);
gv = @(u,v) alpha1par*(beta1par-u.^2.*v);

ustar = alpha2par + beta1par;
vstar = beta1par/((alpha2par+beta1par)^2);
U0 = ustar*ones(n)+1e-5*rand(n);
V0 = vstar*ones(n)+1e-5*rand(n);

nsteps = 25000;
tau = tstar/nsteps;

U = U0;
V = V0;

% Lifting
U = U - ustar;
V = V - vstar;

t = 0;

tic

deltavec = [1;sqrt(cumprod(cvec_rho./bvec_rho))];

S = (Arho.*(deltavec.'))./deltavec;
[Q,lambda] = eig(S,'vector');

Qd = deltavec.*Q;
Qi = (Q./deltavec).'; 
PdirRu = Qd*(phi1(tau*lambda).*Qi);
PdirRv = Qd*(phi1((tau*deltapar)*lambda).*Qi);

[Qth,lambdath] = eig(Atheta,'vector');
d1 = 1./(rho.^(2+lambdapar));
Lmat = d1*lambdath';
pldu = phi1(tau*Lmat);
pldv = phi1((tau*deltapar)*Lmat);

if indic == true
  counter = 1;
  wrho = [rho(1)/2+hrho/2,hrho*ones(1,n(1)-1)];
  ar = pi*rhostar^2;
  meanval(counter) = htheta*(wrho*(rho.*sum(U,2)))/ar + ustar;
end

for jj = 1:nsteps
  gUn = gu(U+ustar,V+vstar);
  gVn = gv(U+ustar,V+vstar);

  Fn = Arho*U + d1.*(U*Atheta')+gUn;
  U = U + tau*PdirRu*(pldu.*(Fn*Qth))*Qth';
  
  Fn = deltapar*(Arho*V + d1.*(V*Atheta'))+gVn;
  V = V + tau*PdirRv*(pldv.*(Fn*Qth))*Qth';
  
  t = t + tau;

  if indic == true
    counter = counter + 1;
    meanval(counter) = htheta*(wrho*(rho.*sum(U,2)))/ar + ustar;
  end
end

toc

% Unlifting
U = U + ustar;
V = V + vstar;

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

