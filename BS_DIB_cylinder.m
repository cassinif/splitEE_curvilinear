clear all
close all

% Example in Section 5.6

rng(2);

indic = true;

n = [160,160,20]; % nrho ntheta nz

rhostar = 25;
zstar = 25;

tstar = 50;

deltapar = 1; 
epsilonpar = 20; 
alpha1par = 1;
alpha2par = 1;
alpha3par = 0.15;
beta1par = 1;
beta2par = 1;
beta3par = 0.15;
zeta1par = 1;
zeta2par = 10;
zeta3par = 1;
zeta4par = 66;
zeta5par = 0.5;
eta1par = 3;
eta2par = 2.5;
eta3par = 0.2;
eta4par = beta2par*(eta1par*(1-zeta5par)*(1-eta3par+eta3par*zeta5par))/(zeta5par*(1+eta3par*zeta5par));
eta5par = 1.5;

f1 = @(u) -alpha1par*(u-alpha2par);
f2 = @(v) -beta1par*(v-beta2par);
f3 = @(u,r,s) zeta1par*(zeta2par*u.*(1-s).*r-zeta3par*r.^3-zeta4par*(s-zeta5par));
f4 = @(v,r,s) zeta1par*(eta1par*v.*(1+eta2par*r).*(1-s).*(1-eta3par*(1-s))-eta4par*(1+eta5par*r).*s.*(1+eta3par*s));

hrho = rhostar/(n(1)-1/2);
rho1 = hrho/2;
rho = linspace(rho1,rhostar,n(1))';

theta = linspace(0,2*pi,n(2)+1)';
theta = theta(1:n(2));
htheta = 2*pi/n(2);

z = linspace(0,zstar,n(3)+1)';
z = z(1:n(3));
hz = zstar/n(3);

[RHO,THETA,Z] = ndgrid(rho,theta,z);
X = RHO.*cos(THETA);
Y = RHO.*sin(THETA);

avec_rho = -2/hrho^2*ones(n(1),1);
bvec_rho = [2/hrho^2;1./hrho^2+1./((2*(2:n(1)-1)'-1)*hrho^2)];
cvec_rho = [2*(1:n(1)-2)'./((2*(1:n(1)-2)'+1)*hrho^2);2/hrho^2];
Arho = diag(avec_rho)+diag(bvec_rho,1)+diag(cvec_rho,-1);

Atheta = toeplitz([-2,1,zeros(1,n(2)-2)]/htheta^2);
Atheta(1,n(2)) = 1/htheta^2;
Atheta(n(2),1) = 1/htheta^2;

Az = toeplitz([-2,1,zeros(1,n(3)-2)]/hz^2);
Az(1,1:2) = [-2,2]/(hz^2);

bcu = @(u,r,s) -(2/hz*alpha3par)*f3(u,r,s);
bcv = @(v,r,s) -(2/hz*beta3par*deltapar)*f4(v,r,s);

ustar = alpha2par;
vstar = beta2par;
rstar = 0;
sstar = zeta5par;

u0 = alpha2par*ones(n);
v0 = beta2par*ones(n);
r0 = rstar*ones(n(1),n(2)) + 1e-2*rand(n(1),n(2));
s0 = sstar*ones(n(1),n(2)) + 1e-2*rand(n(1),n(2));

nsteps = 8000; 
tau = tstar/nsteps;

t = 0;

tic
deltavec_rho = [1;sqrt(cumprod(cvec_rho./bvec_rho))];

S_rho = (Arho.*(deltavec_rho.'))./deltavec_rho;
[Q_rho,lambda_rho] = eig(S_rho,'vector');
Qd_rho = deltavec_rho.*Q_rho;
Qi_rho = (Q_rho./deltavec_rho).'; 
PPrho_u = Qd_rho*(phi1(tau*lambda_rho).*Qi_rho);
PPrho_v = Qd_rho*(phi1((tau*deltapar)*lambda_rho).*Qi_rho);


[Qth,lambdath] = eig(Atheta,'vector');
Qtht = Qth';
drho = 1./(rho.^2);
dp_u = phi1(tensorize(drho,lambdath,tau*ones(n(3),1)));
dp_v = phi1(tensorize(drho,lambdath,(deltapar*tau)*ones(n(3),1)));

deltavec_z = [1;sqrt(2)/2*ones(n(3)-1,1)];
S_z = (Az.*(deltavec_z.'))./deltavec_z;
[Q_z,lambda_z] = eig(S_z,'vector');
Qd_z = deltavec_z.*Q_z;
Qi_z = (Q_z./deltavec_z).'; 
PPz_u = Qd_z*(phi1(tau*lambda_z).*Qi_z);
PPz_v = Qd_z*(phi1((tau*deltapar)*lambda_z).*Qi_z);

Ppr = phi1((tau*drho)*lambdath.');
Pps = phi1((tau*epsilonpar*drho)*lambdath.');
PPrho_r = PPrho_u;
PPrho_s = Qd_rho*(phi1((tau*epsilonpar)*lambda_rho).*Qi_rho);

Lu = @(u) mump(u,Arho,1) + mump(drho.*u,Atheta,2) + mump(u,Az,3);
Lv = @(v) deltapar*Lu(v);
Lr = @(r) Arho*r + (drho.*r)*Atheta.';
Ls = @(s) epsilonpar*Lr(s);

u = u0;
v = v0;
r = r0;
s = s0;

% Lifting
u = u - alpha2par;
v = v - beta2par;

if indic == true
  counter = 1;
  wrho = hrho*[1/4+1/2,ones(1,n(1)-2),1/2];
  ar = pi*rhostar^2;
  meanval(counter) = htheta*(wrho*(rho.*sum(r,2)))/ar;
end

for jj = 1:nsteps
  gun = f1(u+alpha2par);
  gun(:,:,1) = gun(:,:,1) + bcu(u(:,:,1)+alpha2par,r,s);

  gvn = f2(v+beta2par);
  gvn(:,:,1) = gvn(:,:,1) + bcv(v(:,:,1)+beta2par,r,s);

  grn = f3(u(:,:,1)+alpha2par,r,s);

  gsn = f4(v(:,:,1)+beta2par,r,s);

  tmp = Lu(u) + gun;
  u = u + tau*mump(act_2(mump(tmp,PPz_u,3),Qth,Qtht,dp_u),PPrho_u,1);

  tmp = Lv(v) + gvn;
  v = v + tau*mump(act_2(mump(tmp,PPz_v,3),Qth,Qtht,dp_v),PPrho_v,1);
  
  tmp = Lr(r) + grn;
  r = r + tau*(PPrho_r*(Ppr.*(tmp*Qth))*Qth');
  
  tmp = Ls(s) + gsn;
  s = s + tau*(PPrho_s*(Pps.*(tmp*Qth))*Qth');
  
  if indic == true
    counter = counter + 1;
    meanval(counter) = htheta*(wrho*(rho.*sum(r,2)))/ar;
  end

  t = t + tau;

end
toc

% Unlifting
u = u + alpha2par;
v = v + beta2par;

figure
pcolor(X(:,:,1),Y(:,:,1),r)
shading interp
axis equal
colorbar
title(sprintf('r at time t = %.4f',t))
xlabel('x')
ylabel('y')
drawnow

if indic == true
  figure
  plot(0:tau:tstar,meanval,':r')
  xlabel('t')
  ylabel('Integral mean')
  title('<r>')
  drawnow
end

figure;
pcolor(X(:,:,1),Y(:,:,1),u(:,:,1))
shading interp
axis equal
colorbar
title(sprintf('u at z = 0 at time t = %.4f',t))
xlabel('x')
ylabel('y')
drawnow
