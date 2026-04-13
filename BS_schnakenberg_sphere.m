clear all
close all

% Example in Section 5.5

rng(2);

indic = true;

n = [30,50,30]; % nrho ntheta nphi

tstar = 20;

alpha1par = 55;
alpha2par = 0.1;
beta1par = 0.9;
zeta1par = 55;
zeta2par = 5/12;
zeta3par = 5/12;
eta1par = 5;
eta2par = 5;
deltapar = 10;
epsilonpar = 10;
f = @(u,v) alpha2par-u+u.^2.*v;
g = @(u,v) beta1par-u.^2.*v;
h1 = @(u,v,r,s) zeta2par*r-zeta3par*u;
h2 = @(u,v,r,s) eta1par*s-eta2par*v;

rhostar = 1;

hrho = rhostar/n(1);
rho = linspace(0,rhostar,n(1)+1)';
rho = rho(2:n(1)+1);

theta = linspace(0,2*pi,n(2)+1)';
theta = theta(1:n(2));
htheta = 2*pi/n(2);

sigmaphi = fzero(@(x) cot(x*pi/(n(3)-1+2*x))-2*(n(3)-1+2*x)/pi,1/2);
hphi = pi/(n(3)-1+2*sigmaphi);
phi1_gr = sigmaphi*hphi;
phi = linspace(phi1_gr,pi-phi1_gr,n(3))';

[RHO,THETA,PHI] = ndgrid(rho,theta,phi);
X = RHO.*sin(PHI).*cos(THETA);
Y = RHO.*sin(PHI).*sin(THETA);
Z = RHO.*cos(PHI);

avec_rho = -2/hrho^2*ones(n(1),1);
bvec_rho = ((1:n(1)-1)'+1)./((1:n(1)-1)'*hrho^2);
cvec_rho = [(1:n(1)-2)'./(((1:n(1)-2)'+1)*hrho^2);2/hrho^2];
Arho = diag(avec_rho)+diag(bvec_rho,1)+diag(cvec_rho,-1);

Atheta = toeplitz([-2,1,zeros(1,n(2)-2)]/htheta^2);
Atheta(1,n(2)) = 1/htheta^2;
Atheta(n(2),1) = 1/htheta^2;

avec_phi = -2/hphi^2*ones(n(3),1);
bvec_phi = 1./hphi^2+cot(phi1_gr+((1:n(3)-1)'-1)*hphi)/(2*hphi);
cvec_phi = flipud(bvec_phi);
Aphi = diag(avec_phi)+diag(bvec_phi,1)+diag(cvec_phi,-1);

bcu = @(u,v,r,s) (2/hrho+2)*zeta1par*h1(u,v,r,s);
bcv = @(u,v,r,s) (2/hrho+2)*zeta1par*h2(u,v,r,s);

gutmp = @(u,v) alpha1par*f(u,v);
gvtmp = @(u,v) alpha1par*g(u,v);
gr = @(r,s,u,v) zeta1par*(f(r,s)-h1(u,v,r,s));
gs = @(r,s,u,v) zeta1par*(g(r,s)-h2(u,v,r,s));

ustar = alpha2par + beta1par;
vstar = beta1par/((alpha2par + beta1par)^2);
rstar = alpha2par + beta1par;
sstar = beta1par/((alpha2par + beta1par)^2);

u0 = ustar*ones(n) + 1e-3*randn(n);
v0 = vstar*ones(n) + 1e-3*randn(n);
r0 = rstar*ones(n(2),n(3)) + 1e-3*randn(n(2),n(3));
s0 = sstar*ones(n(2),n(3)) + 1e-3*randn(n(2),n(3));

nsteps = 200000;
tau = tstar/nsteps;

t = 0;

tic

deltavec_rho = [1;sqrt(cumprod(cvec_rho./bvec_rho))];
S_rho = (Arho.*(deltavec_rho.'))./deltavec_rho;
[Q_rho,lambda_rho] = eig(S_rho,'vector');
Qd_rho = deltavec_rho.*Q_rho;
Qi_rho = (Q_rho./deltavec_rho).'; 

[Qth,lambdath] = eig(Atheta,'vector');
Qtht = Qth';

deltavec_phi = [1;sqrt(cumprod(cvec_phi./bvec_phi))];
S_phi = (Aphi.*(deltavec_phi.'))./deltavec_phi;
[Q_phi,lambda_phi] = eig(S_phi,'vector');
Qd_phi = deltavec_phi.*Q_phi;
Qi_phi = (Q_phi./deltavec_phi).'; 

drho = 1./rho.^2;
dphi = 1./(sin(phi).^2);
dphi_row = dphi.';
dphi_pag = reshape(dphi,1,1,n(3));

Prhou = Qd_rho*(phi1((tau*lambda_rho)).*Qi_rho);
Pthu = phi1(tensorize(tau*drho,lambdath,dphi));
Pphiu = phi1(tensorize(drho,tau*ones(n(2),1),lambda_phi));

Lu = @(u) mump(u,Arho,1) + dphi_pag.*mump(drho.*u,Atheta,2) + mump(drho.*u,Aphi,3);

Arhov = deltapar*Arho;
drhov = deltapar*drho;

Prhov = Qd_rho*(phi1(((tau*deltapar)*lambda_rho)).*Qi_rho);
Pthv = phi1(tensorize(tau*drhov,lambdath,dphi));
Pphiv = phi1(tensorize(drhov,tau*ones(n(2),1),lambda_phi));

Lv = @(v) mump(v,Arhov,1) + dphi_pag.*mump(drhov.*v,Atheta,2) + mump(drhov.*v,Aphi,3);

Pphir = Qd_phi*(phi1(((tau/rhostar^2)*lambda_phi)).*Qi_phi);
Pthr = phi1(tensorize(tau/rhostar^2*lambdath,dphi));
Lr = @(r) (r*Aphi' + Atheta*(r.*dphi_row))/rhostar^2;

Pphis = Qd_phi*(phi1(((tau*epsilonpar/rhostar^2)*lambda_phi)).*Qi_phi);
Pths = phi1(tensorize(tau*epsilonpar/rhostar^2*lambdath,dphi));
Ls = @(s) (s*Aphi' + Atheta*(s.*dphi_row))*(epsilonpar/rhostar^2);

u = u0;
v = v0;
r = r0;
s = s0;

if indic == true
  counter = 1;
  wrho = hrho*[1,ones(1,n(1)-2),1/2];
  wphi = [phi(1)/2+hphi/2,hphi*ones(1,n(3)-2),hphi/2+phi(1)/2];
  ar = 4/3*pi*rhostar^3;
  ar_r = 4*pi*rhostar^2;
  meanval(counter) = htheta*(wphi*(sin(phi).*(wrho*(rho.^2.*squeeze(sum(u,2)))).'))/ar;
  meanval_r(counter) = rhostar^2*htheta*(wphi*(sin(phi.').*sum(r,1)).')/ar_r;
end

for jj = 1:nsteps
  gun = gutmp(u,v);
  gun(n(1),:,:) = gun(n(1),:,:) + bcu(u(n(1),:,:),v(n(1),:,:),reshape(r,1,n(2),n(3)),reshape(s,1,n(2),n(3)));

  gvn = gvtmp(u,v);
  gvn(n(1),:,:) = gvn(n(1),:,:) + bcv(u(n(1),:,:),v(n(1),:,:),reshape(r,1,n(2),n(3)),reshape(s,1,n(2),n(3)));

  grn = gr(r,s,squeeze(u(n(1),:,:)),squeeze(v(n(1),:,:)));

  gsn = gs(r,s,squeeze(u(n(1),:,:)),squeeze(v(n(1),:,:)));

  tmp = Lu(u) + gun;
  u = u + tau*act_1(act_2(act_3(tmp,Qd_phi,Qi_phi,Pphiu),Qth,Qtht,Pthu),Prhou);

  tmp = Lv(v) + gvn;
  v = v + tau*act_1(act_2(act_3(tmp,Qd_phi,Qi_phi,Pphiv),Qth,Qtht,Pthv),Prhov);

  tmp = Lr(r) + grn;
  r = r + tau*Qth*(Pthr.*(Qth'*(tmp*Pphir')));

  tmp = Ls(s) + gsn;
  s = s + tau*Qth*(Pths.*(Qth'*(tmp*Pphis')));

  t = t + tau;

  if indic == true
    counter = counter + 1;
    meanval(counter) = htheta*(wphi*(sin(phi).*(wrho*(rho.^2.*squeeze(sum(u,2)))).'))/ar;
    meanval_r(counter) = rhostar^2*htheta*(wphi*(sin(phi.').*sum(r,1)).')/ar_r;
  end

end
toc

figure;
surf(squeeze(X(n(1),:,:)),squeeze(Y(n(1),:,:)),squeeze(Z(n(1),:,:)),squeeze(u(n(1),:,:)))
shading interp
axis equal
colorbar
title(sprintf('u at rho = %.2f at time t = %.2f',rhostar,t))
xlabel('x')
ylabel('y')
zlabel('z')
drawnow

figure;
surf(squeeze(X(n(1),:,:)),squeeze(Y(n(1),:,:)),squeeze(Z(n(1),:,:)),r)
shading interp
axis equal
colorbar
title(sprintf('r at time t = %.2f',t))
xlabel('x')
ylabel('y')
zlabel('z')
drawnow

if indic == true
  figure;
  plot(0:tau:tstar,meanval)
  hold on
  plot(0:tau:tstar,meanval_r,':r')
  xlabel('t')
  ylabel('Integral mean')
  title('<u> and <r>')
  legend('<u>','<r>')
  drawnow
end
