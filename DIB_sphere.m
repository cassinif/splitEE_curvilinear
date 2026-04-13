clear all
close all

% Example in Section 5.3

rng(2)

indic = true;

n = [100,50]; % ntheta nphi

tstar = 18;

epsilonpar = 20;
zeta1par = 10;
zeta2par = 10;
zeta3par = 1;
zeta4par = 48;
zeta5par = 0.5;
eta1par = 5;
eta2par = 2.5;
eta3par = 0.2;
eta4par = eta1par*(1-zeta5par)*(1-eta3par+eta3par*zeta5par)/(zeta5par*(1+eta3par*zeta5par));
eta5par = 1.5;
rhostar = 1.1653;

theta = linspace(0,2*pi,n(1)+1)';
theta = theta(1:n(1));
htheta = 2*pi/n(1);

sigmaphi = fzero(@(x) cot(x*pi/(n(2)-1+2*x))-2*(n(2)-1+2*x)/pi,1/2);
hphi = pi/(n(2)-1+2*sigmaphi);
phi1_gr = sigmaphi*hphi;
phi = linspace(phi1_gr,pi-phi1_gr,n(2))';

[THETA,PHI] = ndgrid(theta,phi);
X = rhostar*sin(PHI).*cos(THETA);
Y = rhostar*sin(PHI).*sin(THETA);
Z = rhostar*cos(PHI);

Atheta = toeplitz([-2,1,zeros(1,n(1)-2)]/htheta^2);
Atheta(1,n(1)) = 1/htheta^2;
Atheta(n(1),1) = 1/htheta^2;

avec_phi = -2/hphi^2*ones(n(2),1);
bvec_phi = 1./hphi^2+cot(phi1_gr+((1:n(2)-1)'-1)*hphi)/(2*hphi);
cvec_phi = flipud(bvec_phi);
Aphi = diag(avec_phi)+diag(bvec_phi,1)+diag(cvec_phi,-1);

gu = @(u,v) zeta1par*(zeta2par*(1-v).*u-zeta3par*u.^3-zeta4par*(v-zeta5par));
gv = @(u,v) zeta1par*(eta1par*(1+eta2par*u).*(1-v).*(1-eta3par*(1-v))-eta4par*v.*(1+eta3par*v).*(1+eta5par*u));

U0 = 1e-6*randn(n);
V0 = zeta5par + 1e-6*randn(n);

nsteps = 9000; 
tau = tstar/nsteps;

t = 0;

U = U0;
V = V0;

tic

deltavec = [1;sqrt(cumprod(cvec_phi./bvec_phi))];

S = (Aphi.*(deltavec.'))./deltavec;
[Q,lambda] = eig(S,'vector');

Qd = deltavec.*Q;
Qi = (Q./deltavec).'; 
PdirRu = Qd*(phi1((tau/rhostar^2)*lambda).*Qi);
PdirRv = Qd*(phi1((tau*epsilonpar/rhostar^2)*lambda).*Qi);

[Qth,lambdath] = eig(Atheta,'vector');
d1 = 1./(sin(phi).^2);
Lmat = lambdath*d1';
pldu = phi1((tau/rhostar^2)*Lmat);
pldv = phi1((tau/rhostar^2*epsilonpar)*Lmat);

if indic == true
  counter = 1;
  wphi = [phi(1)/2+hphi/2,hphi*ones(1,n(2)-2),hphi/2+phi(1)/2];
  ar = 4*pi*rhostar^2;  
  meanval(counter) = rhostar^2*htheta*(wphi*(sin(phi.').*sum(U,1)).')/ar;
end

for jj = 1:nsteps
  gUn = gu(U,V);
  gVn = gv(U,V);

  Fn = 1/rhostar^2*(Atheta*(U.*d1') + U*Aphi')+gUn;
  U = U + tau*Qth*(pldu.*(Qth'*(Fn*PdirRu')));
  Fn = epsilonpar/rhostar^2*(Atheta*(V.*d1') + V*Aphi')+gVn;
  V = V + tau*Qth*(pldv.*(Qth'*(Fn*PdirRv')));

  t = t + tau;
  
  if indic == true
    counter = counter + 1;
    meanval(counter) = rhostar^2*htheta*(wphi*(sin(phi.').*sum(U,1)).')/ar;
  end
end

toc

figure;
surf(X,Y,Z,U)
axis equal
colorbar
title(sprintf('u at time t = %.2f',t))
xlabel('x')
ylabel('y')
zlabel('z')
shading interp
drawnow

if indic == true
  figure;
  plot(0:tau:tstar,meanval)
  title('<u>')
  xlabel('t')
  ylabel('Integral mean')
  drawnow
end

