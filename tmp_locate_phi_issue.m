addpath(genpath('irregularWavesMF12'));
g=9.81; h=100; Ux=0; Uy=0; Lx=3000; Ly=3000; Nx=96; Ny=96; t=0; N=60;
rng(1234);

dkx=2*pi/Lx; dky=2*pi/Ly;
kx_idx_all=(-floor(Nx/2)):(ceil(Nx/2)-1);
ky_idx_all=(-floor(Ny/2)):(ceil(Ny/2)-1);
[KXI,KYI]=meshgrid(kx_idx_all,ky_idx_all);
kx_all=KXI(:)*dkx; ky_all=KYI(:)*dky;
kmag_all=hypot(kx_all,ky_all); theta_all=atan2(ky_all,kx_all);
keep=(kx_all>0)|(kx_all==0 & ky_all>0);
kx_all=kx_all(keep); ky_all=ky_all(keep); kmag_all=kmag_all(keep); theta_all=theta_all(keep);
Tp=10; kp=(2*pi/Tp)^2/g; theta0=pi/4; sig_k=0.12*kp; sig_t=deg2rad(12);
W=exp(-0.5*((kmag_all-kp)/sig_k).^2).*exp(-0.5*(angle(exp(1i*(theta_all-theta0)))/sig_t).^2);
[~,idx]=sort(W,'descend'); idx=idx(1:N);
kx=kx_all(idx); ky=ky_all(idx);
amp0=0.04*W(idx)/max(W(idx));
xf=Lx/2; yf=Ly/2;
om=sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
ph=-(kx*xf+ky*yf)+om*t;
a=amp0.*cos(ph); b=amp0.*sin(ph);

x=(0:Nx-1)*(Lx/Nx); y=(0:Ny-1)*(Ly/Ny); [X,Y]=meshgrid(x,y);

c2=coeffsMF12(2,g,h,a,b,kx,ky,Ux,Uy);
c3=coeffsMF12(3,g,h,a,b,kx,ky,Ux,Uy);

[~,p2d]=surfaceMF12_new(2,c2,X,Y,t);
[~,p3d]=surfaceMF12_new(3,c3,X,Y,t);
[~,p2f]=surfaceMF12_new_fixphi(2,c2,X,Y,t);
[~,p3f]=surfaceMF12_new_fixphi(3,c3,X,Y,t);
[~,p2s]=surfaceMF12_spectral(c2,Lx,Ly,Nx,Ny,t);
[~,p3s]=surfaceMF12_spectral(c3,Lx,Ly,Nx,Ny,t);

phi3d = p3d-p2d;
phi3f = p3f-p2f;
phi3s = p3s-p2s;

tmpd=phi3d-phi3s; d_orig=max(abs(tmpd(:)));
r_orig=d_orig/max(abs(phi3d(:)));
tmpf=phi3f-phi3s; d_fix=max(abs(tmpf(:)));
r_fix=d_fix/max(abs(phi3f(:)));

fprintf('orig phi3 rel diff = %.4f%%\n',100*r_orig);
fprintf('fix  phi3 rel diff = %.4f%%\n',100*r_fix);

