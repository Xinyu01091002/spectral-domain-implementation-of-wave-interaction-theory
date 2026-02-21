addpath(genpath(fullfile(pwd,'irregularWavesMF12')));

% Same wave-group setup as plot_phi3_wavegroup_lines.m
g=9.81; h=100; Ux=0; Uy=0; Lx=3000; Ly=3000; Nx=128; Ny=128; t=0; N=120;
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
Sk=exp(-0.5*((kmag_all-kp)/sig_k).^2); St=exp(-0.5*(angle(exp(1i*(theta_all-theta0)))/sig_t).^2);
W=Sk.*St; [~,idxs]=sort(W,'descend'); idx=idxs(1:N);
kx=kx_all(idx); ky=ky_all(idx);
amp=0.04*W(idx)/max(W(idx)); xf=Lx/2; yf=Ly/2;
omega_lin=sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phase=-(kx*xf+ky*yf)+omega_lin*t; a=amp.*cos(phase); b=amp.*sin(phase);

% Diagnostic 1: grid mismatch
mx=max(abs(kx/dkx-round(kx/dkx)));
my=max(abs(ky/dky-round(ky/dky)));
fprintf('grid mismatch max: kx %.3e, ky %.3e\n',mx,my);

% Diagnostic 2: Nyquist overflow for 3rd-order combination wavenumbers
kNx=pi/(Lx/Nx); kNy=pi/(Ly/Ny);
np2m_x=[]; np2m_y=[]; twonpm_x=[]; twonpm_y=[]; nmp_x=[]; nmp_y=[];
for n=1:N
  for m=n+1:N
    np2m_x(end+1,1)=kx(n)+2*kx(m); np2m_y(end+1,1)=ky(n)+2*ky(m);
    twonpm_x(end+1,1)=2*kx(n)+kx(m); twonpm_y(end+1,1)=2*ky(n)+ky(m);
    for p=m+1:N
      nmp_x(end+1,1)=kx(n)+kx(m)+kx(p); nmp_y(end+1,1)=ky(n)+ky(m)+ky(p);
    end
  end
end
over_np2m=mean(abs(np2m_x)>kNx | abs(np2m_y)>kNy);
over_2npm=mean(abs(twonpm_x)>kNx | abs(twonpm_y)>kNy);
over_nmp=mean(abs(nmp_x)>kNx | abs(nmp_y)>kNy);
fprintf('nyquist overflow ratio: np2m %.2f%%, 2npm %.2f%%, nmp %.2f%%\n',100*over_np2m,100*over_2npm,100*over_nmp);

% Compare direct vs spectral phi3 term, before and after snapping k to grid
x=(0:Nx-1)*(Lx/Nx); y=(0:Ny-1)*(Ly/Ny); [X,Y]=meshgrid(x,y);

c2=coeffsMF12(2,g,h,a,b,kx,ky,Ux,Uy); c3=coeffsMF12(3,g,h,a,b,kx,ky,Ux,Uy);
[~,p2d]=surfaceMF12_new(2,c2,X,Y,t); [~,p3d]=surfaceMF12_new(3,c3,X,Y,t);
[~,p2s]=surfaceMF12_spectral(c2,Lx,Ly,Nx,Ny,t); [~,p3s]=surfaceMF12_spectral(c3,Lx,Ly,Nx,Ny,t);
phi3d=p3d-p2d; phi3s=p3s-p2s;
d0=max(abs(phi3d(:)-phi3s(:))); r0=d0/max(abs(phi3d(:)));
fprintf('phi3 term diff (orig k): max %.3e, rel %.3e\n',d0,r0);

% Snap k to grid
kx_sn=round(kx/dkx)*dkx; ky_sn=round(ky/dky)*dky;
c2n=coeffsMF12(2,g,h,a,b,kx_sn,ky_sn,Ux,Uy); c3n=coeffsMF12(3,g,h,a,b,kx_sn,ky_sn,Ux,Uy);
[~,p2d_n]=surfaceMF12_new(2,c2n,X,Y,t); [~,p3d_n]=surfaceMF12_new(3,c3n,X,Y,t);
[~,p2s_n]=surfaceMF12_spectral(c2n,Lx,Ly,Nx,Ny,t); [~,p3s_n]=surfaceMF12_spectral(c3n,Lx,Ly,Nx,Ny,t);
phi3d_n=p3d_n-p2d_n; phi3s_n=p3s_n-p2s_n;
d1=max(abs(phi3d_n(:)-phi3s_n(:))); r1=d1/max(abs(phi3d_n(:)));
fprintf('phi3 term diff (snapped k): max %.3e, rel %.3e\n',d1,r1);

% Also report centerline rel errors
[~,iyc]=min(abs(y-Ly/2));
cl0=max(abs(phi3d(iyc,:)-phi3s(iyc,:)))/max(abs(phi3d(iyc,:)));
cl1=max(abs(phi3d_n(iyc,:)-phi3s_n(iyc,:)))/max(abs(phi3d_n(iyc,:)));
fprintf('centerline rel diff: orig %.3e, snapped %.3e\n',cl0,cl1);
