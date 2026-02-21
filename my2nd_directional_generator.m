%This code is to use fuhrman's function to generate the directional
%subharmonic and compare with the OW3D results
clc;clear;close all;
Akp=0.02;
tic
Tp=12;
kp=0.0279;
Alpha=8;
for kd=[1,2,5]
    
    h=kd/kp;
    % h=150;
    g=9.81;
    omega_p=sqrt(g*kp);
    lambda=225;
    kw=0.004606;
    modk=0:kp/40:10*kp;
    omega=sqrt(g.*modk.*tanh(modk.*h));
    t_init=-15*12;
    dk=modk(2)-modk(1);
    % S=exp(-(modk-kp).^2./(2.*kw.^2));
    for i=1:numel(modk)
        if modk(i)<kp
            S(i)=exp(-(modk(i)-kp).^2/(2*kw.^2));
        else
            %        S(i)=exp(-(k(i)-kp).^2/(2*kw.^2))^(1/alpha);
            kw2=sqrt(kp^2/(2*log(10^Alpha)));
            S(i)=exp(-(modk(i)-kp).^2/(2*kw2.^2));
        end
    end
    theta_p=0;
    %     theta=-pi/2:pi/72:pi/2;
    theta=deg2rad(linspace(-90,90,120));
    % for nz=[9,17,33]
    for spreading_angle=[5,15,25]
        theta_w=deg2rad(spreading_angle);
        D_theta=exp(-(theta-deg2rad(theta_p)).^2./(2*theta_w.^2));
        % D_theta=D_theta+exp(-(theta+deg2rad(theta_p)).^2./(2*theta_w.^2));
        dtheta=theta(2)-theta(1);
        % x=-2000:5:2000;
        % y=-2000:5:2000;
        %         plot(rad2deg(theta),D_theta);
        %%%%%%%%%%%Time discretization
        T_total=30*Tp;
        Nt=30;
        dt=Tp/Nt;
        N_steps=T_total/dt;
        sumsum=sum(D_theta)*sum(S);
        %%%%%Start linear part
        % linear_focus_envelope=Linear_focus_envelope(Akp,9,x);
        for Akp=[0.02]
            clf;
            for ORDER=[2]
                for phi_shift=[0:90:270]
                    %         phi_shift=90;
                    x=linspace(-20*lambda,20*lambda,1025);
                    y=linspace(-15*lambda,15*lambda,257);
                    dx=x(2)-x(1);
                    dy=y(2)-y(1);
                    % x=linspace(0*lambda,40*lambda,512);
                    % y=linspace(-10*lambda,10*lambda,257);
                    [Y,X]=meshgrid(y,x);
                    x_vec=x;
                    y_vec=y;
                    eta=0.*X;
                    phi=0.*X;
                    % phi_shift=deg2rad(90);
                    %                     mask=0.*X;
                    %                     mask(512-300:513+300,256-200:256+200)=1;
                    %     %         mask(128-100:128+100,1025-500:1025+500)=1;
                    tic
                    for i=1:numel(modk)
                        for j=1:numel(D_theta)
                            eta=eta+Akp/kp/sumsum*S(i)*D_theta(j).*cos( modk(i)*cos(theta(j)).*X + modk(i)*sin(theta(j)).*Y+omega(i)*t_init+deg2rad(phi_shift));
                            if omega(i)~=0
                                phi=phi-g./omega(i).*Akp/kp/sumsum*S(i)*D_theta(j).*sin( modk(i)*cos(theta(j)).*X + modk(i)*sin(theta(j)).*Y+omega(i)*t_init+deg2rad(phi_shift));
                            end
                        end
                    end
                    toc
                    %     mesh(X,Y,eta);
                    %                         eta=eta.*mask;
                    %                     test_alpha_central=Linear_focus_envelope(Akp,Alpha,x_vec);
                    %                     hold on
                    %                     plot(test_alpha_central);
                    %                     plot(eta(:,128));
                    Lx    = max(X(:)-min(X(:)));          % 物理域长度 (便于画图，仅用于可视化)
                    Ly    = max(Y(:)-min(Y(:)));
                    [Nx,Ny]=size(eta);
                    dkx   = 2*pi/Lx;
                    dky   = 2*pi/Ly;
                    kx_vec = ( -floor(Nx/2) : +ceil(Nx/2)-1 ) * dkx;   % 长度恰好为 Nx
                    ky_vec = ( -floor(Ny/2) : +ceil(Ny/2)-1 ) * dky;   % 长度恰好为 Ny
                    % kx_vec = ( -floor(Nx/2)+1 : +ceil(Nx/2) ) * dkx;   % 长度恰好为 Nx
                    % ky_vec = ( -floor(Ny/2)+1 : +ceil(Ny/2) ) * dky;   % 长度恰好为 Ny
                    [kyGrid, kxGrid] = meshgrid(ky_vec, kx_vec);

                    Kabs=sqrt(kxGrid.^2+kyGrid.^2);
                    sigma=tanh(Kabs.*h);
                    Omega=sqrt(g.*Kabs.*tanh(Kabs.*h));
                    %%%%%%%%%%Getting velocity at surface with the envelope method
                    S_keta=fft2(eta);
                    %                     mesh(kxGrid,kyGrid,abs(fftshift(S_keta)));
                    %                     keta=real(ifft2(1i.*S_keta.*(ifftshift(kxGrid)))+ifft2(1i.*S_keta.*(ifftshift(kyGrid))));
                    keta=real(ifft2(1i.*S_keta.*(ifftshift(kxGrid))));
                    alpha1=cosh(2.*Kabs.*h);
                    sigma_eta22=(3-sigma.^2)./(4*sigma.^3);
                    sigma_eta22(~isfinite(sigma_eta22))=0;
                    sigma_phi22=sigma_eta22./Omega.*9.81;
                    sigma_phi22(~isfinite(sigma_phi22))=0;
                    sigma_eta33=(27-9.*sigma.^2+9.*sigma.^4-3.*sigma.^6)./(64.*sigma.^6);
                    sigma_eta33(~isfinite(sigma_eta33))=0;
                    sigma_eta33(sigma_eta33>50)=0;
                    sigma_phi33=sigma_eta33./Omega.*9.81;
                    sigma_phi33(~isfinite(sigma_phi33))=0;
                    sigma_eta44 = (24*alpha1.^6 + 116*alpha1.^5 + 214*alpha1.^4 + 188*alpha1.^3 + 133*alpha1.^2 + 101*alpha1 + 34) ./ ...
                        (24*(3*alpha1 + 2).*(alpha1.^4) .* sinh(2.*Kabs.*h));
                    sigma_eta44(~isfinite(sigma_eta44))=0;
                    sigma_eta33(sigma_eta44>50)=0;
                    sigma_phi44=sigma_eta44./Omega.*9.81;
                    sigma_phi44(~isfinite(sigma_phi44))=0;
                    sigma_eta55 = (5*(300*alpha1.^8 + 1579*alpha1.^7 + 3176*alpha1.^6 + 2949*alpha1.^5 + 1188*alpha1.^4 + ...
                        675*alpha1.^3 + 1326*alpha1.^2 + 827*alpha1 + 130)) ./ ...
                        (384*(alpha1 - 1).^6.*(12*alpha1.^2 + 11*alpha1 + 2));
                    sigma_eta55(~isfinite(sigma_eta55))=0;
                    sigma_eta55(sigma_eta55>50)=0;
                    sigma_phi55=sigma_eta55./Omega.*9.81;
                    sigma_phi55(~isfinite(sigma_phi55))=0;

                    test_keta=real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_eta22)))));
                    test_keta_phi=real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_phi22)))));
                    eta_analytic=analytic2D(eta);
                    directional_phase=angle(eta_analytic);
                    eta_22=real(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_eta22)))))).*exp(2i.*directional_phase));
                    phi_22=imag(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_phi22)))))).*exp(2i.*directional_phase));
                    eta_33=real(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_eta33)))))).*envelope(keta).*exp(3i.*directional_phase));
                    phi_33=imag(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_phi33)))))).*envelope(keta).*exp(3i.*directional_phase));
                    eta_44=real(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_eta44)))))).*envelope(keta).^2.*exp(4i.*directional_phase));
                    phi_44=imag(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_phi44)))))).*envelope(keta).^2.*exp(4i.*directional_phase));
                    eta_55=real(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_eta55)))))).*envelope(keta).^3.*exp(5i.*directional_phase));
                    phi_55=imag(envelope(eta).*envelope(real(ifft2(1i.*S_keta.*(ifftshift(kxGrid.*abs(sigma_phi55)))))).*envelope(keta).^3.*exp(5i.*directional_phase));
                    %                     hold on
                    %                     plot(test_env_eta22(:,128));
                    %                     %                                     plot(env_eta22(:,128));
                    %                     plot(eta_p(:,128));
                    %                 hold on
                    %                 plot(test_env_phi22(:,128));
                    %                 plot(env_phi22(:,128));
                    %%%%%%%%%%%%Compare with unidirectional dalzell
                    %                 eta_analytic = analytic2D(eta);
                    %                                 [phi_m, eta_m, phi_p, eta_p] = compute_second_order(eta, dx, dy, h);
                    %                                 %                 phi_m=phi_m*4;
                    %                                 %                 eta_m=eta_m*4;
                    %                                 %                 phi_p=phi_p*4;
                    %                                 %                 eta_p=eta_p*4;
                    %                                 %                 Check with uni-directional results
                    %                                 eta_11_uni=eta(:,128);
                    %                                 [Dalzell_eta22,Dalzell_eta20,Dalzell_phi22,Dalzell_phi20]=dalzell_2d(eta_11_uni,x_vec,h);
                    %                                 eta_p_uni=eta_p(:,128);
                    %                                 eta_m_uni=eta_m(:,128);
                    %                                 etaAm_uni=phi_m(:,128);
                    %                                 etaAp_uni=phi_p(:,128);
                    if ORDER==2
                        tic
                        [phi_m, eta_m, phi_p, eta_p] = compute_second_order(eta, dx, dy, h);
                        %                     phi_m=phi_m*4;
                        %                     eta_m=eta_m*4;
                        %                     phi_p=phi_p*4;
                        %                     eta_p=eta_p*4;
                        
                        toc
                        eta_2nd=eta+eta_p+eta_m;
                        phi_2nd=phi+phi_p+phi_m;
                        %%%%%%%%%%%%%%%%%Test for envelope method on
                        %%%%%%%%%%%%%%%%%surface in finite depth up to 5th
                        %%%%%%%%%%%%%%%%%order
                        eta_2nd=eta_2nd+eta_33+eta_44+eta_55;
                        phi_2nd=phi_2nd+phi_22+phi_33+phi_44+phi_55;

                    end
                    if ORDER==1
                        eta_2nd=eta;
                        phi_2nd=phi;
                    end
                    semilogy(abs(fft(eta_2nd(:,128))));
                    mesh(eta_2nd);
                    drawnow
                    %                 mesh(abs(fftshift(fft2(eta))))
                    %                 set(gca, 'ZScale', 'log');
                    %                     write_path=fullfile(pwd,"test1",sprintf('Order_%d_spread_%d_Akp_%.2f_phi_shift_%d',ORDER,spreading_angle,Akp,phi_shift));
                    write_path=fullfile(pwd,"test1",sprintf('kd%.1f_spread_%d_Akp_%.2f_phi_shift_%d',kd,spreading_angle,Akp,phi_shift));
                    if exist(write_path, 'dir')
                        disp('目录已存在。');
                    else
                        disp('目录不存在，正在创建目录。');
                        % 创建目录
                        mkdir(write_path);
                        disp(['目录已创建：', write_path]);
                    end
                    file_name=fullfile(write_path,'OceanWave3D.init');
                    nx=size(eta_2nd,1);
                    ny=size(eta_2nd,2);
                    % X=x_vec;
                    %                 dx=X(1,2)-X(1,1);
                    % Y=y_vec;
                    %                 dy=Y(2,1)-Y(1,1);
                    Lx=max(X(:))-min(X(:));

                    Ly=max(Y(:))-min(Y(:));

                    fileid=fopen(file_name,'w');

                    fprintf(fileid,' H=%f nx=%d ny=%d dx=%f dy=%f akp=%f shift=%f' ,max(eta_2nd(:)),nx,ny,dx,dy,Akp,phi_shift);

                    fprintf(fileid,'\n%12e %12e %d %d %12e',Lx,Ly,nx,ny,dt);

                    for ry=1:ny
                        for rx=1:nx
                            fprintf(fileid,'\n%12e %12e',eta_2nd(rx,ry),phi_2nd(rx,ry));
                        end
                    end
                    %                 for rx=1:nx
                    %                     for ry=1:ny
                    %                         fprintf(fileid,'\n%12e %12e',eta_2nd(rx,ry),phi_2nd(rx,ry));
                    %                     end
                    %                 end

                    fclose(fileid);

                    %Write inp
                    file_name=fullfile(write_path,'OceanWave3D.inp');
                    fileid=fopen(file_name,'w');

                    fprintf(fileid,'Data for paper directional in %s <-\n',datestr(now,0));

                    fprintf(fileid,'-1 2 <-\n');
                    fprintf(fileid,'%d %d %d %d %d %d 0 0 1 1 1 1 <-\n',Lx,Ly,h,nx,ny,9);
                    %                 fprintf(fileid,'3 0 3 1 1 1 <-\n');
                    fprintf(fileid,'4 4 4 1 1 1 <-\n');
                    fprintf(fileid,'%d %f 1 0. 1 <-\n',round(N_steps)+1,dt);
                    fprintf(fileid,'9.81 <-\n');
                    fprintf(fileid,'1 3 0 55 1e-6 1e-6 1 V 1 1 20 <-\n');
                    fprintf(fileid,'0.05 1.00 1.84 2 0 0 1 6 32 <-\n');
                    fprintf(fileid,'10 1 <-\n');

                    fprintf(fileid,'1 0 <-\n');
                    fprintf(fileid,'0 6 10 0.08 0.08 0.4 <-\n');
                    fprintf(fileid,'0 8. 3 X 0.0 <-\n');
                    fprintf(fileid,'0 0 <-\n');
                    fprintf(fileid,'0 2.0 2 0 0 1 0 <-\n');
                    fprintf(fileid,'0 <-\n');
                    fprintf(fileid,'33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3 <-\n');
                    fclose(fileid);


                    file_name=fullfile(write_path,'OW_readme.txt');
                    fileid=fopen(file_name,'w');
                    fprintf(fileid,'Test 3D T=%1.2fTp H=%f dx=%f, dt=%f\n',0./Tp,max(eta_2nd(:)),dx,dt);
                    fprintf(fileid,'Akp=%f, kp=%f, lambda=%f, Tp=%f\n',Akp,kp,lambda,Tp);
                    % fprintf(fileid,'lambda=%f dx, T=%f dt, cw=%f, CFL=%f, kd=%f\n',lambda/dx,Tp/dt,cw,CFL,kp*h);
                    fprintf(fileid,'Last modified in %s',datestr(now,0));
                    fclose(fileid);
                end
            end
        end
    end
end
% end
