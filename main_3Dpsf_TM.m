% 
% ==============================================================================
% # This script provides the data generation and reconstruction code for self interference digital holography (SIDH)
% imaging of a 3D distributed particles sample.
% #
% # Author:   Tianlong Man   Beijing University of Technology
% # Date:     2026/1/6
% 
%
% #Referene
% 
% # [1] T. Man, Y, Wan, W. Yan, et al., "Adaptive optics via self-interference digital holography
% for non-scanning three-dimensional imaging in biological samples," Biomedical Optics Express, 9(6):2614, 2018.
% # https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix/blob/master/models/networks.py
% #
% # [2] T. Man, W. Zhang, L. Zhang, et al., "Scanning-free three-dimensional fluorescent dipoles imaging by polarization 
% self-interference digital holography (pSIDH)," DOI
% 10.48550/arXiv.2504.10772, 2025
% # https://github.com/csleemooo/Deep_learning_based_on_parameterized_physical_forward_model_for_adaptive_holographic_imaging
% 
% #==============================================================================

% Define parameters %
clc;clear;close all;
um = 1e-6;cm = 1e-2;mm = 1e-3;nm = 1e-9;
lambda = 550*nm;
f_o = 3*mm;    %focal length of objective lens
Pupil_ro = 300;  %back aperture size of the objective in pixel
Pupil_rs = 300;  %aperture size of the SLM in pixel
zs = 2.991*mm;
ds = 0*mm;    %objective to SLM
f_SLM = 250*mm;
M = 1024;N = 1024;
nz = 512;
dhx = 10*um;dhy = 10*um;
dhz = 0.6*um;
k = 2*pi/lambda;

if zs == f_o 
    foci2 = f_SLM;
    z_h = 2*foci2;
    zr = -foci2;
    trans_mag = f_SLM/f_o;
else
    foci1 = ((zs*f_o)/(f_o-zs));
    foci2 = (foci1+ds)*f_SLM/(f_SLM-foci1-ds);
%     z_h = 2*foci1*foci2/(foci1+foci2);
    z_h = 0.5;
    zr = -(foci1+ds+z_h)*(foci2+z_h)/(ds+foci1-foci2);
    trans_mag = (z_h*foci1)/(zs*(foci1+ds));
end


Pupil_o = zeros(N,M); 
for i=1:N
    for j=1:M
        if (i-M/2-1)^2+(j-N/2-1)^2<=(Pupil_ro)^2 
                                    
            Pupil_o(i,j)=1;

        end
    end
end

Pupil_SLM=zeros(N,M);
for i=1:N
    for j=1:M
        if (i-M/2-1)^2+(j-N/2-1)^2<=(Pupil_rs)^2 
            Pupil_SLM(i,j)=1;
        end
    end
end

%========================3D object building======================

o_3d = zeros(M,N,nz);

%random distributed fluorphores

Np = 200;
% xps = randi([390 580],1,Np);
% yps = randi([390 580],1,Np);
xps = randperm(600,Np)+200;
yps = randperm(600,Np)+200;
zps = randperm(200,Np)+200;
% zps = 255;
% o_3d(540,540,245) = 1;
for kk = 1:Np
    o_3d(xps(kk),yps(kk),zps(kk)) = 1;
end

%assign axial coordinates in OBJECT space

for kk = 1:nz
    z(kk) = (kk-(nz/2-1))*dhz+f_o;
end

%================================Done===================================

%=========================Holograms generation==========================
% 
Ip1 = zeros(M,N);
Ip2 = zeros(M,N);
Ip3 = zeros(M,N);
Ip4 = zeros(M,N);



Pupil_ro = 140;
Pupil_o = zeros(N,M); 
for i=1:N
    for j=1:M
        if (i-M/2-1)^2+(j-N/2-1)^2<=(Pupil_ro)^2 
                                    
            Pupil_o(i,j)=1;

        end
    end
end





for kk = 1:nz
    [Ip11,Ip12,Ip13,Ip14]=FINCH_PSF(o_3d(:,:,kk),M,N,dhx,dhy,f_o,z(kk),ds,f_SLM,z_h,lambda,Pupil_SLM,Pupil_o);
    Ip1 = Ip1 + Ip11/10;
    Ip2 = Ip2 + Ip12/10;
    Ip3 = Ip3 + Ip13/10;
    Ip4 = Ip4 + Ip14/10;    
    fprintf(1,'SIDH hologram calculating %4.2f percent\n',(kk/nz*100))
end

% Ip1 = imnoise(Ip1,'poisson');
% Ip2 = imnoise(Ip2,'poisson');
% Ip3 = imnoise(Ip3,'poisson');
% Ip4 = imnoise(Ip4,'poisson');

Ip1 = Ip1/max(max(Ip1));
Ip2 = Ip2/max(max(Ip2));
Ip3 = Ip3/max(max(Ip3));
Ip4 = Ip4/max(max(Ip4));

PSH = (Ip1-Ip3)-1j*(Ip2-Ip4);    %point spread holgoram of the 3D distributed fluorophores


% PSH = imnoise(PSH,'gaussian',0.01,0.005);

% figure,imshow(angle(PSH),[]);
% figure,imshow(Ip1,[]);
imwrite(Ip1,'Hologram1zx.tif');
% imwrite(Ip4,'Hologram4.tif');

% imwrite(angle(PSH)/max(max(angle(PSH))),'Hologram_complexedphase.tif');
%================================Done===================================




%=========================Holograms (Low/High NA) generation==========================

Ipx1 = zeros(M,N);
Ipx2 = zeros(M,N);
Ipx3 = zeros(M,N);
Ipx4 = zeros(M,N);



Pupil_ro = 500;
Pupil_o = zeros(N,M); 
for i=1:N
    for j=1:M
        if (i-M/2-1)^2+(j-N/2-1)^2<=(Pupil_ro)^2 
                                    
            Pupil_o(i,j)=1;

        end
    end
end





for kk = 1:nz
    [Ip11,Ip12,Ip13,Ip14]=FINCH_PSF(o_3d(:,:,kk),M,N,dhx,dhy,f_o,z(kk),ds,f_SLM,z_h,lambda,Pupil_SLM,Pupil_o);
    Ipx1 = Ipx1 + Ip11/10;
    Ipx2 = Ipx2 + Ip12/10;
    Ipx3 = Ipx3 + Ip13/10;
    Ipx4 = Ipx4 + Ip14/10;    
    fprintf(1,'SIDH High NA hologram calculating %4.2f percent\n',(kk/nz*100))
end

% Ip1 = imnoise(Ip1,'poisson');
% Ip2 = imnoise(Ip2,'poisson');
% Ip3 = imnoise(Ip3,'poisson');
% Ip4 = imnoise(Ip4,'poisson');

Ipx1 = Ipx1/max(max(Ipx1));
Ipx2 = Ipx2/max(max(Ipx2));
Ipx3 = Ipx3/max(max(Ipx3));
Ipx4 = Ipx4/max(max(Ipx4));

PSHx = (Ipx1-Ipx3)-1j*(Ipx2-Ipx4);    %point spread holgoram of the 3D distributed fluorophores




%================================Done===================================



% 
% %=============================3D WF images==============================
% Pupil_ro = 800;
% Pupil_o = zeros(N,M); 
% for i=1:N
%     for j=1:M
%         if (i-M/2-1)^2+(j-N/2-1)^2<=(Pupil_ro)^2 
%                                     
%             Pupil_o(i,j)=1;
% 
%         end
%     end
% end
% 
% 
% % Zr_WF = 382*mm:2*mm:682*mm;
% % Zr_WF = 173*27/18*mm:2*27/18*mm:473*27/18*mm;
% Zr_WF = 173*mm:2*mm:473*mm;
% WF = zeros(M,N,size(Zr_WF,2));
% % f_tube = f_o*trans_mag;
% % f_tube = 250*27/18*mm;
% f_tube = 250*mm;
% 
% for jj = 1:nz
%     [Ip]=WF_PSF(o_3d(:,:,jj),M,N,dhx,dhy,f_o,z(jj),ds,f_tube,Zr_WF,lambda,Pupil_SLM,Pupil_o);
%     WF = WF + Ip;  
% %     jj
%     fprintf(1,'WF image calculating %4.2f percent\n',(jj/nz*100))
% end
% WF = WF/max(WF(:));    %global normalization
% WF = WF(M/2-M/4:M/2+M/4-1,M/2-M/4:M/2+M/4-1,:);
% % figure,imshow(WF(:,:,50),[]);
% 
% 
% %================================Done===================================
% 
% 
% %============================Image export===============================
% for kk = 1:size(Zr_WF,2)
%     zr = Zr_WF(kk);
%     filename = ['WF 3D reconstruction\',num2str(zr*1000),'.tif'];
%     imwrite(WF(:,:,kk),filename);
% end
% %================================Done===================================
% 



%=================Fresnel back propagation reconstruction===============

index_z = 0;
Zr = -37*mm:-2*mm:-337*mm;
re_3d_int = zeros(M,N,size(Zr,2));
re_3d_ang = zeros(M,N,size(Zr,2));
for kk = 1:size(Zr,2)
    zr = Zr(kk);
    [re_c,dx2,dy2,x2,y2] = fresnel(PSH,M,N,dhx,dhy,zr,lambda);
    re_int = abs(re_c).^2;
    re_ang = angle(re_c);
    re_ang = re_ang-min(min(re_ang));
    re_ang = re_ang/max(max(re_ang));
    re_3d_int(:,:,kk) = re_int;
    re_3d_ang(:,:,kk) = re_ang;
end

re_3d_int = re_3d_int/max(re_3d_int(:));    %global normalization
re_3d_int = re_3d_int(M/2-M/4:M/2+M/4-1,M/2-M/4:M/2+M/4-1,:);   %crop the image
re_3d_ang = re_3d_ang(M/2-M/4:M/2+M/4-1,M/2-M/4:M/2+M/4-1,:); 
%================================Done===================================


%=================Fresnel back propagation reconstruction===============

index_z = 0;
Zr = -37*mm:-2*mm:-337*mm;
rex_3d_int = zeros(M,N,size(Zr,2));
rex_3d_ang = zeros(M,N,size(Zr,2));
for kk = 1:size(Zr,2)
    zr = Zr(kk);
    [re_c,dx2,dy2,x2,y2] = fresnel(PSHx,M,N,dhx,dhy,zr,lambda);
    re_int = abs(re_c).^2;
    re_ang = angle(re_c);
    re_ang = re_ang-min(min(re_ang));
    re_ang = re_ang/max(max(re_ang));
    rex_3d_int(:,:,kk) = re_int;
    rex_3d_ang(:,:,kk) = re_ang;
end

rex_3d_int = rex_3d_int/max(rex_3d_int(:));    %global normalization
rex_3d_int = rex_3d_int(M/2-M/4:M/2+M/4-1,M/2-M/4:M/2+M/4-1,:);   %crop the image
rex_3d_ang = rex_3d_ang(M/2-M/4:M/2+M/4-1,M/2-M/4:M/2+M/4-1,:); 
%================================Done===================================



%============================Image export===============================
for kk = 1:size(Zr,2)
    zr = Zr(kk);
    filename = ['3D reconstruction intensity\',num2str(zr*1000),'.tif'];
    imwrite(re_3d_int(:,:,kk),filename);
    filename = ['3D reconstruction phase\',num2str(zr*1000),'.tif'];
    imwrite(re_3d_ang(:,:,kk),filename);
end
%================================Done===================================

%============================Image export===============================
for kk = 1:size(Zr,2)
    zr = Zr(kk);
    filename = ['3D reconstruction intensity High NA\',num2str(zr*1000),'.tif'];
    imwrite(rex_3d_int(:,:,kk),filename);
    filename = ['3D reconstruction phase High NA\',num2str(zr*1000),'.tif'];
    imwrite(rex_3d_ang(:,:,kk),filename);
end
%================================Done===================================