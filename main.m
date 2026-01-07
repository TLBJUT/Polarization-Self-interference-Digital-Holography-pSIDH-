% 
% ==============================================================================
% # This script provides the data generation and testing code for
% computational adaptive optics (CAO) aberration correction in self interference digital holography (SIDH).
% #
% # Author:   Tianlong Man   Beijing University of Technology
% # Date:     2026/1/6
% 
%   Important notification: You have to use your own optimizer (for
%   exmample, SPGD, Adam) in this version of code. We are transfering our
%   previous SPGD optimizer to Adam now, and this code will be updated
%   soon!
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
um=1e-6;cm=1e-2;mm=1e-3;
lamda=0.550*um;f1=100*mm;d1=142.5*mm;z1=-0*mm;z2=-1.6*mm;ds=95*mm;f_SLM=2400*mm; 
M=512;N=512;A=1;
hx=M*15*um;hy=N*15*um;
dhx=hx/M;dhy=hy/N;  
% Below two foci after the beam dividing elements
foci1=(d1*f1)/(d1-f1)-ds;
foci2=(foci1*f_SLM)/(f_SLM+foci1);
% Optimizaed recording position of the hologram
zh=2*foci1*foci2/(foci1+foci2);
% Infocused reconstruction distance, pre-calculation
zr=(foci1-zh)*(foci2-zh)/(foci1-foci2);
% =========================Object generation============================
aa=imread('USAF.bmp'); % Use here your own sample image
aa=aa(:,:,1);
aa=double(aa);
a=imresize(aa,[512,512],'bicubic');

ab=imread('11.bmp');
ab=ab(:,:,1);
ab=0;

 b=zeros(M);
  b(256-128:256+127,256-128:256+127)=double(ab);

k=2*pi/lamda;
for ija=1:M
    for ijb=1:N
        if (ija-M/2)^2+(ijb-N/2)^2<=(0/2)^2
            ap(ija,ijb)=1;
        else
            ap(ija,ijb)=0;
        end
    end
end

% Define parameters of CAO %


PARAMS(1) = 512; 	% size of pupil aperture field in pixels (this defines the resolution of the calculation)
PARAMS(8) = 512;%PARAMS 2,3 and 4 are set later in the program because they depend on how the aberrations are read in
PARAMS(5) = 0.550;	% imaging wavelength in microns
PARAMS(6) = 0;		% number of pixels over which PSF is calculated (do not adjust here set in Zernikephase subroutine)
PARAMS(7) = 20; % increase to enhance the display of the wavefront (doesn't affect calculation)
PARAMS(2) = 25.4;   % size of pupil in mm for which PSF and MTF is to be calculated
PARAMS(3) = 25.4;	% size of pupil in mm that Zernike coefficients define. NOTE: You can define the aberrations
					% for any pupil size and calculate their effects for any smaller aperture.
    % Because of the reciprocal relationship between the pupil function and its Fourier transform,
    % it helps to define a pupil that makes up only a small central region of the pupil aperture.
    % Otherwise the point spread function is too small.
PARAMS(4) = 25.4;	param4orig = PARAMS(4); 	% size of pupil field in mm (use a large field to magnify the PSF)

  c(1)=0; 	c(2)=0;%Tilts (Lateral position)
   %defocus and astigmatism
   c(3)=2.0; 	c(4)=0.0;	c(5)=-0.0;
   %coma like
   c(6)=0.0;	c(7)=0.0;	c(8)=0.0;	c(9)=-4;	
   %spherical aberration like
   c(10)=-0.0;	c(11)=0.0;	c(12)=2;	c(13)=0;	c(14)=0.0;	
   %higher order (each row is a new radial order)
   c(15)=0;	c(16)=0;	c(17)=0;	c(18)=0;	c(19)=0;	c(20)=0;	
   c(21)=0;	c(22)=0;	c(23)=0;	c(24)=0;	c(25)=0;	c(26)=0;	c(27)=0;	
   c(28)=0;	c(29)=0;	c(30)=0;	c(31)=0;	c(32)=0;	c(33)=0;	c(34)=0;	c(35)=0;	
   c(36)=0;	c(37)=0;	c(38)=0;	c(39)=0;	c(40)=0;	c(41)=0;	c(42)=0;	c(43)=0;	c(44)=0;	
   c(45)=0;	c(46)=0;	c(47)=0;	c(48)=0;	c(49)=0;	c(50)=0;	c(51)=0;	c(52)=0;	c(53)=0;	c(54)=0;	
   c(55)=0;	c(56)=0;	c(57)=0;	c(58)=0;	c(59)=0;	c(60)=0;	c(61)=0;	c(62)=0;	c(63)=0;	c(64)=0;	c(65)=0;

d=0.1;

[p_ab, pd]= zerfun_PD(PARAMS,c,d);
figure;imshow(angle(p_ab),[]),colormap('jet');
% figure;imshow(fftshift(fft2(p_ab)),[]);
% break;
% =========================Generating the holographic PSF and then holograms========================
[Ip11,Ip12,Ip13,Ip14]=FINCH4_PSF(ap,M,N,dhx,dhy,f1,d1+z1,ds,f_SLM,zh,lamda,p_ab);
maskx=zeros(512,512);
maskx(256-50:256+50,256-50:256+50)=1;
% figure;imshow(Ip11,[]);
% break;
% Aberrated holographic PSF generated %
Ip11=Ip11.*maskx;
Ip12=Ip12.*maskx;
Ip13=Ip13.*maskx;
Ip14=Ip14.*maskx;
% break;
% Aberrated phase-shifted hologram of the extended object generated %
Ip11p=conv2(a,Ip11,'same');
Ip12p=conv2(a,Ip12,'same');
Ip13p=conv2(a,Ip13,'same');
Ip14p=conv2(a,Ip14,'same');
% Aberrated complex-valued hologram of the extended object generated %
f1p=(Ip11p-Ip13p)+1i*(Ip14p-Ip12p);


% Below required for generating the hologram of 3D object (two 2D layers) % 
[Ip21,Ip22,Ip23,Ip24]=FINCH4_PSF(ap,M,N,dhx,dhy,f1,d1+z2,ds,f_SLM,zh,lamda,p_ab);
% break;
Ip21p=conv2(b,Ip21.*maskx,'same');
Ip22p=conv2(b,Ip22.*maskx,'same');
Ip23p=conv2(b,Ip23.*maskx,'same');
Ip24p=conv2(b,Ip24.*maskx,'same');

% f2p=(Ip21p-Ip23p)+1i*(Ip24p-Ip22p);
% Ip21=conv2(b,f2p,'same');
Ip1=Ip11p+Ip21p;
Ip2=Ip12p+Ip22p;
Ip3=Ip13p+Ip23p;
Ip4=Ip14p+Ip24p;
S1=sum(sum(Ip1(:,:)));
S2=sum(sum(Ip2(:,:)));
S3=sum(sum(Ip3(:,:)));
S4=sum(sum(Ip4(:,:)));
Ip1=Ip1+10e4*randn(size(Ip1));
Ip2=Ip2+10e4*randn(size(Ip2));
Ip3=Ip3+10e4*randn(size(Ip3));
Ip4=Ip4+10e4*randn(size(Ip4));
Ip=(Ip1-Ip3)-1j*(Ip2-Ip4);
% Ip=Ip11+Ip21;
figure;imshow(abs(Ip),[]);


% ===========================CAO correction and back-propagation reconstruction===========================
% =====================================Use your own optimizer!============================================

[MM,NN]=size(Ip);
M=MM;N=NN;
dhx=hx/M;dhy=hy/N;
% dhx=4.5*um;dhy=4.5*um;


  c(1)=0; 	c(2)=0;%Tilts (Lateral position)
   %defocus and astigmatism
   c(3)=1.8; 	c(4)=3;	c(5)=-0.0;
   %coma like
   c(6)=0.0;	c(7)=0.0;	c(8)=0;	c(9)=-1.6;	
   %spherical aberration like
   c(10)=-0;	c(11)=0.0;	c(12)=1.5;	c(13)=0;	c(14)=0.0;	
   %higher order (each row is a new radial order)
   c(15)=0;	c(16)=0;	c(17)=0;	c(18)=0;	c(19)=0;	c(20)=0;	
   c(21)=0;	c(22)=0;	c(23)=0;	c(24)=0;	c(25)=0;	c(26)=0;	c(27)=0;	
   c(28)=0;	c(29)=0;	c(30)=0;	c(31)=0;	c(32)=0;	c(33)=0;	c(34)=0;	c(35)=0;	
   c(36)=0;	c(37)=0;	c(38)=0;	c(39)=0;	c(40)=0;	c(41)=0;	c(42)=0;	c(43)=0;	c(44)=0;	
   c(45)=0;	c(46)=0;	c(47)=0;	c(48)=0;	c(49)=0;	c(50)=0;	c(51)=0;	c(52)=0;	c(53)=0;	c(54)=0;	
   c(55)=0;	c(56)=0;	c(57)=0;	c(58)=0;	c(59)=0;	c(60)=0;	c(61)=0;	c(62)=0;	c(63)=0;	c(64)=0;	c(65)=0;

d=0.1;
%Pupil function for CAO correction generated
[p_abx, pd]= zerfun_PD(PARAMS,c,d);
%CAO correction in spatial frequency domain of the reconstruction
F=fftshift(fft2(Ip));
F=F.*p_abx;
figure;imshow(angle(p_abx),[]);
fx=ifft2(F);
figure;imshow(abs(fx),[]);
% imwrite(fx,'re3.5.tif');
% imwrite(f3,'ao0.tif');
% Fk = (fft2(fx));
% Fr = Fk.*conj(p_abx);
% Fr = ifft2(Fr);
% Fr=abs(Fr);
% Fr=Fr/max(max(Fr));
% figure;imshow(Fr,[]);