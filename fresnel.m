% 
% ==============================================================================
% # This script provides the numerical back propagation for
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
% #
% # [2] T. Man, W. Zhang, L. Zhang, et al., "Scanning-free three-dimensional fluorescent dipoles imaging by polarization 
% self-interference digital holography (pSIDH)," DOI
% 10.48550/arXiv.2504.10772, 2025
% # https://github.com/TLBJUT/Polarization-Self-interference-Digital-Holography-pSIDH-.git
% 
% #==============================================================================
function [f1,dx1,dy1,x1,y1] = fresnel(f0,M,N,dx0,dy0,z,lambda)
k=2*pi/lambda;
du=1./(M*dx0);
dv=1./(N*dy0);
u=ones(N,1)*[0:M/2-1 -M/2:-1]*du;                 
v=[0:N/2-1 -N/2:-1]'*ones(1,M)*dv;
H=exp(-i*pi*lambda*z*(u.^2+v.^2));         
% f1=ifftshift(ifft2(fft2(fftshift(f0)).*H));      
% f1=ifft2(fft2(fftshift(f0)).*H);  
f1=ifft2(fft2(f0).*H);  
dx1=dx0;dy1=dy0;
x1=ones(N,1)*[-M/2:M/2-1]*dx1;                         

y1=[-N/2:N/2-1]'*ones(1,M)*dy1;
