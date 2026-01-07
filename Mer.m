% 
% ==============================================================================
% # This script provides the image metric calculation (image sharpness) for
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

clc;clear;close all;
i=double(imread('c12_2.tif'));
i=i/255;
% i=i(444:1809,412:1777);
% i=i(582:921,576:915);
[m n]=size(i);
figure;imshow(i);
I=fftshift(fft2(i));
I=abs(I);
sum_f=sum(sum(I));
Pf=zeros(m,n);
for x=1:n
    for y=1:m
        if (x-n/2)^2+(y-m/2)^2<=(110/2)^2
            Pf(x,y)=1;
        else
            Pf(x,y)=0;
        end
    end
end
I=I.*Pf;
figure;imshow(log((I))/10);
[M, N]=size(I);
S=0;
S_buffer=0;
for n=1:N
    for m=1:M
        S_buffer=I(n,m)*((n-(N-1)/2)^2+(m-(M-1)/2)^2);
        S=S+S_buffer;
    end
end

S=S/sum_f;
