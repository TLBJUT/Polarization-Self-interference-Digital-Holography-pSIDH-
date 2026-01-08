% 
% ==============================================================================
% # This script provides the code for computational adaptive optics (CAO) aberration correction 
% on the phase reconstruction in self interference digital holography (SIDH).
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

clc;close all;clear;
% load('wavefront.mat');   %data loader, import your own data!
% load("f3.mat");     %data loader, import your own data!
% figure;imshow(angle(f3),[]);
% figure;imshow(phase_estimate,[]);
Re = zeros(512,512);
Re(255-180:255+180-1,255-180:255+180-1)=f3;
AO_phase = exp(1j*phase_estimate);
AO_phase = fftshift(AO_phase);
% AO_phase = imresize(AO_phase,[512 512]);
re = (fft2(f3));
re = re.*AO_phase;
% figure;imshow(AO_phase);

% Careful normalization processing is required here on the phase 

re = (ifft2(re));
re = angle(re);
re = re - min(min(re));
re = re/max(max(re));
re = 1-re;
re = re - min(min(re));
re = re/max(max(re));
figure;imshow(re);