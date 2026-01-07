% 
% ==============================================================================
% # This script provides the holographic PSF generation for
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


function [Ip11,Ip12,Ip13,Ip14] = FINCH4_PSF(U1,M,N,dhx,dhy,f1,d1,ds,f_SLM,zh,lamda,p_ab)
sita11=0;sita12=1/2*pi;sita13=pi;sita14=3/2*pi;
k=2*pi/lamda;
x0=U1;
x=ones(N,1)*(-M/2:M/2-1).*dhx;
y=(ones(N,1)*(-M/2:M/2-1).*dhy)';
Ip11=zeros(M,N);
Ip12=zeros(M,N);
Ip13=zeros(M,N);
Ip14=zeros(M,N);
mask=zeros(M);
for xi=1:M
    for yi=1:M
      if (xi-M/2)^2+(yi-M/2)^2<=(M)^2
          mask(xi,yi)=1;
      end
    end
end
for xx=1:N
    for yy=1:N
        if x0(xx,yy)~=0
            xx1=-(yy-M/2-1).*dhx;yy1=-(xx-M/2-1).*dhy; 
            U=exp(1j*k/(2*(d1)).*((x-xx1).^2+(y-yy1).^2));
            x1=ones(N,1)*(-M/2:M/2-1).*dhx;
            y1=(ones(N,1)*(-M/2:M/2-1).*dhy)';
            U1=U.*exp(-1j*k/2/f1.*(x1.^2+y1.^2));
            [f1,dx3,dy3,x2,y2] = fresnel(U1,M,N,dhx,dhy,ds,lamda);    
            x3=ones(N,1)*(-M/2:M/2-1).*dhx;
            y3=(ones(N,1)*(-M/2:M/2-1).*dhy)';
            % The SIDH PSF (i.e. the hologram of a single point target) is
            % generated via coherent superposition (adding the complex
            % light field together) of the two beam after the holographic
            % beam splitter (for example, SLM). We focused on the sample-introduced optical aberrations,
            % so below "p_ab" is applied to both the two beams.
            R1=(0.5*p_ab+0.5*exp((-1i*pi/(lamda*f_SLM))*(x3.^2+y3.^2)+1i*sita11).*p_ab).*mask;
            R2=(0.5*p_ab+0.5*exp((-1i*pi/(lamda*f_SLM))*(x3.^2+y3.^2)+1i*sita12).*p_ab).*mask;
            R3=(0.5*p_ab+0.5*exp((-1i*pi/(lamda*f_SLM))*(x3.^2+y3.^2)+1i*sita13).*p_ab).*mask;
            R4=(0.5*p_ab+0.5*exp((-1i*pi/(lamda*f_SLM))*(x3.^2+y3.^2)+1i*sita14).*p_ab).*mask;
            U31=f1.*R1;
            U32=f1.*R2;
            U33=f1.*R3;
            U34=f1.*R4;
            [f21,dx2,dy2,x2,y2] = fresnel(U31,M,N,dhx,dhy,zh,lamda);
            [f22,dx2,dy2,x2,y2] = fresnel(U32,M,N,dhx,dhy,zh,lamda);
            [f23,dx2,dy2,x2,y2] = fresnel(U33,M,N,dhx,dhy,zh,lamda);
            [f24,dx2,dy2,x2,y2] = fresnel(U34,M,N,dhx,dhy,zh,lamda);
            I1=abs(f21).^2;
            I2=abs(f22).^2;
            I3=abs(f23).^2;
            I4=abs(f24).^2;
            Ip11=Ip11+I1;
            Ip12=Ip12+I2;
            Ip13=Ip13+I3;
            Ip14=Ip14+I4;

        end
    end

end
