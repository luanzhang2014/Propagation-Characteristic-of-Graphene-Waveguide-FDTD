einx=textread('einx.txt');
einy=textread('einy.txt');
einz=textread('einz.txt');

etx3=textread('etx3.txt');
ety3=textread('ety3.txt');
etz3=textread('etz3.txt');
%etx6=textread('etx6.txt');
%ety6=textread('ety6.txt');
%etz6=textread('etz6.txt');
%etx1=textread('etx5.txt');
%ety1=textread('ety5.txt');
%etz1=textread('etz5.txt');
%etx15=textread('etx15.txt');
%ety15=textread('ety15.txt');
%etz15=textread('etz15.txt');

dt=2.145*10^(-15);
nmax=5000;

fs=1/dt;
N=nmax;
n=0:N-1;
einx=fft(einx,N);
einy=fft(einy,N);
einz=fft(einz,N);

etx3=fft(etx3,N);
ety3=fft(ety3,N);
etz3=fft(etz3,N);
%etx6=fft(etx6,N);
%ety6=fft(ety6,N);
%etz6=fft(etz6,N);
%etx1=fft(etx1,N);
%ety1=fft(ety1,N);
%etz1=fft(etz1,N);
%etx15=fft(etx15,N);
%ety15=fft(ety15,N);
%etz15=fft(etz15,N);

f=(0-length(einx)/2:length(einx)/2-1)'*fs/length(einx);

N0=1*10^(13)*(length(einx))/fs;
f0=(0:N0-1)*fs/length(einx);


einx=fftshift(einx);
einy=fftshift(einy);
einz=fftshift(einz);

etx3=fftshift(etx3);
ety3=fftshift(ety3);
etz3=fftshift(etz3);
%etx6=fftshift(etx6);
%ety6=fftshift(ety6);
%etz6=fftshift(etz6);
%etx1=fftshift(etx1);
%ety1=fftshift(ety1);
%etz1=fftshift(etz1);
%etx15=fftshift(etx15);
%ety15=fftshift(ety15);
%etz15=fftshift(etz15);

ein=einx.*conj(einx)+einy.*conj(einy)+einz.*conj(einz);
et3=etx3.*conj(etx3)+ety3.*conj(ety3)+etz3.*conj(etz3);
%et6=etx6.*conj(etx6)+ety6.*conj(ety6)+etz6.*conj(etz6);
%et1=etx1.*conj(etx1)+ety1.*conj(ety1)+etz1.*conj(etz1);
%et15=etx15.*conj(etx15)+ety15.*conj(ety15)+etz15.*conj(etz15);

TT3=(et3./ein).^(1/2);
%TT6=(et6./ein).^(1/2);
%TT1=(et1./ein).^(1/2);
%TT15=et15./ein;

TT03=zeros(N0);
for n=1:N0
    TT03(n)=TT3(length(einx)/2+n);
    
end    
%TT06=zeros(N0);
%for n=1:N0
%    TT06(n)=TT6(length(einx)/2+n);
%end    
%TT01=zeros(N0);
%for n=1:N0
%    TT01(n)=TT1(length(einx)/2+n);
%end    
%TT015=zeros(N0);
%for n=1:N0
%    TT015(n)=TT15(length(einx)/2+n);
%end    

e=1.6*10^(-19);
u=99729*10^(-4);
uc1=0.3*e;
uc2=0.5*e;
uc3=0.6*e;
h=6.626*10^(-34);
hh=h/2/pi;
vf=9.5*10^5;
kb=1.38*10^(-23);
cc=2.99792458e8;  
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

eta1=(muz/epsz)^(1/2);
eta2=(muz/epsz/3)^(1/2);


ns1=uc1^2/pi/(hh*vf)^2;
ns2=uc2^2/pi/(hh*vf)^2;
ns3=uc3^2/pi/(hh*vf)^2;
T1=300;
tao1=u*hh*(ns1*pi)^(1/2)/(e*vf);
tao2=u*hh*(ns2*pi)^(1/2)/(e*vf);
tao3=u*hh*(ns3*pi)^(1/2)/(e*vf);
f=0:10^9:9.4*10^12;
w=2*pi*f;

y1=e^2*uc1*tao1/pi/hh^2./(1+j*w*tao1)+e^2*kb*T1*2*log(exp(-uc1/kb/T1)+1)/pi/hh^2./(1/tao1+j*w);
y2=e^2*uc2*tao2/pi/hh^2./(1+j*w*tao2)+e^2*kb*T1*2*log(exp(-uc2/kb/T1)+1)/pi/hh^2./(1/tao2+j*w);
y3=e^2*uc3*tao3/pi/hh^2./(1+j*w*tao3)+e^2*kb*T1*2*log(exp(-uc3/kb/T1)+1)/pi/hh^2./(1/tao3+j*w);

T1=2*eta2./(eta1+eta2+y1*eta1*eta2);
T2=2*eta2./(eta1+eta2+y2*eta1*eta2);
T3=2*eta2./(eta1+eta2+y3*eta1*eta2);

figure;
%plot(f,real(T1),'r');
%hold on;
%plot(f,real(T2),'g');
%hold on;
%plot(f,real(T3),'b');



plot(f0,TT03,'color','r');
%hold on;
%plot(f0,TT06,'color','b');
%hold on;
%plot(f0,TT01,'color','g');
%hold off;
%legend('0.3eV','0.6eV','1.0eV','0.15eV');
xlabel('ÆµÂÊ/Hz');
ylabel('T');
title('Í¸ÉäÂÊ e0,3e0');





