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
plot(f,real(T1),'r');
hold on;
plot(f,real(T2),'g');
hold on;
plot(f,real(T3),'b');
hold off;

xlabel('频率/Hz');
ylabel('T');
title('透射率(理论值) e0,3e0');