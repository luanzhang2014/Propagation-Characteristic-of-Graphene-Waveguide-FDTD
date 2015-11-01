%***********************************************************************
%     3-D FDTD code with PML absorbing boundary conditions
%***********************************************************************

%    program auther:   Luan Zhang 
%                      Department of Optical Engineering , Zhejiang University
%                      zhangluan676@126.com
%                   
%    Date of this version: April 2014                  
%
%***********************************************************************
clear
clc
%***********************************************************************
%     Fundamental constants
%***********************************************************************
nm=1e-9;
cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space


freq=2*10^12;
omega=2.0*pi*freq;    
lambda=cc/freq;


%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=10;               %number of grid cells in x-direction
je=10;               %number of grid cells in y-direction
ke=116;              %number of grid cells in z-direction

ib=ie+1;
jb=je+1;
kb=ke+1;

is=5;                %location of  hard source in x axis
js=5;                %location of  hard source in y axis
ks=100;              %location of  hard source in z axis

ds=1532*nm;          %space increment of square lattice
dt=2.145*10^(-15);   %time step

nmax=500;            %total number of time steps
 
iebc=0;              %thickness of left and right PML region
jebc=0;              %thickness of front and back PML region
kebc=8;              %thickness of bottom and top PML region
rmax=0.00001;
orderbc=2;
ibbc=iebc+1;
jbbc=jebc+1;
kbbc=kebc+1;

iefbc=ie+2*iebc;
jefbc=je+2*jebc;
kefbc=ke+2*kebc;
ibfbc=iefbc+1;
jbfbc=jefbc+1;
kbfbc=kefbc+1;

%***********************************************************************
%     Material parameters
%***********************************************************************

media=1;

eps=[1.0,3.0];
sig=[0.0,0.0];
mur=[1.0,1.0];
sim=[0.0,0.0];
e=1.6*10^(-19);
u=99729*10^(-4);
uc=0.6*e;
h=6.626*10^(-34);
vf=9.5*10^5;
ns=uc^2/(pi*(vf*h/(2*pi))^2);
tao=u*(h/(2*pi))*(ns*pi)^(1/2)/(e*vf);
A=e^2*uc*tao/(h^2/(4*pi));

%***********************************************************************
%     Wave excitation
%***********************************************************************

for n=1:nmax
source(n)=exp(-4*pi*(n*dt-6.435*10^(-14))^2/(8.04375*10^(-14))^2);
end 

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(ie,jb,kb);           %fields in main grid 
ey=zeros(ib,je,kb);
ez=zeros(ib,jb,ke);

hx=zeros(ib,je,ke);
hy=zeros(ie,jb,ke);
hz=zeros(ie,je,kb);

Etranx=zeros(1,nmax);
Etrany=zeros(1,nmax);
Etranz=zeros(1,nmax);
Etotalx=zeros(1,nmax);
Etotaly=zeros(1,nmax);
Etotalz=zeros(1,nmax);
Esourcex=zeros(1,nmax);
Esourcey=zeros(1,nmax);
Esourcez=zeros(1,nmax);
Htranx=zeros(1,nmax);
Htrany=zeros(1,nmax);
Htranz=zeros(1,nmax);

Jx=zeros(ie,jb,kb);
Jy=zeros(ib,je,kb);


exybcd=zeros(ie,jb,kebc);%%           %fields in bottom PML region
exzbcd=zeros(ie,jb,kebc);%%
eyzbcd=zeros(ib,je,kebc);%%
eyxbcd=zeros(ib,je,kebc);%%
ezxbcd=zeros(ib,jb,kebc);%%
ezybcd=zeros(ib,jb,kebc);%%

hxybcd=zeros(ib,je,kebc);%%
hxzbcd=zeros(ib,je,kebc);%%
hyzbcd=zeros(ie,jb,kebc);%%
hyxbcd=zeros(ie,jb,kebc);%%
hzybcd=zeros(ie,je,kebc);%%
hzxbcd=zeros(ie,je,kebc);%%


exybct=zeros(ie,jb,kbbc);%%          %fields in top PML region
exzbct=zeros(ie,jb,kbbc);%%
eyzbct=zeros(ib,je,kbbc);%% 
eyxbct=zeros(ib,je,kbbc);%%
ezxbct=zeros(ib,jb,kebc);%%
ezybct=zeros(ib,jb,kebc);%%

hxybct=zeros(ib,je,kebc);%%
hxzbct=zeros(ib,je,kebc);%%
hyzbct=zeros(ie,jb,kebc);%%
hyxbct=zeros(ie,jb,kebc);%%
hzxbct=zeros(ie,je,kbbc);%%
hzybct=zeros(ie,je,kbbc);%%


%***********************************************************************
%     Updating coefficients
%***********************************************************************

for i=1:media
  eaf=dt*sig(i)/(2.0*epsz*eps(i));
  ca(i)=(1.0-eaf)/(1.0+eaf);
  cb(i)=dt/epsz/eps(i)/ds/(1.0+eaf);
  haf=dt*sim(i)/(2.0*muz*mur(i));
  da(i)=(1.0-haf)/(1.0+haf);
  db(i)=dt/muz/mur(i)/ds/(1.0+haf);
end

%***********************************************************************
%     main grid
%***********************************************************************

%     Initialize entire main grid to free space

caex(1:ie,1:jb,1:58)=ca(1);     
cbex(1:ie,1:jb,1:58)=cb(1);

caey(1:ib,1:je,1:58)=ca(1);
cbey(1:ib,1:je,1:58)=cb(1);

caez(1:ib,1:jb,1:58)=ca(1);
cbez(1:ib,1:jb,1:58)=cb(1);


dahx(1:ib,1:je,1:58)=da(1);
dbhx(1:ib,1:je,1:58)=db(1);

dahy(1:ie,1:jb,1:58)=da(1);
dbhy(1:ie,1:jb,1:58)=db(1);

dahz(1:ie,1:je,1:58)=da(1);
dbhz(1:ie,1:je,1:58)=db(1);


caex(1:ie,1:jb,59:kb)=ca(1);     
cbex(1:ie,1:jb,59:kb)=cb(1);

caey(1:ib,1:je,59:kb)=ca(1);
cbey(1:ib,1:je,59:kb)=cb(1);

caez(1:ib,1:jb,59:ke)=ca(1);
cbez(1:ib,1:jb,59:ke)=cb(1);


dahx(1:ib,1:je,59:ke)=da(1);
dbhx(1:ib,1:je,59:ke)=db(1);

dahy(1:ie,1:jb,59:ke)=da(1);
dbhy(1:ie,1:jb,59:ke)=db(1);

dahz(1:ie,1:je,59:kb)=da(1);
dbhz(1:ie,1:je,59:kb)=db(1);


%***********************************************************************
%     Fill the PML regions
%***********************************************************************

delbc=kebc*ds;
sigmam=-log(rmax/100.0)*epsz*cc*(orderbc+1)/(2*delbc);
bcfactor=eps(1)*sigmam/(ds*(delbc^orderbc)*(orderbc+1));


% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%     BOTTOM region

caexybcd(1:ie,1:jb,1)=1.0;
cbexybcd(1:ie,1:jb,1)=0.0;
caexzbcd(1:ie,1:jb,1)=1.0;
cbexzbcd(1:ie,1:jb,1)=0.0;
caeyxbcd(1:ib,1:je,1)=1.0;
cbeyxbcd(1:ib,1:je,1)=0.0;
caeyzbcd(1:ib,1:je,1)=1.0;
cbeyzbcd(1:ib,1:je,1)=0.0;

dahzxbcd(1:ie,1:je,1)=1.0;
dbhzxbcd(1:ie,1:je,1)=0.0;
dahzybcd(1:ie,1:je,1)=1.0;
dbhzybcd(1:ie,1:je,1)=0.0;


for k=2:kebc                               % 与sigmaz有关的量
  z1=(kebc-k+1.5)*ds;
  z2=(kebc-k+0.5)*ds;
  sigmaz=bcfactor*(z1^(orderbc+1)-z2^(orderbc+1));
  ca1=exp(-sigmaz*dt/(epsz*eps(1)));
  cb1=(1.0-ca1)/(sigmaz*ds);  
  
  caexzbcd(1:ie,1:jb,k)=ca1;
  cbexzbcd(1:ie,1:jb,k)=cb1;  
  caeyzbcd(1:ib,1:je,k)=ca1;
  cbeyzbcd(1:ib,1:je,k)=cb1;  
  
  caexybcd(1:ie,1:jb,k)=ca(1);
  cbexybcd(1:ie,1:jb,k)=cb(1);  
  caeyxbcd(1:ib,1:je,k)=ca(1);                   
  cbeyxbcd(1:ib,1:je,k)=cb(1);  
  dahzxbcd(1:ie,1:je,k)=da(1);
  dbhzxbcd(1:ie,1:je,k)=db(1);  
  dahzybcd(1:ie,1:je,k)=da(1);
  dbhzybcd(1:ie,1:je,k)=db(1);   
end

sigmaz = bcfactor*(0.5*ds)^(orderbc+1); 
ca1=exp(-sigmaz*dt/(epsz*eps(1)));
cb1=(1.0-ca1)/(sigmaz*ds);
caex(1:ie,1:jb,1)=ca1;
cbex(1:ie,1:jb,1)=cb1;
caey(1:ib,1:je,1)=ca1;
cbey(1:ib,1:je,1)=cb1;

for k=1:kebc                               % 与sigmazs有关的量
  y1=(kebc-k+1)*ds;
  y2=(kebc-k)*ds;
  sigmaz=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  sigmazs=sigmaz*(muz/(epsz*eps(1)));
  da1=exp(-sigmazs*dt/muz);
  db1=(1.0-da1)/(sigmazs*ds);  
  
  dahxzbcd(1:ib,1:je,k)=da1;
  dbhxzbcd(1:ib,1:je,k)=db1;
  dahyzbcd(1:ie,1:jb,k)=da1;
  dbhyzbcd(1:ie,1:jb,k)=db1;   
  
  caezxbcd(1:ib,1:jb,k)=ca(1);              %与sigmaz、sigmazs无关的量      
  cbezxbcd(1:ib,1:jb,k)=cb(1);  
  caezybcd(1:ib,1:jb,k)=ca(1);                   
  cbezybcd(1:ib,1:jb,k)=cb(1);   
  dahxybcd(1:ib,1:je,k)=da(1);  
  dbhxybcd(1:ib,1:je,k)=db(1);  
  dahyxbcd(1:ie,1:jb,k)=da(1);
  dbhyxbcd(1:ie,1:jb,k)=db(1); 
end


%     TOP region

caexybct(1:ie,1:jb,kbbc)=1.0;
cbexybct(1:ie,1:jb,kbbc)=0.0;
caexzbct(1:ie,1:jb,kbbc)=1.0;
cbexzbct(1:ie,1:jb,kbbc)=0.0;
caeyxbct(1:ib,1:je,kbbc)=1.0;
cbeyxbct(1:ib,1:je,kbbc)=0.0;
caeyzbct(1:ib,1:je,kbbc)=1.0;
cbeyzbct(1:ib,1:je,kbbc)=0.0;

dahzxbct(1:ie,1:je,kbbc)=1.0;
dbhzxbct(1:ie,1:je,kbbc)=0.0;
dahzybct(1:ie,1:je,kbbc)=1.0;
dbhzybct(1:ie,1:je,kbbc)=0.0;

for k=1:kebc                               % 与sigmaz有关的量
  z1=(k-0.5)*ds;
  z2=(k-1.5)*ds;  
  sigmaz=bcfactor*(z1^(orderbc+1)-z2^(orderbc+1));
  ca1=exp(-sigmaz*dt/(epsz*eps(1)));
  cb1=(1.0-ca1)/(sigmaz*ds);  
  
  
  caexzbct(1:ie,1:jb,k)=ca1;
  cbexzbct(1:ie,1:jb,k)=cb1;  
  caeyzbct(1:ib,1:je,k)=ca1;
  cbeyzbct(1:ib,1:je,k)=cb1;  
  
end


sigmaz = bcfactor*(0.5*ds)^(orderbc+1); 
ca1=exp(-sigmaz*dt/(epsz*eps(1)));
cb1=(1.0-ca1)/(sigmaz*ds);
caex(1:ie,1:jb,kb)=ca1;
cbex(1:ie,1:jb,kb)=cb1;
caey(1:ib,1:je,kb)=ca1;
cbey(1:ib,1:je,kb)=cb1;



for k=1:kebc                               % 与sigmazs有关的量
  z1=k*ds;
  z2=(k-1)*ds;
  sigmaz=bcfactor*(z1^(orderbc+1)-z2^(orderbc+1));
  sigmazs=sigmaz*(muz/(epsz*eps(1)));
  da1=exp(-sigmazs*dt/muz);
  db1=(1.0-da1)/(sigmazs*ds);  
  
  
  dahxzbct(1:ib,1:je,k)=da1;
  dbhxzbct(1:ib,1:je,k)=db1;
  dahyzbct(1:ie,1:jb,k)=da1;
  dbhyzbct(1:ie,1:jb,k)=db1;   
  
    %与sigmaz、sigmazs无关的量
  caexybct(1:ie,1:jb,k)=ca(1);                   
  cbexybct(1:ie,1:jb,k)=cb(1);  
  caeyxbct(1:ib,1:je,k)=ca(1);                   
  cbeyxbct(1:ib,1:je,k)=cb(1);   
  caezxbct(1:ib,1:jb,k)=ca(1);                   
  cbezxbct(1:ib,1:jb,k)=cb(1);  
  caezybct(1:ib,1:jb,k)=ca(1);                   
  cbezybct(1:ib,1:jb,k)=cb(1);      
  dahxybct(1:ib,1:je,k)=da(1);  
  dbhxybct(1:ib,1:je,k)=db(1);  
  dahyxbct(1:ie,1:jb,k)=da(1);
  dbhyxbct(1:ie,1:jb,k)=db(1);   
  dahzxbct(1:ie,1:je,k)=da(1);
  dbhzxbct(1:ie,1:je,k)=db(1);  
  dahzybct(1:ie,1:je,k)=da(1);
  dbhzybct(1:ie,1:je,k)=db(1); 
end



% ***********************************************************************
%     Movie initialization
% ***********************************************************************

tview(:,:)=real(ex(:,:,80));
sview(:,:)=real(ex(4,:,:));

subplot('position',[0.15 0.45 0.7 0.45]),pcolor(tview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['real Hz(i,j,k=',ks,') ', 'time step = 0']);
xlabel('i coordinate');
ylabel('j coordinate');

subplot('position',[0.15 0.10 0.7 0.25]),pcolor(sview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['real Hz(i=15,j,k), time step = 0']);
xlabel('j coordinate');
ylabel('k coordinate');


% mview(:,:)=real(hz(15,je/2,:));
% tview(:,:)=real(hzxbct(15,je/2,:)+hzybct(15,je/2,:));
% subplot(2,1,1),plot(mview),caxis([-1.0 1.0]);
% ylabel('hz in main')
% subplot(2,1,2),plot(tview),caxis([-1.0 1.0]);
% ylabel('hz in left')

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/2,gcf,rect);

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
%     Update electric fields (Ex 、Ey、Ez) in main grid
%***********************************************************************
    
ex(:,:,100)=ex(:,:,100)+source(n);

ex(1:ie,2:je,58)=caex(1:ie,2:je,58).*ex(1:ie,2:je,58)+...     
           cbex(1:ie,2:je,58).*(hz(1:ie,2:je,58)-hz(1:ie,1:je-1,58)-...
                                    hy(1:ie,2:je,58)+hy(1:ie,2:je,57))-...
                                2*dt*Jx(1:ie,2:je,58)/((epsz+epsz));
ex(1:ie,2:je,2:57)=caex(1:ie,2:je,2:57).*ex(1:ie,2:je,2:57)+...
           cbex(1:ie,2:je,2:57).*(hz(1:ie,2:je,2:57)-hz(1:ie,1:je-1,2:57)-...
                                    hy(1:ie,2:je,2:57)+hy(1:ie,2:je,1:56));
ex(1:ie,2:je,59:ke)=caex(1:ie,2:je,59:ke).*ex(1:ie,2:je,59:ke)+...
           cbex(1:ie,2:je,59:ke).*(hz(1:ie,2:je,59:ke)-hz(1:ie,1:je-1,59:ke)-...
                                    hy(1:ie,2:je,59:ke)+hy(1:ie,2:je,58:ke-1));
                                
ey(2:ie,1:je,58)=caey(2:ie,1:je,58).*ey(2:ie,1:je,58)+... 
           cbey(2:ie,1:je,58).*(hx(2:ie,1:je,58)-hx(2:ie,1:je,57)-...
                              hz(2:ie,1:je,58)+hz(1:ie-1,1:je,58))-...
                                2*dt*Jy(2:ie,1:je,58)/((epsz+epsz));                            
ey(2:ie,1:je,2:57)=caey(2:ie,1:je,2:57).*ey(2:ie,1:je,2:57)+...
           cbey(2:ie,1:je,2:57).*(hx(2:ie,1:je,2:57)-hx(2:ie,1:je,1:56)-...
                              hz(2:ie,1:je,2:57)+hz(1:ie-1,1:je,2:57));
ey(2:ie,1:je,59:ke)=caey(2:ie,1:je,59:ke).*ey(2:ie,1:je,59:ke)+...
           cbey(2:ie,1:je,59:ke).*(hx(2:ie,1:je,59:ke)-hx(2:ie,1:je,58:ke-1)-...
                              hz(2:ie,1:je,59:ke)+hz(1:ie-1,1:je,59:ke));
                          
ez(2:ie,2:je,:)=caez(2:ie,2:je,:).*ez(2:ie,2:je,:)+...
           cbez(2:ie,2:je,:).*(hy(2:ie,2:je,:)-hy(1:ie-1,2:je,:)-...
                               hx(2:ie,2:je,:)+hx(2:ie,1:je-1,:));
 
%***********************************************************************
%     PBC
%***********************************************************************
                                                        
ex(1:ie,1,58)=caex(1:ie,1,58).*ex(1:ie,1,58)+... 
           cbex(1:ie,1,58).*(hz(1:ie,1,58)-hz(1:ie,je,58)-...
                                    hy(1:ie,1,58)+hy(1:ie,1,57))-...
                                2*dt*Jx(1:ie,1,58)/((epsz+epsz));
ex(1:ie,1,2:57)=caex(1:ie,1,2:57).*ex(1:ie,1,2:57)+...
           cbex(1:ie,1,2:57).*(hz(1:ie,1,2:57)-hz(1:ie,je,2:57)-...
                                    hy(1:ie,1,2:57)+hy(1:ie,1,1:56));
ex(1:ie,1,59:ke)=caex(1:ie,1,59:ke).*ex(1:ie,1,59:ke)+...
           cbex(1:ie,1,59:ke).*(hz(1:ie,1,59:ke)-hz(1:ie,je,59:ke)-...
                                    hy(1:ie,1,59:ke)+hy(1:ie,1,58:ke-1));
                                
ex(1:ie,jb,58)=caex(1:ie,jb,58).*ex(1:ie,jb,58)+...
           cbex(1:ie,jb,58).*(hz(1:ie,1,58)-hz(1:ie,je,58)-...   
                                    hy(1:ie,jb,58)+hy(1:ie,jb,57))-...
                                2*dt*Jx(1:ie,jb,58)/((epsz+epsz));
ex(1:ie,jb,2:57)=caex(1:ie,jb,2:57).*ex(1:ie,jb,2:57)+...
           cbex(1:ie,jb,2:57).*(hz(1:ie,1,2:57)-hz(1:ie,je,2:57)-...
                                    hy(1:ie,jb,2:57)+hy(1:ie,jb,1:56));
ex(1:ie,jb,59:ke)=caex(1:ie,jb,59:ke).*ex(1:ie,jb,59:ke)+...
           cbex(1:ie,jb,59:ke).*(hz(1:ie,1,59:ke)-hz(1:ie,je,59:ke)-...
                                    hy(1:ie,jb,59:ke)+hy(1:ie,jb,58:ke-1));

                                
ey(1,1:je,58)=caey(1,1:je,58).*ey(1,1:je,58)+...
           cbey(1,1:je,58).*(hx(1,1:je,58)-hx(1,1:je,57)-...
                              hz(1,1:je,58)+hz(ie,1:je,58))-...
                                2*dt*Jy(1,1:je,58)/((epsz+epsz));
ey(1,1:je,2:57)=caey(1,1:je,2:57).*ey(1,1:je,2:57)+...
           cbey(1,1:je,2:57).*(hx(1,1:je,2:57)-hx(1,1:je,1:56)-...
                              hz(1,1:je,2:57)+hz(ie,1:je,2:57));
ey(1,1:je,59:ke)=caey(1,1:je,59:ke).*ey(1,1:je,59:ke)+...
           cbey(1,1:je,59:ke).*(hx(1,1:je,59:ke)-hx(1,1:je,58:ke-1)-...
                              hz(1,1:je,59:ke)+hz(ie,1:je,59:ke));    
                          
ey(ib,1:je,58)=caey(ib,1:je,58).*ey(ib,1:je,58)+...
           cbey(ib,1:je,58).*(hx(ib,1:je,58)-hx(ib,1:je,57)-...
                              hz(1,1:je,58)+hz(ie,1:je,58))-...
                                2*dt*Jy(ib,1:je,58)/((epsz+epsz));
ey(ib,1:je,2:57)=caey(ib,1:je,2:57).*ey(ib,1:je,2:57)+...
           cbey(ib,1:je,2:57).*(hx(ib,1:je,2:57)-hx(ib,1:je,1:56)-...
                              hz(1,1:je,2:57)+hz(ie,1:je,2:57));
ey(ib,1:je,59:ke)=caey(ib,1:je,59:ke).*ey(ib,1:je,59:ke)+...
           cbey(ib,1:je,59:ke).*(hx(ib,1:je,59:ke)-hx(ib,1:je,58:ke-1)-...
                              hz(1,1:je,59:ke)+hz(ie,1:je,59:ke));    
                                                    
                          
ez(1,2:je,:)=caez(1,2:je,:).*ez(1,2:je,:)+...
           cbez(1,2:je,:).*(hy(1,2:je,:)-hy(ie,2:je,:)-...
                               hx(1,2:je,:)+hx(1,1:je-1,:));
ez(ib,2:je,:)=caez(ib,2:je,:).*ez(ib,2:je,:)+...
           cbez(ib,2:je,:).*(hy(1,2:je,:)-hy(ie,2:je,:)-...
                               hx(ib,2:je,:)+hx(ib,1:je-1,:));                          
ez(2:ie,1,:)=caez(2:ie,1,:).*ez(2:ie,1,:)+...
           cbez(2:ie,1,:).*(hy(2:ie,1,:)-hy(1:ie-1,1,:)-...
                               hx(2:ie,1,:)+hx(2:ie,je,:));                          
ez(2:ie,jb,:)=caez(2:ie,jb,:).*ez(2:ie,jb,:)+...
           cbez(2:ie,jb,:).*(hy(2:ie,jb,:)-hy(1:ie-1,jb,:)-...
                               hx(2:ie,1,:)+hx(2:ie,je,:));

                 
ez(1,1,:)=caez(1,1,:).*ez(1,1,:)+...
           cbez(1,1,:).*(hy(1,1,:)-hy(ie,1,:)-...
                               hx(1,1,:)+hx(1,je,:));
ez(1,jb,:)=caez(1,jb,:).*ez(1,jb,:)+...
           cbez(1,jb,:).*(hy(1,jb,:)-hy(ie,jb,:)-...
                               hx(1,1,:)+hx(1,je,:));                          
ez(ib,1,:)=caez(ib,1,:).*ez(ib,1,:)+...
           cbez(ib,1,:).*(hy(1,1,:)-hy(ie,1,:)-...
                               hx(ib,1,:)+hx(ib,je,:));                          
ez(ib,jb,:)=caez(ib,jb,:).*ez(ib,jb,:)+...
           cbez(ib,jb,:).*(hy(1,jb,:)-hy(ie,jb,:)-...
                               hx(ib,1,:)+hx(ib,je,:));
                           
Jx(1:ie,1:jb,58)=(2*tao-dt)*Jx(1:ie,1:jb,58)/(2*tao+dt)+2*A*dt*ex(1:ie,1:jb,58)/(2*tao+dt)/ds;
Jy(1:ib,1:je,58)=(2*tao-dt)*Jy(1:ib,1:je,58)/(2*tao+dt)+2*A*dt*ey(1:ib,1:je,58)/(2*tao+dt)/ds;


%***********************************************************************
%     Update Exy in PML regions
%***********************************************************************


%   BOTTOM
                              
exybcd(:,2:je,1:kebc)=caexybcd(:,2:je,1:kebc).*exybcd(:,2:je,1:kebc)+...
  cbexybcd(:,2:je,1:kebc).*(hzxbcd(:,2:je,1:kebc)+hzybcd(:,2:je,1:kebc)-...
                                 hzxbcd(:,1:je-1,1:kebc)-hzybcd(:,1:je-1,1:kebc));                           
                             
exybcd(:,1,1:kebc)=caexybcd(:,1,1:kebc).*exybcd(:,1,1:kebc)+...
  cbexybcd(:,1,1:kebc).*(hzxbcd(:,1,1:kebc)+hzybcd(:,1,1:kebc)-...
                                 hzxbcd(:,je,1:kebc)-hzybcd(:,je,1:kebc)); 
                             
exybcd(:,jb,1:kebc)=caexybcd(:,jb,1:kebc).*exybcd(:,jb,1:kebc)+...
  cbexybcd(:,jb,1:kebc).*(hzxbcd(:,1,1:kebc)+hzybcd(:,1,1:kebc)-...
                                 hzxbcd(:,je,1:kebc)-hzybcd(:,je,1:kebc));  
                             
                                                          
                                   
%   TOP

                             
exybct(:,2:je,2:kebc)=caexybct(:,2:je,2:kebc).*exybct(:,2:je,2:kebc)+...
  cbexybct(:,2:je,2:kebc).*(hzxbct(:,2:je,2:kebc)+hzybct(:,2:je,2:kebc)-...
                                 hzxbct(:,1:je-1,2:kebc)-hzybct(:,1:je-1,2:kebc));
                             
exybct(:,1,2:kebc)=caexybct(:,1,2:kebc).*exybct(:,1,2:kebc)+...
  cbexybct(:,1,2:kebc).*(hzxbct(:,1,2:kebc)+hzybct(:,1,2:kebc)-...
                                 hzxbct(:,je,2:kebc)-hzybct(:,je,2:kebc));

exybct(:,jb,2:kebc)=caexybct(:,jb,2:kebc).*exybct(:,jb,2:kebc)+...
  cbexybct(:,jb,2:kebc).*(hzxbct(:,1,2:kebc)+hzybct(:,1,2:kebc)-...
                                 hzxbct(:,je,2:kebc)-hzybct(:,je,2:kebc));
 
 
%***********************************************************************
%     Update Exz in PML regions
%***********************************************************************
 
 
%   BOTTOM

exzbcd(:,1:jb,2:kebc)=caexzbcd(:,1:jb,2:kebc).*exzbcd(:,1:jb,2:kebc)-...
  cbexzbcd(:,1:jb,2:kebc).*(hyzbcd(:,1:jb,2:kebc)+hyxbcd(:,1:jb,2:kebc)-...
                             hyzbcd(:,1:jb,1:kebc-1)-hyxbcd(:,1:jb,1:kebc-1));
                                 
ex(:,1:jb,1)=caex(:,1:jb,1).*ex(:,1:jb,1)-...
  cbex(:,1:jb,1).*(hy(:,1:jb,1)-hyzbcd(:,1:jb,kebc)-hyxbcd(:,1:jb,kebc));


%   TOP

exzbct(:,:,2:kebc)=caexzbct(:,:,2:kebc).*exzbct(:,:,2:kebc)-...
  cbexzbct(:,:,2:kebc).*(hyzbct(:,:,2:kebc)+hyxbct(:,:,2:kebc)-...
                                 hyxbct(:,:,1:kebc-1)-hyzbct(:,:,1:kebc-1));
                             
ex(:,1:jb,kb)=caex(:,1:jb,kb).*ex(:,1:jb,kb)-...
  cbex(:,1:jb,kb).*(hyzbct(:,1:jb,1)+hyxbct(:,1:jb,1)-hy(:,1:jb,ke)); 
                             
                             
                            
%***********************************************************************
%     Update Eyz in PML regions
%***********************************************************************



%   BOTTOM

eyzbcd(1:ib,:,2:kebc)=caeyzbcd(1:ib,:,2:kebc).*eyzbcd(1:ib,:,2:kebc)+...
  cbeyzbcd(1:ib,:,2:kebc).*(hxybcd(1:ib,:,2:kebc)+hxzbcd(1:ib,:,2:kebc)-...
                             hxybcd(1:ib,:,1:kebc-1)-hxzbcd(1:ib,:,1:kebc-1));

ey(1:ib,:,1)=caey(1:ib,:,1).*ey(1:ib,:,1)+...
  cbey(1:ib,:,1).*(hx(1:ib,:,1)-hxybcd(1:ib,:,kebc)-hxzbcd(1:ib,:,kebc)); 


%   TOP

eyzbct(:,:,2:kebc)=caeyzbct(:,:,2:kebc).*eyzbct(:,:,2:kebc)+...
  cbeyzbct(:,:,2:kebc).*(hxybct(:,:,2:kebc)+hxzbct(:,:,2:kebc)-...
                             hxybct(:,:,1:kebc-1)-hxzbct(:,:,1:kebc-1));

ey(1:ib,:,kb)=caey(1:ib,:,kb).*ey(1:ib,:,kb)+...
  cbey(1:ib,:,kb).*(hxybct(1:ib,:,1)+hxzbct(1:ib,:,1)-hx(1:ib,:,ke)); 


%***********************************************************************
%     Update Eyx in PML regions
%***********************************************************************


%   BOTTOM

eyxbcd(2:ie,:,1:kebc)=caeyxbcd(2:ie,:,1:kebc).*eyxbcd(2:ie,:,1:kebc)-...
  cbeyxbcd(2:ie,:,1:kebc).*(hzxbcd(2:ie,:,1:kebc)+hzybcd(2:ie,:,1:kebc)-...
                             hzxbcd(1:ie-1,:,1:kebc)-hzybcd(1:ie-1,:,1:kebc));

eyxbcd(ib,:,1:kebc)=caeyxbcd(ib,:,1:kebc).*eyxbcd(ib,:,1:kebc)-...
  cbeyxbcd(ib,:,1:kebc).*(hzxbcd(1,:,1:kebc)+hzybcd(1,:,1:kebc)-...
                             hzxbcd(ie,:,1:kebc)-hzybcd(ie,:,1:kebc));

eyxbcd(1,:,1:kebc)=caeyxbcd(1,:,1:kebc).*eyxbcd(1,:,1:kebc)-...
  cbeyxbcd(1,:,1:kebc).*(hzxbcd(1,:,1:kebc)+hzybcd(1,:,1:kebc)-...
                             hzxbcd(ie,:,1:kebc)-hzybcd(ie,:,1:kebc));

%   TOP

eyxbct(2:ie,:,2:kebc)=caeyxbct(2:ie,:,2:kebc).*eyxbct(2:ie,:,2:kebc)-...
  cbeyxbct(2:ie,:,2:kebc).*(hzxbct(2:ie,:,2:kebc)+hzybct(2:ie,:,2:kebc)-...
                             hzxbct(1:ie-1,:,2:kebc)-hzybct(1:ie-1,:,2:kebc));

eyxbct(ib,:,2:kebc)=caeyxbct(ib,:,2:kebc).*eyxbct(ib,:,2:kebc)-...
  cbeyxbct(ib,:,2:kebc).*(hzxbct(1,:,2:kebc)+hzybct(1,:,2:kebc)-...
                             hzxbct(ie,:,2:kebc)-hzybct(ie,:,2:kebc));

eyxbct(1,:,2:kebc)=caeyxbct(1,:,2:kebc).*eyxbct(1,:,2:kebc)-...
  cbeyxbct(1,:,2:kebc).*(hzxbct(1,:,2:kebc)+hzybct(1,:,2:kebc)-...
                             hzxbct(ie,:,2:kebc)-hzybct(ie,:,2:kebc));
 
%***********************************************************************
%     Update Ezx in PML regions
%***********************************************************************

                                           
%   BOTTOM

ezxbcd(2:ie,:,1:kebc)=caezxbcd(2:ie,:,1:kebc).*ezxbcd(2:ie,:,1:kebc)+...
  cbezxbcd(2:ie,:,1:kebc).*(hyzbcd(2:ie,:,1:kebc)+hyxbcd(2:ie,:,1:kebc)-...
                       hyzbcd(1:ie-1,:,1:kebc)-hyxbcd(1:ie-1,:,1:kebc));

ezxbcd(ib,:,1:kebc)=caezxbcd(ib,:,1:kebc).*ezxbcd(ib,:,1:kebc)+...
  cbezxbcd(ib,:,1:kebc).*(hyzbcd(1,:,1:kebc)+hyxbcd(1,:,1:kebc)-...
                       hyzbcd(ie,:,1:kebc)-hyxbcd(ie,:,1:kebc));

ezxbcd(1,:,1:kebc)=caezxbcd(1,:,1:kebc).*ezxbcd(1,:,1:kebc)+...
  cbezxbcd(1,:,1:kebc).*(hyzbcd(1,:,1:kebc)+hyxbcd(1,:,1:kebc)-...
                       hyzbcd(ie,:,1:kebc)-hyxbcd(ie,:,1:kebc));

%   TOP

ezxbct(2:ie,:,1:kebc)=caezxbct(2:ie,:,1:kebc).*ezxbct(2:ie,:,1:kebc)+...
  cbezxbct(2:ie,:,1:kebc).*(hyzbct(2:ie,:,1:kebc)+hyxbct(2:ie,:,1:kebc)-...
                       hyzbct(1:ie-1,:,1:kebc)-hyxbct(1:ie-1,:,1:kebc));

ezxbct(ib,:,1:kebc)=caezxbct(ib,:,1:kebc).*ezxbct(ib,:,1:kebc)+...
  cbezxbct(ib,:,1:kebc).*(hyzbct(1,:,1:kebc)+hyxbct(1,:,1:kebc)-...
                       hyzbct(ie,:,1:kebc)-hyxbct(ie,:,1:kebc));

ezxbct(1,:,1:kebc)=caezxbct(1,:,1:kebc).*ezxbct(1,:,1:kebc)+...
  cbezxbct(1,:,1:kebc).*(hyzbct(1,:,1:kebc)+hyxbct(1,:,1:kebc)-...
                       hyzbct(ie,:,1:kebc)-hyxbct(ie,:,1:kebc));

%***********************************************************************
%     Update Ezy in PML regions
%***********************************************************************

               
%   BOTTOM

ezybcd(:,2:je,:)=caezybcd(:,2:je,:).*ezybcd(:,2:je,:)-...
  cbezybcd(:,2:je,:).*(hxybcd(:,2:je,:)+hxzbcd(:,2:je,:)-...
                       hxybcd(:,1:je-1,:)-hxzbcd(:,1:je-1,:));

ezybcd(:,jb,:)=caezybcd(:,jb,:).*ezybcd(:,jb,:)-...
  cbezybcd(:,jb,:).*(hxybcd(:,1,:)+hxzbcd(:,1,:)-...
                       hxybcd(:,je,:)-hxzbcd(:,je,:));

ezybcd(:,1,:)=caezybcd(:,1,:).*ezybcd(:,1,:)-...
  cbezybcd(:,1,:).*(hxybcd(:,1,:)+hxzbcd(:,1,:)-...
                       hxybcd(:,je,:)-hxzbcd(:,je,:));

%   TOP

ezybct(:,2:je,:)=caezybct(:,2:je,:).*ezybct(:,2:je,:)-...
  cbezybct(:,2:je,:).*(hxybct(:,2:je,:)+hxzbct(:,2:je,:)-...
                       hxybct(:,1:je-1,:)-hxzbct(:,1:je-1,:));

ezybct(:,jb,:)=caezybct(:,jb,:).*ezybct(:,jb,:)-...
  cbezybct(:,jb,:).*(hxybct(:,1,:)+hxzbct(:,1,:)-...
                       hxybct(:,je,:)-hxzbct(:,je,:));

ezybct(:,1,:)=caezybct(:,1,:).*ezybct(:,1,:)-...
  cbezybct(:,1,:).*(hxybct(:,1,:)+hxzbct(:,1,:)-...
                       hxybct(:,je,:)-hxzbct(:,je,:));
                   
                   
               
                   
                           
%***********************************************************************
%     Update magnetic fields (Hx、Hy、Hz) in main grid
%***********************************************************************

hx(1:ib,1:je,1:ke)=dahx(1:ib,1:je,1:ke).*hx(1:ib,1:je,1:ke)-... 
              dbhx(1:ib,1:je,1:ke).*(ez(1:ib,2:jb,1:ke)-ez(1:ib,1:je,1:ke)+...
                                ey(1:ib,1:je,1:ke)-ey(1:ib,1:je,2:kb));

hy(1:ie,1:jb,1:ke)=dahy(1:ie,1:jb,1:ke).*hy(1:ie,1:jb,1:ke)-... 
              dbhy(1:ie,1:jb,1:ke).*(ex(1:ie,1:jb,2:kb)-ex(1:ie,1:jb,1:ke)+...
                                ez(1:ie,1:jb,1:ke)-ez(2:ib,1:jb,1:ke));

hz(1:ie,1:je,1:kb)=dahz(1:ie,1:je,1:kb).*hz(1:ie,1:je,1:kb)-... 
              dbhz(1:ie,1:je,1:kb).*(ey(2:ib,1:je,1:kb)-ey(1:ie,1:je,1:kb)+...
                                 ex(1:ie,1:je,1:kb)-ex(1:ie,2:jb,1:kb));


%***********************************************************************
%     Update Hxy in PML regions
%***********************************************************************


%    BOTTOM

hxybcd(:,:,:)=dahxybcd(:,:,:).*hxybcd(:,:,:)-...
  dbhxybcd(:,:,:).*(ezxbcd(:,2:jb,:)+ezybcd(:,2:jb,:)-...
                              ezxbcd(:,1:je,:)-ezybcd(:,1:je,:));

                                                
%    TOP

hxybct(:,:,:)=dahxybct(:,:,:).*hxybct(:,:,:)-...
  dbhxybct(:,:,:).*(ezxbct(:,2:jb,:)+ezybct(:,2:jb,:)-...
                              ezxbct(:,1:je,:)-ezybct(:,1:je,:));

%***********************************************************************
%     Update Hxz in PML regions
%***********************************************************************



%    BOTTOM

                          
hxzbcd(:,:,1:kebc-1)=dahxzbcd(:,:,1:kebc-1).*hxzbcd(:,:,1:kebc-1)+...
  dbhxzbcd(:,:,1:kebc-1).*(eyzbcd(:,:,2:kebc)+eyxbcd(:,:,2:kebc)-...
                              eyzbcd(:,:,1:kebc-1)-eyxbcd(:,:,1:kebc-1));
                          
hxzbcd(:,:,kebc)=dahxzbcd(:,:,kebc).*hxzbcd(:,:,kebc)+...
  dbhxzbcd(:,:,kebc).*(ey(:,:,1)-eyzbcd(:,:,kebc)-eyxbcd(:,:,kebc));                              
                          
                                                                          
%    TOP

hxzbct(:,:,2:kebc)=dahxzbct(:,:,2:kebc).*hxzbct(:,:,2:kebc)+...
  dbhxzbct(:,:,2:kebc).*(eyzbct(:,:,3:kbbc)+eyxbct(:,:,3:kbbc)-...
                              eyzbct(:,:,2:kebc)-eyxbct(:,:,2:kebc));
                          
hxzbct(:,:,1)=dahxzbct(:,:,1).*hxzbct(:,:,1)+...
  dbhxzbct(:,:,1).*(eyzbct(:,:,2)+eyxbct(:,:,2)-ey(:,:,kb));
                          
%***********************************************************************
%     Update Hyz in PML regions
%***********************************************************************
 

%    BOTTOM

hyzbcd(:,:,1:kebc-1)=dahyzbcd(:,:,1:kebc-1).*hyzbcd(:,:,1:kebc-1)-...
  dbhyzbcd(:,:,1:kebc-1).*(exybcd(:,:,2:kebc)+exzbcd(:,:,2:kebc)-...
                              exybcd(:,:,1:kebc-1)-exzbcd(:,:,1:kebc-1));

hyzbcd(:,:,kebc)=dahyzbcd(:,:,kebc).*hyzbcd(:,:,kebc)-...
  dbhyzbcd(:,:,kebc).*(ex(:,:,1)-exybcd(:,:,kebc)-exzbcd(:,:,kebc));                           
                                                
%    TOP

hyzbct(:,:,2:kebc)=dahyzbct(:,:,2:kebc).*hyzbct(:,:,2:kebc)-...
  dbhyzbct(:,:,2:kebc).*(exybct(:,:,3:kbbc)+exzbct(:,:,3:kbbc)-...
                              exybct(:,:,2:kebc)-exzbct(:,:,2:kebc));

hyzbct(:,:,1)=dahyzbct(:,:,1).*hyzbct(:,:,1)-...
  dbhyzbct(:,:,1).*(exybct(:,:,2)+exzbct(:,:,2)-ex(:,:,kb));  


%***********************************************************************
%     Update Hyx in PML regions
%***********************************************************************


%    BOTTOM

hyxbcd(:,:,:)=dahyxbcd(:,:,:).*hyxbcd(:,:,:)+...
  dbhyxbcd(:,:,:).*(ezxbcd(2:ib,:,:)+ezybcd(2:ib,:,:)-...
                              ezxbcd(1:ie,:,:)-ezybcd(1:ie,:,:));                      
                                                
%    TOP

hyxbct(:,:,:)=dahyxbct(:,:,:).*hyxbct(:,:,:)+...
  dbhyxbct(:,:,:).*(ezxbct(2:ib,:,:)+ezybct(2:ib,:,:)-...
                              ezxbct(1:ie,:,:)-ezybct(1:ie,:,:)); 


%***********************************************************************
%     Update Hzx in PML regions
%***********************************************************************


%    BOTTOM

hzxbcd(:,:,:)=dahzxbcd(:,:,:).*hzxbcd(:,:,:)-...
  dbhzxbcd(:,:,:).*(eyzbcd(2:ib,:,:)+eyxbcd(2:ib,:,:)-...
                              eyzbcd(1:ie,:,:)-eyxbcd(1:ie,:,:));                                            
                                                
%    TOP

hzxbct(:,:,2:kebc)=dahzxbct(:,:,2:kebc).*hzxbct(:,:,2:kebc)-...
  dbhzxbct(:,:,2:kebc).*(eyzbct(2:ib,:,2:kebc)+eyxbct(2:ib,:,2:kebc)-...
                              eyzbct(1:ie,:,2:kebc)-eyxbct(1:ie,:,2:kebc));
                          

                          
%***********************************************************************
%     Update Hzy in PML regions
%***********************************************************************

%    BOTTOM

hzybcd(:,:,:)=dahzybcd(:,:,:).*hzybcd(:,:,:)+...
  dbhzybcd(:,:,:).*(exybcd(:,2:jb,:)+exzbcd(:,2:jb,:)-...
                              exybcd(:,1:je,:)-exzbcd(:,1:je,:));
                                                
%    TOP

hzybct(:,:,2:kebc)=dahzybct(:,:,2:kebc).*hzybct(:,:,2:kebc)+...
  dbhzybct(:,:,2:kebc).*(exybct(:,2:jb,2:kebc)+exzbct(:,2:jb,2:kebc)-...
                              exybct(:,1:je,2:kebc)-exzbct(:,1:je,2:kebc));


%***********************************************************************
%     Visualize fields
%***********************************************************************

if mod(n,2)==0;

timestep=int2str(n);
tview(:,:)=real(ex(:,:,110));
sview(:,:)=real(ex(4,:,:));

subplot('position',[0.15 0.45 0.7 0.45]),pcolor(tview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['real Ex(i,j,k=30',ks,') ', 'time step =',timestep]);
xlabel('i coordinate');
ylabel('j coordinate');

subplot('position',[0.15 0.10 0.7 0.25]),pcolor(sview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['real Ex(i=4,j,k), time step =',timestep]);
xlabel('j coordinate');
ylabel('k coordinate');

 mview(:,:)=abs(ex(5,je/2,:));
% mview(:,:)=real(exybct(5,je/2,:)+exzbct(5,je/2,:));
% tview(:,:)=real(hzxbct(15,je/2,:)+hzybct(15,je/2,:));
 subplot(2,1,1),plot(mview),caxis([-1.0 1.0]);
 xlabel('z');
 ylabel('ex');
% subplot(2,1,2),plot(tview),caxis([-1.0 1.0]);
% ylabel('hz in left')

nn=n/2;
M(:,nn)=getframe(gcf,rect);

end;

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************
Etranx(n)=ex(5,5,50);
Etrany(n)=ey(5,5,50);
Etranz(n)=ez(5,5,50);

Etotalx(n)=ex(5,5,70);
Etotaly(n)=ey(5,5,70);
Etotalz(n)=ez(5,5,70);

end

movie(gcf,M,0,10,rect);
 
t=(1*dt:dt:(nmax)*dt);
FEtranx=fft(Etranx,nmax);
FEtrany=fft(Etrany,nmax);
FEtranz=fft(Etranz,nmax);
FEtranx=fftshift(FEtranx);
FEtrany=fftshift(FEtrany);
FEtranz=fftshift(FEtranz);
Etran0=zeros(1,nmax);
for n=1:nmax
Etran0(n)=FEtranx(n)*conj(FEtranx(n))+FEtrany(n)*conj(FEtrany(n))+FEtranz(n)*conj(FEtranz(n));
end

n=1:nmax;

figure;
plot(n,Etran0);
title('Etran0');

figure;
plot(n,abs(Etranx(n)));
title('Etranx');

figure;
plot(n,abs(FEtranx(n)));
title('FEtranx');





