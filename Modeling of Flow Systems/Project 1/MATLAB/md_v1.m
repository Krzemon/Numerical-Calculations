%Symulacja przepływu cieczy metodą dynamiki molekularnej.
%Zajęcia z Modelowania Układow Przepływowych
%(c) M. Zimnoch, v.1.0 30.10.2023
clear;
clc;
%Konfiguracja modelu
npartsx=20;             % Początkowa ilość cząstek w kierunku X
npartsy=20;             % Początkowa ilość cząstek w kierunku X
Lx=10;                  % wymiar domeny w kier x (-Lx - Lx)
Ly=10;                  % wymiar domeny w kier y (-Ly - Ly)
Xh=0.1;                 % czesc obszaru zajety przez obszar przyspieszania
Xc=0.1;                 % czesc obszaru zajety przez termostat
ixh=Xh*2*Lx-Lx;         % pozycja granicy obszaru przeyspieszania
ixc=(1-Xc)*2*Lx-Lx;     % pozycja granicy obszaru chlodzenia
FX=50;                  % Sila przyspieszajaca
Fc=0.99;                % Wspolczynnik termostatowania
nparts=npartsx*npartsy; % Ilość cząstek
%Paramtery potencjału oddziaływania
sigma=2*Lx/npartsx;
epsilon=5;
A=epsilon*sigma^12;
B=0.5*epsilon*sigma^6;
rcut=2.5*sigma;
m=1.0;                  % masa cząstek cieczy
mp=1.0e20;              % masa cząstek przeszkody
nsteps=10000;           % Ilość kroków czasowych
dt=0.005;               % timestep
%Inicjalizacja tablic
xpos=zeros(1,nparts);   % tablica połozen X
ypos=zeros(1,nparts);   % tablica polozen Y
xvel=zeros(1,nparts);   % tablica X skladowych predkosci
yvel=zeros(1,nparts);   % tablica Y skladowych predkosci
xvel_01=zeros(1,nparts);% tablica X skladowych predkosci dla korku t-1/2
yvel_01=zeros(1,nparts);% tablica Y skladowych predkosci dla korku t+1/2
xvel_11=zeros(1,nparts);% tablica X skladowych predkosci dla korku t-1/2
yvel_11=zeros(1,nparts);% tablica Y skladowych predkosci dla korku t+1/2
xfor=zeros(1,nparts);   % tablica X skladowych sily
yfor=zeros(1,nparts);   % tablica Y skladowych sily
mass=ones(1,nparts)*m;  % tablica mas czastek
flag=ones(1,nparts);    % tablica flag czastek

%Rozmieszczenie początkowe cząstek oraz rozlosowanie predkosci i zaznacenie
%czastek ciezkich
for ix=1:npartsx
    for iy=1:npartsy
        %nr czastki
        ipart=(ix-1)*npartsx+iy;
        %polozenia poczatkowe
        xpos(ipart)=(2*Lx/npartsx)*(ix-0.5)-Lx;
        ypos(ipart)=(2*Ly/npartsy)*(iy-0.5)-Ly;
        %Stworzenie scianek rury
        if (iy == 1 || iy == npartsy)
            mass(ipart)=mp;
            xvel(ipart)=0;
            yvel(ipart)=0;
        else
            %predkosci poczatkowe (na razie rozklad normalny)
            xvel(ipart)=randn(1);
            yvel(ipart)=randn(1);       
        end;    
    end;
end; 
xvel_01=xvel;
yvel_01=yvel;

%Petla czasowa
for n=1:nsteps
  %Obliczenie sil
  xfor=zeros(1,nparts);
  yfor=zeros(1,nparts);
  for i=1:nparts-1
      for j=i+1:nparts
          if (xpos(i)<(-Lx+sigma) && (Lx-xpos(j))<sigma)
              dx=2*Lx+xpos(i)-xpos(j);
          elseif ((Lx-xpos(i))<sigma && xpos(j)<(-Lx+sigma))
              dx=-2*Lx+xpos(i)-xpos(j);
          else        
               dx=xpos(i)-xpos(j);
          end;     
          if (ypos(i)<(-Ly+sigma) && (Ly-ypos(j))<sigma)
              dy=2*Ly+ypos(i)-ypos(j);
          elseif ((Ly-ypos(i))<sigma && ypos(j)<(-Ly+sigma))
              dy=-2*Ly+ypos(i)-ypos(j);
          else        
               dy=ypos(i)-ypos(j);
          end;
          d=sqrt(dx*dx+dy*dy);
          if d<rcut
              fx=LJ(dx,d,A, B);
              fy=LJ(dy,d,A, B);
              xfor(i)=xfor(i)+fx;
              xfor(j)=xfor(j)-fx;
              yfor(i)=yfor(i)+fy;
              yfor(j)=yfor(j)-fy;
          end;    
      end;    
   
  end;
  % Dodanie sily przyspieszajacej
  for i=1:nparts
    if (xpos(i) < ixh)  
        xfor(i)=xfor(i)+FX;
    end;        
  end;    
  %Calkowanie rownan ruchu
  xvel_11=xvel_01+dt*xfor./mass;
  xpos=xpos+dt*xvel_11;
  xvel=(xvel_01+xvel_11)./2;
  yvel_11=yvel_01+dt*yfor./mass;
  ypos=ypos+dt*yvel_11;
  yvel=(yvel_01+yvel_11)./2;
  xvel_01=xvel_11;
  yvel_01=yvel_11;
  % Termostat
  for i=1:nparts
      if (xpos(i) > ixc)
          xvel_01(i)=xvel_01(i)*Fc;
          yvel_01(i)=yvel_01(i)*Fc;
      end;
  end;    
  %Periodyczne war brzeg
  ix=fix(xpos./Lx);
  iy=fix(ypos./Ly);
  xpos=xpos-2*Lx*ix;
  ypos=ypos-2*Ly*iy;
  
  %Wizualizacja polozenia 
  
  if (mod(n,1)==0) 
    plot(xpos,ypos,'ob');
    axis([-Lx,Lx,-Ly,Ly]);
    title(['Time step: ' int2str(n)]);
    rectangle('Position',[-Lx -Ly 2*Xh*Lx 2*Ly],'EdgeColor','r','LineWidth',2);
    rectangle('Position',[ixc -Ly 2*Xc*Lx 2*Ly],'EdgeColor','b','LineWidth',2);
    pause(0.01);
  end;  
  
end;