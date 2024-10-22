%function dam_break
%Symulacja 2D pekniecia tamy
%Oparta na Shallow Water Equations
%Calkowanie metoda Lax-Wendroff
%MUP zima2324
%(c) M. Zimnoch 2023
clear;
clc;
global h;
global u;
global v;
%Parametry modelu
Y1_max=100;                 %Szerokosc zbiornika przed tama
Y2_max=40                   %Szerokosc rzeki po tamie
X_max=100                   %dlugosc obszaru symulacji
X_dam=50;                   %polozenie tamy
H1=10;                      %Glebokosc wody przed tama
H2=5;                       %Glebokosc wody za tama
Q_in=0;                     %Wielkosc doplywu do zbiornika
H_out=5;                    %Glebokosc wody przy ujsciu rzeki
dx=2;                      %rozdzielczosc modelu w kier x
dy=2;                      %rozdzielczosc modelu w kier y
dt=0.01;                   %krok czasowy
T_max=120;                 %czas symulacji w sekundach
nplotstep=5;                % Co ile krokow rysowac obraz
g=9.81;                     %przyspieszenie ziemskie
NX=round(X_max/dx);         %Ilosc komorek w kier x
NY=round(Y1_max/dy);        %Ilosc komorek w kier y
ND=round(X_dam/dx);         %Index x polozenia tamy
NY_L=round((Y1_max-Y2_max)/(2*dy));   %Indeks lewego brzegu rzeki za zaporą
NY_R=round((Y1_max+Y2_max)/(2*dy));   %Indeks prawego brzegu rzeki za zapora
NT=round(T_max/dt);
%Przygotowanie tablicy statusow
S=zeros(NX+2,NY+2);             %Tablica statusow wezlow
S(ND+1:NX+2,1:NY_L-1)=-1;     
S(ND+1:NX+2,NY_R+1:NY+2)=-1;
S(1,:) =1;
S(NX+2,NY_L:NY_R)=9; S(ND,1:NY_L-1)=2; S(ND,NY_R+1:NY+2)=2;
S(1:ND+2,1) =3; S(ND+1:NX+2,NY_L) =3;
S(1:ND+2,NY+2)=4; S(ND+1:NX+2,NY_R)=4;
S(1,1) =5;
S(1,NY+2)=6;
S(ND,1)=7;  S(ND,NY_L)=7; S(NX+2,NY_L)=7;
S(ND,NY+2)=8; S(ND,NY_R)=8; S(NX+2,NY_R)=8;

%Inicjalizacja zmiennych
h=zeros(NX+2,NY+2);
u=zeros(NX+2,NY+2);
v=zeros(NX+2,NY+2);
hx=zeros(NX+1,NY+1);
ux=zeros(NX+1,NY+1);
vx=zeros(NX+1,NY+1);
hy=zeros(NX+1,NY+1);
uy=zeros(NX+1,NY+1);
vy=zeros(NX+1,NY+1);
%Warunek poczatkowy
h(1:ND,:)=H1;
h(ND+1:NX+2,NY_L:NY_R)=H2;

[x y]= meshgrid(0:NX+1,0:NY+1);

for n=1:NT      %petla czasowa
    t=n*dt;
    % Warunek brzegowy z uzyciem tablisy S
    for i=1:NX+2
        for j=1:NY+2
            switch S(i,j)
                case 1 % Dopływ
                   h(i,j)=h(i+1,j);   u(i,j)=Q_in;     v(i,j)=v(i+1,j);
                case 2
                   h(i,j)=h(i-1,j);   u(i,j)=0;        v(i,j)=v(i-1,j);   
                case 3    
                   h(i,j)=h(i,j+1);   u(i,j)=u(i,j+1); v(i,j)=0;   
                case 4    
                   h(i,j)=h(i,j-1);   u(i,j)=u(i,j-1); v(i,j)=0;   
                case 5    
                   h(i,j)=h(i+1,j+1); u(i,j)=0;        v(i,j)=0;   
                case 6   
                   h(i,j)=h(i+1,j-1); u(i,j)=0;        v(i,j)=0;   
                case 7    
                   h(i,j)=h(i-1,j+1); u(i,j)=0;        v(i,j)=0;   
                case 8    
                   h(i,j)=h(i-1,j-1); u(i,j)=0;        v(i,j)=0;   
                case 9    
                   h(i,j)=H_out;      u(i,j)=u(i-1,j); v(i,j)=v(i-1,j);   
           end       
        end
    end         
    %Pierwszy krok
    % kierunek x
       i = 1:NX+1;
       j = 1:NY;
   
       % poziom
       hx(i,j) = (h(i+1,j+1)+h(i,j+1))/2 - dt/(2*dx)*(u(i+1,j+1)-u(i,j+1));
   
       % ped x 
       ux(i,j) = (u(i+1,j+1)+u(i,j+1))/2 -  ...
                 dt/(2*dx)*((u(i+1,j+1).^2./h(i+1,j+1) + g/2*h(i+1,j+1).^2) - ...
                            (u(i,j+1).^2./h(i,j+1) + g/2*h(i,j+1).^2));
   
       % ped y 
       vx(i,j) = (v(i+1,j+1)+v(i,j+1))/2 - ...
                 dt/(2*dx)*((u(i+1,j+1).*v(i+1,j+1)./h(i+1,j+1)) - ...
                            (u(i,j+1).*v(i,j+1)./h(i,j+1)));
       
       % kierunek y
       i = 1:NX;
       j = 1:NY+1;
   
       % poziom
       hy(i,j) = (h(i+1,j+1)+h(i+1,j))/2 - dt/(2*dy)*(v(i+1,j+1)-v(i+1,j));
   
       % ped x
       uy(i,j) = (u(i+1,j+1)+u(i+1,j))/2 - ...
                 dt/(2*dy)*((v(i+1,j+1).*u(i+1,j+1)./h(i+1,j+1)) - ...
                            (v(i+1,j).*u(i+1,j)./h(i+1,j)));
       % ped y
       vy(i,j) = (v(i+1,j+1)+v(i+1,j))/2 - ...
                 dt/(2*dy)*((v(i+1,j+1).^2./h(i+1,j+1) + g/2*h(i+1,j+1).^2) - ...
                            (v(i+1,j).^2./h(i+1,j) + g/2*h(i+1,j).^2));
   
       % Drugi krok
       i = 2:NX+1;
       j = 2:NY+1;
   
       % poziom
       h(i,j) = h(i,j) - (dt/dx)*(ux(i,j-1)-ux(i-1,j-1)) - ...
                         (dt/dy)*(vy(i-1,j)-vy(i-1,j-1));
       % ped x
       u(i,j) = u(i,j) - (dt/dx)*((ux(i,j-1).^2./hx(i,j-1) + g/2*hx(i,j-1).^2) - ...
                         (ux(i-1,j-1).^2./hx(i-1,j-1) + g/2*hx(i-1,j-1).^2)) ...
                       - (dt/dy)*((vy(i-1,j).*uy(i-1,j)./hy(i-1,j)) - ...
                         (vy(i-1,j-1).*uy(i-1,j-1)./hy(i-1,j-1)));
       % ped y
       v(i,j) = v(i,j) - (dt/dx)*((ux(i,j-1).*vx(i,j-1)./hx(i,j-1)) - ...
                         (ux(i-1,j-1).*vx(i-1,j-1)./hx(i-1,j-1))) ...
                       - (dt/dy)*((vy(i-1,j).^2./hy(i-1,j) + g/2*hy(i-1,j).^2) - ...
                         (vy(i-1,j-1).^2./hy(i-1,j-1) + g/2*hy(i-1,j-1).^2));

        % Wizualizacja
       if mod(n,nplotstep) == 0
          subplot(1,2,1);
          ts=sprintf('t = %6.2f',t);
          s1=surf(x,y,h);
          %s1.EdgeColor = 'none';
          %s1.FaceColor = 'interp';
          %s1.FaceLighting = 'gouraud';
          view(45,135);
          %colormap('turbo')';
          %colorbar;
          title(ts);
          xlabel('X');
          ylabel('Y');
          zlabel('level (m)');
          axis([0 NX+1 0 NY+1 -1 11]);
          subplot(1,2,2);
          quiver(x,y,v,u);
          xlabel('X');
          ylabel('Y');
          axis([0 NX+1 0 NY+1]);
          drawnow;
       end
end   

