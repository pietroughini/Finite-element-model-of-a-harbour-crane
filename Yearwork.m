clear all
clc
close all
load('structure_fre.mat');
load('structure_frf.mat');
load('structure_mkr.mat');
load('x_Fx.mat');
load('x_Fy.mat');
load('fre_punto8.mat');
load('x_Fx_punto8.mat');


%% Punto 4
omega_vec = 2*pi*[1.34; 1.59; 4.03; 6.53; 7.56; 7.92];
alpha=0.2;
beta=4e-4;

for i = 1:6
    h_vec(i)=0.5*alpha/omega_vec(i,:)+omega_vec(i,:)*0.5*beta;
    omegad_vec(i)=omega_vec(i)*sqrt(1-h_vec(i)^2);
end

%% Punto 5
ndof=68;  %  3 d.o.f x 24 nodes - 4 d.o.c. 

% Partitioning of structural matrices 
MFF=M(1:ndof,1:ndof);
M_modal=phi.'*MFF*phi;
RFF=R(1:ndof,1:ndof);
R_modal=phi.'*RFF*phi;
KFF=K(1:ndof,1:ndof);
K_modal=phi.'*KFF*phi;


%% - FRFa -  
% Input: vertical force at point A
% Output: vertical displacment of point A

F0_42=1;

fmax=10;
df=0.01;
freq=[0:df:fmax];
nf=length(freq);

X0_42a=zeros(1,nf);

for ii=1:4  % we limit the responce to the first 4 modes
    ome=2*pi*freq;
    X0_42_ii=(phi(42,ii)*phi(42,ii)*F0_42)./(-ome.^2*M_modal(ii,ii)+i*ome*R_modal(ii,ii)+K_modal(ii,ii));
    %figure;semilogy(freq,abs(X0_42_ii));grid;title(['FRF mode # ' num2str(ii)]); 
    X0_42a=X0_42a+X0_42_ii;
end

blu         =       [0          0.4470  0.7410];    
giallo      =       [0.9290     0.6940  0.1250];

figure
subplot(2,1,1)
plot(freq, abs(X0_42a), 'Linewidth', 1.5,'LineStyle', '-','Color', giallo) % first 4 modes 
xlabel('[Hz]')
ylabel('[m/N]')
grid on
hold on
plot(freq, abs(x_Fy(42,:)),'Linewidth', 1.5,'LineStyle', '--','Color', blu) % all modes
xlabel('[Hz]')
ylabel('[m/N]')  
title('Y0A/F0yA')
lgd = legend('First 4 modes','All modes');
set(lgd, 'Position', [0.7, 0.8, 0.05, 0.05]);

subplot(2,1,2)
plot(freq, angle(X0_42a),'Linewidth', 1.5,'LineStyle', '-','Color', giallo)  % first 4 modes 
xlabel('[Hz]')
ylabel('[rad]')
grid on
hold on
plot(freq, angle(x_Fy(42,:)),'Linewidth', 1.5,'LineStyle', '--','Color', blu) % all modes
xlabel('[Hz]')
ylabel('[rad]') 
% OSS: we can't rely on modes beyond 10 Hz

%% - FRFb -
% Input: horizontal force at point A
% Output: horizontal displacment of point B

F0_41=1;
fmax=10;
df=0.01;
freq=[0:df:fmax];
nf=length(freq);
X0_23=zeros(1,nf);

for ii=1:4
    ome=2*pi*freq;
    X0_23_ii=(phi(23,ii)*phi(41,ii)*F0_41)./(-ome.^2*M_modal(ii,ii)+i*ome*R_modal(ii,ii)+K_modal(ii,ii));
    %figure;semilogy(freq,abs(X0_23_ii));grid;title(['FRF mode # ' num2str(ii)]);
    X0_23=X0_23+X0_23_ii;
end

blu         =       [0          0.4470  0.7410];    
giallo      =       [0.9290     0.6940  0.1250];

figure
subplot(2,1,1)
plot(freq, abs(X0_23), 'Linewidth', 1.5,'LineStyle', '-','Color', giallo) % first 4 modes 
xlabel('[Hz]')
ylabel('[m/N]')
grid on
hold on
plot(freq, abs(x_Fx(23,:)),'Linewidth', 1.5,'LineStyle', '--','Color', blu) % all modes
xlabel('[Hz]')
ylabel('[m/N]')  
title('X0B/F0xA')
lgd = legend('First 4 modes','All modes');
set(lgd, 'Position', [0.7, 0.8, 0.05, 0.05]);

subplot(2,1,2)
plot(freq, angle(X0_23),'Linewidth', 1.5,'LineStyle', '-','Color', giallo)  % first 4 modes 
xlabel('[Hz]')
ylabel('[rad]')
grid on
hold on
plot(freq, angle(x_Fx(23,:)),'Linewidth', 1.5,'LineStyle', '--','Color', blu) % all modes
xlabel('[Hz]')
ylabel('[rad]') 
% OSS: we can't rely on modes beyond 10 Hz
%% - FRFc -
% Input: vertical force at point A
% Output: vertical acceleration of point A

F0_42=1;
fmax=10;
df=0.01;
freq=[0:df:fmax];
nf=length(freq);
XPP0_42=zeros(1,nf);

for ii=1:4
    ome=2*pi*freq;
    X0_42_ii=(phi(42,ii)*phi(42,ii)*F0_42)./(-ome.^2*M_modal(ii,ii)+i*ome*R_modal(ii,ii)+K_modal(ii,ii));
    XPP0_42_ii=(-ome.^2).*(X0_42_ii);
    %figure;semilogy(freq,abs(XPP0_42_ii));grid;title(['FRF mode # ' num2str(ii)]);
    XPP0_42=XPP0_42+XPP0_42_ii;

end



figure
subplot(2,1,1)
plot(freq, abs(XPP0_42), 'Linewidth', 1.5,'LineStyle', '-','Color', giallo) % first 4 modes 
xlabel('[Hz]')
ylabel('[m^2/N*s]')
grid on
hold on
plot(freq, abs((-ome.^2).*(x_Fy(42,:))),'Linewidth', 1.5,'LineStyle', '--','Color', blu) % all modes
xlabel('[Hz]')
ylabel('[m^2/N*s]') 
title('Y0ppA/F0yA')
lgd = legend('First 4 modes','All modes');
set(lgd, 'Position', [0.7, 0.8, 0.05, 0.05]);

subplot(2,1,2)
plot(freq, angle(XPP0_42),'Linewidth', 1.5,'LineStyle', '-','Color', giallo)  % first 4 modes 
xlabel('[Hz]')
ylabel('[rad]')
grid on
hold on
plot(freq, angle((-ome.^2).*(x_Fy(42,:))),'Linewidth', 1.5,'LineStyle', '--','Color', blu) % all modes
xlabel('[Hz]')
ylabel('[rad]') 
% OSS: we can't rely on modes beyond 10 Hz

%% Punto 6
ndof=68;
ncoordinates=72;

MCF=M(ndof+1:ncoordinates,1:ndof);
RCF=R(ndof+1:ncoordinates,1:ndof);
KCF=K(ndof+1:ncoordinates,1:ndof);

% %% a. 
% Input: vertical force at point A  
% Output: vertical component of constraint force in the hinge O1

fmax=10;
df=0.01;
freq=[0:df:fmax];
nf=length(freq);
i=sqrt(-1);


for kk=1:nf
    ome=2*pi*freq(kk);
    r=(-ome^2*MCF+i*ome*RCF+KCF)*x_Fy(:,kk);
    V1=r(2);
    mod(kk)=abs(V1);
    phase(kk)=angle(V1);
end


figure
subplot(2,1,1)
plot(freq, mod, 'Linewidth', 1.5,'LineStyle', '-','Color', blu) 
xlabel('[Hz]')
ylabel('[N/N]')
title('V01/F0yA')
grid on

subplot(2,1,2)
plot(freq, phase,'Linewidth', 1.5,'LineStyle', '-','Color', blu)   
xlabel('[Hz]')
ylabel('[rad]')
grid on

%% b. 
% Input: horizontal force at point A  
% Output: vertical component of constraint force in the hinge O2

fmax=10;
df=0.01;
freq=[0:df:fmax];
nf=length(freq);
i=sqrt(-1);

for kk=1:nf
    ome=2*pi*freq(kk);
    r=(-ome^2*MCF+i*ome*RCF+KCF)*x_Fx(:,kk);
    V2=r(4);
    mod(kk)=abs(V2);
    phase(kk)=angle(V2);
end

figure
subplot(2,1,1)
plot(freq, mod, 'Linewidth', 1.5,'LineStyle', '-','Color', blu) 
xlabel('[Hz]')
ylabel('[N/N]')
title('V02/F0xA')
grid on

subplot(2,1,2)
plot(freq, phase,'Linewidth', 1.5,'LineStyle', '-','Color', blu)   
xlabel('[Hz]')
ylabel('[rad]')
grid on

%% Punto 7

% Constant
csi=0;
L_2=6;
E=2.06e11;
I_green=6.712e-4;


lamda_2=[0  1   0;  % rotation matrix of node 2 and 3 (in this case are the same since beam 3 has the same direction of beam,
        -1  0   0;  % otherwise orientation of node 3 would be different in general)
         0  0   1];
LAMDA_2=blkdiag(lamda_2, lamda_2);
shape_function_dd=[0;                       % double derivative in space (csi) of shape function
                   12*csi/L_2^3-6/L_2^2;
                   6*csi/L_2^2-4/L_2;
                   0;
                   -12*csi/L_2^3+6/L_2^2;
                    6*csi/L_2^2-2/L_2];
fmax=10;
df=0.01;
freq=[0:df:fmax];
nf=length(freq);

for kk=1:nf
    M_bending=E*I_green*shape_function_dd.'*LAMDA_2*x_Fx(2:7,kk);  
    M_mode(kk)=abs(M_bending);
    M_phase(kk)=angle(M_bending);
end

figure
subplot(3,1,1)
plot(freq, M_mode, '.-')
xlabel('[Hz]')
ylabel('[Nm]')
title('FRF')
grid on
subplot(3,1,2)
plot(freq, M_phase, '.-')
xlabel('[Hz]')
ylabel('[rad]')
grid on
subplot(3,1,3)
loglog(freq, M_mode, '.-') % in log scale
title('LogLog FRF')
xlabel('[Hz]')
ylabel('[Nm]')
grid on

%% Punto 8

% Definition of sampled signal Ypp_rel in [0,20] 

Ma=800;
g=9.81;
delta_t=0.02;    % sampling time (the inverse of sampling frequence fs) 
T=20;            % period of observation of Ypp_rel
N=T/delta_t;     % # of samples
t = 0:delta_t:(N-1)*delta_t; % vector of time


t_up1 = 0.2;   
t_down1 = 0.8;    
t_down2 = 1.2;   
t_up2 = 1.8;    

y_max = 5;  
y_min = -5;

Ypp_rel = zeros(size(t));

% Building the double trapeze
Ypp_rel(t >= 0 & t < t_up1) = linspace(0, y_max, sum(t >= 0 & t < t_up1));
Ypp_rel(t >= t_up1 & t < t_down1) = y_max;
Ypp_rel(t >= t_down1 & t < t_down2) = linspace(y_max, y_min, sum(t >= t_down1 & t < t_down2));
Ypp_rel(t >= t_down2 & t < t_up2) = y_min;
Ypp_rel(t >= t_up2 & t < 2) = linspace(y_min, 0, sum(t >= t_up2 & t < 2));

% Plot of Yrel(t)
figure;
subplot(2,1,1)
plot(t, Ypp_rel, 'LineWidth', 1);
xlabel('t [s]', 'Interpreter', 'latex');
ylabel('$\ddot{Y}_{rel}$ [m/s$^2$]', 'Interpreter', 'latex');
grid on;
title('Time history of the signal');
ylim([-8 8]);

%% ----
% Dft of the input signal
% We are assuming that Yrel is periodical in T, to define its dft
Ypp_rel_dft = fft(Ypp_rel);             
f = (0:N-1) / T ; % vector of frequencies for dft  
Y_magnitude = abs(Ypp_rel_dft);

% Plot of Ypp_rel_dft (Magnitude)
subplot(2,1,2);
plot(f, Y_magnitude, 'LineWidth', 1);
xlabel('Frequency [Hz]','Interpreter', 'latex');
ylabel('$|\ddot{Y}_{rel}|$', 'Interpreter', 'latex');
grid on;
title('Spectrum of the signal (DFT)');
xlim([0 10]); 

% Definition of the new frequency response function G(jW) and so of Ya(W)
ome=2*pi*f;
Gaa=x_Fx_punto8(42,:);  % same # of frequency samples of Ypp_rel_dft --> N

for kk=1:N
    G(kk)=-Gaa(kk)/(1-ome(kk)^2*Ma*Gaa(kk));
end

% Since our input signal is complex, hermitian symmetry isn't satisfied

Ya(1)=G(1)*(Ma*g+Ma*Ypp_rel_dft(1));


for kk=2:N
    Ya(kk)=G(kk)*Ma*Ypp_rel_dft(kk);
end


%  Plot of G(jw)
% figure
% semilogx(f, abs(G), '-')
% xlabel('[Hz]')
% ylabel('[m]')
% title('$|G(\omega)|$', 'Interpreter', 'latex')
figure
plot(f, abs(G), 'Linewidth', 1.5,'LineStyle', '-','Color', blu) 
xlabel('[Hz]')
ylabel('[m]')
title('$|G(\omega)|$', 'Interpreter', 'latex')
grid on

% Plot of Gaa(jw)
figure
plot(f, abs(G), 'Linewidth', 1.5,'LineStyle', '-','Color', blu) 
xlabel('[Hz]')
ylabel('[m]')
title('$|Gaa(\omega)|$', 'Interpreter', 'latex')
grid on
% figure
% semilogx(f, abs(Gaa), '-')
% xlabel('[Hz]')
% ylabel('[m]')
% title('$|Gaa(\omega)|$', 'Interpreter', 'latex')


% Plots of the response
figure
subplot(2,1,1)
plot(f, abs(Ya), '-')
xlabel('[Hz]')
ylabel('[m]')
title('$|Ya(\omega)|$', 'Interpreter', 'latex')
grid on
subplot(2,1,2)
plot(f, angle(Ya), '-')
xlabel('[Hz]')
ylabel('[rad]')
title('$\angle Ya(\omega)$', 'Interpreter', 'latex');
grid on

Ya_time = ifft(Ya, 'symmetric'); % 'symmetric' per forzare valori reali
figure;
plot(t, Ya_time, 'r');
title('Ya(t)');
xlabel('[s]');
ylabel('[t]');


%% Simulation (fictitous case): considering Gaa(jW) as the transfer function

for kk=1:N
    Ya_sim(kk)=Gaa(kk)*(Ma*g+Ma*Ypp_rel_dft(kk));
end

Ya_sim_time = ifft(Ya_sim, 'symmetric'); % 'symmetric' per forzare valori reali
figure;
plot(t, Ya_sim_time, 'r');
title('Ya(t) fictitous');
xlabel('[s]');
ylabel('[t]');

%% Point 9

L=[8,6,6,5.66,5.66,5.66,6,6,6,6,6,6,6,4,4,5,5,5,5.408,5.408,4.5,4.5,5,5,5,6,6,8];
m=[104.83,104.83,104.83,41.97,41.97,41.97,41.97,41.97,41.97,41.97,41.97,41.97,41.97,41.97,41.97,30.51,30.51,30.51,30.51,30.51,30.51,30.51,30.51,30.51,30.51,104.83,104.83,104.83];
M=dot(L,m);
Mb=2000;
Mtot=M+Mb;
Lextra=8.544;  % Length of the two added beams
Mextra=Lextra*30.51;  % Added mass
Increment= Mextra/Mtot*100;  % Percentage of mass increment by adding the two beams --> 2.48%

% By looking at the bode diagrams returned by software:
% Initial maximum amplitude: 0.539
% Final maximum displacement: 0.299
% --> Decrement of 44.5%

















