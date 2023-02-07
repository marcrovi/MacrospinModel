%% Energy contributions
% Constants
% Units
clear ff field freq fields phiHl camps power degH
T = 7.95e5; %A/m
Oe = 1e-4*T;
mu_0 = 4*pi*1e-7;
thetaH = pi/2; %In-plane
phiHl = deg2rad(linspace(0,50,25)); % List of angles
% phiHl = [deg2rad(0) deg2rad(10) deg2rad(20) deg2rad(30) deg2rad(40) deg2rad(45) deg2rad(50) deg2rad(60) deg2rad(70) deg2rad(80) deg2rad(90)];
% phiHl = [deg2rad(45)];
% ----------- Materials ----------------- %
% FeSi
%     Ms = 955e3; % A/m
%     v = 3000; %m/s
%     Ku1 = 0; Ku2=0;
%     Kc1 = -4400; %J/m^3 FeSi
%     Ku4 = 0;
%     Aex = 1e-11;
%     phiA = deg2rad(45);
%     alpha = 0.003;
% Nickel
    Ms = 490e3; % A/m
    v = 3800; %m/s
    u = [1 0 0];
    Ku1 = 500; Ku2=0;
    Kc1 = 0;
    phiA = deg2rad(0); % angle respect to x
    Aex = 1e-11;
    alpha = 0.03;
% Cobalt
%     Ms = 151e4; % A/m
%     v = 3800; %m/s
%     Aex = 1e-11;
%     Kc1 = 000;
%     Ku1 = 0;
%     Ku2 = 0;
%     Ku4 = -3500;
%     phiA = deg2rad(0); % Angle between easy axes and applied field
%     alpha = 0.005;

% from J/m^3 -> mT: Ku4/Ms*1e3
% --------- SAW ---------- %
frequency = 1000e6;
lambda = v/frequency;
k = 2*pi/lambda;
d = 30e-9; %m
% ME constants
b1 = 6e6; % J/m^3 % Ni 3 GHz
% b1 = 6e6; % J/m^3 % Ni 1 GHz
b2 = 5e6; % J/m^3
% --------- Initial Conditions --------- %
% Camp inicial
H0 = 0*T; %T -> A/m
N_fields = 80;
camps = linspace(0,0.025*T,N_fields); % fields in Tesla (*T)
for j=1:length(phiHl)
    phiH = phiHl(j);
    H0 = camps(1); % T -> A/m
for i=1:N_fields
    syms theta phi
    % ---- Magnetization components in spherical coordinates --- %
    m = [sin(theta)*cos(phi); 
        sin(theta)*sin(phi); 
        cos(theta)]; 
    % --- Energy terms --- %
    % They are energy density J/m^3
    % Before input in frequency, we divide by Ms to transform to Tesla
    % ----- Effective fields in energy terms ------ %
    % Zeeman
    H0 = camps(i); % External field
    Ez = -mu_0*Ms.*H0*(sin(theta)*sin(thetaH)*cos(phi-phiH)+cos(theta)*cos(thetaH));

    % --- Uniaxial anisotropy --- %
    % u is the direction of the anisotropy
    % Taken from MuMax3
    Eu = -Ku1*dot(u,m)^2-Ku2*dot(u,m)^4;
    
    % --- Cubic anisotropy --- %
    Ea = Kc1/4*(sin(theta)^4*sin(2*phi)^2+sin(2*theta)^2);

    % --- Demagnetizing fields --- %
    Edem = mu_0*Ms^2/2*cos(theta)^2;
    
    % --- Dipolar energy (-mu_0*h_dip .(dot) m --- %
%     alpha = deg2rad(45); % angle between M and k (SAW at 45 degrees wrt H)
%     P = 1-(1-exp(-k*d))/(k*d);
%     M1 = [P*sin(alpha)^2 0 P*sin(alpha)*cos(alpha); 
%         0           1-P             0         ;
%       P*sin(alpha)*cos(alpha) 0 P*cos(alpha)^2];  
     
%     hdip = M1*m;
%     Edip = -mu_0*dot(hdip,m);
    
    % --- Exchange energy --- %
    hex1 = -2*Aex/(mu_0*Ms^2)*k^2*m;
    Eex1 = mu_0*dot(hex1,m);

    % --- Magnetoelastic waves --- %
    phiSAW = deg2rad(0); % Angle of SAW respect to x
%     eXX = 1.9e-4; eXZ = 3.93e-5; eZZ = 5.64e-5; % 3 GHz Ni
    eXX = 1.43e-4; eXZ = 2.95e-5; eZZ = 4.23e-5; % 1 GHz Ni
%     eXX = 1.66e-4; eXZ = 3.44e-5; eZZ = 4.93e-5; % 1 GHz Co
    mx = m(1); my = m(2); mz = m(3);
    exx = eXX*cos(phiSAW)^2;
    eyy = eXX*sin(phiSAW)^2;
    ezz = eZZ;
    exy = eXX*cos(phiSAW)*sin(phiSAW);
    exz = eXZ*cos(phiSAW);
    eyz = eXZ*sin(phiSAW);
    Eme = b1/Ms^2*(mx^2*exx+my^2*eyy+mz*ezz^2)+2*b2/Ms^2*(exy*mx*my+eyz*my*mz+exz*mx*mz);
    
    % --- TOTAL ENERGY --- %
    E = Ez+Ea+Edem+Eu+Eme+Eex1;
    syms theta
    theta0 = pi/2; % We assume that theta_0 = pi/2, so the magnetization is IN-PLANE
    theta = pi/2;
    
    syms phi theta

    Ep = subs(E);
    dEdp(phi,theta) = diff(subs(Ep,theta,theta0),phi);
    dEdt(phi,theta) = diff(Ep,theta);

    % -- Find angle that minimizes energy -- %
    Ecc = subs(dEdp);
    assume(phi>0)
    phi0 = vpasolve(Ecc==0, phi, [0 pi/2]);
    phi0d = round(rad2deg(phi0),2);
    phi = phi0;

    % Back to syms for double derivatives
    syms phi
    theta = theta0;
    dEdpp = diff(diff(subs(E),phi)) + mu_0*Ms^2*k*d/2*sin(phi)^2;
    syms theta
    phi = phi0;
    dEdtt = diff(diff(subs(E),theta)) - mu_0*Ms^2*k*d/2;

    % Asign the values found of theta_0 and phi_0 to compute
    phi = phi0;
    theta = theta0;

    % --- Frequency --- %
    gamma = 28; %Ghz/T
    % Energy is in terms of A/m, we divide by Ms in A/m %
    ftt = double(subs(dEdtt))/Ms; % Double derivatives evaluated in T
    fpp = double(subs(dEdpp))/Ms; % Double derivatives evaluated in T
    f = gamma/(sin(theta0))*sqrt(abs(fpp)*abs(ftt));
    ff(i) = double(subs(f)); % in GHz
    
    % --- Absorbed power --- %
    V0 = 1e28;%d*1e-3*100e-6; % Volume affected by SAW (I re-scaled)
    w = 2*pi*frequency;
    w0 = 2*pi*ff(i)*1e9; %in Hz
    gam = 28*1e9; % gamma in Hz/(A/m)
    ft = double(subs(dEdt))/Ms; %in T
    fp = double(subs(dEdp))/Ms; % in T
    B = double(fpp*(ft*eXZ)^2+ftt*(fp*eXX)^2-2*w/gam*ft*fp*eXX*eXZ);
    C = double(alpha*w/gam*((ft*eXZ)^2+(fp*eXX)^2));
    kappa(i) = gam*alpha/(1+alpha^2)*(fpp+ftt);
    Pabs(i) = V0/2*gam^2/(1+alpha^2)*w/((w^2-w0^2)^2+(w*kappa(i))^2)*(B*w*kappa(i)+C*(w^2-w0^2));
    field(i) = H0/T;
    fprintf("Step: %i, phi0d: %s, freq: %s \n",i, string(round(rad2deg(phi0),2)), string(round(f,2)))
end
% --- Save variables for each angle --- %
degH(j) = rad2deg(phiH);
freq{j} = double(ff);
fields{j} = field;
power{j} = Pabs;
end

save('Ni1GHzModel_Anis500x_10mT','freq','fields','degH','power')


%% Plot a single angle (last one)
hold on
plot(field(2:end)*1000,ff(2:end),'-b','LineWidth',1.5)
xlabel('$\mu_0 H$ (mT)','Interpreter','latex')
ylabel('$f$ (GHz)','Interpreter','latex')
set(gca,'FontName','CMU Serif','FontSize',20)
% legend('$f_{SAW} = 3$ GHz','$f_{SAW} = 1$ GHz','FMR','Interpreter','latex')
box on

%% Multiple plots
clear Legend
j = 1;
for i=1:4:length(freq)
    hold on
    plot(fields{i}(2:end)*1000,freq{i}(2:end),'-','LineWidth',1.5)
    Legend{j} = [num2str(round(degH(i),0)) '$^\circ$'];
    j = j+1;
end
xlabel('$\mu_0 H$ (mT)','Interpreter','latex')
ylabel('$f$ (GHz)','Interpreter','latex')
set(gca,'FontSize',20,'FontName','CMU Serif','LineWidth',1.2)
leg = legend(Legend,'Interpreter','latex');
title(leg,'$\varphi_H$','Interpreter','latex')
box on
%% Plot Pabs
j = length(power);
fin = 0;
for i=1:length(power)
    Pabb = power{i};
    % --- Build imagesc matrix --- %
    if i == 1
        mat = [double(subs(Pabb(2:end-fin)))];
    else
        mat = [mat; double(subs(Pabb(2:end-fin)))];
    end
   j = j - 1; 
end
h = fields{1}*1e3;
imagesc(h(1:end-fin),degH,abs(mat))
ylabel('$\theta$ (deg)','Interpreter','latex')
xlabel('$\mu_0 H$ (mT)','Interpreter','latex')
set(gca,'YDir','normal','FontSize',20,'FontName','CMU Serif')
% set(gca,'CLim',[0 0.000003])
colorbar
%% Plot 1 Pabs
i = 11;
j = 23;
k = 25;
ini = 2;
fin = 80;
Pabb = power{i};
fieldsP = fields{i}*1e3;
Pabb2 = power{j};
fieldsP2 = fields{j}*1e3;
Pabb3 = power{k};
fieldsP3 = fields{k}*1e3;
hold on
yyaxis left
plot(fieldsP(ini:fin),abs(Pabb(ini:fin)),'LineWidth',3)
ylabel('$P_{\textrm{abs}}$ (arb. units)','Interpreter','latex')
yyaxis right
plot(fieldsP2(ini:fin),abs(Pabb2(ini:fin)),'LineWidth',3)
plot(fieldsP3(ini:fin),abs(Pabb3(ini:fin)),'LineWidth',3)
% ylim([0 max(abs(Pabb))])
xlabel('$\mu_0 H$ (mT)','Interpreter','latex')
ylabel('$P_{\textrm{abs}}$ (arb. units)','Interpreter','latex')
set(gca,'FontSize',20,'FontName','CMU Serif')
legend(['$\theta =$' num2str(round(degH(i),0)-1) '$^\circ$'],['$\theta =$' num2str(round(degH(j),0)-1) '$^\circ$'],...
    ['$\theta =$' num2str(round(degH(k),0)) '$^\circ$'],'Interpreter','latex')
box on
    
    
