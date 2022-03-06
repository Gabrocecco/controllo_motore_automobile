%% 1) Definizione costanti date dal testo, definizioni funzioni f_1(x,u) e f_2(x,u), ricavo coppia di equilibrio ((x_e1, x_e2), u_e)
%

gamma_1 = 0.5;, 
gamma_2 = 0.1;
beta = 1.1;
psi = 0.02;
delta_1 = 5e4;
delta_2 = 0.1;
delta_3 = 0.01;
J = 40;
omega_e = 10;
% definisco funzione x_2' = f_2(x,u) 
x_2 = @(x) delta_1 / J * x - delta_2 / J * omega_e - delta_3 / J * omega_e * omega_e; 

m_e = fzero(x_2, 0)

x_1 = @(x) gamma_1 * (1 - cos(beta * x - psi)) - gamma_2 * omega_e * m_e;   %???

x0=[0 0.018];       %???
x1=[0.018 0.06];

%L'equazione ha due radici, quindi ho due valori di u_e, u_e_0 = 0.0067 e u_e_1 = 0.0297 
u_e_0 = fzero(x_1, x0)  

u_e_1 = fzero(x_1, x1)  %trovo u_e = 0.0297 (valore che scelgo -> da un valore positivo al gudagno nella linearizzazione)

%% Linearisation
    
u_e = u_e_1;   %considero u_e_1 come equilibrio e_e

%scrivo le matrici del sistema calcolate con i jacobiani
A = [(-gamma_2 * omega_e) (-gamma_2 * m_e); (delta_1 / J) (-delta_2 / J - 2 * delta_3 * omega_e / J)];

B = [beta * gamma_1 * sin(beta * u_e - psi); 0];

C = [0 1];

D = 0;

s = tf('s');
%applico la formula della G(s) per trovare la funzione di trasferimento del
%sistema linearizzato che chiamerò GG
GG = zpk(C * inv(s * eye(2) - A) * B + D);

[gnum,  gden] = tfdata(GG);

gnum = cell2mat(gnum);  
gden = cell2mat(gden);

%% Diagramma di bode della G(s) nel range di pulsazioni 10^-2 fino a 10^4

omega_plot_min = 1e-2;
omega_plot_max = 1e4;

figure(1)
h_G = bodeplot(GG,{omega_plot_min,omega_plot_max});
grid on, zoom on;

%% STATIC  REGULATOR %visualizzare zona propibita uscita d(t)

R_s = 1.6; %introduciamo un gudagmno statico pari a 1.6
GGe = R_s*GG;   %calcoliamo la G_e(s) e grafichiamola.

h_GGe = bodeplot(GGe,{omega_plot_min,omega_plot_max});

%% DYNAMIC REGULATOR

Mf_star_a = 83; %deciso precedentemente durante la specifica di sovraelongazione.
omega_c_star_a = 120;   %deciso arbitrariamente tramite tentativi
%mi ricavo |G(jw_c_star)| e arg(G(jw_c_star)) che servono per ricavare
%M_star_a e phi_star_a
[mag_omega_c_star_a, arg_omega_c_star_a, omega_c_star_a] = bode(GGe, omega_c_star_a);

%trasformo in dB |G(jw_c_star)| 
mag_omega_c_star_dB_a = 20*log10(mag_omega_c_star_a);

%calcolo M_star_a e phi_star_a tramite le due formule note.
M_star_a = 10^(-mag_omega_c_star_dB_a/20)
phi_star_a = Mf_star_a - 180 - arg_omega_c_star_a

% ricavo i valori calibrati di tao e alpha con le formule di inversione a
% usando i valori di M_star_a e phi_star_a
tau_b_a = (M_star_a - cos(phi_star_a*pi/180))/omega_c_star_a/sin(phi_star_a*pi/180)
alpha_tau_b_a = (cos(phi_star_a*pi/180) - inv(M_star_a))/omega_c_star_a/sin(phi_star_a*pi/180)
alpha_b_a = alpha_tau_b_a / tau_b_a

% avendo i valori calibrati di tao e alpha posso scrivere diterattemente
% dalla sua formula l'espressione della rete di anticipo.
R_d_a =(1 + tau_b_a*s)/(1 + alpha_tau_b_a * s);

% avendo calcolato sia il regolatore statico e adesso quello dinamico,
% posso scrivere la L(jw)=R_s(jw)*R_d(jw)*G(jw)
LL = R_d_a*GGe;

% controllo che che cos(phi_star_a) > 1/M_star_a
check_flag = cos(phi_star_a*pi/180) - inv(M_star_a)
if check_flag < 0
    return;
end
% regolatore R(jw)
reg = R_d_a * R_s;

[rnum,  rden] = tfdata(reg);

rnum = cell2mat(rnum);
rden = cell2mat(rden);

%%  Visualizzo zona proibita dovuta a disturbo di uscita d(t)

figure(2)
  
% prendo in considerazione il range in cui ha effetto il disturbo di uscita
% basse frequenze (date dal testo)
omega_d_min = 0.0001;
omega_d_MAX = 0.075;
Ad = 45;
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [Ad; Ad; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
    
%% Visualizzo zona proibita dovuta a disturbo di misura n(t)

omega_n_min = 5e3;
omega_n_MAX = 5e6;
An = 45;
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-An; -An; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
    
%% calcolo sovraelongazione e margne di fase.

S_100_spec = 1;
xi = 0.83;  % calclocato nella specifica sulla sovraelongazione
S_100 = 100*exp(-pi*xi/sqrt(1-xi^2))    % applico formula sovraelongazione
Mf_spec = xi*100 % trovo margine di fase (approssimazione di F(s) a poli cc dominanti.
    
%% Disegno zone proibite di d(t) e di margine di fase 

Ta1_spec = 0.05;
omega_Ta_low = 1e-4;
omega_Ta_MAX = 72;  
    
Bnd_Ta_x = [omega_Ta_low; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_low];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
    
h_GGe = bodeplot(LL);   
grid on, zoom on;
    
omega_c_min = 72; % valore minimo di omega_c deciso nella specifica del tempo di assestamento. 
omega_c_MAX = 5e3; % valore massimo di omega_c dovuto alla specifica sull'attenuazione di n(t)  
    
phi_spec = Mf_spec - 180;   %soglia minima della fase per rispettare il vincolo sul margine di fase.
phi_low = -270;
    
Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
    
    
%% Risposta al gradino nel sistema linearizzato
    
FF = LL/(1+LL);
figure(4);
     
WW = 0.5    %altezza gradino
Tfinal = 0.5;
[y_step,t_step]=step(WW*FF,Tfinal);
S_FF = stepinfo(FF,'SettlingTimeThreshold',0.05);
st_FF = S_FF.SettlingTime   %tempo di assestamento effettivo

plot(t_step,y_step)
grid on, zoom on
patch([0,Tfinal,Tfinal,0],[WW*(1+S_100_spec/100),WW*(1+S_100_spec/100),WW+1,WW+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; ylim([0,WW+1]);

patch([Ta1_spec,Tfinal,Tfinal,Ta1_spec],[WW*(1-0.05),WW*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([Ta1_spec,Tfinal,Tfinal,Ta1_spec],[WW*(1+0.05),WW*(1+0.05),WW+1,WW+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
    
figure(5);
bode(FF)
grid on, zoom on
    
%% attenuazione d(t) nel tempo
    
SS = 1/(1+LL);
figure(6);
    
omega_d = 0.01875;
tt=(0:1e-2:1e3)';
DD = 0.3;
dd = DD*(cos(omega_d*tt) + cos(omega_d*tt*2) + cos(omega_d*tt*3) + cos(omega_d*tt*4));
y_d = lsim(SS,dd,tt);   
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('dd','y_d')
    
%% attenuazione n(t) nel tempo
    
figure(7);
    
omega_d = 5e3;
tt=(0:1e-6:4e-3)';
DD = 0.1;
nn = DD*(cos(omega_d*tt) + cos(omega_d*tt*2) + cos(omega_d*tt*3) + cos(omega_d*tt*4));
y_s = lsim(-FF,nn,tt);   
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_s,'b')
grid on
legend('nn','y_s')
    
 %%