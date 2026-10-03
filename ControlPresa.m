clc
clear all
close all
%% Tarea 1: Recoleccion de datos
datos = load('resultados.txt'); 

t = datos(:,1);      
u_ref = datos(:,2);  
x_up = datos(:,3); 
y_dwe = datos(:,4);

Ts = 0.03; 

figure(1);
plot(t, y_dwe, 'LineWidth', 1.2);
title('Respuesta Temporal del Nivel');
xlabel('Tiempo (s)'); 
ylabel('Nivel (mm)');
grid on;

%% Normalización 

norm_u_ref = u_ref - u_ref(1);
norm_x_up = x_up - x_up(1);
norm_y_dwe = y_dwe - y_dwe(1);

%%  Filtrado

y_filt = medfilt1(norm_y_dwe, 25); 
x_filt = medfilt1(norm_x_up, 15);

idx = find(t >= 2.58, 1);
t = t(idx:end) - t(idx); 
du = norm_u_ref(idx:end);
dx = x_filt(idx:end);
dy = y_filt(idx:end);


%%  Análisis Temporal

figure(2);
plot(t, dx, 'LineWidth', 1.2);
title('Respuesta Temporal del Motor');
xlabel('Tiempo (s)'); 
ylabel('Nivel (mm)');
grid on;

figure(3);
plot(t, dy, 'LineWidth', 1.2);
title('Respuesta Temporal del Nivel filtrado y normalizado');
xlabel('Tiempo (s)'); 
ylabel('Nivel (mm)');
grid on;


%%  Datos para identificación

M = iddata(dx, du, Ts);   
G = iddata(dy, dx, Ts);  

%%  Análisis de Frecuencia

m_frec = spa(M);   
g_frec = spa(G);   

% Bode de M(s)
figure(4);
bode(m_frec);
title('Respuesta en Frecuencia M(s) = X_{up} / U_{ref}');
grid on;

% Bode de G(s)
figure(5);
bode(g_frec);
title('Respuesta en Frecuencia G(s) = Y_{dwe} / X_{up}');
grid on;

%% Tarea 2: Identificacion Parametros para PRBS
Ts = 0.030;

omega_f = 10;                   
f = omega_f / (2*pi);           
h_teorico = 1 / (3*f);            
h = max(0.5, h_teorico);          


Tf = 55;                          
N_orden = ceil(log2(Tf/h + 1));   
Mn = 2^N_orden - 1;                

A = 5;                           
n = 2;                            
samples_per_step = round(h / Ts); 

u_prbs_base = idinput(Mn, 'prbs'); 

u_prbs_base = A * u_prbs_base;


u_prbs_expanded = repelem(u_prbs_base, samples_per_step); 
u_prbs = repmat(u_prbs_expanded, n, 1);                   
t_prbs = (0:length(u_prbs)-1)' * Ts;

figure(6);
plot(t_prbs, u_prbs, 'LineWidth', 1.2);
title(sprintf('Señal PRBS Diseñada (h = %.1fs, N = %d, M = %d)', h, N_orden, Mn));
xlabel('Tiempo (s)');
ylabel('Amplitud (mm)');
xlim([0 t_prbs(end)]);
grid on;

%% Tarea 4: Creacion del modelo
% Identificación del Motor M(s) 
np_M = 1; nz_M = 0; 
M_s = tfest(M, np_M, nz_M);


figure(7);
compare(M, M_s);
title('Validación del Motor M(s)');


% Identificación del Canal G(s) 
np_G = 1; nz_G = 0; 
G_s_simple = tfest(G, 1, 0,4.26);

figure(8);
compare(G, G_s_simple);
title('Validación del Canal simple G(s)');
% Identificación del Canal G(z) usando BJ
OBJ = Inf;
orden = [nan nan nan];
for nb = 0 : 4
    for nd = 0 : 4
        for nf = 0 : 2
            sys = bj(G, [nb 0 nd nf 3]);
            obj = sys.Report.Fit.nAIC;
            if obj < OBJ
                OBJ = obj;
                orden = [nb nd nf];
            end
            disp([int2str(nb) ' x ' int2str(nd) ' x ' int2str(nf) ': ' num2str(obj) ' - ' num2str(sys.Report.Fit.FitPercent)])
        end
    end
end

nk_muestras = round(4.26 / Ts); 
ordenes_bj = [3 , 2, 2, 2, nk_muestras]; 
G_z_bj = bj(G, ordenes_bj);
% Paso a modelo continuo
t_sim = (0:Ts:120)';
u_sim = ones(size(t_sim));
G_z_limpio = tf(G_z_bj, 'measured'); 

y_sim_discreta = lsim(G_z_limpio, u_sim, t_sim);
datos_para_s = iddata(y_sim_discreta, u_sim, Ts);
G_s = tfest(datos_para_s, 2, 1, 4.26);

figure(9);
compare(G, G_s);
title('Validación del Canal G(s)');

%% Tarea 5: Analisis dinamico
retardo = 4.26;
[numG,denG] = tfdata(G_s,'v');
[numM,denM] = tfdata(M_s,'v');

[num_pade, den_pade] = pade(retardo, 2);
Pade_tf = tf(num_pade, den_pade);

G_r = tf(numG,denG)*Pade_tf;
[numGr,denGr] = tfdata(G_r,'v');


%Lugar de las raices
La = M_s * G_r; 

figure(10);
rlocus(La);
title('Lugar de las Raíces');
grid on;

%% Tarea 6: PID
G_total = M_s*G_r;

[C_pid, info] = pidtune(G_total, 'PID');

Kp_pid = C_pid.Kp
Ki_pid = C_pid.Ki
Kd_pid = C_pid.Kd

%% Predictor de Smith
polos_G = pole(G_s);
polo_lento = max(polos_G); 
ts_deseado = 20;
tau_cl = ts_deseado / 4; 

G_sin_retardo = tf (G_s.Numerator,G_s.Denominator)
G_total_sin_retardo = M_s * G_sin_retardo;
K_planta = dcgain(G_total_sin_retardo * tf([1 -polo_lento], 1));


Kp_smith = 1 / (K_planta * tau_cl);
Ki_smith = -polo_lento * Kp_smith; 

%Estudio de robustez
retardo_mas50 = 4.26 * 1.5;
[numPade_mas50,denPade_mas50] = pade(retardo_mas50, 2);
Pade_mas50 = tf(numPade_mas50,denPade_mas50);

M_s_mas50 = tf(numM * 1.5, [denM(1), denM(2) * 1.5]);
G_s_mas50 = tf(numG * 1.5, [denG(1), denG(2) * 1.5, denG(3) * 1.5]);

G_total_mas50 =  G_s_mas50 * Pade_mas50;
[numG_mas50, denG_mas50] = tfdata(G_total_mas50, 'v');
[numM_mas50, denM_mas50] = tfdata(M_s_mas50, 'v');

retardo_menos50 = 4.26 * 0.5;
[numPade_menos50,denPade_menos50] = pade(retardo_menos50, 2);
Pade_menos50 = tf(numPade_menos50,denPade_menos50);

M_s_menos50 = tf(numM * 0.5, [denM(1), denM(2) * 0.5]);
G_s_menos50 = tf(numG * 0.5, [denG(1), denG(2) * 0.5, denG(3) * 0.5]);

G_total_menos50 =  G_s_menos50 * Pade_menos50;
[numG_menos50, denG_menos50] = tfdata(G_total_menos50, 'v');
[numM_menos50, denM_menos50] = tfdata(M_s_menos50, 'v');

%% Tarea 7: Control predictivo
s = tf('s');

gp = M_s * G_sin_retardo;
gp.outputdelay = 4.26; 

ftz = c2d(gp, Ts, 'zoh'); 

% Simulación en Espacio de Estados Discreto
sys_ss = ss(ftz);
A_ss = sys_ss.A; B_ss = sys_ss.B; C_ss = sys_ss.C; D_ss = sys_ss.D;
x_planta = zeros(size(A_ss, 1), 1); 

d = ftz.outputdelay; 
P = 350; 
N = 1;  
delta = 3 * eye(P); 
lamda = 15;         

gi = step(ftz, (P + d + 10)*Ts); 
Nm = length(gi) - 1;

% Construcción de la matriz G del DMC
G = zeros(P, N);
G(:,1) = gi(1+d : P+d);
for i = 2:N
    for j = 2:P
        G(j,i) = G(j-1,i-1);
    end
end

% Ganancia óptima del controlador DMC
Mn = inv(G'*delta*G + lamda*eye(N)) * G'*delta;
K1 = Mn(1, :);

nit = 1200; 
u = zeros(1, nit); 
ym = zeros(1, nit); 
r = zeros(1, nit);

% Generación de la rampa 
pendiente = 2 * Ts; 
valor_consigna = 0;
inicio_rampa = 60; 

for k = 1:nit
    if k >= inicio_rampa
        if valor_consigna < 10
            valor_consigna = valor_consigna + pendiente;
            if valor_consigna > 10, valor_consigna = 10; end
        end
        r(k) = valor_consigna;
    else
        r(k) = 0; 
    end
end

duf = zeros(1, Nm);

% Bucle de control
for k = 2:nit
    
    
    ym(k) = C_ss * x_planta + D_ss * u(k-1);
    x_planta = A_ss * x_planta + B_ss * u(k-1);
    
    % Cálculo de la predicción libre f
    f = zeros(1, P);
    for kk = 1:P
        vect_g = zeros(1, Nm); 
        for i = 1:Nm-P
            vect_g(i) = gi(kk+i) - gi(i);
        end
        for i = Nm-P+1:Nm
            vect_g(i) = gi(Nm) - gi(i);
        end
        f(kk) = ym(k) + vect_g*duf';
    end
    
    % Ley de control incremental DMC
    inc_u = K1 * (r(k) - f');
    
    % Apertura acumulada de la compuerta
    u(k) = u(k-1) + inc_u;
    
    % Restricciones físicas de la compuerta (0 a 60 mm)
    if u(k) > 60
        u(k) = 60;
    elseif u(k) < 0
        u(k) = 0;
    end
    
    % Cálculo del incremento real aplicado 
    inc_u_real = u(k) - u(k-1);
    
    % Actualización del histórico de movimientos pasados
    aux_2 = duf(1:Nm-1);
    duf = [inc_u_real, aux_2];
end
% Graficas
t = 0:Ts:(nit-1)*Ts; 

figure(11)
subplot(2,1,1)
plot(t, r, '--b', t, ym, '-r', 'Linewidth', 2.5)
xlabel('Tiempo (s)'); ylabel('Nivel Canal (mm)');
legend('Consigna (y_r)', 'Salida Real (y)', 'Location', 'SouthEast')
grid on; title('Respuesta Suavizada ante Rampa (DMC Optimizado)');

subplot(2,1,2)
plot(t, u, 'b', 'Linewidth', 2.5)
xlabel('Tiempo (s)'); ylabel('Apertura Compuerta (mm)');
legend('Acción de Control (u)')
grid on; title('Señal de Control sin Oscilaciones (0-60 mm)');

%% Tarea 8: Control robusto
%1. Sistema primer orden
K = dcgain(G_total_sin_retardo);     
[y, t] = step(G_total_sin_retardo);
tau = t(find(y >= 0.632 * K, 1));

s = tf('s');
Pade_p = tf([1 -retardo/2 retardo^2/12], [1 retardo/2 retardo^2/12]);
Gp = (K / (tau * s + 1)) * Pade_p;

%Pesos
A_error = 0.001; 
M_pico  = 1.5;   
w_c     = 0.05;  

W1 = (s/M_pico + w_c) / (s + w_c*A_error);
W3 = (s + w_c/10) / (s/10 + w_c*2);

[C_robust, CL, gamma] = mixsyn(Gp, W1, [], W3);

%2. Modelo incierto del 10%
K_inc   = ureal('K_inc', K, 'Percentage', 10);
tau_inc = ureal('tau_inc', tau, 'Percentage', 10);
L_inc   = ureal('L_inc', retardo, 'Percentage', 10);

Pade_incierto = tf([1 -L_inc/2 L_inc^2/12], [1 L_inc/2 L_inc^2/12]);
Gp_incierta = (K_inc / (tau_inc * s + 1)) * Pade_incierto;

Lazo_Abierto_Incierto = C_robust * Gp_incierta;
Sistema_Cerrado_Incierto = feedback(Lazo_Abierto_Incierto, 1);
Fun_Sensibilidad_Incierta = feedback(1, Lazo_Abierto_Incierto);

[stabmarg, wcu] = robuststab(Sistema_Cerrado_Incierto);

[Gm, Pm, Wcg, Wcp] = margin(C_robust * Gp);

MF_min = Pm * stabmarg.UpperBound; 
Ms_max = norm(Fun_Sensibilidad_Incierta, inf);

figure()
hold on;

[G_vertices, info] = usample(Sistema_Cerrado_Incierto, 30); 
t = 0:0.2:100; 

for i = 1:size(G_vertices, 3)
    [y, t_out] = step(G_vertices(:, :, i), t);
    plot(t_out, y, 'Color', [0.75 0.75 0.75], 'LineWidth', 1); 
end

[yp, tp] = step(feedback(C_robust * Gp, 1), t);
plot(tp, yp, 'r', 'LineWidth', 2.5);

grid on;
title('Análisis Robusto H-\infty ante Escalón - Familias de Plantas');
xlabel('Tiempo (segundos)'); 
ylabel('Nivel del Canal (mm)');
legend('Respuestas Inciertas', 'Respuesta', 'Location', 'SouthEast');

