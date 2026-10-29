%% SCRIPT PUNTO 1: VALIDACIÓN ESTOCÁSTICA DE RUIDO (30 CORRIDAS CON SEMILLAS 1..30)
% Responde a las Observaciones 3, 12, 13 y 14 del revisor.
% Evalúa 30 corridas independientes con rng(seed, 'twister') y n(t) = 1.35 * randn:
%   - Modalidad A: Sin amplitud mínima (gaussiano puro)
%   - Modalidad B: Con amplitud mínima |n(t)| >= 0.2
%   - Sensores: Corriente (i_L) y Tensión (v_C)
%   - Umbrales: Fijo (J_th = 0.02) y Dinámico (N = 10000, k = 4, J_min = 0.02)
% Exporta tablas en CSV y gráfica profesional en fondo blanco.

clear; clc; close all;

%% 1. Parámetros Físicos y Numéricos (Idénticos al Modelo de Simulink)
Vin   = 190;
C0    = 0.00285;
L     = 0.005;
R     = 0.082;
Rload = 100;

x_star = [7.6; 380];
P_lyap = [L, 0; 0, C0];

A0 = [-R/L, -1/L; 1/C0, -1/(Rload*C0)];
A1 = [-R/L,  0;   0,    -1/(Rload*C0)];
B  = [1/L; 0];

H = [36.4677, -17.0967; 
     28.5219,  75.6359];

dt        = 2e-6;               % Paso fijo ode4: 2 us (0.002 ms)
t_falla   = 0.2;                % Instante de inyección de falla [s]
t_vec     = (0:dt:0.6)';        % Vector de tiempo de 0 a 0.6 s
num_pasos = length(t_vec);
idx_falla = find(abs(t_vec - t_falla) < 1e-12, 1, 'first'); % Paso 100001
pasos_pos = num_pasos - idx_falla + 1;                      % 200001 muestras con falla

% Parámetros del Umbral Dinámico (Anexo A)
N_ventana  = 10000;             % Ventana móvil: 20 ms
k_sens     = 4;                 % Factor de sensibilidad k = 4
J_min      = 0.02;              % Cota inferior / Umbral fijo = 0.02
sigma_ruido = 1.35;             % Desviación estándar del ruido gaussiano
N_semillas  = 30;               % Número de corridas independientes (seeds 1..30)

%% 2. Simulación Determinista Pre-Falla (de t = 0 a t = 0.2 s)
% Como antes de t = 0.2 s no hay falla (r1 = 0), este tramo es idéntico para todas las corridas.
x          = x_star;
x_hat      = x_star;
x_ctrl_mem = x_star;

for m_step = 1:(idx_falla - 1)
    z = x_ctrl_mem - x_star;
    V0_dot = 2 * z' * P_lyap * (A0 * x_ctrl_mem + B * Vin);
    V1_dot = 2 * z' * P_lyap * (A1 * x_ctrl_mem + B * Vin);
    if V0_dot < V1_dot
        A_s = A0;
    else
        A_s = A1;
    end
    y_med = x;
    
    % Integración ode4 Planta y Observador
    k1_x = A_s * x + B * Vin;
    k2_x = A_s * (x + 0.5*dt*k1_x) + B * Vin;
    k3_x = A_s * (x + 0.5*dt*k2_x) + B * Vin;
    k4_x = A_s * (x + dt*k3_x) + B * Vin;
    x_next = x + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x);
    
    k1_xh = A_s*x_hat + B*Vin + H*(y_med - x_hat);
    k2_xh = A_s*(x_hat + 0.5*dt*k1_xh) + B*Vin + H*(y_med - (x_hat + 0.5*dt*k1_xh));
    k3_xh = A_s*(x_hat + 0.5*dt*k2_xh) + B*Vin + H*(y_med - (x_hat + 0.5*dt*k2_xh));
    k4_xh = A_s*(x_hat + dt*k3_xh) + B*Vin + H*(y_med - (x_hat + dt*k3_xh));
    x_hat_next = x_hat + (dt/6)*(k1_xh + 2*k2_xh + 2*k3_xh + k4_xh);
    
    x          = x_next;
    x_hat      = x_hat_next;
    x_ctrl_mem = y_med;
end

% Guardar estado exacto en t = 0.2 s para arrancar las corridas desde ahí
x_02          = x;
x_hat_02      = x_hat;
x_ctrl_mem_02 = x_ctrl_mem;

%% 3. Bucle Principal: 30 Semillas x 2 Modalidades x 2 Sensores x 2 Umbrales
% Dimensiones: (30 semillas, 2 sensores [1=Corr,2=Tens], 2 modalidades [1=SinMin,2=ConMin], 2 umbrales [1=Fijo,2=Din])
r1_det_mat     = zeros(N_semillas, 2, 2, 2);
detectado_mat  = zeros(N_semillas, 2, 2, 2);
Td_alarma_mat  = zeros(N_semillas, 2, 2, 2); % Retardo en puerto F [ms]
Td_reconf_mat  = zeros(N_semillas, 2, 2, 2); % Retardo hasta reconfiguración [ms]
Jth_max_mat    = zeros(N_semillas, 2, 2, 2); % Máximo de Jth en [0, 0.6] s

% Guardar curvas de Jth_din para graficar (submuestreadas cada 20 pasos)
idx_plot_pos = 1:20:pasos_pos;
t_plot_pos   = t_vec(idx_falla - 1 + idx_plot_pos);
Jth_curvas   = zeros(length(idx_plot_pos), N_semillas, 2, 2); % (tiempo, seed, sensor, modalidad)

fprintf('Ejecutando las 30 corridas estocásticas (seeds 1..30)... Por favor espera unos segundos.\n\n');

for seed = 1:N_semillas
    % Fijar semilla reproducible con generador Mersenne Twister
    rng(seed, 'twister');
    ruido_base = sigma_ruido * randn(pasos_pos, 1);
    
    % Guardar vector de la Semilla 1 en el Workspace para el Caso de Control en Simulink
    if seed == 1
        ruido_seed1_sin_min = ruido_base;
        ruido_seed1_con_min = ruido_base;
        idx_peq = abs(ruido_seed1_con_min) < 0.2;
        signs_p = sign(ruido_seed1_con_min(idx_peq));
        signs_p(signs_p == 0) = 1;
        ruido_seed1_con_min(idx_peq) = 0.2 * signs_p;
    end
    
    for mod_amp = 1:2 % 1 = Sin amplitud mínima | 2 = Con amplitud mínima 0.2
        if mod_amp == 1
            n_vec = ruido_base;
        else
            n_vec = ruido_base;
            idx_bajo = abs(n_vec) < 0.2;
            s_sign   = sign(n_vec(idx_bajo));
            s_sign(s_sign == 0) = 1;
            n_vec(idx_bajo) = 0.2 * s_sign;
        end
        
        for sensor = 1:2 % 1 = Sensor Corriente | 2 = Sensor Tensión
            for tipo_umb = 1:2 % 1 = Umbral Fijo | 2 = Umbral Dinámico
                % Restaurar estado exacto en t = 0.2 s
                x          = x_02;
                x_hat      = x_hat_02;
                x_ctrl_mem = x_ctrl_mem_02;
                F_mem      = 0;
                
                % En t = 0.2 s el buffer ya tiene N_ventana = 10000 muestras nulas (r1 = 0)
                buf    = zeros(N_ventana, 1);
                idx    = 1;
                cnt    = N_ventana;
                sum_r1 = 0;
                sum_r2 = 0;
                
                det_flag   = 0;
                rec_flag   = 0;
                r1_al_det  = NaN;
                td_al_ms   = NaN;
                td_rec_ms  = NaN;
                max_Jth    = J_min;
                
                if tipo_umb == 2
                    Jth_temp_vec = zeros(pasos_pos, 1);
                end
                
                for k_step = 1:pasos_pos
                    % 1. Ley de control conmutado
                    z = x_ctrl_mem - x_star;
                    V0_dot = 2 * z' * P_lyap * (A0 * x_ctrl_mem + B * Vin);
                    V1_dot = 2 * z' * P_lyap * (A1 * x_ctrl_mem + B * Vin);
                    if V0_dot < V1_dot
                        A_s = A0;
                    else
                        A_s = A1;
                    end
                    
                    % 2. Inyección de ruido gaussiano en el sensor correspondiente
                    y_med = x;
                    if sensor == 1
                        y_med(1) = x(1) + n_vec(k_step);
                    else
                        y_med(2) = x(2) + n_vec(k_step);
                    end
                    
                    % 3. Cálculo del residuo cuadrático r1(t)
                    r  = y_med - x_hat;
                    r1 = r(1)^2 + r(2)^2;
                    
                    % 4. Evaluación del umbral (Fijo vs. Dinámico)
                    if tipo_umb == 1
                        Jth = J_min;
                        F   = double(r1 > Jth);
                    else
                        media = sum_r1 / cnt;
                        desv  = sqrt(max((sum_r2 / cnt) - media^2, 0));
                        Jth   = max(J_min, media + k_sens * desv);
                        F     = double(r1 > Jth);
                        
                        if F == 0
                            old_val = buf(idx);
                            sum_r1  = sum_r1 - old_val + r1;
                            sum_r2  = sum_r2 - old_val^2 + r1^2;
                            buf(idx) = r1;
                            idx      = mod(idx, N_ventana) + 1;
                        end
                        if Jth > max_Jth
                            max_Jth = Jth;
                        end
                        Jth_temp_vec(k_step) = Jth;
                    end
                    
                    % Registrar primer instante de detección (F == 1)
                    if det_flag == 0 && F > 0.5
                        det_flag  = 1;
                        r1_al_det = r1;
                        td_al_ms  = (k_step - 1) * dt * 1e3;
                    end
                    % Registrar primer instante de reconfiguración efectiva (F_mem == 1)
                    if rec_flag == 0 && F_mem > 0.5
                        rec_flag  = 1;
                        td_rec_ms = (k_step - 1) * dt * 1e3;
                    end
                    
                    % Si es umbral fijo y ya reconfiguró, podemos cortar la corrida para acelerar
                    if tipo_umb == 1 && rec_flag == 1
                        break;
                    end
                    
                    % 5. Selección de señales de control según F_mem
                    if F_mem > 0.5
                        x_ctrl = x_hat;
                        y_corr = x_star;
                    else
                        x_ctrl = y_med;
                        y_corr = y_med;
                    end
                    
                    % 6. Integración ode4
                    k1_x = A_s * x + B * Vin;
                    k2_x = A_s * (x + 0.5*dt*k1_x) + B * Vin;
                    k3_x = A_s * (x + 0.5*dt*k2_x) + B * Vin;
                    k4_x = A_s * (x + dt*k3_x) + B * Vin;
                    x_next = x + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x);
                    
                    k1_xh = A_s*x_hat + B*Vin + H*(y_corr - x_hat);
                    k2_xh = A_s*(x_hat + 0.5*dt*k1_xh) + B*Vin + H*(y_corr - (x_hat + 0.5*dt*k1_xh));
                    k3_xh = A_s*(x_hat + 0.5*dt*k2_xh) + B*Vin + H*(y_corr - (x_hat + 0.5*dt*k2_xh));
                    k4_xh = A_s*(x_hat + dt*k3_xh) + B*Vin + H*(y_corr - (x_hat + dt*k3_xh));
                    x_hat_next = x_hat + (dt/6)*(k1_xh + 2*k2_xh + 2*k3_xh + k4_xh);
                    
                    x          = x_next;
                    x_hat      = x_hat_next;
                    x_ctrl_mem = x_ctrl;
                    F_mem      = F;
                end
                
                % Guardar métricas de esta corrida
                r1_det_mat(seed, sensor, mod_amp, tipo_umb)    = r1_al_det;
                detectado_mat(seed, sensor, mod_amp, tipo_umb) = det_flag;
                Td_alarma_mat(seed, sensor, mod_amp, tipo_umb) = td_al_ms;
                Td_reconf_mat(seed, sensor, mod_amp, tipo_umb) = td_rec_ms;
                Jth_max_mat(seed, sensor, mod_amp, tipo_umb)   = max_Jth;
                
                if tipo_umb == 2
                    Jth_curvas(:, seed, sensor, mod_amp) = Jth_temp_vec(idx_plot_pos);
                end
            end
        end
    end
end

%% 4. Cálculo de Estadísticos (Media, DE, IC 95% con t_0.975,29 = 2.0452)
t_crit = 2.0452; % Valor crítico t-Student para 95% IC con gl = 29

Escenario_Col = {};
Umbral_Col    = {};
Prob_Det_Porc = [];
Tasa_FA_Porc  = [];
r1_Media      = [];
r1_DE         = [];
r1_IC95_Inf   = [];
r1_IC95_Sup   = [];
Td_Alarma_Med = [];
Td_Alarma_DE  = [];
Td_Reconf_Med = [];
Td_Reconf_DE  = [];
Jth_Max_Med   = [];
Jth_Max_DE    = [];
Jth_Max_Abs   = [];

etiquetas_mod  = {'Sin_Min_0.0', 'Con_Min_0.2'};
etiquetas_sens = {'Corriente', 'Tension'};
etiquetas_umb  = {'Fijo', 'Dinamico'};

for mod_amp = 1:2
    for sensor = 1:2
        for tipo_umb = 1:2
            v_r1   = r1_det_mat(:, sensor, mod_amp, tipo_umb);
            v_det  = detectado_mat(:, sensor, mod_amp, tipo_umb);
            v_tda  = Td_alarma_mat(:, sensor, mod_amp, tipo_umb);
            v_tdr  = Td_reconf_mat(:, sensor, mod_amp, tipo_umb);
            v_jmax = Jth_max_mat(:, sensor, mod_amp, tipo_umb);
            
            m_r1   = mean(v_r1);
            sd_r1  = std(v_r1);
            err_r1 = t_crit * (sd_r1 / sqrt(N_semillas));
            
            Escenario_Col{end+1, 1} = sprintf('Ruido_%s_%s', etiquetas_sens{sensor}, etiquetas_mod{mod_amp});
            Umbral_Col{end+1, 1}    = etiquetas_umb{tipo_umb};
            Prob_Det_Porc(end+1, 1) = mean(v_det) * 100;
            Tasa_FA_Porc(end+1, 1)  = 0.0; % r1 = 0 antes de t = 0.2 s
            r1_Media(end+1, 1)      = m_r1;
            r1_DE(end+1, 1)         = sd_r1;
            r1_IC95_Inf(end+1, 1)   = m_r1 - err_r1;
            r1_IC95_Sup(end+1, 1)   = m_r1 + err_r1;
            Td_Alarma_Med(end+1, 1) = mean(v_tda);
            Td_Alarma_DE(end+1, 1)  = std(v_tda);
            Td_Reconf_Med(end+1, 1) = mean(v_tdr);
            Td_Reconf_DE(end+1, 1)  = std(v_tdr);
            Jth_Max_Med(end+1, 1)   = mean(v_jmax);
            Jth_Max_DE(end+1, 1)    = std(v_jmax);
            Jth_Max_Abs(end+1, 1)   = max(v_jmax);
        end
    end
end

%% 5. Exportar Tablas CSV e Imprimir Resumen Compacto en Command Window
T_Resumen = table(Escenario_Col, Umbral_Col, Prob_Det_Porc, Tasa_FA_Porc, ...
    r1_Media, r1_DE, r1_IC95_Inf, r1_IC95_Sup, ...
    Td_Alarma_Med, Td_Alarma_DE, Td_Reconf_Med, Td_Reconf_DE, ...
    Jth_Max_Med, Jth_Max_DE, Jth_Max_Abs);

writetable(T_Resumen, 'Tabla_Punto1_Resumen_Estadistico_Ruido.csv');

% Tabla detallada de las 30 semillas (Umbral Dinámico)
Semilla = (1:N_semillas)';
r1_Corr_SinMin   = r1_det_mat(:, 1, 1, 2);
TdRec_Corr_SinMin= Td_reconf_mat(:, 1, 1, 2);
Jmax_Corr_SinMin = Jth_max_mat(:, 1, 1, 2);

r1_Tens_SinMin   = r1_det_mat(:, 2, 1, 2);
TdRec_Tens_SinMin= Td_reconf_mat(:, 2, 1, 2);
Jmax_Tens_SinMin = Jth_max_mat(:, 2, 1, 2);

r1_Corr_ConMin   = r1_det_mat(:, 1, 2, 2);
TdRec_Corr_ConMin= Td_reconf_mat(:, 1, 2, 2);
Jmax_Corr_ConMin = Jth_max_mat(:, 1, 2, 2);

r1_Tens_ConMin   = r1_det_mat(:, 2, 2, 2);
TdRec_Tens_ConMin= Td_reconf_mat(:, 2, 2, 2);
Jmax_Tens_ConMin = Jth_max_mat(:, 2, 2, 2);

T_Detalle = table(Semilla, ...
    r1_Corr_SinMin, TdRec_Corr_SinMin, Jmax_Corr_SinMin, ...
    r1_Tens_SinMin, TdRec_Tens_SinMin, Jmax_Tens_SinMin, ...
    r1_Corr_ConMin, TdRec_Corr_ConMin, Jmax_Corr_ConMin, ...
    r1_Tens_ConMin, TdRec_Tens_ConMin, Jmax_Tens_ConMin);

writetable(T_Detalle, 'Tabla_Punto1_Detalle_30Semillas.csv');

fprintf('--- RESUMEN ESTADÍSTICO DE 30 CORRIDAS (SEEDS 1..30, UMBRAL DINÁMICO) ---\n');
fprintf('Escenario                  | Pd(%%) | r1_det (Media +- DE) | IC 95%% r1_det    | Td_rec (ms)       | Jth_max (Media +- DE)\n');
fprintf('--------------------------------------------------------------------------------------------------------------------------\n');
for idx_row = 2:2:height(T_Resumen)
    fprintf('%-26s | %5.1f | %7.4f +- %7.4f | [%6.4f, %6.4f] | %.4f +- %.4f | %.4f +- %.4f (Max: %.4f)\n', ...
        T_Resumen.Escenario_Col{idx_row}, T_Resumen.Prob_Det_Porc(idx_row), ...
        T_Resumen.r1_Media(idx_row), T_Resumen.r1_DE(idx_row), ...
        T_Resumen.r1_IC95_Inf(idx_row), T_Resumen.r1_IC95_Sup(idx_row), ...
        T_Resumen.Td_Reconf_Med(idx_row), T_Resumen.Td_Reconf_DE(idx_row), ...
        T_Resumen.Jth_Max_Med(idx_row), T_Resumen.Jth_Max_DE(idx_row), T_Resumen.Jth_Max_Abs(idx_row));
end
fprintf('--------------------------------------------------------------------------------------------------------------------------\n');
fprintf('VALOR CASO DE CONTROL (Semilla 1, Sin Min): r1_det = %.4f | (Semilla 1, Con Min 0.2): r1_det = %.4f\n\n', ...
    r1_Corr_SinMin(1), r1_Corr_ConMin(1));

%% 6. Generar Gráfica Profesional con Fondo Blanco (30 trayectorias de Jth_din)
figure('Name', 'fig_punto1_ruido_30corridas', 'Color', 'w', 'Position', [80, 80, 1200, 620]);
titulos_sub = {
    '(a) J_{th,din}(t) en 30 corridas - Sensor Corriente (Sin mín.)', ...
    '(b) J_{th,din}(t) en 30 corridas - Sensor Tensión (Sin mín.)', ...
    '(c) J_{th,din}(t) en 30 corridas - Sensor Corriente (Con mín. 0.2)', ...
    '(d) J_{th,din}(t) en 30 corridas - Sensor Tensión (Con mín. 0.2)'
};

sub_idx = 1;
for mod_amp = 1:2
    for sensor = 1:2
        subplot(2, 2, sub_idx);
        hold on;
        % Graficar las 30 semillas en azul claro y el promedio en rojo intenso
        curvas_mat = Jth_curvas(:, :, sensor, mod_amp);
        plot([0; 0.2; t_plot_pos], [J_min; J_min; curvas_mat(:, 1)], '-', 'Color', [0.65 0.80 0.95], 'LineWidth', 0.8);
        for s_i = 2:N_semillas
            plot(t_plot_pos, curvas_mat(:, s_i), '-', 'Color', [0.65 0.80 0.95], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        media_curva = mean(curvas_mat, 2);
        plot([0; 0.2; t_plot_pos], [J_min; J_min; media_curva], '-', 'Color', [0.85 0.15 0.10], 'LineWidth', 1.8);
        yline(0.5776, 'k--', 'LineWidth', 1.2); % Residuo de la falla determinista más pequeña (+10% ganancia i_L)
        
        ax = gca;
        ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
        box on; grid on;
        xlim([0, 0.6]);
        ylim([0, 0.65]);
        xlabel('Tiempo (s)', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
        ylabel('Umbral J_{th,din}(t) [A^2+V^2]', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
        title(titulos_sub{sub_idx}, 'FontSize', 9.5, 'FontWeight', 'bold', 'Color', 'k');
        if sub_idx == 1
            legend('30 corridas (seeds 1..30)', 'Media de las 30 corridas', 'Residuo falla mín. (0.5776)', ...
                'Location', 'northwest', 'FontSize', 8, 'Color', 'w', 'TextColor', 'k');
        end
        sub_idx = sub_idx + 1;
    end
end

exportgraphics(gcf, 'fig_punto1_ruido_30corridas.png', 'Resolution', 300);
disp('¡Tablas CSV y gráfica fig_punto1_ruido_30corridas.png generadas con éxito!');