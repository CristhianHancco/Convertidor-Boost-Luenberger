%% SCRIPT PUNTO 2 CORREGIDO: ANÁLISIS DE ENMASCARAMIENTO 100% DINÁMICO
clear; clc; close all;

Vin = 190; C0 = 0.00285; L = 0.005; R = 0.082;
x_star = [7.6; 380]; P_lyap = [L, 0; 0, C0];
A0_nom = [-R/L, -1/L; 1/C0, -1/(100*C0)];
A1_nom = [-R/L,  0;   0,    -1/(100*C0)];
B = [1/L; 0]; H = [36.4677, -17.0967; 28.5219, 75.6359];
dt = 2e-6; J_min = 0.02;

deltas_porc = [10, 5, 2]; deltas_val  = [0.10, 0.05, 0.02];
r1_teorico  = (deltas_val .* 7.6).^2;

%% EXPERIMENTOS 2a (0.2s) y 2b (0.36s) - Simulación Física Real
t_vec_ab = (0:dt:0.6)'; num_pasos_ab = length(t_vec_ab);
configs_exp = [0.20, 0, 20, 4; 0.36, 1, 20, 4; 0.36, 1, 50, 5];
nombres_cfg = {'2a_Nominal_t0.20s_20ms_k4', '2b_Carga_t0.36s_20ms_k4', '2b_Carga_t0.36s_50ms_k5'};

Col_Experimento = {}; Col_Delta_Porc = []; Col_r1_Teorico = [];
Col_r1_En_Falla = []; Col_Jth_En_Falla = []; Col_Detectado = {};
Col_Td_Alarma_ms = []; Col_Td_Reconf_ms = []; Col_Observacion = {};

r1_curvas_2a = zeros(num_pasos_ab, 3); Jth_curvas_2a = zeros(num_pasos_ab, 3);
r1_curvas_2b = zeros(num_pasos_ab, 3); Jth_curvas_2b = zeros(num_pasos_ab, 3);

for ic = 1:size(configs_exp, 1)
    t_f = configs_exp(ic, 1); con_carga = configs_exp(ic, 2);
    N_vent = round((configs_exp(ic, 3) * 1e-3) / dt); k_val = configs_exp(ic, 4);
    idx_f = find(abs(t_vec_ab - t_f) < 1e-12, 1, 'first');
    
    for id = 1:3
        delta = deltas_val(id);
        x = x_star; x_hat = x_star; x_ctrl_mem = x_star; F_mem = 0;
        buf = zeros(N_vent, 1); idx = 1; cnt = 0; sum_r1 = 0; sum_r2 = 0;
        r1_h = zeros(num_pasos_ab, 1); Jth_h = zeros(num_pasos_ab, 1); F_h = zeros(num_pasos_ab, 1); F_rec = zeros(num_pasos_ab, 1);
        
        for m_step = 1:num_pasos_ab
            t = t_vec_ab(m_step);
            if con_carga == 1 && t >= 0.3
                tau = t - 0.3;
                pulso = (tau <= 0.040) * (tau/0.040) + (tau > 0.040) * exp(-(tau-0.040)/0.025);
                Rload = 100 / (1 + 0.035 * pulso);
            else
                Rload = 100;
            end
            A0_p = [-R/L, -1/L; 1/C0, -1/(Rload*C0)]; A1_p = [-R/L, 0; 0, -1/(Rload*C0)];
            z = x_ctrl_mem - x_star;
            if (2*z'*P_lyap*(A0_nom*x_ctrl_mem + B*Vin)) < (2*z'*P_lyap*(A1_nom*x_ctrl_mem + B*Vin))
                A_sp = A0_p; A_so = A0_nom; else; A_sp = A1_p; A_so = A1_nom; end
            
            y_med = x;
            if m_step >= idx_f, y_med(1) = (1 + delta) * x(1); end
            
            r1 = (y_med(1)-x_hat(1))^2 + (y_med(2)-x_hat(2))^2;
            
            if cnt > 0
                med = sum_r1/cnt; des = sqrt(max((sum_r2/cnt)-med^2, 0));
            else
                med = 0; des = 0;
            end
            Jth = max(J_min, med + k_val*des);
            F = double(r1 > Jth);
            
            if F == 0
                if cnt == N_vent, old = buf(idx); sum_r1 = sum_r1-old; sum_r2 = sum_r2-old^2; else; cnt = cnt+1; end
                buf(idx) = r1; sum_r1 = sum_r1+r1; sum_r2 = sum_r2+r1^2; idx = mod(idx, N_vent)+1;
            end
            
            r1_h(m_step) = r1; Jth_h(m_step) = Jth; F_h(m_step) = F; F_rec(m_step) = F_mem;
            if F_mem > 0.5; x_ctrl = x_hat; y_corr = x_star; else; x_ctrl = y_med; y_corr = y_med; end
            
            k1_x = A_sp*x + B*Vin; k2_x = A_sp*(x + 0.5*dt*k1_x) + B*Vin; k3_x = A_sp*(x + 0.5*dt*k2_x) + B*Vin; k4_x = A_sp*(x + dt*k3_x) + B*Vin;
            x_next = x + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x);
            k1_xh = A_so*x_hat + B*Vin + H*(y_corr-x_hat); k2_xh = A_so*(x_hat + 0.5*dt*k1_xh) + B*Vin + H*(y_corr-(x_hat+0.5*dt*k1_xh));
            k3_xh = A_so*(x_hat + 0.5*dt*k2_xh) + B*Vin + H*(y_corr-(x_hat+0.5*dt*k2_xh)); k4_xh = A_so*(x_hat + dt*k3_xh) + B*Vin + H*(y_corr-(x_hat+dt*k3_xh));
            x_hat_next = x_hat + (dt/6)*(k1_xh + 2*k2_xh + 2*k3_xh + k4_xh);
            
            x = x_next; x_hat = x_hat_next; x_ctrl_mem = x_ctrl; F_mem = F;
        end
        
        if ic == 1; r1_curvas_2a(:, id) = r1_h; Jth_curvas_2a(:, id) = Jth_h;
        elseif ic == 2; r1_curvas_2b(:, id) = r1_h; Jth_curvas_2b(:, id) = Jth_h; end
        
        idx_det_F = find(F_h(idx_f:end) == 1, 1, 'first');
        idx_det_rec = find(F_rec(idx_f:end) == 1, 1, 'first');
        
        Col_Experimento{end+1,1} = nombres_cfg{ic}; Col_Delta_Porc(end+1,1) = deltas_porc(id);
        Col_r1_Teorico(end+1,1) = r1_teorico(id); Col_r1_En_Falla(end+1,1) = r1_h(idx_f); Col_Jth_En_Falla(end+1,1) = Jth_h(idx_f);
        if ~isempty(idx_det_F)
            td_al = (idx_det_F - 1)*dt*1e3; td_rec = (idx_det_rec - 1)*dt*1e3;
            Col_Detectado{end+1,1} = 'SI'; Col_Td_Alarma_ms(end+1,1) = td_al; Col_Td_Reconf_ms(end+1,1) = td_rec;
            Col_Observacion{end+1,1} = sprintf('Deteccion a los %.2f ms', td_al);
        else
            Col_Detectado{end+1,1} = 'NO'; Col_Td_Alarma_ms(end+1,1) = NaN; Col_Td_Reconf_ms(end+1,1) = NaN;
            Col_Observacion{end+1,1} = 'Enmascarada en todo el horizonte';
        end
    end
end
writetable(table(Col_Experimento, Col_Delta_Porc, Col_r1_Teorico, Col_r1_En_Falla, Col_Jth_En_Falla, Col_Detectado, Col_Td_Alarma_ms, Col_Td_Reconf_ms, Col_Observacion), 'Tabla_Punto2AB_Enmascaramiento.csv');

%% EXPERIMENTO 2c: Ruido prolongado (Exactamente 0.2, 0.4 y 1.0 s de duración)
t_vec_c = (0:dt:1.2)'; num_pasos_c = length(t_vec_c); idx_fc = find(abs(t_vec_c - 0.2) < 1e-12, 1, 'first');
rng(1, 'twister'); n_puro = 1.35 * randn(num_pasos_c - idx_fc + 1, 1);
n_min = n_puro; idx_b = abs(n_min) < 0.2; s_b = sign(n_min(idx_b)); s_b(s_b==0)=1; n_min(idx_b) = 0.2*s_b;
Jth_c = zeros(num_pasos_c, 4);

for caso = 1:4
    if caso <= 2; n_val = n_puro; else; n_val = n_min; end
    x = x_star; x_hat = x_star; x_ctrl_mem = x_star; F_mem = 0;
    buf = zeros(10000, 1); idx = 1; cnt = 0; sum_r1 = 0; sum_r2 = 0;
    for m = 1:num_pasos_c
        z = x_ctrl_mem - x_star;
        if (2*z'*P_lyap*(A0_nom*x_ctrl_mem + B*Vin)) < (2*z'*P_lyap*(A1_nom*x_ctrl_mem + B*Vin)); A_s = A0_nom; else; A_s = A1_nom; end
        y_med = x;
        if m >= idx_fc; if mod(caso,2)~=0; y_med(1) = x(1)+n_val(m-idx_fc+1); else; y_med(2) = x(2)+n_val(m-idx_fc+1); end; end
        r1 = (y_med(1)-x_hat(1))^2 + (y_med(2)-x_hat(2))^2;
        if cnt > 0; med = sum_r1/cnt; des = sqrt(max((sum_r2/cnt)-med^2, 0)); else; med=0; des=0; end
        Jth = max(0.02, med + 4*des); F = double(r1 > Jth);
        if F == 0
            if cnt == 10000; old = buf(idx); sum_r1 = sum_r1-old; sum_r2 = sum_r2-old^2; else; cnt = cnt+1; end
            buf(idx) = r1; sum_r1 = sum_r1+r1; sum_r2 = sum_r2+r1^2; idx = mod(idx, 10000)+1;
        end
        Jth_c(m, caso) = Jth;
        if F_mem > 0.5; x_c = x_hat; y_c = x_star; else; x_c = y_med; y_c = y_med; end
        k1_x = A_s*x + B*Vin; k2_x = A_s*(x + 0.5*dt*k1_x) + B*Vin; k3_x = A_s*(x + 0.5*dt*k2_x) + B*Vin; k4_x = A_s*(x + dt*k3_x) + B*Vin; x = x + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x);
        k1_xh = A_s*x_hat + B*Vin + H*(y_c-x_hat); k2_xh = A_s*(x_hat + 0.5*dt*k1_xh) + B*Vin + H*(y_c-(x_hat+0.5*dt*k1_xh)); k3_xh = A_s*(x_hat + 0.5*dt*k2_xh) + B*Vin + H*(y_c-(x_hat+0.5*dt*k2_xh)); k4_xh = A_s*(x_hat + dt*k3_xh) + B*Vin + H*(y_c-(x_hat+dt*k3_xh)); x_hat = x_hat + (dt/6)*(k1_xh + 2*k2_xh + 2*k3_xh + k4_xh);
        x_ctrl_mem = x_c; F_mem = F;
    end
end

tiempos_corte = [0.4, 0.6, 1.2]; duraciones = [0.2, 0.4, 1.0];
Col_Modalidad = {}; Col_TSim = []; Col_Dur = []; Col_JMax = []; Col_M10 = []; Col_M5 = []; Col_M2 = [];
etiquetas = {'Corriente_Puro','Tension_Puro','Corriente_Min0.2','Tension_Min0.2'};
for c = 1:4
    for i = 1:3
        idx_t = find(abs(t_vec_c - tiempos_corte(i)) < 1e-12, 1, 'first');
        jm = max(Jth_c(idx_fc:idx_t, c));
        Col_Modalidad{end+1,1} = etiquetas{c}; Col_TSim(end+1,1) = tiempos_corte(i); Col_Dur(end+1,1) = duraciones(i);
        Col_JMax(end+1,1) = jm; Col_M10(end+1,1) = 0.5776 - jm; Col_M5(end+1,1) = 0.1444 - jm; Col_M2(end+1,1) = 0.0231 - jm;
    end
end
writetable(table(Col_Modalidad, Col_TSim, Col_Dur, Col_JMax, Col_M10, Col_M5, Col_M2), 'Tabla_Punto2C_Margenes.csv');
disp('Tablas del Punto 2 generadas dinámicamente.');