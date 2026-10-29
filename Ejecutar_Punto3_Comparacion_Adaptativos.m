%% SCRIPT PUNTO 3 CORREGIDO: COMPARACIÓN 100% SIMULADA DINÁMICAMENTE
clear; clc;

Vin = 190; C0 = 0.00285; L = 0.005; R = 0.082;
x_star = [7.6; 380]; P_lyap = [L, 0; 0, C0];
A0 = [-R/L, -1/L; 1/C0, -1/(100*C0)]; A1 = [-R/L, 0; 0, -1/(100*C0)]; B = [1/L; 0];
H = [36.4677, -17.0967; 28.5219, 75.6359];
dt = 2e-6; t_vec = (0:dt:0.6)'; num_pasos = length(t_vec);
J_min = 0.02; N_vent = 10000; alpha = 2/(N_vent+1);

%% FUNCIÓN INTERNA PARA SIMULAR CUALQUIER ESCENARIO Y MEDIR DETECCIÓN
% Escenarios: 0=Solo Carga, 1=10% t=0.2, 2=5% t=0.2, 3=2% t=0.2, 4=5% t=0.36, 5=2% t=0.36
run_sim = @(escenario) simular_metodo(escenario, dt, num_pasos, t_vec, Vin, C0, L, R, x_star, P_lyap, A0, A1, B, H, J_min, N_vent, alpha);

[R_Carga] = run_sim(0);
[R_10_20] = run_sim(1);
[R_05_20] = run_sim(2);
[R_02_20] = run_sim(3);
[R_05_36] = run_sim(4);
[R_02_36] = run_sim(5);

% Recuperación tras carga (t > 0.35)
t_rec_din = t_vec(find(t_vec > 0.35 & R_Carga.J_din <= 0.0205, 1, 'first'));
t_rec_ewma = t_vec(find(t_vec > 0.35 & R_Carga.J_ewm <= 0.0205, 1, 'first'));
t_rec_mad = t_vec(find(t_vec > 0.35 & R_Carga.J_mad <= 0.0205, 1, 'first'));

Metricas = {'1. Falsas Alarmas en Carga'; '2. Duracion FA en Carga (ms)'; '3. J_max en Carga'; '4. T. Retorno a 0.02 (s)'; ...
    '5. Td Falla 10% en t=0.2s (ms)'; '6. Td Falla 5% en t=0.2s (ms)'; '7. Td Falla 2% en t=0.2s (ms)'; ...
    '8. Detecta Falla 5% en t=0.36s'; '9. Detecta Falla 2% en t=0.36s'};

fmt_td = @(x) sprintf('%.4f', x);
fmt_det = @(x) repmat(x > 0, 1, 1);

C_Fij = {num2str(sum(R_Carga.F_fij)), num2str(sum(R_Carga.F_fij)*dt*1e3), num2str(max(R_Carga.J_fij)), 'N/A', ...
    fmt_td(R_10_20.Td_fij), fmt_td(R_05_20.Td_fij), fmt_td(R_02_20.Td_fij), num2str(R_05_36.Td_fij>=0), num2str(R_02_36.Td_fij>=0)}';

C_Din = {num2str(sum(R_Carga.F_din)), num2str(sum(R_Carga.F_din)*dt*1e3), num2str(max(R_Carga.J_din)), num2str(t_rec_din), ...
    fmt_td(R_10_20.Td_din), fmt_td(R_05_20.Td_din), fmt_td(R_02_20.Td_din), num2str(R_05_36.Td_din>=0), num2str(R_02_36.Td_din>=0)}';

C_Ewm = {num2str(sum(R_Carga.F_ewm)), num2str(sum(R_Carga.F_ewm)*dt*1e3), num2str(max(R_Carga.J_ewm)), num2str(t_rec_ewma), ...
    fmt_td(R_10_20.Td_ewm), fmt_td(R_05_20.Td_ewm), fmt_td(R_02_20.Td_ewm), num2str(R_05_36.Td_ewm>=0), num2str(R_02_36.Td_ewm>=0)}';

C_Mad = {num2str(sum(R_Carga.F_mad)), num2str(sum(R_Carga.F_mad)*dt*1e3), num2str(max(R_Carga.J_mad)), num2str(t_rec_mad), ...
    fmt_td(R_10_20.Td_mad), fmt_td(R_05_20.Td_mad), fmt_td(R_02_20.Td_mad), num2str(R_05_36.Td_mad>=0), num2str(R_02_36.Td_mad>=0)}';

writetable(table(Metricas, C_Fij, C_Din, C_Ewm, C_Mad), 'Tabla_Punto3_100porc_Dinamica.csv');
disp('Punto 3 calculado con simulación física real de todos los umbrales.');

%% FUNCIÓN DE SIMULACIÓN MULTI-UMBRAL
function Res = simular_metodo(escenario, dt, num_pasos, t_vec, Vin, C0, L, R, x_star, P_lyap, A0_n, A1_n, B, H, J_min, N_vent, alpha)
    x = x_star; x_h = x_star; xc = x_star;
    F_fij = 0; F_din = 0; F_ewm = 0; F_mad = 0;
    
    b_din = zeros(N_vent,1); idin=1; cdin=0; s1=0; s2=0;
    mu_e = 0; var_e = 0;
    b_mad = zeros(N_vent,1); imad=1; cmad=0; m_last=0; mad_last=0;
    
    Res.J_fij = zeros(num_pasos,1); Res.J_din = zeros(num_pasos,1);
    Res.J_ewm = zeros(num_pasos,1); Res.J_mad = zeros(num_pasos,1);
    Res.F_fij = zeros(num_pasos,1); Res.F_din = zeros(num_pasos,1);
    Res.F_ewm = zeros(num_pasos,1); Res.F_mad = zeros(num_pasos,1);
    
    idx_f = -1;
    if escenario==1||escenario==2||escenario==3; idx_f = find(t_vec >= 0.2, 1); end
    if escenario==4||escenario==5; idx_f = find(t_vec >= 0.36, 1); end
    
    for m = 1:num_pasos
        t = t_vec(m);
        if (escenario==0 || escenario==4 || escenario==5) && t >= 0.3
            tau = t - 0.3; pul = (tau<=0.040)*(tau/0.040) + (tau>0.040)*exp(-(tau-0.040)/0.025);
            Rl = 100/(1+0.035*pul);
        else
            Rl = 100;
        end
        A0_p = [-R/L, -1/L; 1/C0, -1/(Rl*C0)]; A1_p = [-R/L, 0; 0, -1/(Rl*C0)];
        z = xc - x_star; if (2*z'*P_lyap*(A0_n*xc+B*Vin)) < (2*z'*P_lyap*(A1_n*xc+B*Vin)); As = A0_p; Aso=A0_n; else; As = A1_p; Aso=A1_n; end
        
        y = x;
        if m >= idx_f
            if escenario==1; y(1) = 1.1*x(1); end
            if escenario==2||escenario==4; y(1) = 1.05*x(1); end
            if escenario==3||escenario==5; y(1) = 1.02*x(1); end
        end
        r1 = (y(1)-x_h(1))^2 + (y(2)-x_h(2))^2;
        
        % Fijo
        Jf = J_min; Ff = r1 > Jf;
        % Dinamico
        if cdin>0; md = s1/cdin; sd = sqrt(max((s2/cdin)-md^2,0)); else; md=0; sd=0; end
        Jd = max(J_min, md + 4*sd); Fd = r1 > Jd;
        if Fd==0; if cdin==N_vent; old=b_din(idin); s1=s1-old; s2=s2-old^2; else; cdin=cdin+1; end; b_din(idin)=r1; s1=s1+r1; s2=s2+r1^2; idin=mod(idin,N_vent)+1; end
        % EWMA
        Je = max(J_min, mu_e + 4*sqrt(max(var_e,0))); Fe = r1 > Je;
        de = r1 - mu_e; mu_e = (1-alpha)*mu_e + alpha*r1; var_e = (1-alpha)*(var_e + alpha*de^2);
        % MAD
        if cmad>0 && s1>1e-10; if mod(m,25)==1; vact=b_mad(1:cmad); m_last=median(vact); mad_last=1.4826*median(abs(vact-m_last)); end; else; m_last=0; mad_last=0; end
        Jm = max(J_min, m_last + 4*mad_last); Fm = r1 > Jm;
        if cmad<N_vent; cmad=cmad+1; end; b_mad(imad)=r1; imad=mod(imad,N_vent)+1;
        
        Res.J_fij(m)=Jf; Res.J_din(m)=Jd; Res.J_ewm(m)=Je; Res.J_mad(m)=Jm;
        Res.F_fij(m)=Ff; Res.F_din(m)=Fd; Res.F_ewm(m)=Fe; Res.F_mad(m)=Fm;
        
        % Integración (planta limpia, solo evaluamos detección, sin realimentar para no alterar el transitorio base de la carga)
        k1_x = As*x + B*Vin; k2_x = As*(x + 0.5*dt*k1_x) + B*Vin; k3_x = As*(x + 0.5*dt*k2_x) + B*Vin; k4_x = As*(x + dt*k3_x) + B*Vin; x = x + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x);
        k1_xh = Aso*x_h + B*Vin + H*(y-x_h); k2_xh = Aso*(x_h + 0.5*dt*k1_xh) + B*Vin + H*(y-(x_h+0.5*dt*k1_xh)); k3_xh = Aso*(x_h + 0.5*dt*k2_xh) + B*Vin + H*(y-(x_h+0.5*dt*k2_xh)); k4_xh = Aso*(x_h + dt*k3_xh) + B*Vin + H*(y-(x_h+dt*k3_xh)); x_h = x_h + (dt/6)*(k1_xh + 2*k2_xh + 2*k3_xh + k4_xh);
        xc = y;
    end
    
    if idx_f > 0
        df = find(Res.F_fij(idx_f:end),1,'first'); if isempty(df); Res.Td_fij=-1; else; Res.Td_fij = (df-1)*dt*1e3; end
        dd = find(Res.F_din(idx_f:end),1,'first'); if isempty(dd); Res.Td_din=-1; else; Res.Td_din = (dd-1)*dt*1e3; end
        de = find(Res.F_ewm(idx_f:end),1,'first'); if isempty(de); Res.Td_ewm=-1; else; Res.Td_ewm = (de-1)*dt*1e3; end
        dm = find(Res.F_mad(idx_f:end),1,'first'); if isempty(dm); Res.Td_mad=-1; else; Res.Td_mad = (dm-1)*dt*1e3; end
    end
end