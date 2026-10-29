%% SCRIPT DE INICIALIZACIÓN: DETECCIÓN DE FALLAS Y RECONFIGURACIÓN (DC-DC BOOST)
% Basado en el artículo: "Sensor Fault Detection and System Reconfiguration 
% for DC-DC Boost Converter" (Sensors, 2018)

clear; clc; close all;

%% 1. Parámetros Físicos del Convertidor Boost (Sección 6)
Vin   = 190;        % Voltaje de entrada [V]
C0    = 0.00285;    % Capacitancia de salida [F]
L     = 0.005;      % Inductancia [H]
R     = 0.082;      % Resistencia parásita del inductor [Ohms]
Rload = 100;        % Resistencia de carga [Ohms]

%% 2. Punto de Equilibrio Deseado y Condiciones Iniciales
iL_ref = 7.6;                   % Corriente deseada del inductor i_L* [A]
vc_ref = 380;                   % Voltaje deseado del capacitor v_c* [V]
x_star = [iL_ref; vc_ref];      % Vector de estados deseado x*

% Condiciones iniciales en estado estacionario (t = 0 s)
x0     = [iL_ref; vc_ref];      % Estado inicial de la planta x(0)
xhat0  = [iL_ref; vc_ref];      % Estado inicial del observador x_hat(0)

%% 3. Matrices del Modelo en Espacio de Estados (Sección 2)
% Modo 1: Interruptor abierto (S = 0)
A0 = [-R/L,   -1/L; 
       1/C0,  -1/(Rload*C0)];

% Modo 2: Interruptor cerrado (S = 1)
A1 = [-R/L,    0; 
       0,     -1/(Rload*C0)];

% Matriz de entrada (B), matriz de salida (C) y acoplamiento de falla (Q)
B  = [1/L; 0];
C  = eye(2);
Q  = eye(2);

%% 4. Parámetros del Observador de Luenberger y Detección (Secciones 3, 5 y 6)
% Constantes de diseño LMI
a11 = 299;
a12 = 152;

% Matriz de ganancia del observador H (obtenida de la LMI del Teorema 1)
H = [36.4677, -17.0967; 
     28.5219,  75.6359];

% Umbral de detección de falla
J_th = 0.02;

%% 5. Parámetros de la Ley de Control Conmutado por Lyapunov (Sección 4)
% Teorema 2 exige: p11 > 0, p22 > 0 y (p11 / L) = (p22 / C0)
p11 = L;
p22 = C0;
P_lyap = [p11, 0; 
          0,   p22];

%% 6. Tiempos de Simulación y Configuración de Fallas (Sección 6)
t_fault = 0.2;      % Instante en el que ocurre la falla [s]
t_sim   = 0.6;      % Tiempo total de simulación [s]

disp('¡Parámetros cargados correctamente en el Workspace!');