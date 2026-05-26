% Procesamiento de datos de telemetría biomédica
clear; clc; close all;

%  DATOS 

f = [10.0; 12.5; 15.0; 17.5; 20.0; 22.5; 25.0; 27.5; 30.0; 32.5; ...
     35.0; 37.5; 40.0; 42.5; 45.0; 47.5; 50.0; 52.5; 55.0; 57.5; ...
     60.0; 62.5; 65.0; 67.5; 70.0; 72.5; 75.0; 77.5; 80.0; 82.5; ...
     85.0; 87.5; 90.0; 92.5; 95.0; 97.5; 100.0; 102.5; 105.0; 107.5];

V = [0.842; 0.911; 0.986; 1.062; 1.143; 1.227; 1.314; 1.401; 1.482; 1.551; ...
     1.216; 1.048; 0.866; 0.689; 0.521; 0.364; 0.223; 0.103; 0.012; -0.041; ...
     -0.057; -0.034; 0.018; 0.096; 0.197; 0.318; 0.452; 0.579; 0.700; ...
     0.809; 0.611; 0.688; 0.756; 0.811; 0.856; 0.894; 0.926; 0.954; 0.980; 1.004];

Z = [182.4; 178.9; 175.1; 171.0; 166.8; 162.7; 158.9; 155.4; 152.0; 149.0; ...
     146.1; 145.2; 145.8; 147.3; 149.9; 153.5; 158.0; 163.2; 168.9; 174.8; ...
     180.5; 186.2; 191.5; 196.2; 200.1; 203.1; 205.2; 206.3; 206.1; ...
     204.7; 198.0; 194.4; 190.9; 187.8; 185.1; 183.0; 181.6; 180.8; 180.6; 180.9];

%% PARTE 1: INTERPOLACIÓN BÁSICA Y SPLINE
fprintf('--- PARTE 1: INTERPOLACIÓN ---\n');

% Función anidada para Lagrange de 2do grado (3 puntos más cercanos)
lagrange2 = @(f_eval, x, y) ...
    y(1)*((f_eval-x(2))*(f_eval-x(3)))/((x(1)-x(2))*(x(1)-x(3))) + ...
    y(2)*((f_eval-x(1))*(f_eval-x(3)))/((x(2)-x(1))*(x(2)-x(3))) + ...
    y(3)*((f_eval-x(1))*(f_eval-x(2)))/((x(3)-x(1))*(x(3)-x(2)));

% --- Evaluación para f = 41.0 kHz ---
idx41 = [12, 13, 14]; % Índices para f = 37.5, 40.0, 42.5
V_Lagrange_41 = lagrange2(41.0, f(idx41), V(idx41));
Z_Lagrange_41 = lagrange2(41.0, f(idx41), Z(idx41));

% --- Evaluación para f = 73.0 kHz ---
idx73 = [25, 26, 27]; % Índices para f = 70.0, 72.5, 75.0
V_Lagrange_73 = lagrange2(73.0, f(idx73), V(idx73));
Z_Lagrange_73 = lagrange2(73.0, f(idx73), Z(idx73));

% --- Spline Cúbico Natural ---
% Nota: csape requiere el Curve Fitting Toolbox para 'variational' (natural).
% Si no se tiene, 'spline' de MATLAB usa condiciones not-a-knot.
try
    spline_V = csape(f, V, 'variational');
    spline_Z = csape(f, Z, 'variational');
    V_Spline_41 = fnval(spline_V, 41.0);
    Z_Spline_41 = fnval(spline_Z, 41.0);
    V_Spline_73 = fnval(spline_V, 73.0);
    Z_Spline_73 = fnval(spline_Z, 73.0);
catch
    
    V_Spline_41 = spline(f, V, 41.0);
    Z_Spline_41 = spline(f, Z, 41.0);
    V_Spline_73 = spline(f, V, 73.0);
    Z_Spline_73 = spline(f, Z, 73.0);
end

fprintf('V(41.0)  -> Lagrange: %.4f V | Spline: %.4f V\n', V_Lagrange_41, V_Spline_41);
fprintf('|Z|(41.0) -> Lagrange: %.3f Ohm | Spline: %.3f Ohm\n', Z_Lagrange_41, Z_Spline_41);
fprintf('V(73.0)  -> Lagrange: %.4f V | Spline: %.4f V\n', V_Lagrange_73, V_Spline_73);
fprintf('|Z|(73.0) -> Lagrange: %.3f Ohm | Spline: %.3f Ohm\n\n', Z_Lagrange_73, Z_Spline_73);


%  PARTE 2: DERIVACIÓN NUMÉRICA

fprintf('--- PARTE 2: DERIVACIÓN NUMÉRICA ---\n');
h = 2.5; % Espaciamiento constante

% --- Frecuencia f = 40.0 kHz (Índice 13) ---
i = 13;
dVdf_C2_40 = (V(i+1) - V(i-1)) / (2*h);
dVdf_C4_40 = (-V(i+2) + 8*V(i+1) - 8*V(i-1) + V(i-2)) / (12*h);

% --- Frecuencia f = 70.0 kHz (Índice 25) ---
i = 25;
dVdf_C2_70 = (V(i+1) - V(i-1)) / (2*h);
dVdf_C4_70 = (-V(i+2) + 8*V(i+1) - 8*V(i-1) + V(i-2)) / (12*h);

% --- Frecuencia f = 100.0 kHz (Índice 37) ---
i = 37;
dVdf_C2_100 = (V(i+1) - V(i-1)) / (2*h);
dVdf_C4_100 = (-V(i+2) + 8*V(i+1) - 8*V(i-1) + V(i-2)) / (12*h);

% --- Extremo inferior f = 10.0 kHz (Progresiva Orden 2) ---
dVdf_Prog_10 = (-3*V(1) + 4*V(2) - V(3)) / (2*h);

fprintf('f = 40.0 kHz  -> Centr. O(2): %.5f | Centr. O(4): %.5f V/kHz\n', dVdf_C2_40, dVdf_C4_40);
fprintf('f = 70.0 kHz  -> Centr. O(2): %.5f | Centr. O(4): %.5f V/kHz\n', dVdf_C2_70, dVdf_C4_70);
fprintf('f = 100.0 kHz -> Centr. O(2): %.5f | Centr. O(4): %.5f V/kHz\n', dVdf_C2_100, dVdf_C4_100);
fprintf('f = 10.0 kHz  -> Progresiva O(2): %.5f V/kHz\n\n', dVdf_Prog_10);


%  PARTE 3: RAÍCES POR BISECCIÓN Y SPLINE

fprintf('--- PARTE 3: LOCALIZACIÓN DE RAÍCES ---\n');

% Función genérica de bisección para este set de datos interpolados linealmente
biseccion_tabla = @(a, b, V_array, f_array, iters) ...
    ejecutar_biseccion(a, b, V_array, f_array, iters);

fprintf('Primer cruce por cero (bisección 3 iteraciones):\n');
raiz1 = biseccion_tabla(55.0, 57.5, V, f, 3);

fprintf('Segundo cruce por cero (bisección 3 iteraciones):\n');
raiz2 = biseccion_tabla(62.5, 65.0, V, f, 3);


%  VISUALIZACIÓN GRÁFICA (Opcional - Recomendada)

figure('Name', 'Caracterización del Front-End Analógico');
subplot(2,1,1);
plot(f, V, 'ro ', 'MarkerFaceColor', 'r'); hold on;
f_smooth = linspace(10, 107.5, 500);
V_smooth = spline(f, V, f_smooth);
plot(f_smooth, V_smooth, 'b-', 'LineWidth', 1.5);
grid on; yline(0, 'k--');
title('Voltaje de Salida V(f)'); xlabel('Frecuencia (kHz)'); ylabel('Voltaje (V)');

subplot(2,1,2);
plot(f, Z, 'mo ', 'MarkerFaceColor', 'm'); hold on;
Z_smooth = spline(f, Z, f_smooth);
plot(f_smooth, Z_smooth, 'g-', 'LineWidth', 1.5);
grid on;
title('Impedancia Equivalente |Z(f)|'); xlabel('Frecuencia (kHz)'); ylabel('Impedancia (\Omega)');

%  FUNCIÓN AUXILIAR PARA BISECCIÓN

function raiz = ejecutar_biseccion(a, b, V_arr, f_arr, max_iter)
    % Buscamos los índices correspondientes en el dominio para evaluar
    interp_V = @(x) interp1(f_arr, V_arr, x, 'linear'); 
    
    for k = 1:max_iter
        m = (a + b) / 2;
        Vm = interp_V(m);
        Va = interp_V(a);
        
        fprintf('  Iter %d: Intervalo = [%.4f, %.4f] -> Pto Medio = %.4f | V(m) = %.5f\n', ...
            k, a, b, m, Vm);
        
        if Va * Vm < 0
            b = m;
        else
            a = m;
        end
    end
    raiz = (a + b) / 2;
    fprintf('  -> Raíz aproximada final: %.4f kHz\n\n', raiz);
end