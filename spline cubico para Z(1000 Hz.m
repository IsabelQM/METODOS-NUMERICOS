clc; clear; close all;

%%  DEFINICIÓN DE LOS DATOS ORIGINALES
datos = [
 100, 152.3; 120, 149.1; 145, 146.8; 170, 144.9; 200, 142.0;
 235, 139.5; 270, 137.9; 310, 136.1; 355, 134.8; 405, 133.6;
 460, 132.7; 520, 131.9; 585, 131.4; 655, 131.1; 730, 130.9;
 810, 131.0; 895, 131.3; 985, 131.9; 1080, 132.7; 1180, 133.8;
 1290, 135.2; 1410, 136.9; 1540, 138.9; 1680, 141.1; 1830, 143.5;
 1990, 146.1; 2160, 149.0; 2340, 152.2; 2530, 155.6; 2730, 159.2
];

x = datos(:,1);
y = datos(:,2);
n = length(x);

f_target = 1000; % Frecuencia objetivo para evaluar


%% CONSTRUCCIÓN Y EVALUACIÓN DEL SPLINE CÚBICO NATURAL

h = diff(x);
a = y;
A = zeros(n, n);
B = zeros(n, 1);

A(1, 1) = 1; 
A(n, n) = 1; 

for i = 2:n-1
    A(i, i-1) = h(i-1);
    A(i, i)   = 2 * (h(i-1) + h(i));
    A(i, i+1) = h(i);
    B(i)      = 3 * ((a(i+1) - a(i)) / h(i) - (a(i) - a(i-1)) / h(i-1));
end

c = A \ B; 
b = zeros(n-1, 1);
d = zeros(n-1, 1);
for i = 1:n-1
    b(i) = (a(i+1) - a(i))/h(i) - h(i)*(2*c(i) + c(i+1))/3;
    d(i) = (c(i+1) - c(i))/(3*h(i));
end

poly_coefficients = [d, c(1:end-1), b, a(1:end-1)];
spline_natural_pp = mkpp(x, poly_coefficients);

% Calcular valor exacto en f = 1000 Hz con el Spline
val_spline_1000 = ppval(spline_natural_pp, f_target);


%% CÁLCULO DEL POLINOMIO SELECCIONADO (Grado 2 con los 3 puntos más cercanos)

% Buscamos los 3 puntos originales más cercanos a 1000 Hz
[~, idx_ordenado] = sort(abs(x - f_target));
idx_3 = sort(idx_ordenado(1:3));

% Ajustamos el polinomio cuadrático local y evaluamos en f = 1000 Hz
p_local = polyfit(x(idx_3), y(idx_3), 2);
val_pol_1000 = polyval(p_local, f_target);


%% IMPRESIÓN DE RESULTADOS EN CONSOLA

fprintf('=========================================================\n');
fprintf(' CÁLCULO E INTERPOLACIÓN EN f = %d Hz\n', f_target);
fprintf('=========================================================\n');
fprintf('Valor usando Spline Cúbico Natural:       %.4f\n', val_spline_1000);
fprintf('Valor usando Polinomio Local (Grado 2):   %.4f\n', val_pol_1000);
fprintf('---------------------------------------------------------\n');
fprintf('Diferencia Absoluta entre modelos:        %.4f\n', abs(val_spline_1000 - val_pol_1000));
fprintf('Diferencia Porcentual relativa:           %.3f %%\n', (abs(val_spline_1000 - val_pol_1000)/val_spline_1000)*100);
fprintf('=========================================================\n\n');



%% EVALUACIÓN EN MALLA FINA Y GRAFICACIÓN

x_fina = linspace(min(x), max(x), 2000);
y_spline_fina = ppval(spline_natural_pp, x_fina);

% Para graficar el polinomio de grado 2 local de forma continua en todo el rango:
y_pol_fina = zeros(size(x_fina));
for k = 1:length(x_fina)
    [~, idx_temp] = sort(abs(x - x_fina(k)));
    idx_temp_3 = sort(idx_temp(1:3));
    p_temp = polyfit(x(idx_temp_3), y(idx_temp_3), 2);
    y_pol_fina(k) = polyval(p_temp, x_fina(k));
end

% Generación de la gráfica
figure('Name', 'Comparativa de Modelos Estables', 'NumberTitle', 'off');
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Datos Originales'); hold on;
plot(x_fina, y_spline_fina, 'b-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico Natural');
plot(x_fina, y_pol_fina, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Polinomio Local Seleccionado (Grado 2)');
plot(f_target, val_spline_1000, 'kx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Evaluación en 1000 Hz');
hold off;

title(sprintf('Interpolación de Impedancia: Evaluación en f = %d Hz', f_target));
xlabel('Frecuencia (Hz)');
ylabel('Impedancia |Z|');
grid on;
legend('Location', 'best');
ylim([min(y)-2, max(y)+2]); % Enfocar en la zona de interés real