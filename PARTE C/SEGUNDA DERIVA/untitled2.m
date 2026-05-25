clc; clear; close all;

%% DEFINICIÓN DE LOS DATOS ORIGINALES
datos = [
 100, 152.3; 120, 149.1; 145, 146.8; 170, 144.9; 200, 142.0;
 235, 139.5; 270, 137.9; 310, 136.1; 355, 134.8; 405, 133.6;
 460, 132.7; 520, 131.9; 585, 131.4; 655, 131.1; 730, 130.9;
 810, 131.0; 895, 131.3; 985, 131.9; 1080, 132.7; 1180, 133.8;
 1290, 135.2; 1410, 136.9; 1540, 138.9; 1680, 141.1; 1830, 143.5;
 1990, 146.1; 2160, 149.0; 2340, 152.2; 2530, 155.6; 2730, 159.2
];

f = datos(:,1);
Z = datos(:,2);
n = length(f);

%% CONSTRUCCIÓN DEL SPLINE CÚBICO NATURAL
h = diff(f);
a = Z;
A = zeros(n, n); B = zeros(n, 1);
A(1, 1) = 1; A(n, n) = 1; 
for i = 2:n-1
    A(i, i-1) = h(i-1);
    A(i, i)   = 2 * (h(i-1) + h(i));
    A(i, i+1) = h(i);
    B(i)      = 3 * ((a(i+1) - a(i)) / h(i) - (a(i) - a(i-1)) / h(i-1));
end
c = A \ B; 
b = zeros(n-1, 1); d = zeros(n-1, 1);
for i = 1:n-1
    b(i) = (a(i+1) - a(i))/h(i) - h(i)*(2*c(i) + c(i+1))/3;
    d(i) = (c(i+1) - c(i))/(3*h(i));
end

%% 3. CREACIÓN DE LAS ESTRUCTURAS DE DERIVADAS ANALÍTICAS (PP)
% Primera Derivada: S'_i = 3*d*(f-f_i)^2 + 2*c*(f-f_i) + b
coef_derivada1 = [3*d, 2*c(1:end-1), b];
derivada1_pp = mkpp(f, coef_derivada1);

% Segunda Derivada: S''_i = 6*d*(f-f_i) + 2*c
coef_derivada2 = [6*d, 2*c(1:end-1)];
derivada2_pp = mkpp(f, coef_derivada2);

%% 4. UBICACIÓN DEL MÍNIMO Y CÁLCULO DE LA SEGUNDA DERIVADA

malla_frecuencias = linspace(min(f), max(f), 200000);
d1_fina = ppval(derivada1_pp, malla_frecuencias);

% Localizar el cruce exacto por cero (- a +)
idx_cambio = find(d1_fina(1:end-1) < 0 & d1_fina(2:end) > 0, 1);
f_minimo = malla_frecuencias(idx_cambio);

% Evaluar numéricamente la segunda derivada en ese punto exacto
d2_en_minimo = ppval(derivada2_pp, f_minimo);


%% 5. DESPLIEGUE DE RESULTADOS Y DISCUSIÓN ANALÍTICA

fprintf('=========================================================\n');
fprintf(' ANÁLISIS DE ESTABILIDAD EN EL MÍNIMO (CRITERIO DE S2)\n');
fprintf('=========================================================\n');
fprintf('Frecuencia del Mínimo (f_min):       %.3f Hz\n', f_minimo);
fprintf('Valor de d^2|Z|/df^2 en el mínimo:   %.6f\n', d2_en_minimo);
fprintf('---------------------------------------------------------\n');
fprintf('DIAGNÓSTICO MATEMÁTICO:\n');

if d2_en_minimo > 0
    fprintf('Signo: POSITIVO (+)\n');
    fprintf('Conclusión: La curva es estrictamente CÓNCAVA HACIA ARRIBA.\n');
    fprintf('El punto crítico f = %.3f Hz es un MÍNIMO LOCAL ESTABLE.\n', f_minimo);
elseif d2_en_minimo < 0
    fprintf('Signo: NEGATIVO (-)\n');
    fprintf('Conclusión: La curva es CÓNCAVA HACIA ABAJO (Máximo).\n');
else
    fprintf('Signo: CERO (0)\n');
    fprintf('Conclusión: Punto de inflexión o silla (Inconcluso).\n');
end
fprintf('=========================================================\n\n');


%% 6. GRAFICACIÓN DE LA SEGUNDA DERIVADA
d2_fina = ppval(derivada2_pp, malla_frecuencias);

figure('Name', 'Segunda Derivada e Inflexión', 'NumberTitle', 'off');
plot(malla_frecuencias, d2_fina, 'b-', 'LineWidth', 2, 'DisplayName', 'd^2|Z|/df^2 Analítica');
hold on;
plot(f_minimo, d2_en_minimo, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9, 'DisplayName', 'Valor en f\_min');
line([min(f) max(f)], [0 0], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off'); % Eje cero

title('Comportamiento de la Segunda Derivada d^2|Z|/df^2');
xlabel('Frecuencia f (Hz)');
ylabel('d^2|Z|/df^2');
grid on;
legend('Location', 'best');
xlim([min(f), max(f)]);