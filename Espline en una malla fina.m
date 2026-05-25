%% DEFINICIÓN DE LOS DATOS ORIGINALES
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

% MALLA FINA DE FRECUENCIA (2000 puntos uniformes)
x_fina = linspace(min(x), max(x), 2000);


%% SPLINE CÚBICO NATURAL
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

% Evaluación en la malla fina
y_spline = ppval(spline_natural_pp, x_fina);



%% EVALUACIÓN DEL POLINOMIO SELECCIONADO 

y_pol_local = zeros(size(x_fina));
grado_modelo = 2;

for k = 1:length(x_fina)
    % Buscamos dinámicamente los 3 puntos más cercanos en la malla original
    [~, idx_ordenado] = sort(abs(x - x_fina(k)));
    idx_3 = sort(idx_ordenado(1:(grado_modelo + 1)));

    % Ajustamos el polinomio cuadrático local
    p_local = polyfit(x(idx_3), y(idx_3), grado_modelo);
    y_pol_local(k) = polyval(p_local, x_fina(k));
end


%% GRAFICACIÓN Y COMPARATIVA

figure('Name', 'Comparativa: Spline Natural vs Polinomio Local Seleccionado', 'NumberTitle', 'off', 'Position', [100, 100, 950, 600]);

% Gráfica del rango completo
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Datos Originales'); hold on;
plot(x_fina, y_spline, 'b-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico Natural');
plot(x_fina, y_pol_local, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Polinomio Local Seleccionado (Grado 2)');
hold off;

title('Evaluación en Malla Fina: Spline Cúbico Natural vs Polinomio Cuadrático Local');
xlabel('Frecuencia (Hz)');
ylabel('Impedancia |Z|');
grid on;
legend('Location', 'best');

% Imprimir comparación específica en el punto f = 1000 Hz
f_test = 1000;
val_spline_1000 = ppval(spline_natural_pp, f_test);
[~, idx_test] = sort(abs(x - f_test));
p_test = polyfit(x(idx_test(1:3)), y(idx_test(1:3)), 2);
val_pol_1000 = polyval(p_test, f_test);

fprintf('=========================================================\n');
fprintf(' COMPARATIVA NUMÉRICA EN f = %d Hz\n', f_test);
fprintf('=========================================================\n');
fprintf('Valor estimado por Spline Cúbico Natural: %.5f\n', val_spline_1000);
fprintf('Valor estimado por Polinomio Local Seleccionado: %.5f\n', val_pol_1000);
fprintf('Diferencia absoluta: %.5f\n', abs(val_spline_1000 - val_pol_1000));
fprintf('=========================================================\n');