% 1. DEFINICION DE DATOS
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
N = length(x); % Número de puntos (30)
grado = N - 1; % Grado del polinomio (29)

x_eval = linspace(min(x), max(x), 1000);


% MÉTODO 1: MÉTODO MATRICIAL (Matriz de Vandermonde)
V = zeros(N, N);
for i = 1:N
    V(:, i) = x.^(N - i); 
end

coef_matricial = V \ y;

% Evaluamos el polinomio matricial en los puntos de interés
y_matricial = polyval(coef_matricial, x_eval);


% METODO 2: INTERPOLACION DE LAGRANGE
y_lagrange = zeros(size(x_eval));

for k = 1:length(x_eval)
    sumatoria = 0;
    for i = 1:N
        % Calcular el polinomio productor L_i(x)
        L_i = 1;
        for j = 1:N
            if i ~= j
                L_i = L_i * (x_eval(k) - x(j)) / (x(i) - x(j));
            end
        end
        sumatoria = sumatoria + y(i) * L_i;
    end
    y_lagrange(k) = sumatoria;
end


% VISUALIZACION Y GRAFICA
figure('Name', 'Interpolación Polinómica de Grado 29', 'NumberTitle', 'off');
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Datos Originales');
hold on;
plot(x_eval, y_matricial, 'b-', 'LineWidth', 2, 'DisplayName', 'Método Matricial (Vandermonde)');
plot(x_eval, y_lagrange, 'g--', 'LineWidth', 2, 'DisplayName', 'Método de Lagrange');

title(sprintf('Interpolación Polinómica Global (Grado %d)', grado));
xlabel('Eje X');
ylabel('Eje Y');
grid on;
legend('Location', 'best');

% Ajustar límites verticales porque el fenómeno de Runge disparará los extremos
ylim([min(y)-10, max(y)+10]); 
hold off;

% Mostrar mensaje 
fprintf('Procesamiento completado para %d puntos (Polinomio de grado %d).\n', N, grado);
if cond(V) > 1e12
    fprintf('¡ADVERTENCIA!: El número de condición de la matriz de Vandermonde es %.2e.\n', cond(V));
    fprintf('Esto causa alta inestabilidad numérica en el Método Matricial.\n');
end
