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
f_data = datos(:,1); Z_data = datos(:,2); n = length(f_data);
Z_th = 150; % Umbral crítico

%% CONSTRUCCIÓN DEL SPLINE CÚBICO NATURAL Y SU DERIVADA
h = diff(f_data); a = Z_data; A = zeros(n, n); B = zeros(n, 1);
A(1,1)=1; A(n,n)=1;
for i=2:n-1
    A(i,i-1)=h(i-1); A(i,i)=2*(h(i-1)+h(i)); A(i,i+1)=h(i);
    B(i)=3*((a(i+1)-a(i))/h(i) - (a(i)-a(i-1))/h(i-1));
end
c = A\B; b = zeros(n-1,1); d = zeros(n-1,1);
for i=1:n-1
    b(i)=(a(i+1)-a(i))/h(i) - h(i)*(2*c(i)+c(i+1))/3;
    d(i)=(c(i+1)-c(i))/(3*h(i));
end
spline_pp = mkpp(f_data, [d, c(1:end-1), b, a(1:end-1)]);
d1_pp = mkpp(f_data, [3*d, 2*c(1:end-1), b]); % Primera derivada analítica

% Definición de funciones auxiliares usando estructuras de evaluación nativas
F_objetivo = @(freq) ppval(spline_pp, freq) - Z_th;
Derivada_F = @(freq) ppval(d1_pp, freq);

tol = 1e-5; % Tolerancia para garantizar al menos 4 cifras significativas


%% MÉTODO DE BISECCIÓN

% Raíz 1 (Intervalo 100 - 150)
[r1_bis, iter1_bis] = biseccion(F_objetivo, 100, 150, tol);
% Raíz 2 (Intervalo 2100 - 2300)
[r2_bis, iter2_bis] = biseccion(F_objetivo, 2100, 2300, tol);


%% 4. MÉTODO DE NEWTON-RAPHSON

% Raíz 1 (Aproximación inicial f0 = 110)
[r1_nr, iter1_nr] = newton_raphson(F_objetivo, Derivada_F, 110, tol);
% Raíz 2 (Aproximación inicial f0 = 2200)
[r2_nr, iter2_nr] = newton_raphson(F_objetivo, Derivada_F, 2200, tol);


%% ANÁLISIS DE SENSIBILIDAD EN LA SEGUNDA RAÍZ (f ≈ 2218 Hz)

dZ_df = ppval(d1_pp, r2_nr); % Tasa de cambio instantánea d|Z|/df
df_dZ = 1 / dZ_df;           % Sensibilidad inversa analítica


%% DESPLIEGUE DE REPORTES EN CONSOLA

fprintf('===================================================================\n');
fprintf('        COMPARATIVA DE CONVERGENCIA Y RESULTADOS DE RAÍCES\n');
fprintf('===================================================================\n');
fprintf('MÉTODO DE BISECCIÓN:\n');
fprintf('  -> Raíz 1 (Límite Inferior): %.4f Hz (%d iteraciones)\n', r1_bis, iter1_bis);
fprintf('  -> Raíz 2 (Límite Superior): %.4f Hz (%d iteraciones)\n\n', r2_bis, iter2_bis);
fprintf('MÉTODO DE NEWTON-RAPHSON:\n');
fprintf('  -> Raíz 1 (Límite Inferior): %.4f Hz (%d iteraciones)\n', r1_nr, iter1_nr);
fprintf('  -> Raíz 2 (Límite Superior): %.4f Hz (%d iteraciones)\n', r2_nr, iter2_nr);
fprintf('-------------------------------------------------------------------\n');
fprintf('ANÁLISIS DE SENSIBILIDAD EN f = %.2f Hz:\n', r2_nr);
fprintf('  -> Derivada d|Z|/df:          %.5f Ohms/Hz\n', dZ_df);
fprintf('  -> Sensibilidad df/d|Z|:      %.4f Hz/Ohm\n', df_dZ);
fprintf('  -> Impacto: Un error de +1 Ohm desplaza la raíz por %.2f Hz.\n', df_dZ);
fprintf('===================================================================\n');

%% VISUALIZACIÓN GRÁFICA DE LA BANDA SEGURA
malla = linspace(min(f_data), max(f_data), 2000);
Z_malla = ppval(spline_pp, malla);

figure('Name', 'Caracterización de Banda Segura', 'NumberTitle', 'off');
plot(f_data, Z_data, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Mediciones'); hold on;
plot(malla, Z_malla, 'b-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico |Z|(f)');
line([min(f_data) max(f_data)], [Z_th Z_th], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Umbral Z_{th} = 150 \Omega');

% Destacar raíces encontradas
plot([r1_nr, r2_nr], [Z_th, Z_th], 'gX', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Límites de Banda');

title('Espectro de Operación Segura del Módulo de Transmisión');
xlabel('Frecuencia f (Hz)'); ylabel('Impedancia |Z| (\Omega)');
grid on; legend('Location', 'best');
xlim([min(f_data), max(f_data)]); ylim([min(Z_malla)-5, max(Z_malla)+5]);
hold off;

%% FUNCIONES INTERNAS DE LOS ALGORITMOS NUMÉRICOS
function [raiz, iter] = biseccion(F, a, b, tol)
    iter = 0;
    while (b - a)/2 > tol
        iter = iter + 1;
        c = (a + b)/2;
        if F(c) == 0
            break;
        elseif F(a)*F(c) < 0
            b = c;
        else
            a = c;
        end
    end
    raiz = (a + b)/2;
end

function [raiz, iter] = newton_raphson(F, dF, x0, tol)
    iter = 0; error_est = 1;
    while error_est > tol
        iter = iter + 1;
        x1 = x0 - F(x0)/dF(x0);
        error_est = abs(x1 - x0);
        x0 = x1;
    end
    raiz = x1;
end