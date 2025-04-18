clear all; close all; clc;
load("full_1000_3months.mat")
vector = zeros(1,n_debris);
for u = 1:n_debris
    time2 = TimeObserved{u};
    if ~isempty(TimeObserved{u})
        vector(u) = TimeObserved{u}(end);
    else
        vector(u) = 0;
    end
end
% Eliminar los ceros (objetos no detectados)
detected_times = vector(vector > 0);

% Ordenar los tiempos de detección
sorted_times = sort(detected_times);

% Calcular porcentaje acumulado de objetos detectados
cumulative_percentage = arrayfun(@(t) sum(sorted_times <= t) / sum(visited>0) * 100, sorted_times);

f = fit(sorted_times',cumulative_percentage', 'a*log(x)+b', 'StartPoint', [0,0] );

t_target = exp((80-f.b)/f.a);
p_nuevo = f.a*log(5*365*24*3600)+f.b;

figure;
plot(sorted_times, cumulative_percentage, 'o', 'MarkerSize', 8, 'LineWidth', 2); % Datos
hold on;
plot(f, sorted_times, cumulative_percentage); % Curva ajustada
xlabel('Tiempo (días)');
ylabel('Porcentaje (%)');
title('Ajuste Logarítmico de los Datos');
legend('Datos', 'Curva Ajustada');
grid on;

t_original = sorted_times;
p_original = cumulative_percentage;
log_func = @(params, t) params(1) * log(t + 1e-6) + params(2);

params = nlinfit(t_original, p_original, log_func, [1, 1]);  % Ajuste no lineal
a0 = params(1);  % Parámetro 'a' ajustado de la curva original
b0 = params(2);  % Parámetro 'b' ajustado de la curva original

% Datos del punto de ajuste
x_punto = 6*30*3600*24;  % El punto de tiempo donde queremos que pase la curva
y_punto = 13.85;  % El porcentaje en ese tiempo

% Ajuste logarítmico con los puntos (0, 0) y (x, y)
epsilon = 1e-6;  % Valor pequeño para evitar log(0)

% Calcular el parámetro a ajustado para que pase por (x_punto, y_punto)
a = y_punto / log(x_punto + epsilon);

% Ecuación logarítmica ajustada
log_func = @(t) a * log(t + epsilon);

% Graficar los datos y la curva ajustada
figure;
plot(t_original, p_original, 'o', 'MarkerSize', 8, 'LineWidth', 2); % Datos originales
hold on;

% Graficar la curva logarítmica ajustada
t_fine = linspace(min(t_original), max(t_original), 100);  % Más puntos para suavizar la curva
p_fine = log_func(t_fine);  % Curva ajustada
plot(t_fine, p_fine, '-r', 'LineWidth', 2); % Curva ajustada
xlabel('Tiempo (días)');
ylabel('Porcentaje (%)');
title('Ajuste Logarítmico con los puntos [0, 0] y [x, y]');
legend('Datos', 'Curva Ajustada');
grid on;









% Datos del punto de ajuste
x_punto = 6*30*3600*24;  % El punto de tiempo donde queremos que pase la curva
y_punto = 15.39;  % El porcentaje en ese tiempo

%x_punto = 2*24*3600;
%y_punto = 9.68;

% Función logarítmica que queremos ajustar a los datos
% La función es de la forma: f(t) = a * log(t + epsilon) + b
log_func = @(params, t) params(1) * log(t + 1e-6) + params(2);  % a*log(t) + b

% Estimación inicial de los parámetros a y b
initial_params = [1, 1];  % Estimación inicial para [a, b]

% Ajuste de los parámetros usando nlinfit
params = nlinfit(t_original, p_original, log_func, initial_params);

a_ajustado = (y_punto-params(2)) / log(x_punto + 1e-6);

% Calculamos el parámetro 'b' usando la curva ajustada
b_ajustado = params(2);  % Usamos el valor de 'b' obtenido del ajuste original

% Definir la nueva función logarítmica con los parámetros ajustados
log_func_ajustada = @(t) a_ajustado * log(t + 1e-6) + b_ajustado;


% Graficar los datos y la curva ajustada
figure;
plot(t_original, p_original, 'o', 'MarkerSize', 8, 'LineWidth', 2); % Datos originales
hold on;

% Graficar la curva logarítmica ajustada
t_fine = linspace(min(t_original), max(t_original), 100);  % Más puntos para suavizar la curva
p_fine = log_func_ajustada(t_fine);  % Curva ajustada
plot(t_fine, p_fine, '-r', 'LineWidth', 2); % Curva ajustada
xlabel('Tiempo (días)');
ylabel('Porcentaje (%)');
title('Ajuste Logarítmico a los Datos');
legend('Datos', 'Curva Ajustada');
grid on;


porcentaje_nuevo = log_func(params,5*365*24*3600)


