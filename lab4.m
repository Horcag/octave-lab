clear all;
clc;
format long;

% Задача Коши: y' + xy = (1 + x)e^(-x)y^2, y(0) = 1, [0, 1]
% Точное решение y = exp(x)
f1 = @(x) exp(x); % Точное решение
f = @(x, y) (1 + x) * exp(-x) * y^2 - x * y; % Правая часть ДУ
a = 0; % Начало отрезка
b = 1; % Конец отрезка
y0 = 1; % Начальное условие

% Пункт 1: Нахождение шага интегрирования для метода Рунге-Кутты 4 порядка с точностью 10^(-4)
required_accuracy = 1e-4; % Требуемая точность
p = 4; % Порядок точности метода Рунге-Кутты
N_init = 10; % Начальное количество интервалов
N = N_init;
error_rk4 = Inf;

printf('Поиск оптимального шага для метода Рунге-Кутты 4-го порядка\n');
printf('с заданной точностью %e\n\n', required_accuracy);
printf('   N   |    h    | Погрешность | Достаточно?\n');
printf('---------------------------------------\n');

while error_rk4 > required_accuracy && N <= 10000
    % Решение с шагом h
    h = (b - a) / N;
    x_h = a:h:b;
    y_h = zeros(1, N+1);
    y_h(1) = y0;
    
    % Метод Рунге-Кутты с шагом h
    for n = 1:N
        xn = x_h(n);
        yn = y_h(n);
        
        k1 = h * f(xn, yn);
        k2 = h * f(xn + h/2, yn + k1/2);
        k3 = h * f(xn + h/2, yn + k2/2);
        k4 = h * f(xn + h, yn + k3);
        
        y_h(n+1) = yn + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    
    % Решение с шагом h/2
    h_half = h / 2;
    N_half = N * 2;
    x_half = a:h_half:b;
    y_half = zeros(1, N_half+1);
    y_half(1) = y0;
    
    % Метод Рунге-Кутты с шагом h/2
    for n = 1:N_half
        xn = x_half(n);
        yn = y_half(n);
        
        k1 = h_half * f(xn, yn);
        k2 = h_half * f(xn + h_half/2, yn + k1/2);
        k3 = h_half * f(xn + h_half/2, yn + k2/2);
        k4 = h_half * f(xn + h_half, yn + k3);
        
        y_half(n+1) = yn + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    
    % Сравнение решений в общих точках (оценка по правилу Рунге)
    y_exact = f1(x_h);
    error_h = max(abs(y_exact - y_h));
    error_half = max(abs(f1(x_half) - y_half));
    
    % Оценка погрешности по правилу Рунге
    % |y_h - y_exact| ≈ |y_h - y_h/2| / (2^p - 1)
    error_runge = max(abs(y_h - y_half(1:2:end))) / (2^p - 1);
    error_rk4 = max(error_half, error_runge);
    
    if (error_rk4 <= required_accuracy)
      is_sufficient = "Да";
    else
      is_sufficient = "Нет";
    endif
    printf('%5d | %7.5f | %e | %5s\n', N, h, error_rk4, is_sufficient);
    
    N *= 2; % Удвоение количества интервалов для следующей итерации
end

if (error_rk4 <= required_accuracy)
    printf('\nНайденный шаг h = %f обеспечивает точность %e\n\n', h, error_rk4);
    N = N / 2; % Возвращаемся к последнему удовлетворяющему значению
    h = (b - a) / N;
else
    printf('\nНе удалось достичь требуемой точности при N = %d. Используем максимальное N.\n\n', N/2);
    N = 10000;
    h = (b - a) / N;
endif
% Пункт 2: Решение задачи методами Рунге-Кутты 4-го порядка и Эйлера
% Инициализация переменных
x = a:h:b;
y_rk4 = zeros(1, N+1);
y_euler = zeros(1, N+1);
y_rk4(1) = y0;
y_euler(1) = y0;
y_an = f1(x); % Точное решение

printf('Решение задачи Коши методом Рунге-Кутты 4-го порядка и методом Эйлера\n');
printf('Количество шагов N = %d, шаг h = %.6f\n\n', N, h);

% Метод Рунге-Кутты 4-го порядка
for n = 1:N
    xn = x(n);
    yn = y_rk4(n);

    k1 = h * f(xn, yn);
    k2 = h * f(xn + h/2, yn + k1/2);
    k3 = h * f(xn + h/2, yn + k2/2);
    k4 = h * f(xn + h, yn + k3);

    y_rk4(n+1) = yn + (k1 + 2*k2 + 2*k3 + k4) / 6;
end

% Метод Эйлера
for n = 1:N
    xn = x(n);
    yn = y_euler(n);

    y_euler(n+1) = yn + h * f(xn, yn);
end
% Пункт 3: Сравнение с точным решением
delta_rk4 = abs(y_an - y_rk4);
delta_euler = abs(y_an - y_euler);

norm_rk4 = norm(delta_rk4);
norm_euler = norm(delta_euler);

max_abs_rk4 = max(delta_rk4);
max_abs_euler = max(delta_euler);

printf('Сравнение приближенных решений с точным:\n');
printf('Метод Рунге-Кутты 4-го порядка:\n');
printf('  - Норма невязки: %e\n', norm_rk4);
printf('  - Максимальное отклонение: %e\n', max_abs_rk4);
printf('Метод Эйлера:\n');
printf('  - Норма невязки: %e\n', norm_euler);
printf('  - Максимальное отклонение: %e\n\n', max_abs_euler);

% Пункт 4: Построение графика
figure;
plot(x, y_an, 'g-', 'LineWidth', 2, 'DisplayName', 'Точное решение (e^x)');
hold on;
plot(x, y_rk4, 'b-', 'LineWidth', 2, 'DisplayName', 'Метод Рунге-Кутта 4-го порядка');
plot(x, y_euler, 'r--', 'LineWidth', 2, 'DisplayName', 'Метод Эйлера');
xlabel('x');
ylabel('y');
title('Сравнение решений задачи Коши');
legend('show');
grid on;
hold off;

% Пункт 5: Создание таблицы результатов
% Выберем некоторое количество точек для отображения в таблице
num_table_points = min(10, N+1);
table_indices = round(linspace(1, N+1, num_table_points));

% Форматированный вывод таблицы
printf('\nТаблица 4.1. Результаты расчетов\n\n');
printf('%-10s %-15s %-15s %-15s %-15s %-15s\n', 'x', 'y', 'y_Эйлера', 'Δy_Эйлера', 'y_Рунге-Кутты', 'Δy_Рунге-Кутты');
printf('%-10s %-15s %-15s %-15s %-15s %-15s\n', '---------', '--------------', '--------------', '--------------', '--------------', '--------------');

for i = table_indices
    xi = x(i);
    y_exact = y_an(i);
    y_euler_i = y_euler(i);
    y_rk4_i = y_rk4(i);
    delta_euler_i = delta_euler(i);
    delta_rk4_i = delta_rk4(i);
    
    printf('%-10.6f %-15.8f %-15.8f %-15.8e %-15.8f %-15.8e\n', ...
        xi, y_exact, y_euler_i, delta_euler_i, y_rk4_i, delta_rk4_i);
end

% Сохранение результатов в CSV файл с
results_table = [x', y_an', y_euler', delta_euler', y_rk4', delta_rk4'];
headers = {'x', 'y_точное', 'y_Эйлера', 'Δy_Эйлера', 'y_Рунге-Кутты', 'Δy_Рунге-Кутты'};

% Создаем файл и записываем заголовки
fid = fopen('results_table.csv', 'w');
% Записываем строку заголовков
fprintf(fid, '%s,%s,%s,%s,%s,%s\n', headers{:});
fclose(fid);

% Дописываем данные в файл
dlmwrite('results_table.csv', results_table, '-append');
printf('\nРезультаты также сохранены в файл results_table.csv\n');
