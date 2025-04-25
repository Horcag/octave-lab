clear all;
clc;
format long; % Устанавливаем формат вывода для высокой точности

% --- Определения функций ---
F1 = @(x, y) sin(x + 1) - y - 1.2;
F2 = @(x, y) 2 * x + cos(y) - 2;
F = @(x, y) [F1(x, y); F2(x, y)]; % Исходная система

% Функции для метода простых итераций (для сравнения/проверки)
phi1_y = @(x) sin(x + 1) - 1.2;  % y_next = phi1(x_current)
phi2_x = @(y) 1 - 0.5 * cos(y); % x_next = phi2(y_current)

% Якобиан
J = @(x, y) [cos(x + 1), -1;
             2,         -sin(y)];

% --- Графический анализ ---
figure;
hold on; grid on; xlabel('x'); ylabel('y'); title('Графики системы');
x_range = linspace(-2, 2, 200);
y1_plot = sin(x_range + 1) - 1.2;
plot(x_range, y1_plot, 'r', 'LineWidth', 1.5);
y_range_plot = linspace(-4, 2, 200);
x2_plot = (2 - cos(y_range_plot)) / 2;
plot(x2_plot, y_range_plot, 'b', 'LineWidth', 1.5);
legend('sin(x+1) - y = 1.2', '2x + cos(y) = 2');
axis equal;
printf('Начальное приближение (из графика): x=0.5, y=-1.5\n\n');

% --- Параметры точности и итераций ---
tol = 1e-4;      % Точность для критерия остановки ||Xk - Xk-1|| < tol
max_iter = 500;  % Максимальное число итераций
% --- Стандартный Метод Ньютона ---
printf('--- Стандартный Метод Ньютона ---\n');
xy_newton = [0; 0]; % Начальное приближение
iter_newton = 0;
error_newton = inf;

printf('Итерация |     x_k    |     y_k    |   ||Xk-Xk-1|| \n');
printf('-----------------------------------------------------\n');
printf('%8d | %10.6f | %10.6f |      -       \n', iter_newton, xy_newton(1), xy_newton(2));

while error_newton > tol && iter_newton < max_iter
    iter_newton = iter_newton + 1;
    xy_prev = xy_newton; % Сохраняем предыдущее

    J_val = J(xy_prev(1), xy_prev(2));
    F_val = F(xy_prev(1), xy_prev(2));

    if abs(det(J_val)) < 1e-10
        warning('Ньютон: Якобиан близок к вырожденному на итерации %d.', iter_newton);
        iter_newton = max_iter + 1; % Прерываем цикл с ошибкой
        break;
    endif

    delta = J_val \ (-F_val);
    xy_newton = xy_prev + delta;
    error_newton = norm(xy_newton - xy_prev);

    printf('%8d | %10.6f | %10.6f |  %e\n', iter_newton, xy_newton(1), xy_newton(2), error_newton);
endwhile

printf('-----------------------------------------------------\n');
if iter_newton > max_iter
    fprintf('Стандартный метод Ньютона НЕ сошелся за %d итераций.\n', max_iter);
else
    fprintf('Стандартный метод Ньютона сошелся за %d итераций.\n', iter_newton);
    fprintf('Решение (Ньютон): x = %.15f, y = %.15f\n', xy_newton(1), xy_newton(2));
endif
printf('\n');


% --- Модифицированный Метод Ньютона ---
printf('--- Модифицированный Метод Ньютона ---\n');
xy_mod_newton = [0; 0]; % Начальное приближение
X0_vec_mod = xy_mod_newton;
iter_mod_newton = 0;
error_mod_newton = inf;

% Вычисляем Якобиан один раз
J0_val_mod = J(X0_vec_mod(1), X0_vec_mod(2));
printf('Фиксированный Якобиан J(X0):\n');
disp(J0_val_mod);
if abs(det(J0_val_mod)) < 1e-10
    error('Мод. Ньютон: Начальный Якобиан близок к вырожденному.');
endif

printf('Итерация |     x_k    |     y_k    |   ||Xk-Xk-1|| \n');
printf('-----------------------------------------------------\n');
printf('%8d | %10.6f | %10.6f |      -       \n', iter_mod_newton, xy_mod_newton(1), xy_mod_newton(2));

while error_mod_newton > tol && iter_mod_newton < max_iter
    iter_mod_newton = iter_mod_newton + 1;
    xy_prev_mod = xy_mod_newton; % Сохраняем предыдущее

    F_val_mod = F(xy_prev_mod(1), xy_prev_mod(2));
    delta_mod = J0_val_mod \ (-F_val_mod); % Используем J0
    xy_mod_newton = xy_prev_mod + delta_mod;
    error_mod_newton = norm(xy_mod_newton - xy_prev_mod);

    printf('%8d | %10.6f | %10.6f |  %e\n', iter_mod_newton, xy_mod_newton(1), xy_mod_newton(2), error_mod_newton);
endwhile

printf('-----------------------------------------------------\n');
if iter_mod_newton > max_iter
    fprintf('Модифицированный метод Ньютона НЕ сошелся за %d итераций.\n', max_iter);
else
    fprintf('Модифицированный метод Ньютона сошелся за %d итераций.\n', iter_mod_newton);
    fprintf('Решение (Мод. Ньютон): x = %.15f, y = %.15f\n', xy_mod_newton(1), xy_mod_newton(2));
endif
printf('\n');

% --- Метод Простых Итераций (для сравнения) ---
printf('--- Метод Простых Итераций ---\n');
xy_mpi = [0; 0]; % Начальное приближение
iter_mpi = 0;
error_mpi = inf;
max_iter_mpi = 500; % МПИ может требовать больше итераций

printf('Итерация |     x_k    |     y_k    |   ||Xk-Xk-1|| \n');
printf('-----------------------------------------------------\n');
printf('%8d | %10.6f | %10.6f |      -       \n', iter_mpi, xy_mpi(1), xy_mpi(2));

while error_mpi > tol && iter_mpi < max_iter_mpi
    iter_mpi = iter_mpi + 1;
    xy_prev_mpi = xy_mpi;

    x_next_mpi = phi2_x(xy_prev_mpi(2));
    y_next_mpi = phi1_y(xy_prev_mpi(1));
    xy_mpi = [x_next_mpi; y_next_mpi];
    error_mpi = norm(xy_mpi - xy_prev_mpi);

    printf('%8d | %10.6f | %10.6f |  %e\n', iter_mpi, xy_mpi(1), xy_mpi(2), error_mpi);
     if any(isinf(xy_mpi)) || any(isnan(xy_mpi)) || norm(xy_mpi) > 1e6
       warning('МПИ: Возможная расходимость на итерации %d.', iter_mpi);
       iter_mpi = max_iter_mpi + 1; break;
    endif
endwhile
printf('-----------------------------------------------------\n');
if iter_mpi > max_iter_mpi
    fprintf('Метод простых итераций НЕ сошелся за %d итераций.\n', max_iter_mpi);
else
     fprintf('Метод простых итераций сошелся за %d итераций.\n', iter_mpi);
     fprintf('Решение (МПИ): x = %.15f, y = %.15f\n', xy_mpi(1), xy_mpi(2));
endif
printf('\n');


% --- Итоговая таблица сравнения итераций ---
printf('--- Сравнение количества итераций для точности %e ---\n', tol);
printf('----------------------------------------------------------\n');
printf('| Метод                       | Итераций | Сошелся? |\n');
printf('----------------------------------------------------------\n');

% Результат для стандартного Ньютона
if iter_newton <= max_iter
    printf('| Стандартный метод Ньютона  | %8d |    Да    |\n', iter_newton);
else
    printf('| Стандартный метод Ньютона  | %8d |    Нет   |\n', max_iter);
endif

% Результат для модифицированного Ньютона
if iter_mod_newton <= max_iter
    printf('| Модифиц. метод Ньютона   | %8d |    Да    |\n', iter_mod_newton);
else
    printf('| Модифиц. метод Ньютона   | %8d |    Нет   |\n', max_iter);
endif

% Результат для МПИ
if iter_mpi <= max_iter_mpi
    printf('| Метод простых итераций   | %8d |    Да    |\n', iter_mpi);
else
    printf('| Метод простых итераций   | %8d |    Нет   |\n', max_iter_mpi);
endif
printf('----------------------------------------------------------\n');

% Финальная проверка разницы решений (если все сошлись)
if iter_newton <= max_iter && iter_mod_newton <= max_iter && iter_mpi <= max_iter_mpi
    norm_diff_N_ModN = norm(xy_newton - xy_mod_newton);
    norm_diff_N_MPI = norm(xy_newton - xy_mpi);
    printf('Норма разницы ||X_Newton - X_ModNewton|| = %e\n', norm_diff_N_ModN);
    printf('Норма разницы ||X_Newton - X_MPI||      = %e\n', norm_diff_N_MPI);
endif