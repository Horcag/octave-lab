clear all;

% Функция f1: реализует разложение Холецкого матрицы A матричным методом
function G = f1 (A)
    n = size(A, 1);
    m = size(A, 2);

    if n != m
        error("Matrix A must be square.");
    endif

    G = zeros(n, m);

    for i = 1:n
        G(i, i) = sqrt(A(i, i) - G(i, 1:i - 1) * G(i, 1:i - 1)');
        G(i + 1:n, i) = (A(i + 1:n, i) - G(i + 1:n, 1:i - 1) * G(i, 1:i - 1)') / G(i, i);

    endfor

endfunction

% Функция f2: реализует разложение Холецкого матрицы A векторным методом
function G = f2 (A)
    n = size(A, 1);
    m = size(A, 2);

    if n != m
        error("Matrix A must be square.");
    endif

    % Инициализируем массив для хранения столбцов матрицы G
    G_cols = cell(n, 1);

    for j = 1:n
        % Создаем j-й столбец G
        g_col = zeros(n, 1);

        % Вычисление диагонального элемента (G[j,j])
        sum_sq = 0;

        for k = 1:j - 1
            % Берем j-й элемент из k-го столбца
            g_kj = G_cols{k}(j);
            sum_sq = sum_sq + g_kj^2;
        endfor

        g_col(j) = sqrt(A(j, j) - sum_sq);

        % Вычисление элементов ниже диагонали в столбце j
        for i = j + 1:n
            sum_prod = 0;

            for k = 1:j - 1
                % Скалярное произведение соответствующих элементов векторов
                g_ki = G_cols{k}(i);
                g_kj = G_cols{k}(j);
                sum_prod = sum_prod + g_ki * g_kj;
            endfor

            g_col(i) = (A(i, j) - sum_prod) / g_col(j);
        endfor

        % Сохраняем вычисленный столбец
        G_cols{j} = g_col;
    endfor

    % Собираем финальную матрицу G из столбцов
    G = zeros(n, n);

    for j = 1:n
        G(:, j) = G_cols{j};
    endfor

endfunction

% Функция f3: реализует разложение Холецкого матрицы A скалярным (поэлементным) методом
function G = f3 (A)
    n = size(A, 1);
    m = size(A, 2);

    if n != m
        error("Matrix A must be square.");
    endif

    % Инициализируем нижнюю треугольную матрицу G нулями
    G = zeros(n, n);

    % Проходим по столбцам
    for j = 1:n
        % Сначала вычисляем диагональный элемент G(j, j)
        sum_sq = 0;
        % Суммируем квадраты элементов строки j левее диагонали
        for k = 1:j-1
            sum_sq = sum_sq + G(j, k)^2;
        endfor

        % Вычисляем элемент под корнем
        diag_val_sq = A(j, j) - sum_sq;

        % Находим диагональный элемент
        G(j, j) = sqrt(diag_val_sq);

        % Затем вычисляем поддиагональные элементы в столбце j (строки i > j)
        for i = j+1:n
            sum_prod = 0;
            % Суммируем произведения элементов из строк i и j левее столбца j
            for k = 1:j-1
                sum_prod = sum_prod + G(i, k) * G(j, k);
            endfor
            % Находим поддиагональный элемент
            % Используем A(i, j), т.к. для симметричной A(i,j) == A(j,i)
            G(i, j) = (A(i, j) - sum_prod) / G(j, j); % Делим на уже вычисленный G(j,j)
        endfor
    endfor

endfunction

n = 10;
A = tril(rand(n, n));
A = A + A' + n * eye(n);
G = f3(A);
G_t = G';

disp(norm(A - G * G_t))
