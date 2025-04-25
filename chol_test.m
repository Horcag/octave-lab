% --- НАЧАЛО КОДА ---
clear all;

% Функция f1: Матричный метод разложения Холецкого
function G = f1 (A)
    n = size(A, 1);
    m = size(A, 2);
    G = zeros(n, m);
    for i=1:n
        G(i,1:i-1) = A(i,1:i-1) * inv(transpose(G(1:i-1,1:i-1)));
        G(i,i) = sqrt(A(i,i) - G(i,1:i-1) * transpose(G(i,1:i-1)));
      endfor
endfunction

% Функция f2: Векторный метод разложения Холецкого
function G = f2 (A)
    n = size(A, 1);
    m = size(A, 2);

    G = zeros(n, n);
    for i = 1:n 
        for j = 1:i-1
            sum_k = G(i, 1:j-1) * G(j, 1:j-1)'; %
            G(i, j) = (A(i, j) - sum_k) / G(j, j);
        endfor
        sum_k_diag = G(i, 1:i-1) * G(i, 1:i-1)';
        val_sqrt = A(i, i) - sum_k_diag;
        G(i, i) = sqrt(val_sqrt);
    endfor
endfunction

% Функция f3: Скалярный (поэлементный) метод разложения Холецкого
function G = f3 (A)
    n = size(A, 1);
    m = size(A, 2);

    G = zeros(n, n); % Инициализация G

    for j = 1:n % Цикл по столбцам j
        sum_sq = 0;
        % Цикл по k для диагонального элемента G(j,j)
        for k = 1:j-1
            sum_sq = sum_sq + G(j, k)^2;
        endfor
        diag_val_sq = A(j, j) - sum_sq;
        G(j, j) = sqrt(diag_val_sq); % Диагональный G(j,j)

        % Цикл по строкам i ниже диагонали для столбца j
        for i = j+1:n
            sum_prod = 0;
            % Цикл по k для поддиагонального элемента G(i,j)
            for k = 1:j-1
                sum_prod = sum_prod + G(i, k) * G(j, k);
            endfor
            % Используем A(i,j) из-за симметрии A(i,j)=A(j,i)
            G(i, j) = (A(i, j) - sum_prod) / G(j, j); % Поддиагональный G(i,j)
        endfor
    endfor
endfunction

% --- Пример использования и проверка ---
n = 20; % Размер матрицы

% Создаем симметричную, положительно определенную матрицу A
A_rand = rand(n, n); % Случайная матрица
A = A_rand * A_rand' + n * eye(n); % A = M*M' + диаг. преобладание

fprintf('Тестирование разложения Холецкого для матрицы %dx%d\n', n, n);

% Тест f1 (матричный)
try
    tic; G1 = f1(A); time1 = toc;
    err1 = norm(A - G1*G1', 'fro');
    fprintf('Метод f1 (матричный): Время = %.6f c, Ошибка ||A-GG''|| = %e\n', time1, err1);
catch ME
    fprintf('Метод f1 (матричный): ОШИБКА - %s\n', ME.message);
end_try_catch

% Тест f2 (векторный)
try
    tic; G2 = f2(A); time2 = toc;
    err2 = norm(A - G2*G2', 'fro');
    fprintf('Метод f2 (векторный): Время = %.6f c, Ошибка ||A-GG''|| = %e\n', time2, err2);
catch ME
    fprintf('Метод f2 (векторный): ОШИБКА - %s\n', ME.message);
end_try_catch

% Тест f3 (скалярный)
try
    tic; G3 = f3(A); time3 = toc;
    err3 = norm(A - G3*G3', 'fro');
    fprintf('Метод f3 (скалярный): Время = %.6f c, Ошибка ||A-GG''|| = %e\n', time3, err3);
 catch ME
    fprintf('Метод f3 (скалярный): ОШИБКА - %s\n', ME.message);
end_try_catch

% --- КОНЕЦ КОДА ---
