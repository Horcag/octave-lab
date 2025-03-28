clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% БЛОЧНЫЙ АЛГОРИТМ РЕШЕНИЯ НИЖНЕТРЕУГОЛЬНОЙ СИСТЕМЫ УРАВНЕНИЙ
%
% ОПИСАНИЕ:
% Эта функция реализует блочный вариант алгоритма прямой подстановки
% для решения системы линейных уравнений LX = B, где:
% - L: нижнетреугольная матрица размера n×n
% - B: матрица правых частей размера n×m
% - X: искомая матрица решений размера n×m
%
% ПРЕИМУЩЕСТВА БЛОЧНОГО ПОДХОДА:
% 1. Лучшая локальность данных, что снижает количество кэш-промахов
% 2. Возможность использования оптимизированных BLAS операций над блоками
% 3. Более эффективное использование параллелизма
%
% ПАРАМЕТРЫ:
% L - нижнетреугольная матрица коэффициентов
% B - матрица правых частей
% p - размер блока (n должно делиться на p без остатка)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = block_forward_substitution(L, B, p)
    % Определяем размеры входных матриц:
    % n - количество строк/столбцов в L (квадратная матрица)
    % m - количество столбцов в B (количество систем для решения)
    n = size(L, 1); % Размер матрицы L (количество строк)
    m = size(B, 2); % Количество столбцов в B

    % Вычисляем количество блоков r = n/p
    % Весь алгоритм основан на разделении матриц L, B и X на блоки размера p×p
    r = n / p; % Количество блоков (предполагается, что n делится на p)

    % Создаем матрицу результата, заполненную нулями
    % В процессе работы алгоритма она будет заполняться решением
    X = zeros(n, m); % Инициализация результата

    % Основной цикл по блокам строк матрицы L
    % i изменяется от 1 до r, где r - количество блоков
    printf("L =\n");
    disp(L)

    printf("B =\n");
    disp(B)
    for i = 1:r
        printf("i = %d\n", i);
        % Вычисляем индексы строк текущего диагонального блока L_ii
        % Например, если p=2 и i=2, то rows_Lii = [3, 4]
        rows_Lii = (i - 1) * p + 1:i * p; % Индексы строк для текущего блока Lii
        printf("rows_Lii =\n");
        disp(rows_Lii);
        printf("\n");
        % Извлекаем текущий диагональный блок L_ii размера p×p
        % Этот блок будет использован для решения системы уравнений
        Lii = L(rows_Lii, rows_Lii); % Диагональный блок L
        printf("Lii =\n");
        disp(Lii) 
        printf("\n")

        % Индексы для соответствующих строк в матрице B
        % Извлекаем подблок B_i из матрицы B
        rows_Bii = (i - 1) * p + 1:i * p; % Индексы строк для текущего блока B
        printf("rows_Bii =\n");
        disp(rows_Bii) 
        printf("\n")
        Bi = B(rows_Bii, :); % Текущий блок B
        printf("Bi =\n");
        disp(Bi) 
        printf("\n")
        % Инициализируем матрицу S для накопления суммы L_ij * X_j
        % Эта сумма представляет влияние уже найденных блоков X_j на текущий блок
        S = zeros(p, m); % Инициализация суммы для коррекции Bi
        printf("S =\n");
        disp(S) 
        printf("\n")

        % Внутренний цикл для вычисления суммы L_ij * X_j для всех j < i
        for j = 1:i - 1
            printf("j = %d\n", j);
            % Вычисляем индексы столбцов для блока L_ij
            % Например, если p=2 и j=1, то cols_Lij = [1, 2]
            cols_Lij = (j - 1) * p + 1:j * p; % Индексы столбцов для блока Lij
            printf("cols_Lij =\n");
            disp(cols_Lij) 
            printf("\n")
            
            % Извлекаем недиагональный блок L_ij размера p×p
            Lij = L(rows_Lii, cols_Lij); % Блок Lij
            printf("Lij =\n");
            disp(Lij) 
            printf("\n")
            
            % Извлекаем ранее вычисленный блок решения X_j
            Xj = X(cols_Lij, :); % Уже вычисленный блок Xj
            printf("Xj =\n");
            disp(Xj) 
            printf("\n")

            % Обновляем сумму: S = S + L_ij * X_j
            % Это матричное произведение размера p×m
            S = S + Lij * Xj; % Обновление суммы
            printf("S =\n");
            disp(S) 
            printf("\n")
        endfor

        % Вычисляем скорректированную правую часть:
        % B_i' = B_i - Σ(L_ij * X_j) для всех j < i
        Bii = Bi - S; % Коррекция блока B
        printf("Bii =\n");
        disp(Bii) 
        printf("\n")
        
        % Решаем систему L_ii * X_i = B_i'
        % Оператор "\" в Octave решает систему линейных уравнений
        Xi = Lii \ Bii; % Решение для текущего блока
        printf("Xi =\n");
        disp(Xi) 
        printf("\n")
        
        % Сохраняем найденный блок решения X_i в соответствующие позиции матрицы X
        X(rows_Lii, :) = Xi; % Сохранение результата
        printf("X =\n");
        disp(X) 
        printf("\n")
    endfor
endfunction

% ТЕСТИРОВАНИЕ ФУНКЦИИ:
% Создаем тестовые матрицы и проверяем правильность работы алгоритма

# n = 100; m = 50; % Размеры матриц
# r = 50; % Размер блока
# L = tril(rand(n, n)); % Создаем случайную нижнетреугольную матрицу
# L = L + diag(rand(n, 1) + 1); % Диагональные элементы от 1 до 2 (невырожденная L)
# B = rand(n, m); % Создаем случайную матрицу правых частей

# % Применяем нашу блочную функцию для решения системы
# X = block_forward_substitution(L, B, r);

# % Используем встроенную функцию для проверки (эталонное решение)
# X_expected = L \ B; % Стандартное решение

# % Сравниваем результаты: вычисляем норму разности между нашим решением и эталонным
# printf("Норма разности X и X_expected = %d\n", norm(X - X_expected));
# % Если норма близка к нулю (порядка 1e-15), значит алгоритм работает корректно
L = [2 0 0 0; 
     2 4 0 0; 
     4 2 6 0; 
     2 4 3 4]; 
B = [10; 22; 38; 80]; 
p = 2;
X = block_forward_substitution(L, B, p); 
  
X