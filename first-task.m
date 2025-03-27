clear all;

## Функция, реализующая умножение матриц A и B скалярным способом
function C = f1(A, B, n)
    C = zeros(n, n); % Матрица, в которую будет записан результат умножения (заполнена нулями)

    for i = 1:n

        for j = 1:n

            for k = 1:n
                C(i, j) = C(i, j) + A(i, k) * B(k, j); % Скалярное произведение строки i матрицы A на столбец j матрицы B
            endfor

        endfor

    endfor

    printf("f1. Норма = %d\n", norm(A * B - C)); % Норма разности матриц AB и C должна быть близка к нулю (порядка 1e-15)

endfunction

## Функция, реализующая умножение матриц A и B с помощью операции saxpy
function C = f2(A, B, n)
    %% SAXPY - Одинарной точности A*X плюс Y
    %
    % Выполняет векторную операцию: y = a*x + y, где 'a' - скаляр, а 'x' и 'y' - векторы.
    % Это фундаментальная операция в линейной алгебре и ключевой компонент в BLAS
    % (Базовые подпрограммы линейной алгебры) уровня 1.
    %
    % SAXPY отличается от других базовых операций тем, что:
    % - Ограничена пропускной способностью памяти, а не вычислительной мощностью
    % - Служит распространенным эталоном для параллельных вычислительных систем
    % - Имеет отличные свойства для распараллеливания
    % - Является "s" (одинарной точности) вариантом семейства функций AXPY
    %
    % Название SAXPY означает "Single-precision A·X Plus Y" (Одинарной точности A·X плюс Y).
    C = zeros(n, n);

    for k = 1:n

        for i = 1:n

            C(i, :) = C(i, :) + A(i, k) * B(k, :); % saxpy-операция

        endfor

    endfor

    printf("f2. Норма = %d\n", norm(A * B - C));
endfunction

## Функция, реализующая умножение матриц A и B с помощью внешнего произведения
function C = f3 (A, B, n)
    C = zeros(n, n);

    for k = 1:n

        C= C + A(:, k) * B(k, :); % mm - внешнее произведение (матричная операция)

    endfor

    printf("f3. Норма = %d\n", norm(A * B - C));

endfunction

## Функция, реализующая блочное умножение матриц A и B
function f4 (A, B, n, N, alpha)
    
    C = zeros(n, n); 

    for kb = 1:N
        k = (kb - 1) * alpha + 1:kb * alpha; % Индексы блока k

        for ib = 1:N
            i = (ib - 1) * alpha + 1:ib * alpha; % Индексы блока i
            for jb = 1:N
                j = (jb - 1) * alpha + 1:jb * alpha; % Индексы блока j
                C(i, j) = C(i, j) + A(i, k) * B(k, j); 
            endfor

        endfor

    endfor

    # printf("f4. Норма = %d\n", norm(A * B - C)); 
  
endfunction


n = 400;
A = rand(n, n); B = rand(n, n); % Матрицы заполняются случайными числами от 0 до 1
alpha = 200; % Размер блока
N = n / alpha; % Количество блоков


# tic;
# f1(A, B, n);
# time_f1 = toc;
# printf("Время выполнения f1: %f секунд\n", time_f1);

# tic;
# f2(A, B, n);
# time_f2 = toc;
# printf("Время выполнения f2: %f секунд\n", time_f2);

# tic;
# f3(A, B, n);
# time_f3 = toc;
# printf("Время выполнения f3: %f секунд\n", time_f3);

# tic;
# f4(A, B, n, N, alpha);
# time_f4 = toc;
# printf("Время выполнения f4: %f секунд\n", time_f4);


iterations = 50;
time_f4 = 0;
for i = 1:iterations
    tic;
    f4(A, B, n, N, alpha);
    time_f4 = time_f4 + toc;
endfor
time_f4 = time_f4 / iterations;
printf("Среднее время f4: %f секунд\n", time_f4);


