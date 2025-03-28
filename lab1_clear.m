clear all;

% БЛОЧНЫЙ АЛГОРИТМ РЕШЕНИЯ НИЖНЕТРЕУГОЛЬНОЙ СИСТЕМЫ УРАВНЕНИЙ
% Функция реализует блочный вариант алгоритма прямой подстановки
% для решения системы линейных уравнений LX = B.
% Параметры:
% L - нижнетреугольная матрица коэффициентов
% B - матрица правых частей
% p - размер блока (n должно делиться на p без остатка)
function X = block_forward_substitution(L, B, p)
    n = size(L, 1);
    m = size(B, 2);
    r = n / p;
    X = zeros(n, m);

    for i = 1:r
        rows_Lii = (i - 1) * p + 1:i * p;
        Lii = L(rows_Lii, rows_Lii);
        
        rows_Bii = (i - 1) * p + 1:i * p;
        Bi = B(rows_Bii, :);
        
        S = zeros(p, m);
        
        for j = 1:i - 1
            cols_Lij = (j - 1) * p + 1:j * p;
            Lij = L(rows_Lii, cols_Lij);
            Xj = X(cols_Lij, :);
            S = S + Lij * Xj;
        endfor

        Bii = Bi - S;
        Xi = Lii \ Bii;
        X(rows_Lii, :) = Xi;
    endfor
endfunction

% Тестовый пример
L = [2 0 0 0; 
     2 4 0 0; 
     4 2 6 0; 
     2 4 3 4]; 
B = [10; 22; 38; 80]; 
p = 2;
X = block_forward_substitution(L, B, p);

X