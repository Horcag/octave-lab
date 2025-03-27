clear all;
n = 12; r = 6; q = 12;

% Генерация нижнетреугольной матрицы L и правой части B
L = tril(rand(n, n));
L = L + diag(rand(n, 1) + 1); % Диагональные элементы увеличены для устойчивости
B = rand(n, q);
X = zeros(n, q);

% Разделение на блоки
N = ceil(n / r);
block_indices = cell(N, 1);

for j = 1:N - 1
    block_indices{j} = ((j - 1) * r + 1):(j * r);
endfor
block_indices{N} = ((N - 1) * r + 1):n;


% Решение блочного уравнения
for j = 1:N
    j_indices = block_indices{j};
    X(j_indices, :) = L(j_indices, j_indices) \ B(j_indices, :);

    for i = j + 1:N
        i_indices = block_indices{i};
        B(i_indices, :) = B(i_indices, :) - L(i_indices, j_indices) * X(j_indices, :);
    end
endfor




