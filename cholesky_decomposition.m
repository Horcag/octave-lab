clear all;

function G = f1 (A)
    n = size(A, 1);
    m = size(A, 2);

    if n != m
        error("Matrix A must be square.");
    endif

    G = zeros(n, m);

    for i = 1:n
        G(i, i) = sqrt(A(i, i) - G(i, 1:i - 1) * G(i, 1:i - 1)');

        for j = (i + 1):n
            G(j, i) = (A(j, i) - G(j, 1:i - 1) * G(i, 1:i - 1)') / G(i, i);
        endfor

    endfor

endfunction

A = [4 2 2;
    2 5 3;
    2 3 6];
G = f1(A);

disp(G);
