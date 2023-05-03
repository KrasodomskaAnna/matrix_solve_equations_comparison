clc
clear all
close all

% odpowiednie fragmenty kodu mozna wykonac poprzez zaznaczenie i wcisniecie F9 w Matlabie
% komentowanie/odkomentowywanie: ctrl+r / ctrl+t

% Zadanie A
%------------------
N = 10;
density = 3; % parametr decydujacy o gestosci polaczen miedzy stronami
[Edges] = generate_network(N, density);
%-----------------

% Zadanie B
%------------------
% generacja macierzy I, A, B i wektora b
% macierze A, B i I musza byc przechowywane w formacie sparse (rzadkim)
B = sparse(Edges(2,:), Edges(1,:), ones(), N, N);
L = sum(B);
A = spdiags(transpose(1./L), 0, N, N);
I = speye(N);
d = 0.85;
M = I - d*B*A;
b(1:N, 1) = (1-d)/N;

if ~issparse(M)
    disp("M is not sparse");
    return;
end
%-----------------
r = M\b;
r_lower = find(r < 0);
if ~isempty(r_lower)
    disp("Not every element of r is grater or equal than zero");
    return;
end

% Zadanie D
%------------------
clc
clear all
close all

N = [500, 1000, 3000, 6000, 12000];

density = 10; % parametr decydujacy o gestosci polaczen miedzy stronami
d = 0.85;
for i = 1:5
    [Edges] = generate_network(N(i), density);
    B = sparse(Edges(2,:), Edges(1,:), ones(), N(i), N(i));
    L = sum(B);
    A = spdiags(transpose(1./L), 0, N(i), N(i));
    I = speye(N(i));
    M = I - d*B*A;
    b(1:N(i), 1) = (1-d)/N(i);
    tic
    % obliczenia start
    r = M\b;
    % obliczenia stop
    czas_Gauss(i) = toc;
end
plot(N, czas_Gauss);
title(["Wykres zależności czasu obliczeń", "od liczby stron metodą bezpośrednią Gaussa"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieD;
%------------------



% Zadanie E
%------------------
clc
clear all
close all

% sprawdz przykladowe dzialanie funkcji tril, triu, diag:
% Z = rand(4,4)
% tril(Z,-1)
% triu(Z,1)
% diag(diag(Z))

N = [500, 1000, 3000, 6000, 12000];
density = 10;
d = 0.85;
% r(k+1) =−D−1(L+U)r(k) +D−1b
for i = 1:5
    [Edges] = generate_network(N(i), density);
    B = sparse(Edges(2,:), Edges(1,:), ones(), N(i), N(i));
    L = sum(B);
    A = spdiags(transpose(1./L), 0, N(i), N(i));
    I = speye(N(i));
    M = I - d*B*A;
    b(1:N(i), 1) = (1-d)/N(i);
    r = ones(N(i), 1);
    L = tril(M,-1);
    U = triu(M,1);
    D = diag(diag(M));
    res_norm = norm((M*r)-b);
    iteration = 0;
    czas_Jacobi(i) = 0;
    while res_norm > 10^-14
        tic
        r = -(D\((L+U)*r))+D\b;
        czas_Jacobi(i) = czas_Jacobi(i) + toc;
        iteration = iteration+1;
        if N(i) == 1000
            residual_err(iteration) = res_norm;
        end
        res_norm = norm((M*r)-b);
    end
    % obliczenia stop
    iterations(i) = iteration;
end

figure(1);
plot(N, czas_Jacobi);
title(["Wykres zależności czasu wyznaczenia rozwiązania", "od liczby stron metodą iteracyjną Jacobiego"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieE_czas;

figure(2);
plot(N, iterations);
title(["Wykres zależności liczby iteracji wymaganej do osiągnięcia rozwiązania", "od liczby stron metodą iteracyjną Jacobiego"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieE_iteracje;

figure(3);
semilogy(1:length(residual_err), residual_err);
title(["Wykres normy błędu rezydualnego w kolejnych iteracjach badanego algorytmu", "dla N = 1000 metodą iteracyjną Jacobiego"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieE_norma;

%------------------


% Zadanie F
%------------------
clc
clear all
close all

N = [500, 1000, 3000, 6000, 12000];
density = 10;
d = 0.85;
% r(k+1) =−(D+L)−1Ur(k)+(D+L)−1b
for i = 1:5
    [Edges] = generate_network(N(i), density);
    B = sparse(Edges(2,:), Edges(1,:), ones(), N(i), N(i));
    L = sum(B);
    A = spdiags(transpose(1./L), 0, N(i), N(i));
    I = speye(N(i));
    M = I - d*B*A;
    b(1:N(i), 1) = (1-d)/N(i);
    r = ones(N(i), 1);
    L = tril(M,-1);
    U = triu(M,1);
    D = diag(diag(M));
    res_norm = norm((M*r)-b);
    iteration = 0;
    czas_Gauss_Seidl(i) = 0;
    while res_norm > 10^-14
        tic
        r = -((D+L)\(U*r)) +(D+L)\b;
        czas_Gauss_Seidl(i) = czas_Gauss_Seidl(i) + toc;
        iteration = iteration+1;
        if N(i) == 1000
            residual_err(iteration) = res_norm;
        end
        res_norm = norm((M*r)-b);
    end
    iterations(i) = iteration;
end

figure(1);
plot(N, czas_Gauss_Seidl);
title(["Wykres zależności czasu wyznaczenia rozwiązania", "od liczby stron metodą iteracyjną Gaussa-Seidla"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieF_czas;

figure(2);
plot(N, iterations);
title(["Wykres zależności liczby iteracji wymaganej do osiągnięcia rozwiązania", "od liczby stron metodą iteracyjną Gaussa-Seidla"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieF_iteracje;

figure(3);
semilogy(1:length(residual_err), residual_err);
title(["Wykres normy błędu rezydualnego w kolejnych iteracjach badanego algorytmu", "dla N = 1000 metodą iteracyjną Gaussa-Seidla"]);
xlabel("liczba stron N");
ylabel("czas obliczeń [s]");
print -dpng zadanieF_norma;


% Zadanie G
%------------------
clc
clear all
close all
load Dane_Filtr_Dielektryczny_lab3_MN.mat

tic
r = M\b;
czas_Gauss = toc;
norm_Gauss = norm((M*r)-b);
disp("norma metody bezpośredniej: " + norm_Gauss);


r = ones(length(M), 1);
L = tril(M,-1);
U = triu(M,1);
D = diag(diag(M));
res_norm = norm((M*r)-b);
i = 0;
czas_Jacobi = 0;
min_x = 0;
min_y = 10000000000;
while i <= 1000 && res_norm > 10^-14 && ~isnan(res_norm) && res_norm ~= Inf
    tic
    r = -(D\((L+U)*r))+D\b;
    czas_Jacobi = czas_Jacobi + toc;
    i = i+1;
    if min_y > res_norm
        min_y = res_norm;
        min_x = i;
    end
    norm_Jacobi(i) = res_norm;
    res_norm = norm((M*r)-b);
end
disp("norma metoda Jacobiego: " + norm_Jacobi(i));
disp("minimum lokalne: " + min_y + ", dla liczby stron wynoszącej: " + min_x);

r = ones(length(M), 1);
res_norm = norm((M*r)-b);
i = 0;
czas_Gauss_Seidl = 0;
min_x = 0;
min_y = 10000000000;
while i <= 1000 && res_norm > 10^-14 && ~isnan(res_norm) && res_norm ~= Inf
    tic
    r = -((D+L)\(U*r)) +(D+L)\b;
    czas_Gauss_Seidl = czas_Gauss_Seidl + toc;
    i = i+1;
    if min_y > res_norm
        min_y = res_norm;
        min_x = i;
    end
    norm_Gauss_Seidl(i) = res_norm;
    res_norm = norm((M*r)-b);
end
disp("norma metoda Gaussa-Seidla: " + norm_Gauss_Seidl(i));
disp("minimum lokalne: " + min_y + ", dla liczby stron wynoszącej: " + min_x);

figure(1);
semilogy(1:length(norm_Jacobi), norm_Jacobi);
title("Wykres normy błędu rezydualnego w kolejnych iteracjach metody Jacobiego dla N = 1000");
xlabel("wartość N");
ylabel("liczba iteracji");

figure(2);
semilogy(1:length(norm_Gauss_Seidl), norm_Gauss_Seidl);
title("Wykres normy błędu rezydualnego w kolejnych iteracjach metody Gaussa-Seidla dla N = 1000");
xlabel("wartość N");
ylabel("liczba iteracji");