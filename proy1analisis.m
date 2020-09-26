clear; clc; format long;
%Nombre de la función
f = 'rosenbrock';
%Punto inicial
x0 = [3.5 4.5]';
%Llamamos al método de región de confianza y medimos el tiempo de máquina
tic;
[xf, iter, deltas, puntos] = miregion(f, x0);
toc;
%Vectores de interés
delta=deltas;
x1=puntos(:,1);
x2=puntos(:,2);
iteracion=[1:iter]';


%Tabla de resultados
Resultados=table(iteracion, delta, x1, x2);
disp(Resultados);

%Valores en el óptimo
disp('Resultados en el óptimo:')
normGrad=norm(gradiente(f,xf));
optimo=table(iter, xf(1), xf(2), normGrad);
disp(optimo);
