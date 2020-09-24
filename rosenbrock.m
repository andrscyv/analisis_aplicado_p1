function fx = rosenbrock(x)
% funci�n de Rosenbrock: f: R^2 --> R
% cuyo m�nimo local es muy dif�cil de alcanzar por medio
% de los m�todos de optimizaci�n.
%In
% x.- vector de longitud 2
% fx.- n�mero real con el valor de la funci�n.
%  
% ITAM
% An�lisis Aplicado
% 12 de agosto de 2020
%---------------------------------------------

fx = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;