function [xf, iter] = miregion(f, x0)
% Método de región de confianza para aproximar un mínimo
% de f: Rn->R dos veces continuamente diferenciable
% In
% f cadena de caracteres con la codigo en Matlab de la función a minimizar.
% x0 vector columna de dimensión n con el punto inicial.
% Out
% xf vector columna de dimension n con la aproximación final.
% iter número de iteraciones que se realizaron.
%
% Andres Cruz y Vera C.U.155899
% Javier Montiel González C.U.159216
%-------------------------------------------------------------------------
% parámetros
deltamin = 1.e-04;
deltamax = 5; 
eta = 0.25;
tol = 1.e-06; % tolerancia para la norma del gradiente.
maxiter = 100;    % número máximo de iteraciones externas permitidas
maxregion = 20;    %es el numero máximo que permanece region de confianza 
%en un solo punto
tol = 1.e-06; % tolerancia para la norma del gradiente.

% valores iniciales
iter = 0;        % contador para las iteraciones externas
jregion = 0;       % contador interno

n = length(x0);
g = gradiente(fname,x);
ng = norm(g);
H = hessian(fname,x0);
delta = deltamax;
xk = x0;

while(ng>tol && iter<maxiter && jregion<maxregion)
    pk = doblez(B, g, delta);
    redact = feval(f,xk)-feval(f,xk+pk);
    redpre = -((1/2)*pk'*H*pk+g'*pk);
end


end

