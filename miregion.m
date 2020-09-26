function [xf, iter, delta] = miregion(f, x0)
% Método de región de confianza para aproximar un mínimo
% de f: Rn->R dos veces continuamente diferenciable
% In
% f cadena de caracteres con el código en Matlab de la función a minimizar.
% x0 vector columna de dimensión n con el punto inicial.
% Out
% xf vector columna de dimension n con la aproximación final.
% iter número de iteraciones que se realizaron.
%
% Andres Cruz y Vera C.U.155899
% Javier Montiel González C.U.159216
%--------------------------------------------------------------------------
% Parámetros

deltamin = 1.e-04;
deltamax = 5; 
delta= 1;
eta = 0.25;
maxiter = 100;    % número máximo de iteraciones externas permitidas
maxregion = 20;    %es el numero máximo que permanece region de confianza 
%en un solo punto
tol = 1.e-06; % tolerancia para la norma del gradiente.

% valores iniciales
iter = 0;        % contador para las iteraciones externas
jregion = 0;     % contador interno un punto xk
iDelta = 0;      % contador interno para un valor delta

g = gradiente(f,x0);
ng = norm(g);
B = hessian(f,x0);
xk = x0;

while(ng>tol && iter<maxiter && jregion<maxregion && iDelta<20)
    pk = doblez(B, g, delta);
    redact = feval(f,xk)-feval(f,xk+pk);
    redpre = -((1/2)*pk'*B*pk+g'*pk);
    rho = redact/redpre;
    
    if rho >= (1 - eta)
        nuevoDelta=2*delta;
        if(nuevoDelta>deltamin && nuevoDelta<deltamax)
            delta = nuevoDelta;
            iDelta=0;
        else
            iDelta= iDelta+1;
        end
        xk = xk + pk;
        jregion = 0;
    elseif eta <= rho && rho < (1-eta)
        xk = xk + pk;
        jregion = 0;
        iDelta=0;
    else
        nuevoDelta=(1/2)*delta;
        if(nuevoDelta>deltamin && nuevoDelta<deltamax)
            delta = nuevoDelta;
            iDelta=0;
        else
            iDelta= 20;
        end
        jregion = jregion + 1;
    end
    
    B = hessian(f, xk);
    g = gradiente(f, xk);
    ng = norm(g);
    iter = iter + 1;
end

xf = xk;
end

