function [xf, iter, deltas, puntos] = miregion(f, x0)
% M�todo de regi�n de confianza para aproximar un m�nimo
% de f: Rn->R dos veces continuamente diferenciable

% Entradas
% f cadena de caracteres con con el nombre del c�digo en Matlab de 
% la funci�n a minimizar.
% x0 vector columna de dimensi�n n con el punto inicial.

% Salidas
% xf vector columna de dimension n con la aproximaci�n final.
% iter n�mero de iteraciones que se realizaron.
%
% Andres Cruz y Vera C.U.155899
% Javier Montiel Gonz�lez C.U.159216
%--------------------------------------------------------------------------
% Par�metros

deltamin = 1.e-04;
deltamax = 5; 
delta= 1;
eta = 0.25;        % Constate que determina si se acepta el paso
maxiter = 100;     % n�mero m�ximo de iteraciones externas permitidas
maxregion = 20;    %es el numero m�ximo que permanece region de confianza 
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

%Vectores auxiliares para guardar el resultado de cada iteraci�n
deltas=[];
puntos=[];

%Definimos los criterios de parada del m�todo
while(ng>tol && iter<maxiter && jregion<maxregion && iDelta<20)
    %Almacenamos el valor de delta y xk en cada iteraci�n
    deltas=[deltas; delta];
    puntos=[puntos; xk'];
    %El paso que se obtiene por el m�todo del doblez
    pk = doblez(B, g, delta);
    %Calculamos la recduci�n actual 
    redact = feval(f,xk)-feval(f,xk+pk);
    %Calculamos la reducci�n que se predice
    redpre = -((1/2)*pk'*B*pk+g'*pk);
    %Definimos el coeficiente de la reducci�n actual entre la reducci�n 
    %que se predice
    rho = redact/redpre;
    
    %Evaluamos los tres posibles casos del valor de rho
    if rho >= (1 - eta)
        %Se acepta el paso y puede aumentar el tama�o de delta
        nuevoDelta=2*delta;
        %Verificamos que delta no rebase deltamax
        if(nuevoDelta<deltamax)
            delta = nuevoDelta;
            iDelta=0;
        else
        %En caso que rebase deltamax, el valor de delta no cambia
            iDelta= iDelta+1;
        end
        xk = xk + pk;
        jregion = 0;
    elseif eta <= rho && rho < (1-eta)
        %Se acepta el paso pero no cambia delta
        xk = xk + pk;
        jregion = 0;
        iDelta=0;
    else
        %Se rechaza el paso y se puede reducir delta
        nuevoDelta=(1/2)*delta;
        if(nuevoDelta>deltamin)
            delta = nuevoDelta;
            iDelta=0;
        else
            iDelta= 20;
        end
        jregion = jregion + 1;
    end
    
    %Calculamos los valores de la matriz hessiana y el vector gradiente en 
    %el nuevo punto
    B = hessian(f, xk);
    g = gradiente(f, xk);
    ng = norm(g);
    iter = iter + 1;
end

xf = xk;
end

