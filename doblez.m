function [ps] = doblez(B, g, Delta)
%T�cnica del doblez para el problema cuadr�tico de regi�n de confianza
% Min (1/2)*p'*B*p + g'*p +f(x)
% s.a. || p || <= Delta

pN = -B\g; %direcci�n de Newton
pC = -((g'*g)/(g'*B*g))*g; %punto de Cauchy


if norm(pN) <= Delta
    %Caso i)
    ps = pN;
else
    %Caso ii)
    %Verificamos si la norma de pC es mayor o igual a delta
    if norm(pC) >= Delta
        ps = (Delta/norm(pC))*pC;
    else
        % Ecuaci�n de 2 grado en la variable t
        % || pC + t (pN - pC)||^2 - (Delta)^2 = 0
        % (pC + t*u)'(pC + t*u) - (Delta)^2 = 0 / u= pN - pC.
        % a*(t*t) + b*t + c = 0
        u = pN - pC;
        a = u'*u; b= 2*pC'*u; c= pC'*pC-Delta^2;
        t = roots([a b c]);
        %Tomamos la ra�z positiva
        if t(1)>0
            ts = t(1);
        else
            ts = t(2);
        end
        ps = pC + ts*u;
    end
end

end
