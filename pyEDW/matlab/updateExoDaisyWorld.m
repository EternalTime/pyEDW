function xNew = updateExoDaisyWorld(x,dt,theta)
%UPDATEEXODAISYWORLD timestep for the exo-Daisy World model
%
%   x (DoF Vector)
%    1:     f_B
%    2:     f_W
%    3:     T
%    4:     L
%   
%   theta (Parameter Vector)
%    1:     f
%    2:     1 - A_G
%    3:     A_B - A_G
%    4:     A_W - A_G
%    5:     Q
%    6:     gamma_D / gamma_G
%    7:     1 / gamma_G tau_E
%    8:     1 / gamma_G tau_S
%    9:     8 (T_opt / Delta T)^4
%    10:    delta
%    11:    lambda
%

    Z = randn();
    S = randsample([-1 1],1);
    
    [a,b] = EoM(x   , theta);
    k1 = real(a*dt + (Z-S)*b*sqrt(dt));
    
    [a,b] = EoM(x+k1, theta);
    k2 = real(a*dt + (Z+S)*b*sqrt(dt));

    xNew = x + 0.5*(k1 + k2);
    xNew(1:2) = max(xNew(1:2),0);


    function [a,b] = EoM(x,theta)
        dAf  =  theta(3)*x(1) + theta(4)*x(2);
        Tg4  =  x(3)^4 + theta(5)*dAf;
        df   =  theta(1) - x(1) - x(2);

        a = [  ...
                df*W((Tg4-theta(5)*theta(3))^0.25 - 1, theta(9))*x(1) - theta(6)*x(1) , ...
                df*W((Tg4-theta(5)*theta(4))^0.25 - 1, theta(9))*x(2) - theta(6)*x(2) , ...
                theta(7)*((1 - dAf/theta(2))*x(4) - x(3)^4) , ...
                theta(8)*(1 + theta(11) - x(4)) ...
            ];
        b = [  ...
                0 , ...
                0 , ...
                0 , ...
                sqrt(2*theta(8))*theta(10)*(1 + theta(11)) ...
            ];
    end

    function g = W(T,alpha)
        g = exp(-alpha*T^4);
    end


end