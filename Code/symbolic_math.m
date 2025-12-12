clear
clc
close all

% Symbolic sanity check


x0 = [0, 10];
equations(x0)
% x = fsolve(@equations, x0)
% disp(['Omega = ' num2str(x(1)) ' rad'])
% disp(['  Tau = ' num2str(x(2)) ' hrs'])

x = [1 1];
equations(x)

function roots = equations(x)
    syms k N_m d k_e d_e v_a d_a N_e A_e w t real

    J0 = [k - 2*N_e/N_m, -d*k_e/d_e*N_e;
          v_a, -d_a];
    
    Jt = [-d*k_e/d_e*A_e 0; 0 0];

    determinant = det(J0 + Jt*(cos(w*t) - 1i*sin(w*t)) - eye(2)*(1i*w));
    
    eq1 = real(determinant);
    eq2 = imag(determinant);
    
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = 0.5;
    Nm = 1.2; 
    k = 0.9;
    
    %calculate equilibrium points
    Ne = k / (k/Nm + d*ke/de*va/da);
    Ae = va/da*Ne;
    
    eq1_subs = subs(eq1, {d, k_e, d_e, v_a, d_a, N_m, k, N_e, A_e}, {d, ke, de, va, da, Nm, k, Ne, Ae});
    eq2_subs = subs(eq2, {d, k_e, d_e, v_a, d_a, N_m, k, N_e, A_e}, {d, ke, de, va, da, Nm, k, Ne, Ae});
    
    disp(eq1_subs)
    % roots(1) = double(eq1_subs);
    % roots(2) = double(eq2_subs);
    roots = zeros(2,1)
end


