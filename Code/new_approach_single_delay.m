clear
clc
close all

fzero(@single_delay2, 1e10)
% single_delay(1e100)

function root = single_delay(omega)
    %initialize model parameters:
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = mean([0.274, 1.17]);%0.5;
    Nm = mean([1.17, 1.25]); %1.2; 
    k = mean([0.885 0.936]); %0.9;

    Ne = k / (k/Nm + d*ke/de*va/da);
    Ae = va/da*Ne;

    C1 = -k*da + 2*k*Ne / Nm*da + Ne*d*ke/de*va;
    alpha = d*ke*Ae/de;
    C2 = da - k + 2*k*Ne/Nm;
    
    cos_num = -C1 - omega^2*C2 / da;
    cos_denom = alpha*da + alpha*omega^2/da;
    cos_function = cos_num / cos_denom;
    
    sin_function = 1/alpha * (C2*omega + alpha * omega*cos_function);
    root = cos_function^2 + sin_function^2 - 1;
end

function root = single_delay2(omega)
    %initialize model parameters:
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = mean([0.274, 1.17]);%0.5;
    Nm = mean([1.17, 1.25]); %1.2; 
    k = mean([0.885 0.936]); %0.9;

    Ne = k / (k/Nm + d*ke/de*va/da);
    Ae = va/da*Ne;

    C1 = -k*da + 2*k*Ne / Nm*da + Ne*d*ke/de*va;
    alpha = d*ke*Ae/de;
    C2 = da - k + 2*k*Ne/Nm;
    
    root = (da^2 + omega^2)^2 * (-C1*da - omega^2*C2)^2 + 1/(alpha^2 * da^2) * (alpha*da^2 + alpha*omega^2)^2 * (C2*omega*(da^2 + omega^2) - (C1*da*omega - omega^3*C2))^2 - (da^2 + omega^2)^2 * (alpha*da^2 + alpha*omega^2)^2;  
end

