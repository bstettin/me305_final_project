d = 4e-3;
ke = 5;
de = 2;
va = 4.8e-7;
da = mean([0.274, 1.17]);%0.5;
Nm = mean([1.17, 1.25]) *1e9; %1.2; 
k = mean([0.885 0.936]); %0.9;

Ae = k*de / (d*ke);
Ne = da*Ae / va;

Ne 

% Ne / Nm