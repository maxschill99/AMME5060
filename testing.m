clc
clear

il = 1
ih = 10
dx = 1/9;

for i = il:ih
    x(i) = (i-1)*dx;
    T(1,i) = sin(pi*x(i));
end 

% fprintf('%.8f', pi)




