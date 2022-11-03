clc
clear

nx = 8;
ny = nx;
ndx = 1/(nx-1);
ndy = ndx;
nalpha = 1e-4;
dtold = 25.5;
ndt = 1/ceil(1/dtold)
nCFL = 0.25;

syms CFL dx alpha dt

eqn = CFL == (1/(dx^2) + 1/(dx^2))*alpha*dt;

% THIS IS THE VALUE THAT YOU ARE CALCULATING FROM THE SOLVER
% NEED TO UPDATE THE CORRESPNDING VALUES ABOVE
%%%%%%%%%%%%%%%%%%
% input = dt;
input = CFL;
%%%%%%%%%%%%%%%%%%

if input== dt
    eqn = subs(eqn,CFL,nCFL);
    eqn = subs(eqn,dx,ndx);
    eqn = subs(eqn,alpha,nalpha);
elseif input==CFL
    eqn = subs(eqn,dx,ndx);
    eqn = subs(eqn,alpha,nalpha);
    eqn = subs(eqn,dt,ndt);
end

S = solve(eqn,input);
sol = vpa(S,3)