
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% For further referecens see the main file runGRWunconf.m



function [x]=solveEQ_CPU(a,b1,b2,b3,b4,x1,x2,x3,x4,rhs)

% Solve the equation:
% a.*x + b1.*x1 + b2.*x2 + b3.*x3 + b4.*x4 = rhs
% for x
x=(rhs-b1.*x1-b2.*x2-b3.*x3-b4.*x4)./a;


