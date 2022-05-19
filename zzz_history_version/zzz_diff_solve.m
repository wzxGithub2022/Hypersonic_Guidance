clear;clc;

syms x
dx(x) = diff(x*x,x);

value = double(dx(3));

