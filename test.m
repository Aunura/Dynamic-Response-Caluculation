clear all;clc;close all;
delta_t=1;
gamma = 0.5;
beta = 0.25;

am1 = 1/beta/delta_t^2;
am2 = 1/beta/delta_t;
am3 = 1/2/beta;

ac1 = gamma/beta/delta_t;
ac2 = gamma/beta;
ac3 = (gamma/2/beta-1)*delta_t;