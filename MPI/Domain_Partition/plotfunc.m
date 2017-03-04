clear all; close all; clc;
m1 = load('mandelbrotpart12.dat');
m2 = load('mandelbrotpart24.dat');
%m3 = load('mandelbrot36.dat');
m3 = load('mandelbrotpart48.dat');
m4 = load('mandelbrotpart60.dat');
m5 = load('mandelbrotpart96.dat');
m6 = load('mandelbrotpart120.dat');

m = m1;
m(2,:) = m2;
m(3,:) = m3;
m(4,:) = m4;
m(5,:) = m5;
m(6,:) = m6;
%m(7,:) = m7;

%Time required to run the serial code
t = 317.86;

for i = 1:6
    m(i,4) = t/m(i,2);
end

%Calculating the Karp-Flatt Metric
for i = 1:6
    m(i,5) = ((1/m(i,4))-(1/m(i,1)))/(1-(1/m(i,1)));
end

figure(1);
plot(m(:,1), m(:,4), '-b', 'LineWidth', 2);
title('Speedup Factor vs Number of Processors, strong scaling for Mandelbrot Set');
xlabel('Number of processors');
ylabel('Speedup Factor');
axis square;

figure(2);
plot(m(:,1), m(:,5), '-g', 'LineWidth',2);
title('Karp-Flatt Serial Fraction plot - Mandelbrot set with partition');
xlabel('Number of processors');
ylabel('Karp-Flatt serial metric');
axis square;