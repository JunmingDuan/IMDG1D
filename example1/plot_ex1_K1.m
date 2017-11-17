clear all;
exact = load('ex1_exact.dat');
numer1 = load('example1_Nx20_K1.dat');
numer2 = load('example1_Nx40_K1.dat');
numer3 = load('example1_Nx80_K1.dat');
numer4 = load('example1_Nx160_K1.dat');
numer5 = load('example1_Nx320_K1.dat');
x1 = numer1(:,1);
x2 = numer2(:,1);
x3 = numer3(:,1);
x4 = numer4(:,1);
x5 = numer5(:,1);
y1 = numer1(:,2);
y2 = numer2(:,2);
y3 = numer3(:,2);
y4 = numer4(:,2);
y5 = numer5(:,2);
plot(exact(:,1), exact(:,2), '-k', x5, y5, 'ro');
ex1 = (sin(4*x1)-8*sin(2*x1)+12*x1)/32;
ex2 = (sin(4*x2)-8*sin(2*x2)+12*x2)/32;
ex3 = (sin(4*x3)-8*sin(2*x3)+12*x3)/32;
ex4 = (sin(4*x4)-8*sin(2*x4)+12*x4)/32;
ex5 = (sin(4*x5)-8*sin(2*x5)+12*x5)/32;
sqrt(norm(ex1-y1,2)^2/size(ex1,1))
sqrt(norm(ex2-y2,2)^2/size(ex2,1))
sqrt(norm(ex3-y3,2)^2/size(ex3,1))
sqrt(norm(ex4-y4,2)^2/size(ex4,1))
sqrt(norm(ex5-y5,2)^2/size(ex5,1))
