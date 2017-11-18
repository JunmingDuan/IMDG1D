function plot_ex1(K, n);
% para: n, P_n polynomial

f = @(x) (sin(4*x)-8*sin(2*x)+12*x)/32;
exact = load('ex1_exact.dat');
err = zeros(1,5);

switch n;
case 20;
  numer1 = load(['example1_Nx20_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,2);
  plot(exact(:,1), exact(:,2), '-k', x1, y1, 'ro');
  ex1 = f(x1);
  err(1) = sqrt(norm(ex1-y1,2)^2/size(ex1,1));
  disp(err);
case 40;
  numer2 = load(['example1_Nx40_K',num2str(K),'.dat']);
  x2 = numer2(:,1); y2 = numer2(:,2);
  plot(exact(:,1), exact(:,2), '-k', x2, y2, 'ro');
  ex2 = f(x2);
  err(2) = sqrt(norm(ex2-y2,2)^2/size(ex2,1));
  fprintf('%.4e\n',err);
case 80;
  numer3 = load(['example1_Nx80_K',num2str(K),'.dat']);
  x3 = numer3(:,1); y3 = numer3(:,2);
  plot(exact(:,1), exact(:,2), '-k', x3, y3, 'ro');
  ex3 = f(x3);
  err(3) = sqrt(norm(ex3-y3,2)^2/size(ex3,1));
  fprintf('%.4e\n',err);
case 160;
  numer4 = load(['example1_Nx160_K',num2str(K),'.dat']);
  x4 = numer4(:,1); y4 = numer4(:,2);
  plot(exact(:,1), exact(:,2), '-k', x4, y4, 'ro');
  ex4 = f(x4);
  err(4) = sqrt(norm(ex4-y4,2)^2/size(ex4,1));
  fprintf('%.4e\n',err);
case 320;
  numer5 = load(['example1_Nx320_K',num2str(K),'.dat']);
  x5 = numer5(:,1); y5 = numer5(:,2);
  plot(exact(:,1), exact(:,2), '-k', x5, y5, 'ro');
  ex5 = f(x5);
  err(5) = sqrt(norm(ex5-y5,2)^2/size(ex5,1));
  fprintf('%.4e\n',err);
case 320;
case 0;
  numer1 = load(['example1_Nx20_K',num2str(K),'.dat']);
  numer2 = load(['example1_Nx40_K',num2str(K),'.dat']);
  numer3 = load(['example1_Nx80_K',num2str(K),'.dat']);
  numer4 = load(['example1_Nx160_K',num2str(K),'.dat']);
  numer5 = load(['example1_Nx320_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,2);
  x2 = numer2(:,1); y2 = numer2(:,2);
  x3 = numer3(:,1); y3 = numer3(:,2);
  x4 = numer4(:,1); y4 = numer4(:,2);
  x5 = numer5(:,1); y5 = numer5(:,2);
  ex1 = f(x1);
  ex2 = f(x2);
  ex3 = f(x3);
  ex4 = f(x4);
  ex5 = f(x5);
  err(1) = sqrt(norm(ex1-y1,2)^2/size(ex1,1));
  err(2) = sqrt(norm(ex2-y2,2)^2/size(ex2,1));
  err(3) = sqrt(norm(ex3-y3,2)^2/size(ex3,1));
  err(4) = sqrt(norm(ex4-y4,2)^2/size(ex4,1));
  err(5) = sqrt(norm(ex5-y5,2)^2/size(ex5,1));
  plot(exact(:,1), exact(:,2), '-k', x1, y1, 'o', x2, y2, '*', x3, y3, '--', x4, ...
  y4, '^', x5, y5, 'v');
  order = [
  log(err(1)/err(2))/log(2);
  log(err(2)/err(3))/log(2);
  log(err(3)/err(4))/log(2);
  log(err(4)/err(5))/log(2)];
  fprintf('%.4e\n',err,order);
  %min(y1)
  %min(y2)
  %min(y3)
  %min(y4)
  %min(y5)
otherwise
  disp('Wrong choice of plot');
end

