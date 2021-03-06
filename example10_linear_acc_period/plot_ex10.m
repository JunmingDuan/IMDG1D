function plot_ex9(K, n, limiter);
% para: n, P_n polynomial;

addpath('../src/');
format long ;
err = zeros(5,3);
t_end = 0.5;

switch n;
case 20;
  numer1 = load(['example10_Nx20_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  ex1 = exact(x1, t_end);
  plot(x1, ex1, '-k', x1, y1, 'ro');
  nx = 20; h = 1/nx;
  err(1,1) = cal_norm(ex1-y1, numer1(:,2), h, 2);
  err(1,2) = cal_norm(ex1-y1, numer1(:,2), h, 1);
  err(1,3) = cal_norm(ex1-y1, numer1(:,2), h, 'inf');
  disp(err);
case 40;
  numer2 = load(['example10_Nx40_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  x2 = numer2(:,1); y2 = numer2(:,3);
  ex2 = exact(x2, t_end);
  plot(x2, ex2, '-k', x2, y2, 'ro');
  nx = 40; h = 1/nx;
  err(2,1) = cal_norm(ex2-y2, numer2(:,2), h, 2);
  err(2,2) = cal_norm(ex2-y2, numer2(:,2), h, 1);
  err(2,3) = cal_norm(ex2-y2, numer2(:,2), h, 'inf');
  disp(err);
case 80;
  numer3 = load(['example10_Nx80_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  x3 = numer3(:,1); y3 = numer3(:,3);
  ex3 = exact(x3, t_end);
  plot(x3, ex3, '-k', x3, y3, 'ro');
  nx = 80; h = 1/nx;
  err(3,1) = cal_norm(ex3-y3, numer3(:,2), h, 2);
  err(3,2) = cal_norm(ex3-y3, numer3(:,2), h, 1);
  err(3,3) = cal_norm(ex3-y3, numer3(:,2), h, 'inf');
  disp(err);
case 160;
  numer4 = load(['example10_Nx160_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  x4 = numer4(:,1); y4 = numer4(:,3);
  ex4 = exact(x4, t_end);
  plot(x4, ex4, '-k', x4, y4, 'ro');
  nx = 160; h = 1/nx;
  err(4,1) = cal_norm(ex4-y4, numer4(:,2), h, 2);
  err(4,2) = cal_norm(ex4-y4, numer4(:,2), h, 1);
  err(4,3) = cal_norm(ex4-y4, numer4(:,2), h, 'inf');
  disp(err);
case 320;
  numer5 = load(['example10_Nx320_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  x5 = numer5(:,1); y5 = numer5(:,3);
  ex5 = exact(x5, t_end);
  plot(x5, ex5, '-k', x5, y5, 'ro');
  nx = 320; h = 1/nx;
  err(5,1) = cal_norm(ex5-y5, numer5(:,2), h, 2);
  err(5,2) = cal_norm(ex5-y5, numer5(:,2), h, 1);
  err(5,3) = cal_norm(ex5-y5, numer5(:,2), h, 'inf');
  disp(err);
case 0;
  numer1 = load(['example10_Nx20_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  numer2 = load(['example10_Nx40_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  numer3 = load(['example10_Nx80_K',num2str(K), '_PP',num2str(limiter),'.dat']);
  numer4 = load(['example10_Nx160_K',num2str(K),'_PP',num2str(limiter),'.dat']);
  numer5 = load(['example10_Nx320_K',num2str(K),'_PP',num2str(limiter),'.dat']);
  %numer6 = load(['example10_Nx640_K',num2str(K),'_PP',num2str(limiter),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  x2 = numer2(:,1); y2 = numer2(:,3);
  x3 = numer3(:,1); y3 = numer3(:,3);
  x4 = numer4(:,1); y4 = numer4(:,3);
  x5 = numer5(:,1); y5 = numer5(:,3);
  %x6 = numer6(:,1); y6 = numer6(:,3);
  ex1 = exact(x1, t_end);
  ex2 = exact(x2, t_end);
  ex3 = exact(x3, t_end);
  ex4 = exact(x4, t_end);
  ex5 = exact(x5, t_end);
  plot(x4, ex4, '-k', x1, y1, 'o', x2, y2, '*', x3, y3, '--', x4, ...
  y4, '^', x5, y5, 'v', 'MarkerSize', 1);%, x6, y6, '<');
  %cal_order
  nx = 20; h = 1/nx;
  err(1,1) = cal_norm(ex1-y1, numer1(:,2), h, 2);
  err(1,2) = cal_norm(ex1-y1, numer1(:,2), h, 1);
  err(1,3) = cal_norm(ex1-y1, numer1(:,2), h, 'inf');
  nx = 40; h = 1/nx;
  err(2,1) = cal_norm(ex2-y2, numer2(:,2), h, 2);
  err(2,2) = cal_norm(ex2-y2, numer2(:,2), h, 1);
  err(2,3) = cal_norm(ex2-y2, numer2(:,2), h, 'inf');
  nx = 80; h = 1/nx;
  err(3,1) = cal_norm(ex3-y3, numer3(:,2), h, 2);
  err(3,2) = cal_norm(ex3-y3, numer3(:,2), h, 1);
  err(3,3) = cal_norm(ex3-y3, numer3(:,2), h, 'inf');
  nx = 160; h = 1/nx;
  err(4,1) = cal_norm(ex4-y4, numer4(:,2), h, 2);
  err(4,2) = cal_norm(ex4-y4, numer4(:,2), h, 1);
  err(4,3) = cal_norm(ex4-y4, numer4(:,2), h, 'inf');
  nx = 320; h = 1/nx;
  err(5,1) = cal_norm(ex5-y5, numer5(:,2), h, 2);
  err(5,2) = cal_norm(ex5-y5, numer5(:,2), h, 1);
  err(5,3) = cal_norm(ex5-y5, numer5(:,2), h, 'inf');
  min_y = [min(y1);min(y2);min(y3);min(y4);min(y5)];
  %min_y = [min(y1);min(y2);min(y3);min(y4);];%min(y5)];
  order = [
  zeros(1,3);
  log2(err(1,:)./err(2,:));
  log2(err(2,:)./err(3,:));
  log2(err(3,:)./err(4,:));
  log2(err(4,:)./err(5,:))];
  N = [20;40;80;160;320];
  diary table2.dat
  diary on;
  for n = 1:5
    fprintf('%3d ', N(n));
    for i = 1:3
      fprintf('%.3e %.2f ', err(n,i), order(n,i));
    end
    fprintf('%.3e\n', min_y(n));
  end
  axis([0,2*pi,-0.1,2]);
  diary off;

otherwise
  disp('Wrong choice of plot');
end
end

function u = exact(x, t)
u = 1+sin(x-t);
end

