function plot_ex1(K, n);
% para: n, P_n polynomial;

format long ;
f = @(x) (sin(4*x)-8*sin(2*x)+12*x)/32;
exact = load('ex1_exact.dat');
err = zeros(5,3);

switch n;
case 20;
  numer1 = load(['example1_Nx20_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  plot(exact(:,1), exact(:,2), '-k', x1, y1, 'ro');
  ex1 = f(x1);
  err(1,1) = sqrt(sum((ex1-y1).^2.*numer1(:,2)*1/20/2));
  err(1,2) = sum(abs(ex1-y1).*numer1(:,2)*1/20/2);
  for i = 1:size(numer1,1)
    if abs(ex1(i)-y1(i)) > err(1,3)
      err(1,3) = abs(ex1(i)-y1(i));
    end
  end
  disp(err);
case 40;
  numer2 = load(['example1_Nx40_K',num2str(K),'.dat']);
  x2 = numer2(:,1); y2 = numer2(:,3);
  plot(exact(:,1), exact(:,2), '-k', x2, y2, 'ro');
  ex2 = f(x2);
  err(2,1) = sqrt(sum((ex2-y2).^2.*numer2(:,2)*1/40/2));
  err(2,2) = sum(abs(ex2-y2).*numer2(:,2)*1/40/2);
  for i = 1:size(numer2,1)
    if abs(ex2(i)-y2(i)) > err(2,3)
      err(2,3) = abs(ex2(i)-y2(i));
    end
  end
  disp(err);
case 80;
  numer3 = load(['example1_Nx80_K',num2str(K),'.dat']);
  x3 = numer3(:,1); y3 = numer3(:,3);
  plot(exact(:,1), exact(:,2), '-k', x3, y3, 'ro');
  ex3 = f(x3);
  err(3,1) = sqrt(sum((ex3-y3).^2.*numer3(:,2)*1/80/2));
  err(3,2) = sum(abs(ex3-y3).*numer3(:,2)*1/80/2);
  for i = 1:size(numer3,1)
    if abs(ex3(i)-y3(i)) > err(3,3)
      err(3,3) = abs(ex3(i)-y3(i));
    end
  end
  disp(err);
case 160;
  numer4 = load(['example1_Nx160_K',num2str(K),'.dat']);
  x4 = numer4(:,1); y4 = numer4(:,3);
  plot(exact(:,1), exact(:,2), '-k', x4, y4, 'ro');
  ex4 = f(x4);
  err(4,1) = sqrt(sum((ex4-y4).^2.*numer4(:,2)*1/160/2));
  err(4,2) = sum(abs(ex4-y4).*numer4(:,2)*1/160/2);
  for i = 1:size(numer4,1)
    if abs(ex4(i)-y4(i)) > err(4,3)
      err(4,3) = abs(ex4(i)-y4(i));
    end
  end
  disp(err);
case 320;
  numer5 = load(['example1_Nx320_K',num2str(K),'.dat']);
  x5 = numer5(:,1); y5 = numer5(:,3);
  plot(exact(:,1), exact(:,2), '-k', x5, y5, 'ro');
  ex5 = f(x5);
  err(5,1) = sqrt(sum((ex5-y5).^2.*numer5(:,2)*1/320/2));
  err(5,2) = sum(abs(ex5-y5).*numer5(:,2)*1/320/2);
  for i = 1:size(numer5,1)
    if abs(ex5(i)-y5(i)) > err(5,3)
      err(5,3) = abs(ex5(i)-y5(i));
    end
  end
  disp(err);
case 0;
  numer1 = load(['example1_Nx20_K',num2str(K),'.dat']);
  numer2 = load(['example1_Nx40_K',num2str(K),'.dat']);
  numer3 = load(['example1_Nx80_K',num2str(K),'.dat']);
  numer4 = load(['example1_Nx160_K',num2str(K),'.dat']);
  numer5 = load(['example1_Nx320_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  x2 = numer2(:,1); y2 = numer2(:,3);
  x3 = numer3(:,1); y3 = numer3(:,3);
  x4 = numer4(:,1); y4 = numer4(:,3);
  x5 = numer5(:,1); y5 = numer5(:,3);
  ex1 = f(x1);
  ex2 = f(x2);
  ex3 = f(x3);
  ex4 = f(x4);
  ex5 = f(x5);
  err(1,1) = sqrt(sum((ex1-y1).^2.*numer1(:,2)*1/20/2));
  err(1,2) = sum(abs(ex1-y1).*numer1(:,2)*1/20/2);
  for i = 1:size(numer1,1)
    if abs(ex1(i)-y1(i)) > err(1,3)
      err(1,3) = abs(ex1(i)-y1(i));
    end
  end
  err(2,1) = sqrt(sum((ex2-y2).^2.*numer2(:,2)*1/40/2));
  err(2,2) = sum(abs(ex2-y2).*numer2(:,2)*1/40/2);
  for i = 1:size(numer2,1)
    if abs(ex2(i)-y2(i)) > err(2,3)
      err(2,3) = abs(ex2(i)-y2(i));
    end
  end
  err(3,1) = sqrt(sum((ex3-y3).^2.*numer3(:,2)*1/80/2));
  err(3,2) = sum(abs(ex3-y3).*numer3(:,2)*1/80/2);
  for i = 1:size(numer3,1)
    if abs(ex3(i)-y3(i)) > err(3,3)
      err(3,3) = abs(ex3(i)-y3(i));
    end
  end
  err(4,1) = sqrt(sum((ex4-y4).^2.*numer4(:,2)*1/160/2));
  err(4,2) = sum(abs(ex4-y4).*numer4(:,2)*1/160/2);
  for i = 1:size(numer4,1)
    if abs(ex4(i)-y4(i)) > err(4,3)
      err(4,3) = abs(ex4(i)-y4(i));
    end
  end
  err(5,1) = sqrt(sum((ex5-y5).^2.*numer5(:,2)*1/320/2));
  err(5,2) = sum(abs(ex5-y5).*numer5(:,2)*1/320/2);
  for i = 1:size(numer5,1)
    if abs(ex5(i)-y5(i)) > err(5,3)
      err(5,3) = abs(ex5(i)-y5(i));
    end
  end
  plot(exact(:,1), exact(:,2), '-k', x1, y1, 'o', x2, y2, '*', x3, y3, '--', x4, ...
  y4, '^', x5, y5, 'v');
  %cal_order
  min_y = [min(y1);min(y2);min(y3);min(y4);min(y5)];
  order = [
  zeros(1,3);
  log2(err(1,:)./err(2,:));
  log2(err(2,:)./err(3,:));
  log2(err(3,:)./err(4,:));
  log2(err(4,:)./err(5,:))];
  N = [20;40;80;160;320];
  for n = 1:5
    fprintf('%3d ', N(n));
    for i = 1:3
      fprintf('%.3e %.2f ', err(n,i), order(n,i));
    end
    fprintf('%.3e\n', min_y(n));
  end

otherwise
  disp('Wrong choice of plot');
end

