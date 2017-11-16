exact = load('ex1_exact.dat');
numer = load('example1.dat');
plot(exact(:,1), exact(:,2), '-k', numer(:,1), numer(:,2), 'ro');
