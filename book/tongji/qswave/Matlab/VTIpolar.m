close;figure;
N=8;

hold on;

polarization_compare(N, 2, 1, 0.3, 0.1, 0, 'b');

axis([-3.14 3.14 -3.14 3.14])

set(gca,'xtick',[-3 -2 -1 0 1 2 3])


print -depsc junk_ml.eps;
%quit
