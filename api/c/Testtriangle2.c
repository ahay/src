#include <stdio.h>
#include <stdlib.h>

#include "triangle2.h"
#include "dottest.h"

int main(void) {
    double dot1[2], dot2[2];
    static int n1=100, n2=100, f1=10, f2=5;

    sf_triangle2_init(f1, f2, n1, n2, 3);
    sf_dot_test(sf_triangle2_lop, n1*n2, n1*n2, dot1, dot2);
    sf_triangle2_close();

    printf ("%12.3f ? %12.3f\n",dot1[0],dot1[1]);
    printf ("%12.3f ? %12.3f\n",dot2[0],dot2[1]);

    exit(0);
}
