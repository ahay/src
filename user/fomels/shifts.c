void shifts2(int s1, int s2, int n1, int n2, float **inp, float ***sft)
/*< generate shifts >*/
{
    int i1, i2, j1, j2, i, k1, k2, l1, l2;
    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i=0;
	    for (j2=0; j2 < s2; j2++) {
		for (j1=0; j1 < s1; j1++) {
		    if (i==0) {
			sft[i][i2][i1] = inp[i2][i1];
		    } else {			
			k2 = i2+j2;
			k1 = i1+j1;
			l2 = i2-j2;
			l1 = i1-j1;
			if (l1 < 0 || l2 < 0) {
			    sft[i][i2][i1] = inp[k2][k1];
			} else if (k1 >= n1 || k2 >= n2) {
			    sft[i][i2][i1] = inp[l2][l1];
			} else {
			    sft[i][i2][i1] = inp[k2][k1] + inp[l2][l1];
			}
		    }
		    i++;
		}
	    }
	}
    }
}
