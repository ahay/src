class VAI {
    int ii;
public:
    // constructor
    VAI(int m0, int m1) {
	nd=2;
	n= new int(nd);
	n[0]=m0;
	n[1]=m1;
    }
   // destructor
    ~VAI(){};

    int & operator() (int i0,int i1) {
	ii = (i1-1)*n[0] + i0 - 1;
	return(ii);
    }
    
//    int & operator() (int, int);

private:
    int *n;
    int nd;
};

