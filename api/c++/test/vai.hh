class VAI {
//    int ii;
public:
    // constructor
    VAI(int,int);
   // destructor
    ~VAI(){};

    int operator() (int,int);
    
private:
    int n[2];
};

