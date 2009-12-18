#include <petscksp.h>

int main(int argc,char* argv[]) {
    KSP Solver; Mat A;
    MPI_Comm comm = MPI_COMM_WORLD;
    PetscInitialize (&argc, &argv, 0, 0);
    MatCreate (comm, &A);
    KSPCreate (comm, &Solver);
    KSPSetType (Solver, KSPGMRES);
    PetscFinalize ();
    return 0;
}
