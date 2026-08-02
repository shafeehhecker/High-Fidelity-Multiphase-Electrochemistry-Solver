// petsc_solver.hpp
class PETScNonlinearSolver {
private:
    MPI_Comm comm;
    DM dm;                  // DMPlex mesh manager
    SNES snes;              // Nonlinear solver
    Vec u;                  // Solution vector
    Mat J;                  // Jacobian matrix
    Vec F;                  // Residual vector
    
    // Physical parameters
    double kappa;           // Electrolyte conductivity
    double j0;              // Exchange current density
    double alpha_a, alpha_c;// Transfer coefficients
    double phi_s;           // Electrode potential
    double E_eq;            // Equilibrium potential
    
public:
    PETScNonlinearSolver(MPI_Comm comm);
    ~PETScNonlinearSolver();
    
    // Core solver methods
    void setMesh(DM dm_in);
    void setParameters(const SolverParams& params);
    void setBoundaryConditions(const BCInfo& bc);
    void solve(Vec& solution);
    
    // PETSc callback functions (static)
    static PetscErrorCode FormFunction(SNES snes, Vec u, Vec F, void* ctx);
    static PetscErrorCode FormJacobian(SNES snes, Vec u, Mat J, Mat P, void* ctx);
    
    // Assembly helpers
    void assembleResidual(Vec u, Vec F);
    void assembleJacobian(Vec u, Mat J);
};
