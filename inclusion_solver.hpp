#ifndef INCLUSION_SOLVER
#define INCLUSION_SOLVER

/* #include "../common/mesh_extras.hpp" */
/* #include "electromagnetics.hpp" */

#include "pfem_extras.hpp"

#include <string>
#include <map>


using namespace mfem;
using namespace std;

class H1_ParFESpace; // * H1FESpace_;    // Continuous space for phi
class ND_ParFESpace; // * HCurlFESpace_; // Tangentially continuous space for E
class RT_ParFESpace; // * HDivFESpace_;  // Normally continuous space for D
class L2_ParFESpace; // * L2FESpace_;    // Discontinuous space for rho

class InclusionSolver
{
public:
   InclusionSolver(ParMesh & pmesh, int order,
               Array<int> & dbcs,
               Coefficient & epsCoef,
               double (*phi_bc )(const Vector&));
   ~InclusionSolver();

   HYPRE_BigInt GetProblemSize();

   void PrintSizes();

   void Assemble();

   void Update();

   void Solve();

   void GetErrorEstimates(Vector & errors);

   void RegisterVisItFields(VisItDataCollection & visit_dc);

   void WriteVisItFields(int it = 0);

   void InitializeGLVis();

   void DisplayToGLVis();

   const ParGridFunction & GetVectorPotential() { return *phi_; }

   void WriteToVtk(const char* name, int res=1);

   ParGridFunction * e_;         // Electric Field
private:

   int myid_;      // Local processor rank
   int num_procs_; // Number of processors
   int order_;     // Basis function order

   ParMesh * pmesh_;

   Array<int> * dbcs_; // Dirichlet BC Surface Attribute IDs

   VisItDataCollection * visit_dc_; // To prepare fields for VisIt viewing

   H1_ParFESpace * H1FESpace_;    // Continuous space for phi
   ND_ParFESpace * HCurlFESpace_; // Tangentially continuous space for E
   RT_ParFESpace * HDivFESpace_;  // Normally continuous space for D
   L2_ParFESpace * L2FESpace_;    // Discontinuous space for rho

   ParBilinearForm * divEpsGrad_; // Laplacian operator
   ParBilinearForm * h1Mass_;     // For Volumetric Charge Density Source
   ParBilinearForm * h1SurfMass_; // For Surface Charge Density Source
   ParBilinearForm * hDivMass_;   // For Computing D from E

   ParMixedBilinearForm * hCurlHDivEps_; // For computing D from E
   ParMixedBilinearForm * hCurlHDiv_;    // For computing D from E and P
   ParMixedBilinearForm * weakDiv_;      // For computing the source term from P

   ParLinearForm * rhod_; // Dual of Volumetric Charge Density Source

   ParLinearForm * l2_vol_int_;  // Integral of L2 field
   ParLinearForm * rt_surf_int_; // Integral of H(Div) field over boundary

   ParDiscreteGradOperator * grad_; // For Computing E from phi
   ParDiscreteDivOperator  * div_;  // For Computing rho from D

   ParGridFunction * phi_;       // Electric Scalar Potential
   ParGridFunction * rho_;       // Volumetric Charge Density (Div(D))
   //ParGridFunction * e_;         // Electric Field
   ParGridFunction * d_;         // Electric Flux Density (aka Dielectric Flux)

   ConstantCoefficient oneCoef_;   // Coefficient equal to 1
   Coefficient       * epsCoef_;   // Dielectric Permittivity Coefficient
   Coefficient       * phiBCCoef_; // Scalar Potential Boundary Condition
   Coefficient       * rhoCoef_;   // Charge Density Coefficient
   VectorCoefficient * pCoef_;     // Polarization Vector Field Coefficient

   // Source functions
   double (*phi_bc_func_ )(const Vector&);          // Scalar Potential BC


   std::map<std::string,socketstream*> socks_; // Visualization sockets

   Array<int> ess_bdr_, ess_bdr_tdofs_; // Essential Boundary Condition DoFs
};


#endif // INCLUSION_SOLVER
