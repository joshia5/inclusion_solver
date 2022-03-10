
#include "inclusion_solver.hpp"
#include <fstream>
#include <iostream>

// Permittivity Functions
Coefficient * SetupPermittivityCoefficient();

static Vector pw_eps_(0);     // Piecewise permittivity values
static Vector ds_params_(0);  // Center, Radius, and Permittivity
//                               of dielectric sphere
double dielectric_sphere(const Vector &);

// Charge Density Function
static Vector cs_params_(0);  // Center, Radius, and Total Charge
//                               of charged sphere
double charged_sphere(const Vector &);

// Point Charges
static Vector pc_params_(0); // Point charge locations and charges

// Polarization
static Vector vp_params_(0);  // Axis Start, Axis End, Cylinder Radius,
//                               and Polarization Magnitude
void voltaic_pile(const Vector &, Vector &);

// Phi Boundary Condition
static Vector e_uniform_(0);
double phi_bc_uniform(const Vector &);

// Prints the program's logo to the given output stream
void display_banner(ostream & os);

int main(int argc, char *argv[])
{
   MPI_Session mpi(argc, argv);


   // Parse command-line options.
   const char *mesh_file = "";
   int order = 1;
   int maxit = 100;
   int parallel_ref_levels = 0;
   bool visualization = true;
   bool visit = true;

   Array<int> dbcs;
   Array<int> nbcs;

   Vector dbcv;
   Vector nbcv;

   bool dbcg = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&parallel_ref_levels, "-rp", "--parallel-ref-levels",
                  "Number of parallel refinement levels.");
   args.AddOption(&e_uniform_, "-uebc", "--uniform-e-bc",
                  "Specify if the three components of the constant "
                  "electric field");
   args.AddOption(&pw_eps_, "-pwe", "--piecewise-eps",
                  "Piecewise values of Permittivity");
   args.AddOption(&ds_params_, "-ds", "--dielectric-sphere-params",
                  "Center, Radius, and Permittivity of Dielectric Sphere");
   args.AddOption(&cs_params_, "-cs", "--charged-sphere-params",
                  "Center, Radius, and Total Charge of Charged Sphere");
   args.AddOption(&pc_params_, "-pc", "--point-charge-params",
                  "Charges and locations of Point Charges");
   args.AddOption(&vp_params_, "-vp", "--voltaic-pile-params",
                  "Axis End Points, Radius, and "
                  "Polarization of Cylindrical Voltaic Pile");
   args.AddOption(&dbcs, "-dbcs", "--dirichlet-bc-surf",
                  "Dirichlet Boundary Condition Surfaces");
   args.AddOption(&dbcv, "-dbcv", "--dirichlet-bc-vals",
                  "Dirichlet Boundary Condition Values");
   args.AddOption(&dbcg, "-dbcg", "--dirichlet-bc-gradient",
                  "-no-dbcg", "--no-dirichlet-bc-gradient",
                  "Dirichlet Boundary Condition Gradient (phi = -z)");
   args.AddOption(&nbcs, "-nbcs", "--neumann-bc-surf",
                  "Neumann Boundary Condition Surfaces");
   args.AddOption(&nbcv, "-nbcv", "--neumann-bc-vals",
                  "Neumann Boundary Condition Values");
   args.AddOption(&maxit, "-maxit", "--max-amr-iterations",
                  "Max number of iterations in the main AMR loop.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
                  "Enable or disable VisIt visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root())
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (mpi.Root())
   {
      args.PrintOptions(cout);
   }

   // Read the (serial) mesh from the given mesh file on all processors.  We
   // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   // and volume meshes with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int sdim = mesh->SpaceDimension();

   if (mpi.Root())
   {
      cout << "Starting initialization." << endl;
   }


   // Define a parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh pmesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   // Refine this mesh in parallel to increase the resolution.
   int par_ref_levels = parallel_ref_levels;
   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh.UniformRefinement();
   }
   // Make sure tet-only meshes are marked for local refinement.
   pmesh.Finalize(true);

   // If the gradient bc was selected but the E field was not specified
   // set a default vector value.
   if ( dbcg && e_uniform_.Size() != sdim )
   {
      e_uniform_.SetSize(sdim);
      e_uniform_ = 0.0;
      e_uniform_(sdim-1) = 1.0;
   }

   // If values for Dirichlet BCs were not set assume they are zero
   if ( dbcv.Size() < dbcs.Size() && !dbcg )
   {
      dbcv.SetSize(dbcs.Size());
      dbcv = 0.0;
   }

   // If values for Neumann BCs were not set assume they are zero
   if ( nbcv.Size() < nbcs.Size() )
   {
      nbcv.SetSize(nbcs.Size());
      nbcv = 0.0;
   }

   // Create a coefficient describing the dielectric permittivity
   Coefficient * epsCoef = SetupPermittivityCoefficient();

   // Create the Electrostatic solver
   InclusionSolver Volta(pmesh, order, dbcs, dbcv, nbcs, nbcv, *epsCoef,
                     ( e_uniform_.Size() > 0 ) ? phi_bc_uniform    : NULL,
                     ( cs_params_.Size() > 0 ) ? charged_sphere    : NULL,
                     ( vp_params_.Size() > 0 ) ? voltaic_pile      : NULL,
                     pc_params_);

   // Initialize GLVis visualization
   if (visualization)
   {
      Volta.InitializeGLVis();
   }

   // Initialize VisIt visualization
   VisItDataCollection visit_dc("Volta-AMR-Parallel", &pmesh);

   if ( visit )
   {
      Volta.RegisterVisItFields(visit_dc);
   }
   if (mpi.Root()) { cout << "Initialization done." << endl; }

   // The main AMR loop. In each iteration we solve the problem on the current
   // mesh, visualize the solution, estimate the error on all elements, refine
   // the worst elements and update all objects to work with the new mesh.  We
   // refine until the maximum number of dofs in the nodal finite element space
   // reaches 10 million.
   const int max_dofs = 10000000;
   for (int it = 1; it <= maxit; it++)
   {
      if (mpi.Root())
      {
         cout << "\nAMR Iteration " << it << endl;
      }

      // Display the current number of DoFs in each finite element space
      Volta.PrintSizes();

      // Assemble all forms
      Volta.Assemble();

      // Solve the system and compute any auxiliary fields
      Volta.Solve();

      // Determine the current size of the linear system
      int prob_size = Volta.GetProblemSize();

      // Write fields to disk for VisIt
      if ( visit )
      {
         Volta.WriteVisItFields(it);
      }

      // Send the solution by socket to a GLVis server.
      if (visualization)
      {
         Volta.DisplayToGLVis();
      }

      if (mpi.Root())
      {
         cout << "AMR iteration " << it << " complete." << endl;
      }

      // Check stopping criteria
      if (prob_size > max_dofs)
      {
         if (mpi.Root())
         {
            cout << "Reached maximum number of dofs, exiting..." << endl;
         }
         break;
      }
      if (it == maxit)
      {
         break;
      }

      // Wait for user input. Ask every 10th iteration.
      char c = 'c';
      if (mpi.Root() && (it % 10 == 0))
      {
         cout << "press (q)uit or (c)ontinue --> " << flush;
         cin >> c;
      }
      MPI_Bcast(&c, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

      if (c != 'c')
      {
         break;
      }

      // Estimate element errors using the Zienkiewicz-Zhu error estimator.
      Vector errors(pmesh.GetNE());
      Volta.GetErrorEstimates(errors);

      double local_max_err = errors.Max();
      double global_max_err;
      MPI_Allreduce(&local_max_err, &global_max_err, 1,
                    MPI_DOUBLE, MPI_MAX, pmesh.GetComm());

      // Refine the elements whose error is larger than a fraction of the
      // maximum element error.
      const double frac = 0.7;
      double threshold = frac * global_max_err;
      if (mpi.Root()) { cout << "Refining ..." << endl; }
      pmesh.RefineByError(errors, threshold);

      // Update the electrostatic solver to reflect the new state of the mesh.
      Volta.Update();

      if (pmesh.Nonconforming() && mpi.WorldSize() > 1)
      {
         if (mpi.Root()) { cout << "Rebalancing ..." << endl; }
         pmesh.Rebalance();

         // Update again after rebalancing
         Volta.Update();
      }
   }

   delete epsCoef;

   return 0;
}

// The Permittivity is a required coefficient which may be defined in
// various ways so we'll determine the appropriate coefficient type here.
Coefficient *
SetupPermittivityCoefficient()
{
   Coefficient * coef = NULL;

   if ( ds_params_.Size() > 0 )
   {
      coef = new FunctionCoefficient(dielectric_sphere);
   }
   else if ( pw_eps_.Size() > 0 )
   {
      coef = new PWConstCoefficient(pw_eps_);
   }
   else
   {
      coef = new ConstantCoefficient(epsilon0_);
   }

   return coef;
}

// A sphere with constant permittivity.  The sphere has a radius,
// center, and permittivity specified on the command line and stored
// in ds_params_.
double dielectric_sphere(const Vector &x)
{
   double r2 = 0.0;

   for (int i=0; i<x.Size(); i++)
   {
      r2 += (x(i)-ds_params_(i))*(x(i)-ds_params_(i));
   }

   if ( sqrt(r2) <= ds_params_(x.Size()) )
   {
      return ds_params_(x.Size()+1) * epsilon0_;
   }
   return epsilon0_;
}

// A sphere with constant charge density.  The sphere has a radius,
// center, and total charge specified on the command line and stored
// in cs_params_.
double charged_sphere(const Vector &x)
{
   double r2 = 0.0;
   double rho = 0.0;

   if ( cs_params_(x.Size()) > 0.0 )
   {
      switch ( x.Size() )
      {
         case 2:
            rho = cs_params_(x.Size()+1) /
                  (M_PI * pow(cs_params_(x.Size()), 2));
            break;
         case 3:
            rho = 0.75 * cs_params_(x.Size()+1) /
                  (M_PI * pow(cs_params_(x.Size()), 3));
            break;
         default:
            rho = 0.0;
      }
   }

   for (int i=0; i<x.Size(); i++)
   {
      r2 += (x(i) - cs_params_(i)) * (x(i) - cs_params_(i));
   }

   if ( sqrt(r2) <= cs_params_(x.Size()) )
   {
      return rho;
   }
   return 0.0;
}

// A Cylindrical Rod of constant polarization.  The cylinder has two
// axis end points, a radius, and a constant electric polarization oriented
// along the axis.
void voltaic_pile(const Vector &x, Vector &p)
{
   p.SetSize(x.Size());
   p = 0.0;

   Vector  a(x.Size());  // Normalized Axis vector
   Vector xu(x.Size());  // x vector relative to the axis end-point

   xu = x;

   for (int i=0; i<x.Size(); i++)
   {
      xu[i] -= vp_params_[i];
      a[i]   = vp_params_[x.Size()+i] - vp_params_[i];
   }

   double h = a.Norml2();

   if ( h == 0.0 )
   {
      return;
   }

   double  r = vp_params_[2 * x.Size()];
   double xa = xu * a;

   if ( h > 0.0 )
   {
      xu.Add(-xa / (h * h), a);
   }

   double xp = xu.Norml2();

   if ( xa >= 0.0 && xa <= h*h && xp <= r )
   {
      p.Add(vp_params_[2 * x.Size() + 1] / h, a);
   }
}

// To produce a uniform electric field the potential can be set
// to (- Ex x - Ey y - Ez z).
double phi_bc_uniform(const Vector &x)
{
   double phi = 0.0;

   for (int i=0; i<x.Size(); i++)
   {
      phi -= x(i) * e_uniform_(i);
   }

   return phi;
}
