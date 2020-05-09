// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <array>
#include <bitset>
#include <map>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <cmath>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/localfunctions/lagrange/pk.hh>

template <class GV>
int izracunajProfilMatrice(GV const &gv, std::vector<std::set<int>> &profilMatrice) {
      int m = 0;
      for (auto const &element : elements(gv)) {
          for (auto const & is : intersections(gv,element)){
              if (is.boundary()){
                  ++m;
              }
          }
      }
      const int dim = GV::dimensionworld;
      // Broj vrhova mreže
      const auto N = gv.size(dim);
      int M = N - m;
      profilMatrice.resize(M);

      // Mapper koji indeksira sve vrhove mreže
      auto const &idset = gv.indexSet();

      using vektor = Dune::FieldVector<double, dim>;
      std::unordered_set<vektor> a;
      
      // Po svim elementima mreže
      for (auto const &element : elements(gv)) {
          for (auto const & is : intersections(gv,element)){
              if (is.boundary()){
                  a.insert(is.geometry().corner(0));
              }
          }
      }
      for (auto const &element : elements(gv)) {
          int m2 = 0;
          for (auto const & is : intersections(gv,element)){
              if (is.boundary()){
                  ++m2;
              }
          }
        if (m2 > 0) continue;

        int m3 = 0;
        for (auto const &vertex : vertices(gv)){
            for (auto k : a){
                if (k == vertex.geometry().corner(0)) {
                    ++m3;
                }
                if (m3 > 0) break;
            }
            if (m3 > 0) break;
        }
        if (m3 > 0) continue;

        // Broj vrhova u elementu
        int n_v = element.geometry().corners();

        // A_{i,j} je različito od nule samo ako vrhovi s indeksima i,j pripadaju
        // istom elementu. Stioga je dovoljno na svakom element skupiti sve parove
        // indeksa vrhova.

        // VAŠ KOD DOLAZI OVDJE!
        for (int i = 0; i < n_v; ++i){
            auto index_i = idset.subIndex(element, i, dim);
            for (int j = 0; j < n_v; ++j){
                auto index_j = idset.subIndex(element, j, dim);
                profilMatrice[index_i].insert(index_j);
            }
        }
      }
      std::cout << " Profil matrice je izračunat.\n";
      return M;
}

template <typename GV, typename Mat, typename Vec>
void init(GV const &gv, Mat &A, Vec &b)
{
  const int dim = GV::dimensionworld;
  std::vector<std::set<int>> profilMatrice;

  // Broj vrhova mreže
  int N = izracunajProfilMatrice(gv, profilMatrice);

  // dimenzije matrice A i vektora desne strane b
  A.setSize(N, N);
  A.setBuildMode(Mat::random);
  b.resize(N, false);

  // Postavimo broj ne-nul elemenata u svakom retku
  for (int i = 0; i < N; ++i){
      A.setrowsize(i, profilMatrice[i].size());
  }

  A.endrowsizes();  // dimenzioniranje redaka kompletirano

  // Definiraj indekse stupaca u svakom retki
  for (int i = 0; i < N; ++i){
      for (auto j : profilMatrice[i])
        A.addindex(i, j);
  }
  A.endindices();  // definiranje stupaca završeno
  // Profil matrice je time fiksiran
  // inicijalizacija nulama
  A = 0.0;
  b = 0.0;
  std::cout << " Matrica i vektor desne strane su inicijalizirani.\n";
}

template <typename GV, typename FEM, typename Mat, typename Vec>
void assemble(GV &gv, FEM const &fem, Mat &A, Vec &b) {
  // A = matrica krutosti, b = vektor desne strane
  const int dim = GV::dimensionworld;

  using FEMGradient = typename FEM::Traits::LocalBasisType::Traits::JacobianType;
  using FEMRange    = typename FEM::Traits::LocalBasisType::Traits::RangeType;

  init(gv, A, b);

  auto const &mapper = gv.indexSet();

  using vektor = Dune::FieldVector<double, dim>;
  std::unordered_set<vektor> a1;

  // Po svim elementima mreže
  for (auto const &element : elements(gv)) {
      for (auto const & is : intersections(gv,element)){
          if (is.boundary()){
              a1.insert(is.geometry().corner(0));
          }
      }
  }

  for (auto const &element : elements(gv)) {

    int m = 0;
    for (auto const & is : intersections(gv,element)){
          if (is.boundary()) ++m;
    }
    if (m > 0) continue;

    int m3 = 0;
    for (auto const &vertex : vertices(gv)){
        for (auto k : a1){
            if (k == vertex.geometry().corner(0)){
                ++m3;
            }
            if (m3 > 0) break;
        }
        if (m3 > 0) break;
    }
    if (m3 > 0) continue;

    auto gt = element.type();
    auto vertexsize = element.geometry().corners();
    auto basis_size = fem.localBasis().size();
    // provjera da baznih funkcija ima koliko i vrhova
    assert(vertexsize == basis_size);

    const auto &rule = Dune::QuadratureRules<double, dim>::rule(gt, 2);

    for (auto const &qpoint : rule)
    {
        // VAŠ KOD DOLAZI OVDJE!!!
        auto const weight = qpoint.weight();
        auto const & point = qpoint.position();
        auto dx = element.geometry().integrationElement(point);
        auto jac = element.geometry().jacobianInverseTransposed(point);
        auto acoeff = 0.0;
        auto fcoeff = 0.0;
        
        auto kcoeff = 1.0;

        std::vector<FEMRange> phihat(basis_size);
        fem.localBasis().evaluateFunction(point, phihat);

        std::vector<FEMGradient> gradphihat(basis_size);
        fem.localBasis().evaluateJacobian(point, gradphihat);

        //  Po svim baznim funkcijama
        for (unsigned int i = 0; i < basis_size; i++) {

          auto index_i = mapper.subIndex(element, i, dim);

          b[index_i] += fcoeff * phihat[i] * dx * weight;  // volumni dio desne strane

          Dune::FieldVector<double, dim> grad1;
          jac.mv(gradphihat[i][0], grad1);  // grad1 = jac *gradphihat[i][0]

          // Po svim baznim funkcijama
          for (int j = 0; j < basis_size; j++) {
          
            auto index_j = mapper.subIndex(element, j, dim);
            Dune::FieldVector<double, dim> grad2;
            jac.mv(gradphihat[j][0], grad2);  // grad2 = jac *gradphihat[j][0]
            
            A[index_i][index_j] += ( kcoeff * (grad1 * grad2) + (acoeff) * phihat[i] * phihat[j] ) * dx * weight;
          }
        }
    }  
  }
}  

int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is krilo." << std::endl;
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
        <<" processes!"<<std::endl;

    constexpr int dim = 2;
    typedef Dune::UGGrid<dim> GridType;
    using GridView = GridType::LeafGridView;
    GridType *pgrid = Dune::GmshReader<GridType>::read("src_dir/grids/airfoil_exterior.msh", true, false);

    typedef GridType::LeafGridView GV;
    const GV &gv = pgrid->leafGridView();


    using Matrix = Dune::BCRSMatrix< Dune::FieldMatrix<double, 1, 1> >;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>    >;

    Matrix A;     // Matrica krutosti
    Vector x, b;  // rješenje i desna strana

    using FEM = Dune::PkLocalFiniteElement<double, double, dim, 1>;
    FEM fem;

    assemble(gv, fem, A, b);

    // formiranje rješavača
    Dune::MatrixAdapter<Matrix, Vector, Vector> op(A);
    Dune::SeqILU<Matrix, Vector, Vector> ilu1(A, 1, 0.92);
    Dune::BiCGSTABSolver<Vector> bcgs(op, ilu1, 1e-15, 5000, 0);
    Dune::InverseOperatorResult r;

    x.resize(b.N());
    x = 0.0;

    // rješavanje sustava
    bcgs.apply(x, b, r);
    std::cout << " Sustav je riješen.\n";

    if(r.converged){
        std::cout << "Broj iteracija = "        << r.iterations
                  << ", redukcija residuala = " << r.reduction << "\n";
    }
    else
        std::cout << "Solver nije konvergirao.\n";

    Dune::VTKWriter<GridView> vtkwriter(gv);
    vtkwriter.addVertexData(x, "u");
    vtkwriter.write("fem", Dune::VTK::OutputType::ascii);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
