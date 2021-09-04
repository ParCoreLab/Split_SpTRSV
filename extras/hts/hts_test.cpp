// SILU driver
// Created: 03-12-2019
// Author: Najeeb Ahmad


//#include "common.h"
#include "matrix.h"
#include "util.h"
//#include "SILU.h"
//#include "quicksort.h"
//#include "policy.h"
#include "dummyFactorize.h"
#include "mmio.h"
//#include <cuda_runtime.h>
#include "hts_test_util.hpp"
#include "hts.hpp"
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace Experimental;

#define ITERATIONS 30

namespace Experimental {
namespace htsimpl {
template <typename T> bool is_complex () { return false; }
template <> bool is_complex<std::complex<double> >() { return true; }
template <> bool is_complex<std::complex<float> >() { return true; }

template<typename Int, typename Size, typename Sclr> class Tester {
  typedef HTS<Int, Size, Sclr> ihts;
  typedef typename ihts::Real Real;
  typedef util<Int, Size, Sclr> ut;
  typedef typename ut::TestOptions TestOptions;
  typedef typename ut::Data Data;

  static int are_same (const Data& A, const Data& B) {
#define same(expr, ret) if ( ! (expr)) return ret;
    same(A.m == B.m, 1);
    same(A.ir[A.m] == B.ir[B.m], 3);
    const Size nnz = A.ir[A.m];
    for (Int k = 0; k <= A.m; ++k)
      same(B.ir[k] == A.ir[k], 4);
    for (Size k = 0; k < nnz; ++k) {
      same(B.jc[k] == A.jc[k], 5);
      same(B.v[k] == A.v[k], 6);
    }
    return 0;
#undef same
  }

  static int test_transpose (const int verbose, const Int n) {
    Data ao, at, a;
  
    TestOptions to;
    to.n = n;
    to.matrix_type = TestOptions::sparse;
    ut::gen_tri_matrix(to, ao);
  
    ut::transpose(ao, at);
    ut::transpose(at, a);

    const int nerr = are_same(ao, a);
    if (nerr && verbose)
      std::cout << "failed: test_transpose\n";
    return nerr;
  }

static int test_lt (const TestOptions& to, csr_matrix *sp_mat, const bool print_options=false,
                   const bool exception_expected=false) {
    int nerr = 0;

    const Int max_nrhs = 1;
    const Real tol = std::numeric_limits<Real>::epsilon()*1e6;
    float times[100];
    int my_time = 0;
    high_resolution_clock::time_point t1;
    Data d;
    int nnz;
    d.clear();
    {
      ut::populate_sparse_matrix(d, sp_mat);
      nnz = d.ir.back();
      ut::gen_rand_perm(d.m, d.p);   // row permutation vector
      ut::gen_rand_perm(d.m, d.q);   // col permutation vector
      ut::gen_rand_vector(d.m, d.r); // row scaling

    }

    const Int ldb = d.m + 3, ldx = d.m + 4;
    std::vector<Sclr> b(ldb*max_nrhs), xt(d.m*max_nrhs), x(ldx*max_nrhs);
    {
      ut::gen_rand_vector(xt.size(), xt);
      // Generate the rhs b.
      std::vector<Sclr> y(d.m);

      for (Int irhs = 0; irhs < max_nrhs; ++irhs) {
        const Sclr* const xtp = xt.data() + irhs*d.m;
        Sclr* const bp = b.data() + irhs*ldb;
        for (Int i = 0; i < d.m; ++i)
          x[i] = xtp[d.q[i]];
        ut::mvp(d, to.transpose, to.conjugate, x.data(), y.data());
        if (to.has_unit_diag())
          for (Int i = 0; i < d.m; ++i)
            y[i] += x[i];
        for (Int i = 0; i < d.m; ++i)
          bp[d.p[i]] = y[i];
        for (Int i = 0; i < d.m; ++i)
          bp[i] *= d.r[i];
      }
    }


    std::vector<Sclr> bo(b);

    typename ihts::CrsMatrix* T;
    typename ihts::Options opts;

    {
      T = ihts::make_CrsMatrix(d.m, d.ir.data(), d.jc.data(), d.v.data(),
                               to.transpose, to.conjugate);
      if ( ! to.reprocess && to.nthreads > 1)
        ihts::register_Deallocator(T, &d);
      if (to.solve_type == TestOptions::ls_only)
        ihts::set_level_schedule_only(opts);
      else if (to.solve_type == TestOptions::rb_only)
        opts.min_lset_size = d.m + 1;
      // To really test things well, choose very small block sizes. (This is bad
      // for performance.) These parameters are not meant to be set by the user
      // ordinarily, but they are exposed in case an expert is tuning performance
      // or, as here, for testing.
      
      opts.min_block_size = 6;
      opts.min_parallel_rows = 2;
      opts.pp_min_block_size = 12;
      if (to.matrix_type == TestOptions::block_sparse)
        opts.levelset_block_size = to.block_size;
    }

    {
      typename ihts::Impl* impl;
      try {
        std::vector<Sclr> Td;
        if (to.reprocess) {
          // Save the true values.
          Td = d.v;
          // Zero the matrix values to simulate reprocessing.
          d.v.assign(d.v.size(), 1);
        }
        impl = ihts::preprocess(T, max_nrhs - 1 /* For testing; see below. */,
                                to.nthreads, to.reprocess, d.p.data(), d.q.data(),
                                d.r.data(), &opts);
        if (to.reprocess) {
          // Restore the values.
          d.v = Td;
        }
      } catch (...) {
        if ( ! exception_expected) {
          std::cerr << "Unexpected exception on ";
          to.print(std::cerr);
          std::cerr << "\n";
          ut::write_matrixmarket(d, "unexpected_exception.mm");
        }
        ihts::delete_CrsMatrix(T);
        throw;
      }
      if (print_options)
        ihts::print_options(impl, std::cout);
      if (to.reprocess) {
        // Run 2 times to test idempotency.
        for (int rep = 0; rep < 2; ++rep)
          ihts::reprocess_numeric(impl, T, d.r.data());
      }
      // Exercise reset_max_nrhs.
      ihts::reset_max_nrhs(impl, max_nrhs);
      if (ihts::is_lower_tri(impl) &&
          ((to.upper && ! to.transpose) || ( ! to.upper && to.transpose)) &&
          d.m > 1 && nnz > static_cast<Size>(d.m) /* not diag */)
        ++nerr;
      for (int slv = 0; slv <= 0; ++slv) {
        // Check each solve interface.
        switch (slv) {
        case 0:
          t1 = start_profile_CPU();
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), ldb, ldx);
          times[my_time] = end_profile_CPU(t1);
          printf("Time: %f msec\n", times[my_time]);
          my_time++;
          break;
        case 1:
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), 0.0,  1.0,
                          ldb, ldx);
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), 2.0, -1.0,
                          ldb, ldx);
          break;
        case 2:
          ihts::solve_omp(impl, b.data(), to.nrhs, ldb);
          break;
        }
        double rd = 0;
        for (Int i = 0; i < to.nrhs; ++i)
          rd = std::max(
            rd, ut::reldif(xt.data() + i*d.m,
                           (slv == 2 ? b.data() + i*ldb : x.data() + i*ldx),
                           d.m));
        if (slv == 2) b = bo;
        if (rd >= tol) {
          ++nerr;
          if (to.verbose) std::cout << "rd " << slv << ": " << rd << "\n";
        }
      }
      ihts::delete_Impl(impl);
    }

    if (to.nthreads == 1) {
      std::vector<Sclr> xb(d.m*to.nrhs), w(d.m);
      for (Int irhs = 0; irhs < to.nrhs; ++irhs) {
        const Sclr* const bp = b.data() + irhs*ldb;
        Sclr* const xbp = xb.data() + irhs*d.m;
        for (Int i = 0; i < d.m; ++i)
          xbp[i] = bp[i];
      }
      ihts::solve_serial(T, ! to.upper, to.has_unit_diag(), xb.data(), to.nrhs,
                         d.p.data(), d.q.data(), d.r.data(), w.data());
      const double rd = ut::reldif(xt.data(), xb.data(), d.m*to.nrhs);
      if (rd >= tol) {
        ++nerr;
        if (to.verbose) std::cout << "serial rd: " << rd << "\n";
      }
    }

    ihts::delete_CrsMatrix(T);

    if (to.verbose) {
      const bool print = to.verbose == 2 || (to.verbose == 1 && nerr);
      if (print) {
        std::cout << (nerr ? "fail" : "pass") << "ed: ";
        to.print(std::cout);
        std::cout << "\n";
      }
    }        

    return nerr;
  }

  // Run one valid triangle. Solve for x in
  //     P' T Q' x = R \ b,
  // where P = diag(p) and similarly for Q and R, and T is a triangle.
  static int test (const TestOptions& to, csr_matrix *sp_mat, const bool print_options=false,
                   const bool exception_expected=false) {
    int nerr = 0;

    const Int max_nrhs = 1;
    const Real tol = std::numeric_limits<Real>::epsilon()*1e6;

    // Generate matrix data.
    Data d;
    Size nnz;
    {
      ut::gen_tri_matrix(to, d);     // triangle
      nnz = d.ir.back();
      ut::gen_rand_perm(d.m, d.p);   // row permutation vector
      ut::gen_rand_perm(d.m, d.q);   // col permutation vector
      ut::gen_rand_vector(d.m, d.r); // row scaling
    }

    const Int ldb = d.m + 3, ldx = d.m + 4;
    std::vector<Sclr> b(ldb*max_nrhs), xt(d.m*max_nrhs), x(ldx*max_nrhs);
    {
      // True x.
      ut::gen_rand_vector(xt.size(), xt);
      // Generate the rhs b.
      std::vector<Sclr> y(d.m);
      for (Int irhs = 0; irhs < max_nrhs; ++irhs) {
        const Sclr* const xtp = xt.data() + irhs*d.m;
        Sclr* const bp = b.data() + irhs*ldb;
        for (Int i = 0; i < d.m; ++i)
          x[i] = xtp[d.q[i]];
        ut::mvp(d, to.transpose, to.conjugate, x.data(), y.data());
        if (to.has_unit_diag())
          for (Int i = 0; i < d.m; ++i)
            y[i] += x[i];
        for (Int i = 0; i < d.m; ++i)
          bp[d.p[i]] = y[i];
        for (Int i = 0; i < d.m; ++i)
          bp[i] *= d.r[i];
      }
    }
    std::vector<Sclr> bo(b);

    typename ihts::CrsMatrix* T;
    typename ihts::Options opts;
    {
      T = ihts::make_CrsMatrix(d.m, d.ir.data(), d.jc.data(), d.v.data(),
                               to.transpose, to.conjugate);
      if ( ! to.reprocess && to.nthreads > 1)
        ihts::register_Deallocator(T, &d);
      if (to.solve_type == TestOptions::ls_only)
        ihts::set_level_schedule_only(opts);
      else if (to.solve_type == TestOptions::rb_only)
        opts.min_lset_size = d.m + 1;
      // To really test things well, choose very small block sizes. (This is bad
      // for performance.) These parameters are not meant to be set by the user
      // ordinarily, but they are exposed in case an expert is tuning performance
      // or, as here, for testing.
      opts.min_block_size = 6;
      opts.min_parallel_rows = 2;
      opts.pp_min_block_size = 12;
      if (to.matrix_type == TestOptions::block_sparse)
        opts.levelset_block_size = to.block_size;
    }

    {
      typename ihts::Impl* impl;
      try {
        std::vector<Sclr> Td;
        if (to.reprocess) {
          // Save the true values.
          Td = d.v;
          // Zero the matrix values to simulate reprocessing.
          d.v.assign(d.v.size(), 1);
        }
        impl = ihts::preprocess(T, max_nrhs - 1 /* For testing; see below. */,
                                to.nthreads, to.reprocess, d.p.data(), d.q.data(),
                                d.r.data(), &opts);
        if (to.reprocess) {
          // Restore the values.
          d.v = Td;
        }
      } catch (...) {
        if ( ! exception_expected) {
          std::cerr << "Unexpected exception on ";
          to.print(std::cerr);
          std::cerr << "\n";
          ut::write_matrixmarket(d, "unexpected_exception.mm");
        }
        ihts::delete_CrsMatrix(T);
        throw;
      }
      if (print_options)
        ihts::print_options(impl, std::cout);
      if (to.reprocess) {
        // Run 2 times to test idempotency.
        for (int rep = 0; rep < 2; ++rep)
          ihts::reprocess_numeric(impl, T, d.r.data());
      }
      // Exercise reset_max_nrhs.
      ihts::reset_max_nrhs(impl, max_nrhs);
      if (ihts::is_lower_tri(impl) &&
          ((to.upper && ! to.transpose) || ( ! to.upper && to.transpose)) &&
          d.m > 1 && nnz > static_cast<Size>(d.m) /* not diag */)
        ++nerr;
      for (int slv = 0; slv <= 2; ++slv) {
        // Check each solve interface.
        switch (slv) {
        case 0:
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), ldb, ldx);
          break;
        case 1:
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), 0.0,  1.0,
                          ldb, ldx);
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), 2.0, -1.0,
                          ldb, ldx);
          break;
        case 2:
          ihts::solve_omp(impl, b.data(), to.nrhs, ldb);
          break;
        }
        double rd = 0;
        for (Int i = 0; i < to.nrhs; ++i)
          rd = std::max(
            rd, ut::reldif(xt.data() + i*d.m,
                           (slv == 2 ? b.data() + i*ldb : x.data() + i*ldx),
                           d.m));
        if (slv == 2) b = bo;
        if (rd >= tol) {
          ++nerr;
          if (to.verbose) std::cout << "rd " << slv << ": " << rd << "\n";
        }
      }
      ihts::delete_Impl(impl);
    }

    if (to.nthreads == 1) {
      std::vector<Sclr> xb(d.m*to.nrhs), w(d.m);
      for (Int irhs = 0; irhs < to.nrhs; ++irhs) {
        const Sclr* const bp = b.data() + irhs*ldb;
        Sclr* const xbp = xb.data() + irhs*d.m;
        for (Int i = 0; i < d.m; ++i)
          xbp[i] = bp[i];
      }
      ihts::solve_serial(T, ! to.upper, to.has_unit_diag(), xb.data(), to.nrhs,
                         d.p.data(), d.q.data(), d.r.data(), w.data());
      const double rd = ut::reldif(xt.data(), xb.data(), d.m*to.nrhs);
      if (rd >= tol) {
        ++nerr;
        if (to.verbose) std::cout << "serial rd: " << rd << "\n";
      }
    }

    ihts::delete_CrsMatrix(T);

    if (to.verbose) {
      const bool print = to.verbose == 2 || (to.verbose == 1 && nerr);
      if (print) {
        std::cout << (nerr ? "fail" : "pass") << "ed: ";
        to.print(std::cout);
        std::cout << "\n";
      }
    }

    return nerr;
  }

  // Run one invalid triangle and make sure an exception is thrown.
  static int test_for_exception (const TestOptions& to) {
    bool correct_exception_thrown = false;
    try {
      test(to, false, true);
    } catch (const hts::NotFullDiagonalException&) {
      correct_exception_thrown = to.matrix_type == TestOptions::missing_some_diag;
    } catch (const hts::NotTriangularException&) {
      correct_exception_thrown = ut::is_not_tri(to.matrix_type);
    }
    const int nerr = correct_exception_thrown ? 0 : 1;
    if (nerr && to.verbose) {
      std::cout << "test_for_exception failed for ";
      to.print(std::cout);
      std::cout << "\n";
    }
    return nerr;
  }

public:

  static high_resolution_clock::time_point start_profile_CPU()
  {
      high_resolution_clock::time_point hrt1 = high_resolution_clock::now();
      return hrt1;
  }

  static float end_profile_CPU(high_resolution_clock::time_point hrt1)
  {
        high_resolution_clock::time_point hrt2 = high_resolution_clock::now();
        duration<float> t_diff = hrt2 - hrt1;
        nanoseconds t_diff_ns = duration_cast<nanoseconds>(t_diff);
        float elapsed_time = t_diff_ns.count()/1e6;
    
        return elapsed_time;
  }

  // Run all the tests.
  static int test (const int verbose, csr_matrix *sp_mat) {
    int nerr = 0;

    // Test that we throw on an unsigned Int.
    {
      bool caught = false;
      try {
        HTS<unsigned int, int, double>::make_CrsMatrix(1, 0, 0, 0);
      } catch (const hts::Exception& e) {
        caught = true;
      }
      if ( ! caught) ++nerr;
    }
        // Test our own transpose to make sure it's OK for subsequent use.
    nerr += test_transpose(verbose, 277);

    const int ns[] = {1, 2, 3, 21, 300};
    const int max_nthreads = omp_get_max_threads();
    const int nthreads_step = max_nthreads > 40 ? 11 : 3;

    TestOptions to;
    //to.block_size = 3; // only for block_sparse
    to.verbose = verbose;
    bool print_options = to.verbose;
    for (int nthreads = 1; ; ) {
      to.nthreads = nthreads;
      if (to.verbose >= 1) {
        std::cout << "Threads: " << nthreads << "    | ";
        std::cout.flush();
      }

      for (int ti = 0; ti < 1; ++ti) {
        to.transpose = ti % 2;
        if ( ! is_complex<Sclr>() && ti >= 2) break;
        if (ti >= 2) to.conjugate = true;
        for (size_t ni = 0; ni < 1; ++ni) {
          to.n = ns[ni];
          for (size_t si = 0; si < 1; ++si) {
            to.solve_type = static_cast<typename ut::TestOptions::SolveType>(si);
            for (int ui = 0; ui < 1; ++ui) {
              to.upper = ui;
              for (int ri = 0; ri < 1; ++ri) {
                for (int nrhs = 1; nrhs < 2; nrhs += 2) {
                  to.nrhs = nrhs;
                  to.reprocess = ri;
                  //to.matrix_type = TestOptions::diag;
                  //nerr += test(to, print_options);
                  print_options = false;
                  //to.matrix_type = TestOptions::dense; nerr += test(to);
                  to.matrix_type = TestOptions::sparse; nerr += test_lt(to, sp_mat);
                  //to.matrix_type = TestOptions::block_sparse; nerr += test(to);
                  //to.matrix_type = TestOptions::implicit_unit_diag; nerr += test(to);
                  //to.matrix_type = TestOptions::block_sparse_implicit_unit_diag; nerr += test(to);
                  if (to.n > 2) {
                    //printf("Was here\n");
                    to.matrix_type = TestOptions::not_tri;
                    nerr += test_for_exception(to);
                    to.matrix_type = TestOptions::missing_some_diag;
                    nerr += test_for_exception(to);
                    to.matrix_type = TestOptions::not_tri_almost_diag;
                    nerr += test_for_exception(to);
                  }
                }
              }
            }
          }
        }
      }
      if (nthreads == max_nthreads) break;
      nthreads = std::min(max_nthreads, nthreads + nthreads_step);
    }

    return nerr;
  }
};
} // namespace htsimpl
} // namespace Experimental

static high_resolution_clock::time_point start_profile_CPU()
  {
      high_resolution_clock::time_point hrt1 = high_resolution_clock::now();
      return hrt1;
  }

  static float end_profile_CPU(high_resolution_clock::time_point hrt1)
  {
        high_resolution_clock::time_point hrt2 = high_resolution_clock::now();
        duration<float> t_diff = hrt2 - hrt1;
        nanoseconds t_diff_ns = duration_cast<nanoseconds>(t_diff);
        float elapsed_time = t_diff_ns.count()/1e6;
    
        return elapsed_time;
  }

  static double urand () {
#if 0 // Not all compilers have this, it seems.
  static std::default_random_engine generator;
  static std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
#else
  return rand() / ((double) RAND_MAX + 1.0);
#endif
}

template<typename T> static inline double square (const T& x) { return x*x; }

 static double
  reldif (const double* const a, const double* const b, const int n) {
    double num = 0, den = 0;
    for (int i = 0; i < n; ++i) {
      num += square(a[i] - b[i]);
      den += square(a[i]);
    }
    return std::sqrt(num/den);
  }



int main(int argc, char **argv)
{
    if(argc < 3)
    {
    	std::cout << "Error: Missing arguments\n";
    	std::cout << "Usage: " << argv[0] << " matrix_id output_folder\n";
    	return EXIT_FAILURE;
    }

    csr_matrix matrix;
    string matrix_name;
    ofstream outfile("results.txt", ios::app);
    if(outfile.is_open())
    {
    }
    else
    {
        cout << "Failed to open file!" << endl;
        return EXIT_FAILURE;
    }
    if (argc == 4)
    {
        printf("Reading .mtx file\n");
        int retCode = 0;
	   matrix_name = argv[4];
        retCode = mm_read_unsymmetric_sparse(argv[3], &matrix.m, &matrix.n, &matrix.nnz,
                &matrix.csrVal, &matrix.csrRowPtr, &matrix.csrColIdx);
        coo2csr_in(matrix.m, matrix.nnz, matrix.csrVal, matrix.csrRowPtr, matrix.csrColIdx);
        printf("Done reading file\n");
        printf("Matrix Rows: %d\n", matrix.m);
        printf("Matrix Cols: %d\n", matrix.n);

        if(retCode == -1)
        {
            cout << "Error reading input .mtx file\n";
            return EXIT_FAILURE;
        }
    }
    else
    {
        printf("Downloading matrix from Suitsparse\n");
        string matrix_id(argv[1]), output_folder(argv[2]);
        int mat_id = stoi(matrix_id);

        uf_matrix_reader reader;
        matrix_name = reader.get_matrix_name(mat_id);
        if (does_file_exist(path_join(output_folder, matrix_name))) 
        {
            std::cout << "Already calculated for " << matrix_name << "\n";
            return EXIT_SUCCESS;
        }
        matrix = reader.get_matrix_csr_libufget(mat_id);

    }
    csr_matrix L, U;
    high_resolution_clock::time_point t1;
    dummyFactorize(&matrix, &L, &U, SUBSTITUTION_FORWARD);
    typedef HTS<int, int, double> ihts;
    typedef typename ihts::Real Real;
    typename ihts::CrsMatrix* T;
    std::vector<int> ir;
    std::vector<int> jc, p, q;
    std::vector<double> v;
    std::vector<double> r;
    const double tol = std::numeric_limits<Real>::epsilon()*1e6;
    float times[ITERATIONS];

    // Copy lower-triangular matrix to CrsMatrix
    for(int i = 0; i < L.nnz; i++)
    {
      jc.push_back(L.csrColIdx[i]);
      v.push_back(L.csrVal[i]);
    }

    for(int i = 0; i <= L.n; i++)
    {
      ir.push_back(L.csrRowPtr[i]);      
    }

    // Create p and q matrices
    p.resize(L.m);
    q.resize(L.m);
    r.resize(L.m);
    for (size_t i = 0; i < L.m; ++i)
    {
      p[i] = i;
      q[i] = i;
      r[i] = 1.0;
    }

    const int ldb = L.m + 3, ldx = L.m + 4;
    std::vector<double> b(ldb), xt(L.m), x(ldx);  

    // Create CSR matrix for HTS library
    T = ihts::make_CrsMatrix(L.m, ir.data(), jc.data(), v.data());
    const int max_nthreads = omp_get_max_threads();
    // printf("OMP_NUM_THREADS=%d\n", max_nthreads);
    int nrhs = 1;

    // Generate right-hand-side
    xt.resize(L.m);
    for (size_t i = 0; i < L.m; ++i)
      xt[i] = urand() - 0.5;

    std::vector<double> y(L.m);

    const double* const xtp = xt.data();
    double* const bp = b.data();
    for (int i = 0; i < L.m; ++i)
        x[i] = xtp[q[i]];
    
    // mv product
    for (int i = 0; i < L.m; ++i) {
        double acc = 0;
        const int iri = ir[i], irip1 = ir[i+1];
        for (int j = iri; j < irip1; ++j) {
          acc += v[j]*x[jc[j]];
        }
        y[i] = acc;
      }

    for (int i = 0; i < L.m; ++i)
       bp[p[i]] = y[i];
    for (int i = 0; i < L.m; ++i)
       bp[i] *= r[i];

    std::vector<double> bo(b);
    typename ihts::Options opts;
    //ihts::set_level_schedule_only(opts);
    float preprocess_time = 0.0;
    t1 = start_profile_CPU();
    ihts::Impl* Limpl = ihts::preprocess(T, nrhs, max_nthreads, 0, p.data(), q.data(), r.data());
    preprocess_time = end_profile_CPU(t1);

    float my_time = 0.0;
    for(int i = 0; i < ITERATIONS; i++)
    {
      t1 = start_profile_CPU();
      ihts::solve_omp(Limpl, b.data(), 1, x.data(), ldb, ldx);
      times[i] = end_profile_CPU(t1);
      //printf("%f\n",times[i]);      
      my_time += times[i];  
    }

    std::vector<double> xb(L.m), w(L.m);
    const double* const bp1 = b.data();
    double* const xbp = xb.data();
    for (int i = 0; i < L.m; ++i)
      xbp[i] = bp1[i];

    ihts::solve_serial(T, 1, 0, xb.data(), 1,
                         p.data(), q.data(), r.data(), w.data());

    // for(int i = 0; i < L.m; i++)
    // {
    //   printf("%f, %f, %f\n", xt[i], xb[i],x[i]);
    // }

    const double rd = reldif(xt.data(), x.data(), L.m);
    if (rd >= tol) {
        printf("HTS: Test failed\n");
      }
      else
      {
        printf("HTS: Test passed\n");
        printf("%s %0.6f\n", "SpTRSV", my_time/ITERATIONS);
        printf("%s %0.6f\n", "Preprocess", preprocess_time);
      }
    
    //int nerr = Experimental::htsimpl::Tester<int, int, double>::test(1, &L);
    //std::cout << "HTS Test: " << (nerr ? "Failed" : "Passed") << "\n";
    
	return EXIT_SUCCESS;
}
