#include "./alglib/dataanalysis.h"
#include "./alglib/ap.h"
#include "./alglib/stdafx.h"
#define ae_bool bool
#define ae_true true
#define ae_false false

#ifndef kmeans
#define kmeans

//this should be the function called from clustering.cpp
void custom_kmeans(const alglib::clusterizerstate &s, const alglib_impl::ae_int_t k, alglib::kmeansreport &rep, const alglib::xparams _xparams = alglib::xdefault);

void custom_kmeans_inner(alglib_impl::clusterizerstate* s, alglib_impl::ae_int_t k, alglib_impl::kmeansreport* rep, alglib_impl::ae_state *_state);

double average(std::vector<double> x);

void kmeans_internal(alglib_impl::ae_matrix* xy, alglib_impl::ae_int_t npoints, alglib_impl::ae_int_t nvars, alglib_impl::ae_int_t k, alglib_impl::ae_int_t initalgo, alglib_impl::ae_int_t seed, alglib_impl::ae_int_t maxits, alglib_impl::ae_int_t restarts, ae_bool kmeansdbgnoits, alglib_impl::ae_int_t* info, alglib_impl::ae_int_t* iterationscount, alglib_impl::ae_matrix* ccol, ae_bool needccol, alglib_impl::ae_matrix* crow, ae_bool needcrow, alglib_impl::ae_vector* xyc, double* energy, alglib_impl::kmeansbuffers* buf, alglib_impl::ae_state *_state);

#endif
