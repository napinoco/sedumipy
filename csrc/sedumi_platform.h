/* ************************************************************
   HEADER sedumi_platform.h

   Platform shim that lets the SeDuMi C kernels build either as before
   (MATLAB/Octave MEX files, via mex.h) or as part of a standalone,
   MATLAB/Octave-free build (define SEDUMI_STANDALONE).

   This is part of the ongoing MATLAB-free port (see python_port/README.md).
   The MEX build path is untouched -- this file only adds a second,
   parallel path that avoids mex.h entirely.
   ************************************************************ */

#if !defined(SEDUMI_PLATFORM_H)
#define SEDUMI_PLATFORM_H

#if defined(SEDUMI_STANDALONE)

#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>

/* Plain replacements for the MEX API's index/size types. mwIndex and
   mwSize are both unsigned and the same width on 32-bit-index MATLAB/
   Octave builds; size_t matches that here. */
typedef size_t mwIndex;
typedef size_t mwSize;
typedef ptrdiff_t mwSignedIndex;

/* Reference BLAS / OpenBLAS default to 32-bit Fortran integers (LP64,
   not ILP64) on Linux/macOS/Windows, unlike MATLAB's mwSignedIndex-sized
   blasint. Matches the FORT()-wrapped calls in sdmauxScalarmul.c,
   sdmauxRdot.c and blkchol2.c (dcopy/dscal/daxpy/ddot/idamax).

   SEDUMI_BLAS_ILP64 switches this to a 64-bit Fortran integer instead,
   for linking against an ILP64 BLAS build -- e.g. scipy-openblas64 (see
   tools/build_libsedumi.sh), which ships prebuilt, pip-installable
   OpenBLAS binaries for Linux/macOS/Windows and is preferred there over
   requiring a system BLAS install or building one. */
#if defined(SEDUMI_BLAS_ILP64)
#include <stdint.h>
typedef int64_t blasint;
#else
typedef int blasint;
#endif

/* scipy-openblas64 (and other co-installable OpenBLAS variants) rename
   every exported symbol with a prefix and/or suffix, so several such
   builds can coexist in one process without clashing -- e.g. "dscal_"
   becomes "scipy_dscal_64_" (prefix "scipy_", suffix "64_"). The build
   passes BLAS_SYMBOL_PREFIX/BLAS_SYMBOL_SUFFIX as bare-identifier macros
   (matching scipy_openblas64's own pkg-config Cflags) when linking such
   a build; default them to nothing otherwise so FORT(x) reduces to the
   plain x_/x it always was. The double indirection through
   SEDUMI_BLAS_CAT is required so the *values* of BLAS_SYMBOL_PREFIX/
   BLAS_SYMBOL_SUFFIX get substituted before token-pasting, not their
   literal names. */
#ifndef BLAS_SYMBOL_PREFIX
#define BLAS_SYMBOL_PREFIX
#endif
#ifndef BLAS_SYMBOL_SUFFIX
#define BLAS_SYMBOL_SUFFIX
#endif

#define SEDUMI_BLAS_CAT_(a, b) a##b
#define SEDUMI_BLAS_CAT(a, b) SEDUMI_BLAS_CAT_(a, b)

#if defined(SEDUMI_BLAS_NO_UNDERSCORE)
#define SEDUMI_BLAS_MANGLED(x) x
#else
#define SEDUMI_BLAS_MANGLED(x) x##_
#endif

#define FORT(x) \
    SEDUMI_BLAS_CAT(SEDUMI_BLAS_CAT(BLAS_SYMBOL_PREFIX, SEDUMI_BLAS_MANGLED(x)), BLAS_SYMBOL_SUFFIX)

extern void FORT(dcopy)(const blasint *n, const double *x, const blasint *incx,
                         double *y, const blasint *incy);
extern void FORT(dscal)(const blasint *n, const double *alpha, double *x,
                         const blasint *incx);
extern void FORT(daxpy)(const blasint *n, const double *alpha, const double *x,
                         const blasint *incx, double *y, const blasint *incy);
extern double FORT(ddot)(const blasint *n, const double *x, const blasint *incx,
                          const double *y, const blasint *incy);
extern blasint FORT(idamax)(const blasint *n, const double *x, const blasint *incx);

/* mexErrMsgTxt/mexPrintf/mxAssert equivalents: a standalone build can't pop
   a MATLAB error dialog, so report to stderr and abort instead. Kernels
   ported to be called from Python should prefer returning an error code
   over reaching this, but the legacy call sites (ported largely unchanged
   from the MEX sources) still expect a "this never returns" error path. */
void sedumi_fatal(const char *msg);
#define mexErrMsgTxt(msg) sedumi_fatal(msg)
#define mexWarnMsgTxt(msg) sedumi_warn(msg)
#define mexPrintf printf
#define mxAssert(cond, msg) do { if (!(cond)) sedumi_fatal(msg); } while (0)

/* A few of the "computational core" functions (e.g. symbfwmat() in
   symbfwblk.c) grow a work array with mxRealloc() directly, rather than
   only inside their file's mexFunction gateway. Map the MEX allocator
   family onto the standard one so those call sites need no change. */
#define mxCalloc(n, sz) calloc((n), (sz))
#define mxMalloc(sz) malloc(sz)
#define mxRealloc(p, sz) realloc((p), (sz))
#define mxFree(p) free(p)

void sedumi_warn(const char *msg);

#define SEDUMI_ASSERT(cond, msg) do { if (!(cond)) sedumi_fatal(msg); } while (0)

#else /* !SEDUMI_STANDALONE: original MEX build, unchanged */

#include "mex.h"
#ifdef OCTAVE
#include "f77blas.h"  /* defines "blasint" data type */
#define FORT(x) BLASFUNC(x)
#else /* Matlab */
#include "blas.h"
typedef ptrdiff_t blasint;
/**
 * For Matlab R2019a (probably before) and newer, when including
 * "blas.h" the respective BLAS identifiers are already defined,
 * e.g. for "dcopy":
 *
 *   #define dcopy FORTRAN_WRAPPER(dcopy)
 *
 * thus calling FORT(dcopy) == FORTRAN_WRAPPER(dcopy) inside SeDuMi
 * would result in "dcopy__" and already resulted in some bug reports.
 *
 * Compiling with -DFWRAPPER restores the previous behavior.
 */
#ifdef FWRAPPER
#define FORT(x) FORTRAN_WRAPPER(x)
#else
#define FORT(x) x
#endif
#endif

#define SEDUMI_ASSERT(cond, msg) mxAssert(cond, msg)

#endif /* SEDUMI_STANDALONE */

/* Guards the mwIndex -> blasint narrowing that every FORT()-wrapped BLAS
   call below does (`blasint one=1,nn=n;`). The BLAS interface integer can
   be narrower than this code's own index type in the standalone build: a
   reference-BLAS/OpenBLAS/Accelerate LP64 interface takes a 32-bit
   Fortran INTEGER, while mwIndex is size_t. Exceeding it needs a single
   contiguous vector longer than 2^31 doubles (~17 GB), which no problem
   this solver can set up will reach -- but the conversion is implicit
   and would silently truncate rather than fail, so trap it instead of
   computing a wrong answer. (An ILP64 build -- SEDUMI_BLAS_ILP64, e.g.
   scipy-openblas64 -- has no such 32-bit limit, same as MATLAB below;
   this can only actually trigger with a 32-bit blasint.)

   Written as a round-trip rather than a comparison against INT_MAX so it
   stays correct for every build: blasint is ptrdiff_t under MATLAB (an
   ILP64 BLAS, where this can never trigger) and whatever OpenBLAS's
   f77blas.h says under Octave. Casting back through mwIndex also catches
   the sign flip when bit 31 is set, since blasint is signed and mwIndex
   is not. */
#define SEDUMI_ASSERT_BLASINT_FITS(nn, n)                                \
    SEDUMI_ASSERT((mwIndex)(nn) == (mwIndex)(n),                         \
                  "vector length exceeds this build's BLAS integer "     \
                  "width -- see sedumi_platform.h")

#endif /* SEDUMI_PLATFORM_H */
