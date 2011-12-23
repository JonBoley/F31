/* rng/gsl_rng.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef __GSL_RNG_H__
#define __GSL_RNG_H__
#include <stdlib.h>
/* #include <gsl/gsl_errno.h> */
#include "gsl_errno.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct
  {
    const char *name;
    unsigned long int max;
    unsigned long int min;
    size_t size;
    void (*set) (void *state, unsigned long int seed);
    unsigned long int (*get) (void *state);
    double (*get_double) (void *state);
  }
gsl_rng_type;

typedef struct
  {
    const gsl_rng_type * type;
    void *state;
  }
gsl_rng;


/* These structs also need to appear in default.c so you can select
   them via the environment variable GSL_RNG_TYPE */

#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_borosh13;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_borosh13;
#else
extern const gsl_rng_type *gsl_rng_borosh13;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_coveyou;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_coveyou;
#else
extern const gsl_rng_type *gsl_rng_coveyou;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_cmrg;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_cmrg;
#else
extern const gsl_rng_type *gsl_rng_cmrg;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_fishman18;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_fishman18;
#else
extern const gsl_rng_type *gsl_rng_fishman18;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_fishman20;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_fishman20;
#else
extern const gsl_rng_type *gsl_rng_fishman20;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_fishman2x;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_fishman2x;
#else
extern const gsl_rng_type *gsl_rng_fishman2x;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_gfsr4;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_gfsr4;
#else
extern const gsl_rng_type *gsl_rng_gfsr4;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_knuthran;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_knuthran;
#else
extern const gsl_rng_type *gsl_rng_knuthran;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_knuthran2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_knuthran2;
#else
extern const gsl_rng_type *gsl_rng_knuthran2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_lecuyer21;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_lecuyer21;
#else
extern const gsl_rng_type *gsl_rng_lecuyer21;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_minstd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_minstd;
#else
extern const gsl_rng_type *gsl_rng_minstd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_mrg;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_mrg;
#else
extern const gsl_rng_type *gsl_rng_mrg;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_mt19937;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_mt19937;
#else
extern const gsl_rng_type *gsl_rng_mt19937;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_r250;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_r250;
#else
extern const gsl_rng_type *gsl_rng_r250;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ran0;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ran0;
#else
extern const gsl_rng_type *gsl_rng_ran0;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ran1;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ran1;
#else
extern const gsl_rng_type *gsl_rng_ran1;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ran2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ran2;
#else
extern const gsl_rng_type *gsl_rng_ran2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ran3;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ran3;
#else
extern const gsl_rng_type *gsl_rng_ran3;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_rand;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_rand;
#else
extern const gsl_rng_type *gsl_rng_rand;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_rand48;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_rand48;
#else
extern const gsl_rng_type *gsl_rng_rand48;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random128_bsd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random128_bsd;
#else
extern const gsl_rng_type *gsl_rng_random128_bsd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random128_glibc2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random128_glibc2;
#else
extern const gsl_rng_type *gsl_rng_random128_glibc2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random128_libc5;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random128_libc5;
#else
extern const gsl_rng_type *gsl_rng_random128_libc5;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random256_bsd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random256_bsd;
#else
extern const gsl_rng_type *gsl_rng_random256_bsd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random256_glibc2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random256_glibc2;
#else
extern const gsl_rng_type *gsl_rng_random256_glibc2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random256_libc5;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random256_libc5;
#else
extern const gsl_rng_type *gsl_rng_random256_libc5;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random32_bsd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random32_bsd;
#else
extern const gsl_rng_type *gsl_rng_random32_bsd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random32_glibc2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random32_glibc2;
#else
extern const gsl_rng_type *gsl_rng_random32_glibc2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random32_libc5;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random32_libc5;
#else
extern const gsl_rng_type *gsl_rng_random32_libc5;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random64_bsd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random64_bsd;
#else
extern const gsl_rng_type *gsl_rng_random64_bsd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random64_glibc2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random64_glibc2;
#else
extern const gsl_rng_type *gsl_rng_random64_glibc2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random64_libc5;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random64_libc5;
#else
extern const gsl_rng_type *gsl_rng_random64_libc5;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random8_bsd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random8_bsd;
#else
extern const gsl_rng_type *gsl_rng_random8_bsd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random8_glibc2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random8_glibc2;
#else
extern const gsl_rng_type *gsl_rng_random8_glibc2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random8_libc5;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random8_libc5;
#else
extern const gsl_rng_type *gsl_rng_random8_libc5;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random_bsd;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random_bsd;
#else
extern const gsl_rng_type *gsl_rng_random_bsd;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random_glibc2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random_glibc2;
#else
extern const gsl_rng_type *gsl_rng_random_glibc2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_random_libc5;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_random_libc5;
#else
extern const gsl_rng_type *gsl_rng_random_libc5;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_randu;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_randu;
#else
extern const gsl_rng_type *gsl_rng_randu;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranf;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranf;
#else
extern const gsl_rng_type *gsl_rng_ranf;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlux;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlux;
#else
extern const gsl_rng_type *gsl_rng_ranlux;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlux389;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlux389;
#else
extern const gsl_rng_type *gsl_rng_ranlux389;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlxd1;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlxd1;
#else
extern const gsl_rng_type *gsl_rng_ranlxd1;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlxd2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlxd2;
#else
extern const gsl_rng_type *gsl_rng_ranlxd2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlxs0;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlxs0;
#else
extern const gsl_rng_type *gsl_rng_ranlxs0;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlxs1;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlxs1;
#else
extern const gsl_rng_type *gsl_rng_ranlxs1;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranlxs2;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranlxs2;
#else
extern const gsl_rng_type *gsl_rng_ranlxs2;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_ranmar;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_ranmar;
#else
extern const gsl_rng_type *gsl_rng_ranmar;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_slatec;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_slatec;
#else
extern const gsl_rng_type *gsl_rng_slatec;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_taus;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_taus;
#else
extern const gsl_rng_type *gsl_rng_taus;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_transputer;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_transputer;
#else
extern const gsl_rng_type *gsl_rng_transputer;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_tt800;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_tt800;
#else
extern const gsl_rng_type *gsl_rng_tt800;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_uni;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_uni;
#else
extern const gsl_rng_type *gsl_rng_uni;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_uni32;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_uni32;
#else
extern const gsl_rng_type *gsl_rng_uni32;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_vax;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_vax;
#else
extern const gsl_rng_type *gsl_rng_vax;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_waterman14;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_waterman14;
#else
extern const gsl_rng_type *gsl_rng_waterman14;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_zuf;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_zuf;
#else
extern const gsl_rng_type *gsl_rng_zuf;
#endif

const gsl_rng_type ** gsl_rng_types_setup(void);

#ifdef GSL_EXPORTS
__declspec(dllexport) const gsl_rng_type *gsl_rng_default;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) const gsl_rng_type *gsl_rng_default;
#else
extern const gsl_rng_type *gsl_rng_default;
#endif
#ifdef GSL_EXPORTS
__declspec(dllexport) unsigned long int gsl_rng_default_seed;
#elif defined(GSL_IMPORTS)
__declspec(dllimport) unsigned long int gsl_rng_default_seed;
#else
extern unsigned long int gsl_rng_default_seed;
#endif

gsl_rng *gsl_rng_alloc (const gsl_rng_type * T);
int gsl_rng_memcpy (gsl_rng * dest, const gsl_rng * src);
gsl_rng *gsl_rng_clone (const gsl_rng * r);

void gsl_rng_free (gsl_rng * r);

void gsl_rng_set (const gsl_rng * r, unsigned long int seed);
unsigned long int gsl_rng_max (const gsl_rng * r);
unsigned long int gsl_rng_min (const gsl_rng * r);
const char *gsl_rng_name (const gsl_rng * r);
size_t gsl_rng_size (const gsl_rng * r);
void * gsl_rng_state (const gsl_rng * r);

void gsl_rng_print_state (const gsl_rng * r);

const gsl_rng_type * gsl_rng_env_setup (void);

unsigned long int gsl_rng_get (const gsl_rng * r);
double gsl_rng_uniform (const gsl_rng * r);
double gsl_rng_uniform_pos (const gsl_rng * r);
unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n);


#ifdef HAVE_INLINE
extern inline unsigned long int gsl_rng_get (const gsl_rng * r);

extern inline unsigned long int
gsl_rng_get (const gsl_rng * r)
{
  return (r->type->get) (r->state);
}

extern inline double gsl_rng_uniform (const gsl_rng * r);

extern inline double
gsl_rng_uniform (const gsl_rng * r)
{
  return (r->type->get_double) (r->state);
}

extern inline double gsl_rng_uniform_pos (const gsl_rng * r);

extern inline double
gsl_rng_uniform_pos (const gsl_rng * r)
{
  double x ;
  do
    {
      x = (r->type->get_double) (r->state) ;
    }
  while (x == 0) ;

  return x ;
}

extern inline unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n);

extern inline unsigned long int
gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
{
  unsigned long int offset = r->type->min;
  unsigned long int range = r->type->max - offset;
  unsigned long int scale = range / n;
  unsigned long int k;

  if (n > range) 
    {
      GSL_ERROR_VAL ("n exceeds maximum value of generator",
			GSL_EINVAL, 0) ;
    }

  do
    {
      k = (((r->type->get) (r->state)) - offset) / scale;
    }
  while (k >= n);

  return k;
}
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_RNG_H__ */
