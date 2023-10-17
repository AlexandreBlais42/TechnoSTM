/*
 *  $Id: gwyomp.h 22348 2019-07-31 11:04:26Z yeti-dn $
 *  Copyright (C) 2018-2019 David Necas (Yeti).
 *  E-mail: yeti@gwyddion.net.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
 */

/*< private_header >*/

#ifndef __GWY_OMP_H__
#define __GWY_OMP_H__

#include "config.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <glib.h>
#include <libgwyddion/gwymacros.h>
#include <libgwyddion/gwythreads.h>
#include <libgwyddion/gwyutils.h>

G_BEGIN_DECLS

#ifdef _OPENMP
#include <omp.h>
#define gwy_omp_max_threads() \
    (gwy_threads_are_enabled() ? omp_get_max_threads() : 1)
#define gwy_omp_num_threads() \
    (gwy_threads_are_enabled() ? omp_get_num_threads() : 1)
#define gwy_omp_thread_num omp_get_thread_num
#define gwy_omp_chunk_start(size) \
    (omp_get_thread_num()*(size)/gwy_omp_num_threads())
#define gwy_omp_chunk_end(size) \
    ((omp_get_thread_num() + 1)*(size)/gwy_omp_num_threads())
#else
#define gwy_omp_max_threads() 1
#define gwy_omp_num_threads() 1
#define gwy_omp_thread_num() 0
#define gwy_omp_chunk_start(size) 0
#define gwy_omp_chunk_end(size) (size)
#endif

/*
 * Helpers for common array reductions -- which we do not use from OpenMP
 * because, first, they require OpenMP 4.5 and, second, I observed mixed
 * performance results when I tried them, so I do not trust the compiler
 * support much so far.
 *
 * The allocator reduces to no-op in the single-thread case.  We might want to
 * use an aligned memory allocator to avoid false sharing (cache line
 * cross-talk)...
 *
 * The reduction helpers should be paired with it and free the array the
 * allocator created, making the code transparent.
 */
#define gwy_omp_if_threads_new0(a,n) \
    (gwy_omp_num_threads() > 1 ? g_malloc0(sizeof((a)[0])*(n)) : (a))

#define gwy_omp_if_threads_new(a,n) \
    (gwy_omp_num_threads() > 1 ? g_malloc(sizeof((a)[0])*(n)) : (a))

#define gwy_omp_if_threads_dup(a,n) \
    (gwy_omp_num_threads() > 1 ? g_memdup((a), sizeof((a)[0])*(n)) : (a))

#define gwy_omp_if_threads_free(a,ta) \
    do { \
        if ((ta) != (a)) \
            g_free(ta); \
    } while (0)

static inline void
gwy_omp_if_threads_sum_double(G_GNUC_UNUSED gdouble *a,
                              G_GNUC_UNUSED gdouble *ta,
                              G_GNUC_UNUSED guint n)
{
#ifdef _OPENMP
    guint i;

    if (ta == a)
       return;

#pragma omp critical
    for (i = 0; i < n; i++)
        a[i] += ta[i];

    g_free(ta);
#endif
}

static inline void
gwy_omp_if_threads_sum_uint(G_GNUC_UNUSED guint *a,
                            G_GNUC_UNUSED guint *ta,
                            G_GNUC_UNUSED guint n)
{
#ifdef _OPENMP
    guint i;

    if (ta == a)
       return;

#pragma omp critical
    for (i = 0; i < n; i++)
        a[i] += ta[i];

    g_free(ta);
#endif
}

/* Atomic value helpers, used namely for cancellation flags. */

static inline gboolean
gwy_omp_atomic_read_boolean(gboolean *ptr)
{
    gboolean value;

#ifdef _OPENMP
#pragma omp atomic read
    value = *ptr;
#else
    value = *ptr;
#endif
    return value;
}

/* XXX XXX XXX: GCC prints
 * warning: parameter ‘value’ set but not used
 * without the G_GNUC_UNUSED silencer.  No idea why.  Is this the correct
 * way to do this?  But I essentially replicated atomic2.c from OpenMP
 * examples. */
static inline void
gwy_omp_atomic_write_boolean(gboolean *ptr, G_GNUC_UNUSED gboolean value)
{
#ifdef _OPENMP
#pragma omp atomic write
    *ptr = value;
#else
    *ptr = value;
#endif
}

static inline gint
gwy_omp_atomic_increment_int(gint *ptr)
{
    gint value;

#ifdef _OPENMP
#pragma omp atomic capture
    value = (*ptr)++;
#else
    value = (*ptr)++;
#endif
    return value;
}

static inline guint
gwy_omp_atomic_increment_uint(guint *ptr)
{
    guint value;

#ifdef _OPENMP
#pragma omp atomic capture
    value = (*ptr)++;
#else
    value = (*ptr)++;
#endif
    return value;
}

static inline gboolean
gwy_omp_set_fraction_check_cancel(GwySetFractionFunc set_fraction,
                                  gint i, gint ifrom, gint ito,
                                  gboolean *pcancelled)
{
    /* Cannot cancel without set_fraction() function. */
    if (!set_fraction)
        return FALSE;

    /* Update progress and check cancellation.  Do this only in the master
     * thread.  We do not use pragma master because of nesting limitations.
     * Thread numbers can be used freely. */
    if (!gwy_omp_thread_num() && !set_fraction((i-ifrom + 1.0)/(ito - ifrom)))
        gwy_omp_atomic_write_boolean(pcancelled, TRUE);

    /* In all threads, report the cancellation flag. */
    return gwy_omp_atomic_read_boolean(pcancelled);
}

G_END_DECLS

#endif /* __GWY_OMP_H__ */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
