/*
 *  $Id: gwymath.c 25415 2023-06-06 14:11:08Z yeti-dn $
 *  Copyright (C) 2003-2022 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
 *
 *  The quicksort algorithm was copied from GNU C library, Copyright (C) 1991, 1992, 1996, 1997, 1999 Free Software
 *  Foundation, Inc.  See below.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 *  License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any
 *  later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "config.h"
#include <string.h>
#include <stdlib.h>
#include <libgwyddion/gwymacros.h>
#include <libgwyddion/gwymath.h>
#include <libgwyddion/gwythreads.h>
#include "gwyomp.h"

/* Lower symmetric part indexing */
/* i MUST be greater or equal than j */
#define SLi(a, i, j) a[(i)*((i) + 1)/2 + (j)]

#define DSWAP(x, y) GWY_SWAP(gdouble, x, y)
#define ISWAP(x, y) GWY_SWAP(guint, x, y)

#define GWY_TWO_PI 6.28318530717958647692528676655900576839433879875016

GType
gwy_xy_get_type(void)
{
    /* Threads: type registered from gwy_types_init(). */
    static GType xy_type = 0;

    if (G_UNLIKELY(!xy_type)) {
        xy_type = g_boxed_type_register_static("GwyXY",
                                               (GBoxedCopyFunc)gwy_xy_copy,
                                               (GBoxedFreeFunc)gwy_xy_free);
    }

    return xy_type;
}

/**
 * gwy_xy_new:
 * @x: X-coordinate.
 * @y: Y-coordinate.
 *
 * Creates Cartesian coordinates in plane.
 *
 * This is mostly useful for language bindings.
 *
 * Returns: New XY structure.  The result must be freed using gwy_xy_free(), not g_free().
 *
 * Since: 2.47
 **/
GwyXY*
gwy_xy_new(gdouble x, gdouble y)
{
    GwyXY *xy = g_slice_new(GwyXY);
    xy->x = x;
    xy->y = y;
    return xy;
}

/**
 * gwy_xy_copy:
 * @xy: Cartesian coordinates in plane.
 *
 * Copies Cartesian coordinates in plane.
 *
 * Returns: A copy of @xy. The result must be freed using gwy_xy_free(), not g_free().
 *
 * Since: 2.45
 **/
GwyXY*
gwy_xy_copy(const GwyXY *xy)
{
    g_return_val_if_fail(xy, NULL);
    return g_slice_copy(sizeof(GwyXY), xy);
}

/**
 * gwy_xy_free:
 * @xy: Cartesian coordinates in plane.
 *
 * Frees Cartesian coordinates in plane created with gwy_xy_copy().
 *
 * Since: 2.45
 **/
void
gwy_xy_free(GwyXY *xy)
{
    g_slice_free1(sizeof(GwyXY), xy);
}

GType
gwy_xyz_get_type(void)
{
    /* Threads: type registered from gwy_types_init(). */
    static GType xyz_type = 0;

    if (G_UNLIKELY(!xyz_type)) {
        xyz_type = g_boxed_type_register_static("GwyXYZ",
                                                (GBoxedCopyFunc)gwy_xyz_copy,
                                                (GBoxedFreeFunc)gwy_xyz_free);
    }

    return xyz_type;
}

/**
 * gwy_xyz_new:
 * @x: X-coordinate.
 * @y: Y-coordinate.
 * @z: Z-coordinate.
 *
 * Creates Cartesian coordinates in space.
 *
 * This is mostly useful for language bindings.
 *
 * Returns: New XYZ structure.  The result must be freed using gwy_xyz_free(), not g_free().
 *
 * Since: 2.47
 **/
GwyXYZ*
gwy_xyz_new(gdouble x, gdouble y, gdouble z)
{
    GwyXYZ *xyz = g_slice_new(GwyXYZ);
    xyz->x = x;
    xyz->y = y;
    xyz->z = z;
    return xyz;
}

/**
 * gwy_xyz_copy:
 * @xyz: Cartesian coordinates in space.
 *
 * Copies Cartesian coordinates in space.
 *
 * Returns: A copy of @xyz. The result must be freed using gwy_xyz_free(), not g_free().
 *
 * Since: 2.45
 **/
GwyXYZ*
gwy_xyz_copy(const GwyXYZ *xyz)
{
    g_return_val_if_fail(xyz, NULL);
    return g_slice_copy(sizeof(GwyXYZ), xyz);
}

/**
 * gwy_xyz_free:
 * @xyz: Cartesian coordinates in space.
 *
 * Frees Cartesian coordinates in space created with gwy_xyz_copy().
 *
 * Since: 2.45
 **/
void
gwy_xyz_free(GwyXYZ *xyz)
{
    g_slice_free1(sizeof(GwyXYZ), xyz);
}

/**
 * gwy_math_humanize_numbers:
 * @unit: The smallest possible step.
 * @maximum: The maximum possible value.
 * @precision: A location to store printf() precession, if not %NULL.
 *
 * Finds a human-friendly representation for a range of numbers.
 *
 * Returns: The magnitude i.e., a power of 1000.
 **/
gdouble
gwy_math_humanize_numbers(gdouble unit,
                          gdouble maximum,
                          gint *precision)
{
    gdouble lm, lu, mag, q;

    g_return_val_if_fail(unit >= 0.0, 0.0);
    g_return_val_if_fail(maximum >= 0.0, 0.0);

    if (G_UNLIKELY(unit == 0.0 || maximum == 0.0)) {
        if (unit > 0.0)
            maximum = unit;
        else if (maximum > 0.0)
            unit = maximum;
        else {
            if (precision)
                *precision = 1;
            return 1.0;
        }
    }

    lm = log10(maximum) + 1e-12;
    lu = log10(unit) - 1e-12;
    mag = 3.0*floor(lm/3.0);
    q = 3.0*ceil(lu/3.0);
    if (q > mag)
        q = 3.0*ceil((lu - 1.0)/3.0);
    if (lu > -0.5 && lm < 3.1) {
        while (lu > mag+2)
            mag += 3.0;
    }
    else if (lm <= 0.5 && lm > -1.5) {
        mag = 0.0;
    }
    else {
        while (q > mag)
            mag += 3.0;
    }

    if (precision) {
        *precision = MAX(0, ceil(mag - lu));
        *precision = MIN(*precision, 16);
    }

    return pow10(mag);
}

/**
 * gwy_math_is_in_polygon:
 * @x: The x coordinate of the test point.
 * @y: The y coordinate of the test point.
 * @poly: An array of coordinate pairs (points) that define a polygon.
 * @n: The number of corners of the polygon.
 *
 * Establishes wether the test point @x, @y is inside the polygon @poly. The polygon can be defined either clockwise
 * or anti-clockwise and can be a concave, convex or self-intersecting polygon.
 *
 * <warning> Result can be either TRUE or FALSE if the test point is *exactly* on an edge. </warning>
 *
 * Returns: %TRUE if the test point is inside poly and %FALSE otherwise.
 *
 * Since: 2.7
 **/
/* This neat little check algorithm  was found at http://alienryderflex.com/polygon and has been adapted */
gboolean
gwy_math_is_in_polygon(gdouble x,
                       gdouble y,
                       const gdouble *poly,
                       guint n)
{
    guint i, j = 0;
    gboolean inside = FALSE;
    gdouble xx, yy;

    for (i = 0; i < n; i++) {
        j++;
        if (j == n)
            j = 0;
        if ((poly[2*i + 1] < y && poly[2*j + 1] >= y) || (poly[2*j + 1] < y && poly[2*i + 1] >= y)) {
            xx = poly[2*j] - poly[2*i];
            yy = poly[2*j + 1] - poly[2*i + 1];
            if (poly[2*i] + ((y - poly[2*i + 1])/yy)*xx < x)
                inside = !inside;
        }
    }

    return inside;
}

/**
 * gwy_math_find_nearest_line:
 * @x: X-coordinate of the point to search.
 * @y: Y-coordinate of the point to search.
 * @d2min: Where to store the squared minimal distance, or %NULL.
 * @n: The number of lines (i.e. @coords has 4@n items).
 * @coords: Line coordinates stored as x00, y00, x01, y01, x10, y10, etc.
 * @metric: Metric matrix (2x2, but stored sequentially by rows: m11, m12, m21, m22), it must be positive definite.
 *          Vector norm is then calculated as m11*x*x + (m12 + m21)*x*y + m22*y*y. It can be %NULL, standard Euclidean
 *          metric is then used.
 *
 * Finds the line from @coords nearest to the point (@x, @y).
 *
 * Returns: The line number. It may return -1 if (@x, @y) doesn't lie in the orthogonal stripe of any of the lines.
 **/
gint
gwy_math_find_nearest_line(gdouble x, gdouble y,
                           gdouble *d2min,
                           gint n, const gdouble *coords,
                           const gdouble *metric)
{
    gint i, m;
    gdouble vx, vy, d, d2m = G_MAXDOUBLE;

    g_return_val_if_fail(n > 0, -1);
    g_return_val_if_fail(coords, -1);

    m = -1;
    if (metric) {
        for (i = 0; i < n; i++) {
            gdouble xx = x - (coords[0] + coords[2])/2;
            gdouble yy = y - (coords[1] + coords[3])/2;

            vx = (coords[2] - coords[0])/2;
            vy = (coords[3] - coords[1])/2;
            coords += 4;
            if (vx == 0.0 && vy == 0.0)
                continue;
            d = metric[0]*vx*vx + (metric[1] + metric[2])*vx*vy + metric[3]*vy*vy;
            if (d <= 0.0) {
                g_warning("Metric does not evaluate as positive definite");
                continue;
            }
            d = -(metric[0]*vx*xx + (metric[1] + metric[2])*(vx*yy + vy*xx)/2 + metric[3]*vy*yy)/d;
            /* Out of orthogonal stripe */
            if (d < -1.0 || d > 1.0)
                continue;
            d = metric[0]*(xx + vx*d)*(xx + vx*d) + (metric[1] + metric[2])*(xx + vx*d)*(yy + vy*d)
                + metric[3]*(yy + vy*d)*(yy + vy*d);
            if (d < d2m) {
                d2m = d;
                m = i;
            }
        }
    }
    else {
        for (i = 0; i < n; i++) {
            gdouble xl0 = *(coords++);
            gdouble yl0 = *(coords++);
            gdouble xl1 = *(coords++);
            gdouble yl1 = *(coords++);

            vx = yl1 - yl0;
            vy = xl0 - xl1;
            if (vx == 0.0 && vy == 0.0)
                continue;
            if (vx*(y - yl0) < vy*(x - xl0))
                continue;
            if (vx*(yl1 - y) < vy*(xl1 - x))
                continue;
            d = vx*(x - xl0) + vy*(y - yl0);
            d *= d/(vx*vx + vy*vy);
            if (d < d2m) {
                d2m = d;
                m = i;
            }
        }
    }
    if (d2min)
      *d2min = d2m;

    return m;
}

/**
 * gwy_math_find_nearest_point:
 * @x: X-coordinate of the point to search.
 * @y: Y-coordinate of the point to search.
 * @d2min: Location to store the squared minimal distance to, or %NULL.
 * @n: The number of points (i.e. @coords has 2@n items).
 * @coords: Point coordinates stored as x0, y0, x1, y1, x2, y2, etc.
 * @metric: Metric matrix (2x2, but stored sequentially by rows: m11, m12, m21, m22).  Vector norm is then calculated
 *          as m11*x*x + (m12 + m21)*x*y + m22*y*y. It can be %NULL, standard Euclidean metric is then used.
 *
 * Finds the point from @coords nearest to the point (@x, @y).
 *
 * Returns: The point number.
 **/
gint
gwy_math_find_nearest_point(gdouble x, gdouble y,
                            gdouble *d2min,
                            gint n, const gdouble *coords,
                            const gdouble *metric)
{
    gint i, m;
    gdouble d, xd, yd, d2m = G_MAXDOUBLE;

    g_return_val_if_fail(n > 0, -1);
    g_return_val_if_fail(coords, -1);

    m = 0;
    if (metric) {
        for (i = 0; i < n; i++) {
            xd = *(coords++) - x;
            yd = *(coords++) - y;
            d = metric[0]*xd*xd + (metric[1] + metric[2])*xd*yd + metric[3]*yd*yd;
            if (d < d2m) {
                d2m = d;
                m = i;
            }
        }
    }
    else {
        for (i = 0; i < n; i++) {
            xd = *(coords++) - x;
            yd = *(coords++) - y;
            d = xd*xd + yd*yd;
            if (d < d2m) {
                d2m = d;
                m = i;
            }
        }
    }
    if (d2min)
      *d2min = d2m;

    return m;
}

/**
 * gwy_math_lin_solve:
 * @n: The size of the system.
 * @matrix: The matrix of the system (@n times @n), ordered by row, then column.
 * @rhs: The right hand side of the sytem.
 * @result: Where the result should be stored.  May be %NULL to allocate a fresh array for the result.
 *
 * Solve a regular system of linear equations.
 *
 * Returns: The solution (@result if it wasn't %NULL), may be %NULL if the matrix is singular.
 **/
gdouble*
gwy_math_lin_solve(gint n, const gdouble *matrix,
                   const gdouble *rhs,
                   gdouble *result)
{
    gdouble *m, *r;

    g_return_val_if_fail(n > 0, NULL);
    g_return_val_if_fail(matrix && rhs, NULL);

    m = (gdouble*)g_memdup(matrix, n*n*sizeof(gdouble));
    r = (gdouble*)g_memdup(rhs, n*sizeof(gdouble));
    result = gwy_math_lin_solve_rewrite(n, m, r, result);
    g_free(r);
    g_free(m);

    return result;
}

/**
 * gwy_math_lin_solve_rewrite:
 * @n: The size of the system.
 * @matrix: The matrix of the system (@n times @n), ordered by row, then column.
 * @rhs: The right hand side of the sytem.
 * @result: Where the result should be stored.  May be %NULL to allocate a fresh array for the result.
 *
 * Solves a regular system of linear equations.
 *
 * This is a memory-conservative version of gwy_math_lin_solve() overwriting @matrix and @rhs with intermediate
 * results.
 *
 * Returns: The solution (@result if it wasn't %NULL), may be %NULL if the matrix is singular.
 **/
gdouble*
gwy_math_lin_solve_rewrite(gint n, gdouble *matrix,
                           gdouble *rhs,
                           gdouble *result)
{
    gint *perm;
    gint i, j, jj;

    g_return_val_if_fail(n > 0, NULL);
    g_return_val_if_fail(matrix && rhs, NULL);

    perm = (n <= 12 ? g_newa(gint, n) : g_new(gint, n));

    /* elimination */
    for (i = 0; i < n; i++) {
        gdouble *row = matrix + i*n;
        gdouble piv = 0;
        gint pivj = 0;

        /* find pivot */
        for (j = 0; j < n; j++) {
            if (fabs(row[j]) > piv) {
                pivj = j;
                piv = fabs(row[j]);
            }
        }
        if (piv == 0.0) {
            g_warning("Singluar matrix");
            if (n > 12)
                g_free(perm);
            return NULL;
        }
        piv = row[pivj];
        perm[i] = pivj;

        /* subtract */
        for (j = i+1; j < n; j++) {
            gdouble *jrow = matrix + j*n;
            gdouble q = jrow[pivj]/piv;

            for (jj = 0; jj < n; jj++)
                jrow[jj] -= q*row[jj];

            jrow[pivj] = 0.0;
            rhs[j] -= q*rhs[i];
        }
    }

    /* back substitute */
    if (!result)
        result = g_new(gdouble, n);
    for (i = n-1; i >= 0; i--) {
        gdouble *row = matrix + i*n;
        gdouble x = rhs[i];

        for (j = n-1; j > i; j--)
            x -= result[perm[j]]*row[perm[j]];

        result[perm[i]] = x/row[perm[i]];
    }
    if (n > 12)
        g_free(perm);

    return result;
}

/**
 * gwy_math_fit_polynom:
 * @ndata: The number of items in @xdata, @ydata.
 * @xdata: Independent variable data (of size @ndata).
 * @ydata: Dependent variable data (of size @ndata).
 * @n: The degree of polynom to fit.
 * @coeffs: An array of size @n+1 to store the coefficients to, or %NULL (a fresh array is allocated then).
 *
 * Fits a polynom through a general (x, y) data set.
 *
 * Returns: The coefficients of the polynom (@coeffs when it was not %NULL, otherwise a newly allocated array).
 **/
gdouble*
gwy_math_fit_polynom(gint ndata,
                     const gdouble *xdata,
                     const gdouble *ydata,
                     gint n,
                     gdouble *coeffs)
{
    gdouble *sumx, *m;
    gint i, j;

    g_return_val_if_fail(ndata >= 0, NULL);
    g_return_val_if_fail(n >= 0, NULL);

    sumx = g_new0(gdouble, 2*n+1);

    if (!coeffs)
        coeffs = g_new0(gdouble, n+1);
    else
        gwy_clear(coeffs, n+1);

    for (i = 0; i < ndata; i++) {
        gdouble x = xdata[i];
        gdouble y = ydata[i];
        gdouble xp;

        xp = 1.0;
        for (j = 0; j <= n; j++) {
            sumx[j] += xp;
            coeffs[j] += xp*y;
            xp *= x;
        }
        for (j = n+1; j <= 2*n; j++) {
            sumx[j] += xp;
            xp *= x;
        }
    }

    m = g_new(gdouble, (n+1)*(n+2)/2);
    for (i = 0; i <= n; i++) {
        gdouble *row = m + i*(i+1)/2;

        for (j = 0; j <= i; j++)
            row[j] = sumx[i+j];
    }
    if (!gwy_math_choleski_decompose(n+1, m))
        gwy_clear(coeffs, n+1);
    else
        gwy_math_choleski_solve(n+1, m, coeffs);

    g_free(m);
    g_free(sumx);

    return coeffs;
}

/**
 * gwy_math_choleski_decompose:
 * @n: The dimension of @a.
 * @matrix: Lower triangular part of a symmetric matrix, stored by rows, i.e., matrix = [a_00 a_10 a_11 a_20 a_21 a_22
 *          a_30 ...].
 *
 * Decomposes a symmetric positive definite matrix in place.
 *
 * Returns: Whether the matrix was really positive definite.  If %FALSE, the decomposition failed and @a does not
 *          contain any meaningful values.
 **/
gboolean
gwy_math_choleski_decompose(gint dim, gdouble *a)
{
    gint i, j, k;
    gdouble s, r;

    for (k = 0; k < dim; k++) {
        /* diagonal element */
        s = SLi(a, k, k);
        for (i = 0; i < k; i++)
            s -= SLi(a, k, i) * SLi(a, k, i);
        if (s <= 0.0)
            return FALSE;
        SLi(a, k, k) = s = sqrt(s);

        /* nondiagonal elements */
        for (j = k+1; j < dim; j++) {
            r = SLi(a, j, k);
            for (i = 0; i < k; i++)
                r -= SLi(a, k, i) * SLi(a, j, i);
            SLi(a, j, k) = r/s;
        }
    }

    return TRUE;
}

/**
 * gwy_math_choleski_invert:
 * @n: Matrix size.
 * @matrix: Lower triangular part of a symmetric matrix, stored by rows, i.e., matrix = [a_00 a_10 a_11 a_20 a_21 a_22
 *          a_30 ...].
 *
 * Inverts a symmetric positive definite matrix in place.
 *
 * Returns: Whether the matrix was really positive definite.  If %FALSE, the inversion failed and @matrix does not
 *          contain any meaningful values.
 *
 * Since: 2.46
 **/
gboolean
gwy_math_choleski_invert(gint n, gdouble *a)
{

    gint q = 0, m;
    gdouble s, t;
    gdouble *x;
    gint k, i, j;

    x = g_newa(gdouble, n);
    for (k = n-1; k >= 0; k--) {
        s = a[0];
        if (s <= 0)
            return FALSE;
        m = 0;
        for (i = 0; i < n-1; i++) {
            q = m+1;
            m += i+2;
            t = a[q];
            x[i] = -t/s;      /* note use temporary x */
            if (i >= k)
                x[i] = -x[i];
            for (j = q; j < m; j++)
                a[j - (i+1)] = a[j+1] + t * x[j - q];
        }
        a[m] = 1.0/s;
        for (i = 0; i < n-1; i++)
            a[q + i] = x[i];
    }

    return TRUE;
}

/**
 * gwy_math_choleski_solve:
 * @n: The dimension of @a.
 * @decomp: Lower triangular part of Choleski decomposition as computed by gwy_math_choleski_decompose().
 * @rhs: Right hand side vector.  Is is modified in place, on return it contains the solution.
 *
 * Solves a system of linear equations with predecomposed symmetric positive definite matrix @a and right hand side
 * @b.
 **/
void
gwy_math_choleski_solve(gint dim, const gdouble *a, gdouble *b)
{
    gint i, j;

    /* back-substitution with the lower triangular matrix */
    for (j = 0; j < dim; j++) {
        for (i = 0; i < j; i++)
            b[j] -= SLi(a, j, i)*b[i];
        b[j] /= SLi(a, j, j);
    }

    /* back-substitution with the upper triangular matrix */
    for (j = dim-1; j >= 0; j--) {
        for (i = j+1; i < dim; i++)
            b[j] -= SLi(a, i, j)*b[i];
        b[j] /= SLi(a, j, j);
    }
}

/**
 * gwy_math_tridiag_solve_rewrite:
 * @n: The dimension of @d.
 * @d: The diagonal of a tridiagonal matrix, its contents will be overwritten.
 * @a: The above-diagonal stripe (it has @n-1 elements).
 * @b: The below-diagonal stripe (it has @n-1 elements).
 * @rhs: The right hand side of the system, upon return it will contain the solution.
 *
 * Solves a tridiagonal system of linear equations.
 *
 * Returns: %TRUE if the elimination suceeded, %FALSE if the system is (numerically) singular.  The contents of @d and
 *          @rhs may be overwritten in the case of failure too, but not to any meaningful values.
 **/
gboolean
gwy_math_tridiag_solve_rewrite(gint n,
                               gdouble *d,
                               const gdouble *a,
                               const gdouble *b,
                               gdouble *rhs)
{
    gint i;

    g_return_val_if_fail(n > 0, FALSE);

    /* Eliminate b[elow diagonal] */
    for (i = 0; i < n-1; i++) {
        /* If d[i] is zero, elimination fails (now or later) */
        if (!d[i])
            return FALSE;
        d[i+1] -= b[i]/d[i]*a[i];
        rhs[i+1] -= b[i]/d[i]*rhs[i];
    }
    if (!d[n-1])
        return FALSE;

    /* Eliminate a[bove diagonal], calculating the solution */
    for (i = n-1; i > 0; i--) {
        rhs[i] /= d[i];
        rhs[i-1] -= a[i-1]*rhs[i];
    }
    rhs[0] /= d[0];

    return TRUE;
}

static inline void
order_3(gdouble *array)
{
    GWY_ORDER(gdouble, array[0], array[1]);
    if (array[2] < array[1]) {
        DSWAP(array[1], array[2]);
        GWY_ORDER(gdouble, array[0], array[1]);
    }
}

static gdouble
kth_rank_simple(gsize n, gdouble *array, gsize k)
{
    gdouble a, b, c, d;
    gsize i;

    if (n == 1)
        return array[0];

    if (n == 2) {
        GWY_ORDER(gdouble, array[0], array[1]);
        return array[k];
    }

    if (n == 3 && k == 1) {
        order_3(array);
        return array[1];
    }

    if (k == 0) {
        a = array[0];
        for (i = 1; i < n; i++) {
            c = array[i];
            if (c < a) {
                array[i] = a;
                array[0] = a = c;
            }
        }
        return a;
    }

    if (k == n-1) {
        a = array[n-1];
        for (i = 0; i < n-1; i++) {
            c = array[i];
            if (c > a) {
                array[i] = a;
                array[n-1] = a = c;
            }
        }
        return a;
    }

    if (k == 1) {
        GWY_ORDER(gdouble, array[0], array[1]);
        a = array[0];
        b = array[1];

        for (i = 2; i < n; i++) {
            c = array[i];
            if (c < b) {
                if (c < a) {
                    array[i] = b;
                    array[1] = b = a;
                    array[0] = a = c;
                }
                else {
                    array[i] = b;
                    array[1] = b = c;
                }
            }
        }
        return b;
    }

    if (k == n-2) {
        GWY_ORDER(gdouble, array[n-2], array[n-1]);
        a = array[n-1];
        b = array[n-2];

        for (i = 0; i < n-2; i++) {
            c = array[i];
            if (c > b) {
                if (c > a) {
                    array[i] = b;
                    array[n-2] = b = a;
                    array[n-1] = a = c;
                }
                else {
                    array[i] = b;
                    array[n-2] = b = c;
                }
            }
        }
        return b;
    }

    if (k == 2) {
        order_3(array);
        a = array[0];
        b = array[1];
        c = array[2];

        for (i = 3; i < n; i++) {
            d = array[i];
            if (d < c) {
                if (d < b) {
                    if (d < a) {
                        array[i] = c;
                        array[2] = c = b;
                        array[1] = b = a;
                        array[0] = a = d;
                    }
                    else {
                        array[i] = c;
                        array[2] = c = b;
                        array[1] = b = d;
                    }
                }
                else {
                    array[i] = c;
                    array[2] = c = d;
                }
            }
        }
        return c;
    }

    if (k == n-3) {
        order_3(array + n-3);
        a = array[n-1];
        b = array[n-2];
        c = array[n-3];

        for (i = 0; i < n-3; i++) {
            d = array[i];
            if (d > c) {
                if (d > b) {
                    if (d > a) {
                        array[i] = c;
                        array[n-3] = c = b;
                        array[n-2] = b = a;
                        array[n-1] = a = d;
                    }
                    else {
                        array[i] = c;
                        array[n-3] = c = b;
                        array[n-2] = b = d;
                    }
                }
                else {
                    array[i] = c;
                    array[n-3] = c = d;
                }
            }
        }
        return c;
    }

    g_assert_not_reached();
}

/**
 * gwy_math_kth_rank:
 * @n: Number of items in @array.
 * @array: Array of doubles.  It is shuffled by this function.
 * @k: Rank of the value to find (from lowest to highest).
 *
 * Finds k-th item of an array of values using Quick select algorithm.
 *
 * The value positions change as follows.  The returned value is guaranteed to be at @k-th position in the array (i.e.
 * correctly ranked).  All other values are correctly ordered with respect to this value: preceeding values are
 * smaller (or equal) and following values are larger (or equal).
 *
 * Returns: The @k-th value of @array if it was sorted.
 *
 * Since: 2.50
 **/
gdouble
gwy_math_kth_rank(gsize n, gdouble *array, gsize k)
{
    gsize lo, hi;
    gsize middle, ll, hh;
    gdouble m;

    g_return_val_if_fail(k < n, 0.0);

    lo = 0;
    hi = n-1;
    while (TRUE) {
        if (hi <= lo+2 || k <= lo+2 || k+2 >= hi)
            return kth_rank_simple(hi+1 - lo, array + lo, k - lo);

        /* Find median of lo, middle and hi items; swap into position lo */
        middle = (lo + hi)/2;

        GWY_ORDER(gdouble, array[middle], array[hi]);
        GWY_ORDER(gdouble, array[lo], array[hi]);
        GWY_ORDER(gdouble, array[middle], array[lo]);

        /* Swap low item (now in position middle) into position (lo+1) */
        DSWAP(array[middle], array[lo+1]);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = lo+1;
        hh = hi;
        m = array[lo];
        while (TRUE) {
            do {
                ll++;
            } while (m > array[ll]);
            do {
                hh--;
            } while (array[hh] > m);

            if (hh < ll)
                break;

            DSWAP(array[ll], array[hh]);
        }

        /* Swap middle item (in position lo) back into correct position */
        array[lo] = array[hh];
        array[hh] = m;

        /* Re-set active partition */
        if (hh <= k)
            lo = hh;
        if (hh >= k)
            hi = hh-1;
    }
}

/**
 * gwy_math_median:
 * @n: Number of items in @array.
 * @array: Array of doubles.  It is shuffled by this function.
 *
 * Finds median of an array of values using Quick select algorithm.
 *
 * See gwy_math_kth_rank() for details of how the values are shuffled.
 *
 * Returns: The median value of @array.
 **/
gdouble
gwy_math_median(gsize n, gdouble *array)
{
    return gwy_math_kth_rank(n, array, n/2);
}

/* When there are many values to find, just sort the entire thing and read the values at the corresponding ranks. */
static void
kth_ranks_brute(gsize n, gdouble *array,
                guint nk, const guint *k, gdouble *values)
{
    guint j;

    gwy_math_sort(n, array);
    for (j = 0; j < nk; j++)
        values[j] = array[k[j]];
}

/* We assume k[] is uniq-sorted which eliminates a bunch of cases. */
static void
kth_ranks_small(gsize n, gdouble *array,
                guint nk, const guint *k, gdouble *values)
{
    guint k0, k1, d0, d1;

    if (!nk)
        return;

    k0 = k[0];
    if (nk == 1) {
        values[0] = gwy_math_kth_rank(n, array, k0);
        return;
    }

    k1 = k[1];
    d0 = (k0 <= n/2) ? n/2 - k0 : k0 - n/2;
    d1 = (k1 <= n/2) ? n/2 - k1 : k1 - n/2;
    if (d0 <= d1) {
        values[0] = gwy_math_kth_rank(n, array, k0);
        k0++;
        values[1] = gwy_math_kth_rank(n-k0, array+k0, k1-k0);
    }
    else {
        values[1] = gwy_math_kth_rank(n, array, k1);
        values[0] = gwy_math_kth_rank(k1, array, k0);
    }
}

static void
kth_ranks_recurse(gsize n, gdouble *array,
                  guint nk, guint *k, gdouble *values)
{
    guint jmid, kmid, j;

    if (nk <= 2) {
        kth_ranks_small(n, array, nk, k, values);
        return;
    }

    jmid = nk/2;
    kmid = k[jmid];
    values[jmid] = gwy_math_kth_rank(n, array, kmid);

    /* Now recurse into the halfs.  Both are non-empty because nk >= 3. */
    kth_ranks_recurse(kmid, array, jmid, k, values);

    jmid++;
    for (j = jmid; j < nk; j++)
        k[j] -= kmid+1;
    kth_ranks_recurse(n-1-kmid, array+kmid+1, nk-jmid, k+jmid, values+jmid);
    for (j = jmid; j < nk; j++)
        k[j] += kmid+1;
}

static guint
bisect_lower_guint(const guint *a, guint n, guint x)
{
    guint lo = 0, hi = n-1;

    if (G_UNLIKELY(x < a[lo]))
        return 0;
    if (G_UNLIKELY(x >= a[hi]))
        return n-1;

    while (hi - lo > 1) {
        guint mid = (hi + lo)/2;

        if (x < a[mid])
            hi = mid;
        else
            lo = mid;
    }

    return lo;
}

static void
kth_ranks_fastpath(gsize n, gdouble *array,
                   guint nk, const guint *k, gdouble *values)
{
    if (nk < 2 || k[0] < k[1])
        kth_ranks_small(n, array, nk, k, values);
    else if (k[0] == k[1])
        values[0] = values[1] = gwy_math_kth_rank(n, array, k[0]);
    else {
        guint ksorted[2];

        ksorted[0] = k[1];
        ksorted[1] = k[0];
        kth_ranks_small(n, array, nk, ksorted, values);
        DSWAP(values[0], values[1]);
    }
}

/**
 * gwy_math_kth_ranks:
 * @n: Number of items in @array.
 * @array: Array of doubles.  It is shuffled by this function.
 * @nk: Number of ranked values to find (sizes of arrays @k and @values).
 * @k: Ranks of the value to find.
 * @values: Array where to store values with ranks (from smallest to highest) given in @k.
 *
 * Finds simultaneously several k-th items of an array of values.
 *
 * The values are shuffled similarly to gwy_math_kth_rank(), except that the guarantee holds for all given ranks
 * simultaneously.  All values with explicitly requested ranks are at their correct positions and all values lying
 * between them in the array are also between them numerically.
 *
 * Since: 2.50
 **/
void
gwy_math_kth_ranks(gsize n, gdouble *array,
                   guint nk, const guint *k, gdouble *values)
{
    guint t, j, nkred;
    gdouble logn;
    guint *ksorted;
    gdouble *valsorted;

    for (j = 0; j < nk; j++) {
        g_return_if_fail(k[j] < n);
    }

    if (nk <= 2) {
        kth_ranks_fastpath(n, array, nk, k, values);
        return;
    }
    if (n < 30) {
        kth_ranks_brute(n, array, nk, k, values);
        return;
    }

    logn = log(n);
    if (nk > 0.12*exp(0.3*logn)*logn*logn) {
        kth_ranks_brute(n, array, nk, k, values);
        return;
    }

    if (nk <= 64) {
        ksorted = g_newa(guint, nk);
        valsorted = g_newa(gdouble, nk);
    }
    else {
        ksorted = g_new(guint, nk);
        valsorted = g_new(gdouble, nk);
    }

    /* Do uniq-sort on the k values.  The uniq is more to avoid some odd cases in the recursion than for efficiency. */
    gwy_assign(ksorted, k, nk);
    gwy_guint_sort(nk, ksorted);
    t = 0;
    for (j = 0; j < nk; j++) {
        if (ksorted[j] != ksorted[t]) {
            t++;
            ksorted[t] = ksorted[j];
        }
    }
    nkred = t+1;
    /* The recursion can be at most log2(nkred) deep. */
    kth_ranks_recurse(n, array, nkred, ksorted, valsorted);
    /* Assign the values to the original array. */
    for (j = 0; j < nk; j++)
        values[j] = valsorted[bisect_lower_guint(ksorted, nkred, k[j])];

    if (nk > 64) {
        g_free(ksorted);
        g_free(valsorted);
    }
}

/**
 * gwy_math_trimmed_mean:
 * @n: Number of items in @array.
 * @array: Array of doubles.  It is shuffled by this function.
 * @nlowest: The number of lowest values to discard.
 * @nhighest: The number of highest values to discard.
 *
 * Finds trimmed mean of an array of values.
 *
 * At least one value must remain after the trimming, i.e. @nlowest + @nhighest must be smaller than @n.  Usually one
 * passes the same number as both @nlowest and @nhighest, but it is not a requirement.
 *
 * The function can be also used to calculate normal mean values as it implements efficiently the cases when no
 * trimming is done at either end.
 *
 * Returns: The trimmed mean.
 *
 * Since: 2.50
 **/
gdouble
gwy_math_trimmed_mean(gsize n, gdouble *array, guint nlowest, guint nhighest)
{
    gsize i, nred;
    gdouble s;

    g_return_val_if_fail(nlowest < n, 0.0);
    g_return_val_if_fail(nhighest < n - nlowest, 0.0);

    /* Note that when using the k-th rank functions, we ask for the first discarded item, and hence ignore it the
     * returned values. */
    if (!nlowest) {
        if (!nhighest)
            nred = n;
        else {
            nred = n-nhighest;
            gwy_math_kth_rank(n, array, nred);
        }
    }
    else if (!nhighest) {
        nred = n-nlowest;
        gwy_math_kth_rank(n, array, nlowest-1);
        array += nlowest;
    }
    else {
        guint k[2];
        gdouble v[2];

        k[0] = nlowest-1;
        k[1] = n-nhighest;
        nred = n - (nlowest + nhighest);
        gwy_math_kth_ranks(n, array, 2, k, v);
        array += nlowest;
    }

    /* Now the reduced array block is guaranteed to contain only non-discarded values. So just calculate the average. */
    s = 0.0;
    for (i = nred; i; i--, array++)
        s += *array;
    return s/nred;
}

static inline guint
percentile_to_rank(gsize n, gdouble p, GwyPercentileInterpolationType interp)
{
    gdouble kreal, x, eps;
    guint k;

    g_return_val_if_fail(p >= 0.0, 0);
    g_return_val_if_fail(p <= 100.0, n-1);

    kreal = p/100.0*(n - 1);
    k = (guint)floor(kreal);
    x = kreal - k;
    eps = n*3e-16;
    if (interp == GWY_PERCENTILE_INTERPOLATION_LOWER || x <= eps)
        return k;
    if (interp == GWY_PERCENTILE_INTERPOLATION_HIGHER || x >= 1.0 - eps)
        return k+1;
    /* Nearest. */
    return x < 0.5 ? k : k+1;
}

/**
 * gwy_math_percentiles:
 * @n: Number of items in @array.
 * @array: Array of doubles.  It is shuffled by this function.
 * @interp: Interpolation method to use for percentiles that do not correspond exactly to an integer rank.
 * @np: Number of percentiles to find.
 * @p: Array of size @np with the percentiles to compute.  The values are in percents, i.e. from the range [0,100].
 * @values: Array where to store values with percentiles given in @p.
 *
 * Finds simultaneously several percentiles of an array of values.
 *
 * The values in @array are shuffled similarly to gwy_math_kth_ranks(). However, it is difficult to state how exactly
 * @p translates to the values that become correctly ranked (and it depends on @interp).  Hence you can only assume
 * the set of values is preserved.
 *
 * Since: 2.50
 **/
void
gwy_math_percentiles(gsize n, gdouble *array,
                     GwyPercentileInterpolationType interp,
                     guint np, const gdouble *p, gdouble *values)
{
    gdouble *v2;
    guint *k;
    guint j;

    if (interp == GWY_PERCENTILE_INTERPOLATION_LOWER
        || interp == GWY_PERCENTILE_INTERPOLATION_HIGHER
        || interp == GWY_PERCENTILE_INTERPOLATION_NEAREST) {
        k = g_new(guint, np);
        for (j = 0; j < np; j++)
            k[j] = percentile_to_rank(n, p[j], interp);

        gwy_math_kth_ranks(n, array, np, k, values);
        g_free(k);
        return;
    }

    k = g_new(guint, 2*np);
    v2 = g_new(gdouble, 2*np);

    /* When p corresponds exactly to an integer rank we store the same k twice.
     * But gwy_math_kth_ranks() does a uniq-sort of k[] anyway. */
    for (j = 0; j < np; j++) {
        k[2*j] = percentile_to_rank(n, p[j], GWY_PERCENTILE_INTERPOLATION_LOWER);
        k[2*j + 1] = percentile_to_rank(n, p[j], GWY_PERCENTILE_INTERPOLATION_HIGHER);
    }
    gwy_math_kth_ranks(n, array, 2*np, k, v2);
    for (j = 0; j < np; j++) {
        gdouble vlower = v2[2*j], vupper = v2[2*j + 1], x;

        if (k[2*j + 1] == k[2*j])
            values[j] = vlower;
        else if (interp == GWY_PERCENTILE_INTERPOLATION_MIDPOINT)
            values[j] = 0.5*(vlower + vupper);
        else {
            g_assert(k[2*j + 1] == k[2*j] + 1);
            x = p[j]/100.0*(n - 1) - k[2*j];
            values[j] = x*vupper + (1.0 - x)*vlower;
        }
    }

    g_free(k);
    g_free(v2);
}

/**
 * gwy_math_curvature_at_apex:
 * @coeffs: Array of the six polynomial coefficients of a quadratic surface in the following order: 1, x, y, x², xy, y².
 * @kappa1: Location to store the smaller curvature to.
 * @kappa2: Location to store the larger curvature to.
 * @phi1: Location to store the direction of the smaller curvature to.
 * @phi2: Location to store the direction of the larger curvature to.
 * @xc: Location to store x-coordinate of the centre of the quadratic surface.
 * @yc: Location to store y-coordinate of the centre of the quadratic surface.
 * @zc: Location to store value at the centre of the quadratic surface.
 *
 * Calculates curvature parameters at the apex from two-dimensional quadratic polynomial coefficients.
 *
 * See also gwy_math_curvature_at_origin() which computes the local surface curvature at @x=0 and @y=0.
 *
 * If the quadratic surface was obtained by fitting the dimensions of the fitted area should not differ, in the
 * lateral coordinates used, by many orders from 1.  Otherwise the recognition of flat surfaces might not work.
 *
 * Curvatures have signs, positive mean a concave (cup-like) surface, negative mean a convex (cap-like) surface.  They
 * are ordered including the sign.
 *
 * Directions are angles from the interval (-π/2, π/2].
 *
 * If the quadratic surface is degenerate, i.e. flat in at least one direction, the centre is undefined.  The centre
 * is then chosen as the closest point the origin of coordinates.  For flat surfaces this means the origin is simply
 * returned as the centre position.  Consequently, you should use Cartesian coordinates with origin in a natural
 * centre, for instance centre of image or grain.
 *
 * Returns: The number of curved dimensions (0 to 2).
 *
 * Since: 2.61
 **/
guint
gwy_math_curvature_at_apex(const gdouble *coeffs,
                           gdouble *pkappa1,
                           gdouble *pkappa2,
                           gdouble *pphi1,
                           gdouble *pphi2,
                           gdouble *pxc,
                           gdouble *pyc,
                           gdouble *pzc)
{
    gdouble a = coeffs[0], bx = coeffs[1], by = coeffs[2], cxx = coeffs[3], cxy = coeffs[4], cyy = coeffs[5];
    gdouble cm, cp, bx1, by1, xc, yc, phi, cx, cy;
    guint degree;

    /* Check if any sane quadratic term is present. */
    if (fabs(cxx) + fabs(cxy) + fabs(cyy) <= 1e-14*(fabs(bx) + fabs(by))) {
        /* Linear gradient */
        if (pkappa1)
            *pkappa1 = 0.0;
        if (pkappa2)
            *pkappa2 = 0.0;
        if (pxc)
            *pxc = 0.0;
        if (pyc)
            *pyc = 0.0;
        if (pzc)
            *pzc = a;
        if (pphi1)
            *pphi1 = 0.0;
        if (pphi2)
            *pphi2 = G_PI_2;
        return 0;
    }

    /* At least one quadratic term.  Eliminate the mixed term, if any. */
    cm = cxx - cyy;
    cp = cxx + cyy;
    phi = 0.5*atan2(cxy, cm);
    cx = cp + hypot(cm, cxy);
    cy = cp - hypot(cm, cxy);
    bx1 = bx*cos(phi) + by*sin(phi);
    by1 = -bx*sin(phi) + by*cos(phi);

    /* Eliminate linear terms */
    if (fabs(cx) < 1e-14*fabs(cy)) {
        /* Only y quadratic term */
        xc = 0.0;
        yc = -by1/cy;
        degree = 1;
    }
    else if (fabs(cy) < 1e-14*fabs(cx)) {
        /* Only x quadratic term */
        xc = -bx1/cx;
        yc = 0.0;
        degree = 1;
    }
    else {
        /* Two quadratic terms */
        xc = -bx1/cx;
        yc = -by1/cy;
        degree = 2;
    }

    if (pxc)
        *pxc = xc*cos(phi) - yc*sin(phi);
    if (pyc)
        *pyc = xc*sin(phi) + yc*cos(phi);
    if (pzc)
        *pzc = a + xc*bx1 + yc*by1 + xc*xc*cx + yc*yc*cy;

    if (cx > cy) {
        GWY_SWAP(gdouble, cx, cy);
        phi += G_PI_2;
    }
    /* Compenstate the left-handed coordinate system */
    phi = -phi;

    if (pkappa1)
        *pkappa1 = cx;
    if (pkappa2)
        *pkappa2 = cy;

    if (pphi1)
        *pphi1 = gwy_canonicalize_angle(phi, FALSE, FALSE);
    if (pphi2)
        *pphi2 = gwy_canonicalize_angle(phi + G_PI_2, FALSE, FALSE);

    return degree;
}

/**
 * gwy_math_curvature:
 * @coeffs: Array of the six polynomial coefficients of a quadratic surface in the following order: 1, x, y, x², xy, y².
 * @kappa1: Location to store the smaller curvature to.
 * @kappa2: Location to store the larger curvature to.
 * @phi1: Location to store the direction of the smaller curvature to.
 * @phi2: Location to store the direction of the larger curvature to.
 * @xc: Location to store x-coordinate of the centre of the quadratic surface.
 * @yc: Location to store y-coordinate of the centre of the quadratic surface.
 * @zc: Location to store value at the centre of the quadratic surface.
 *
 * Calculates curvature parameters at the apex from two-dimensional quadratic polynomial coefficients.
 *
 * This is an old name for gwy_math_curvature_at_apex().  See the description there.
 *
 * Returns: The number of curved dimensions (0 to 2).
 *
 * Since: 2.22
 **/
guint
gwy_math_curvature(const gdouble *coeffs,
                   gdouble *pkappa1,
                   gdouble *pkappa2,
                   gdouble *pphi1,
                   gdouble *pphi2,
                   gdouble *pxc,
                   gdouble *pyc,
                   gdouble *pzc)
{
    return gwy_math_curvature_at_apex(coeffs, pkappa1, pkappa2, pphi1, pphi2, pxc, pyc, pzc);
}

/**
 * gwy_math_curvature_at_origin:
 * @coeffs: Array of the six polynomial coefficients of a quadratic surface in the following order: 1, x, y, x², xy, y².
 * @kappa1: Location to store the smaller curvature to.
 * @kappa2: Location to store the larger curvature to.
 * @phi1: Location to store the direction of the smaller curvature to.
 * @phi2: Location to store the direction of the larger curvature to.
 *
 * Calculates curvature parameters at origin from two-dimensional quadratic polynomial coefficients.
 *
 * See gwy_math_curvature() for discussion of scaling and sign convenrions.  This function function differs from it
 * by computing the local surface curvature at @x=0 and @y=0, whereas gwy_math_curvature() computes the curvature at
 * the apex of the parabolic surface.
 *
 * The array @coeffs is consistent with gwy_math_curvature_at_apex(), even though here the constant term is not used.
 *
 * Returns: The number of curved dimensions (0 to 2).
 *
 * Since: 2.61
 **/
guint
gwy_math_curvature_at_origin(const gdouble *coeffs,
                             gdouble *pkappa1,
                             gdouble *pkappa2,
                             gdouble *pphi1,
                             gdouble *pphi2)
{
    gdouble bx = coeffs[1], by = coeffs[2], cxx = coeffs[3], cxy = coeffs[4], cyy = coeffs[5];
    gdouble b2, tpcoeffs[6];

    /* Transform to a coordinate system in the local tangent plane to the surface (x points along the gradient, y is
     * horizontal, z is tilted to be normal to the surface). We only care about curvatures so we set all the
     * uninteresting coefficients to zeros and use gwy_math_curvature_at_apex(), knowing the apex and origin now
     * coincide. */
    tpcoeffs[0] = tpcoeffs[1] = tpcoeffs[2] = 0.0;
    b2 = bx*bx + by*by;
    if (b2 == 0.0)
        gwy_assign(tpcoeffs + 3, coeffs + 3, 3);
    else {
        gdouble beta2 = 1.0 + b2, beta = sqrt(beta2);
        tpcoeffs[3] = (bx*bx*cxx + bx*by*cxy + by*by*cyy)/(beta2*beta*b2);
        tpcoeffs[4] = (2.0*bx*by*(cyy - cxx) + (bx*bx - by*by)*cxy)/(beta2*b2);
        tpcoeffs[5] = (by*by*cxx - bx*by*cxy + bx*bx*cyy)/(beta*b2);
    }

    return gwy_math_curvature_at_apex(tpcoeffs, pkappa1, pkappa2, pphi1, pphi2, NULL, NULL, NULL);
}

/**
 * gwy_math_refine_maximum_1d:
 * @y: Array of length 3, containing the neighbourhood values with the maximum in the centre.
 * @x: Location to store the refined @x-coordinate.
 *
 * Performs subpixel refinement of parabolic a one-dimensional maximum.
 *
 * The central value corresponds to x-coordinate 0, distances between values are unity.  The refinement is based by
 * fitting a parabola through the maximum.  If it fails or the calculated maximum lies farther than the surrounding
 * values the function sets the refined maximum to the origin and returns %FALSE.
 *
 * Returns: %TRUE if the refinement succeeded, %FALSE if it failed.  The value of @x is usable regardless of the
 *          return value.
 *
 * Since: 2.49
 **/
gboolean
gwy_math_refine_maximum_1d(const gdouble *y, gdouble *x)
{
    gdouble b, D;

    *x = 0.0;
    D = y[2] + y[0] - 2.0*y[1];
    b = 0.5*(y[0] - y[2]);
    if (D == 0.0 || fabs(D) < fabs(b))
        return FALSE;

    *x = b/D;
    return TRUE;
}

/**
 * gwy_math_refine_maximum_2d:
 * @z: Array of length 9, containing the square 3x3 neighbourhood values in matrix order and with the maximum in the
 *     centre.
 * @x: Location to store the refined @x-coordinate.
 * @y: Location to store the refined @y-coordinate.
 *
 * Performs subpixel refinement of parabolic a two-dimensional maximum.
 *
 * The central value corresponds to coordinates (0,0), distances between values are unity.  The refinement is based by
 * fitting a two-dimensional parabola through the maximum.  If it fails or the calculated maximum lies farther than
 * the surrounding values the function sets the refined maximum to the origin and returns %FALSE.
 *
 * Returns: %TRUE if the refinement succeeded, %FALSE if it failed.  The values of @x and @y are usable regardless of
 *          the return value.
 *
 * Since: 2.49
 **/
gboolean
gwy_math_refine_maximum_2d(const gdouble *z,
                           gdouble *x, gdouble *y)
{
    gdouble sz, szx, szy, szxx, szxy, szyy;
    gdouble bx, by, cxx, cxy, cyy, D, sx, sy;
    gdouble m[6], rhs[3];

    *x = *y = 0;

    sz = z[0] + z[1] + z[2] + z[3] + z[4] + z[5] + z[6] + z[7] + z[8];
    szx = -z[0] + z[2] - z[3] + z[5] - z[6] + z[8];
    szy = -z[0] - z[1] - z[2] + z[6] + z[7] + z[8];
    szxx = z[0] + z[2] + z[3] + z[5] + z[6] + z[8];
    szxy = z[0] - z[2] - z[6] + z[8];
    szyy = z[0] + z[1] + z[2] + z[6] + z[7] + z[8];

    m[0] = 9.0;
    m[1] = m[2] = m[3] = m[5] = 6.0;
    m[4] = 4.0;
    gwy_math_choleski_decompose(3, m);

    rhs[0] = sz;
    rhs[1] = szxx;
    rhs[2] = szyy;
    gwy_math_choleski_solve(3, m, rhs);

    bx = szx/6.0;
    by = szy/6.0;
    cxx = rhs[1];
    cxy = szxy/4.0;
    cyy = rhs[2];

    D = 4.0*cxx*cyy - cxy*cxy;
    /* Don't try the sub-pixel refinement if bad cancellation occurs.  Zero D can means a line-like maximum that we
     * could still refine in the orthogonal direction but that seems a fringe case. */
    if (D == 0.0 || fabs(D) < 1e-8*MAX(fabs(4.0*cxx*cyy), fabs(cxy*cxy)))
        return FALSE;

    sx = (by*cxy - 2.0*bx*cyy)/D;
    sy = (bx*cxy - 2.0*by*cxx)/D;

    /* Don't trust the sub-pixel refinement if it moves the maximum too far from the centre. */
    if (sx*sx + sy*sy > 2.0)
        return FALSE;

    *x = sx;
    *y = sy;
    return TRUE;
}

/**
 * gwy_math_refine_maximum:
 * @z: Array of length 9, containing the square 3x3 neighbourhood values in matrix order and with the maximum in the
 *     centre.
 * @x: Location to store the refined @x-coordinate.
 * @y: Location to store the refined @y-coordinate.
 *
 * Performs subpixel refinement of parabolic a two-dimensional maximum.
 *
 * An alias for gwy_math_refine_maximum_2d().
 *
 * Returns: %TRUE if the refinement succeeded, %FALSE if it failed.  The values of @x and @y are usable regardless of
 *          the return value.
 *
 * Since: 2.42
 **/
gboolean
gwy_math_refine_maximum(const gdouble *z, gdouble *x, gdouble *y)
{
    return gwy_math_refine_maximum_2d(z, x, y);
}

static gboolean
interpolate_parabolic(const GwyXY *xy, gdouble *x)
{
    gdouble u1 = (xy[1].x - xy[0].x)*(xy[2].y - xy[1].y);
    gdouble u2 = (xy[2].x - xy[1].x)*(xy[0].y - xy[1].y);
    gdouble tx;

    if (fabs(u2 + u1) <= 1e-12*(fabs(u1) + fabs(u2)))
        return FALSE;

    tx = 0.5*(xy[1].x + (u2*xy[2].x + u1*xy[0].x)/(u2 + u1));
    if (tx <= xy[0].x || tx >= xy[2].x)
        return FALSE;

    *x = tx;
    return TRUE;
}

/**
 * gwy_compare_double:
 * @a: Pointer to a double.
 * @b: Pointer to a double.
 *
 * Compares two double values, given as pointers.
 *
 * This function is suitable as #GCompareFunc and can be also used with plain qsort(). The typical usage is sorting
 * of arrays containing structs where the first item is a floating point value/coordinate, followed by additional
 * data. For sorting of plain arrays of doubles use gwy_math_sort().
 *
 * It should only be used to sort normal numbers. The behaviour for NaNs is undefined.
 *
 * Since: 2.62
 **/
gint
gwy_compare_double(gconstpointer a, gconstpointer b)
{
    const gdouble *da = (const gdouble*)a;
    const gdouble *db = (const gdouble*)b;

    if (*da < *db)
        return -1;
    if (*da > *db)
        return 1;
    return 0;
}

static gboolean
find_min_in_array(const GwyXY *xy, guint n, guint *pimin)
{
    gboolean any_variation = FALSE;
    gdouble y, yy;
    guint imin, i;

    imin = n/2;
    for (i = 1; i <= n; i++) {
        if (xy[i].y < xy[imin].y)
            imin = i;
    }

    y = xy[imin].y;
    if (imin > 0) {
       yy = xy[imin-1].y;
       if (yy - y > 6e-16*(fabs(y) + fabs(yy)))
           any_variation = TRUE;
    }
    if (imin+1 < n) {
       yy = xy[imin+1].y;
       if (yy - y > 6e-16*(fabs(y) + fabs(yy)))
           any_variation = TRUE;
    }

    *pimin = imin;
    return any_variation;
}

/**
 * gwy_math_find_minimum_1d:
 * @function: Function to minimize.
 * @a: First interval endpoint.
 * @b: Second interval endpoint.
 * @user_data: User data passed to @function.
 *
 * Finds a minimum of a real function in a finite interval.
 *
 * The function simply does what it says on the tin.  If there are multiple minima in [a,b] any of them can be
 * returned, even though some effort to scan the interval is made.  There is no requiement for the minimum to lie
 * inside [a,b]; if it occurrs at one of the endpoints, the endpoint is returned.
 *
 * Since: 2.51
 **/
gdouble
gwy_math_find_minimum_1d(GwyRealFunc function,
                         gdouble a, gdouble b,
                         gpointer user_data)
{
    enum { initial_n = 12 };
    GwyXY xy[initial_n+1];
    gdouble x, y, xeps;
    guint i, imin, n, iter;
    gboolean at_left_edge = FALSE, at_right_edge = FALSE;

    GWY_ORDER(gdouble, a, b);
    if (b-a < 1.2e-16*(fabs(a) + fabs(b)))
        return 0.5*(a + b);

    /* Initial scan of the interval. */
    imin = 0;
    xy[0].x = a;
    xy[0].y = function(a, user_data);
    for (i = 1; i <= initial_n; i++) {
        x = (i == initial_n) ? b : b/initial_n*i + a/initial_n*(initial_n-i);
        xy[i].x = x;
        xy[i].y = y = function(x, user_data);
    }

    if (!find_min_in_array(xy, initial_n, &imin))
        return 0.5*(a + b);

    gwy_debug("initial minimum at %g, point #%u", xy[imin].x, imin);
    /* Use the first 4-5 values to keep points while iterating. */
    if (imin == 0)
        at_left_edge = TRUE;
    else if (imin == initial_n) {
        at_right_edge = TRUE;
        memmove(xy, xy + initial_n-2, 3*sizeof(GwyXY));
    }
    else
        memmove(xy, xy + imin-1, 3*sizeof(GwyXY));

    gwy_debug("initial subinterval [%g..%g]", xy[0].x, xy[2].x);

    iter = 0;
    while (xy[2].x - xy[0].x > 1.2e-15*(fabs(xy[0].x) + fabs(xy[2].x))) {
        gwy_debug("new iter, interval [%.16g..%.16g] %g, edges: %d %d",
                  xy[0].x, xy[2].x, xy[2].x - xy[0].x,
                  at_left_edge, at_right_edge);
        n = 3;
        /* Just split the interval closer to edge when the minimum seems at the edge. */
        if (at_left_edge)
            xy[n++].x = 0.8*xy[0].x + 0.2*xy[1].x;
        else if (at_right_edge)
            xy[n++].x = 0.2*xy[1].x + 0.8*xy[2].x;
        else {
            /* Optimistic bisection of the larger interval, always try this point.
             * XXX: This is not very efficient, we can end up only improving the interval from the bisection side. */
            if (xy[1].x - xy[0].x >= xy[2].x - xy[1].x) {
                xy[n++].x = 0.2*xy[0].x + 0.8*xy[1].x;
                gwy_debug("bisect-left %.16g (0.2)", xy[n-1].x);
            }
            else {
                xy[n++].x = 0.8*xy[1].x + 0.2*xy[2].x;
                gwy_debug("bisect-right %.16g (0.8)", xy[n-1].x);
            }

            /* Parabolic interpolation, use if it yields distinct point inside the interval. */
            if (interpolate_parabolic(xy, &x)) {
                xeps = 1.2e-15*(fabs(x));
                if (x - xy[0].x > xeps && xy[2].x - x > xeps && fabs(x - xy[3].x) > xeps) {
                    xy[n++].x = x;
                    gwy_debug("parabolic %.16g (%g)", x, (x - xy[0].x)/(xy[2].x - xy[0].x));
                }
            }
        }

        /* Find the new three points bracketing the minimum.  We should not change state from non-edge to edge, but we
         * can change state from edge to non-edge. */
        for (i = 3; i < n; i++) {
            xy[i].y = function(xy[i].x, user_data);
            gwy_debug("point %.16g, value %.16g", xy[i].x, xy[i].y);
        }
        qsort(xy, n, sizeof(GwyXY), gwy_compare_double);

        if (!find_min_in_array(xy, n, &imin))
            return xy[imin].x;

        gwy_debug("minimum at %g, point #%u (%g)", xy[imin].x, imin, (xy[imin].x - xy[0].x)/(xy[n-1].x - xy[0].x));
        if (imin == 0)
            at_left_edge = TRUE;
        else if (imin == n-1) {
            at_right_edge = TRUE;
            memmove(xy, xy + n-3, 3*sizeof(GwyXY));
        }
        else
            memmove(xy, xy + imin-1, 3*sizeof(GwyXY));

        if (iter++ == 50)
            break;
    }

    return xy[1].x;
}

static guint
estimate_regular_res(gdouble *pos, gint n, gdouble *minpos, gdouble *maxpos)
{
    gdouble maxstep = 0.0;
    guint k, res;

    gwy_math_sort(n, pos);
    *minpos = pos[0];
    *maxpos = pos[n-1];
    gwy_debug("range [%g,%g]", *minpos, *maxpos);
    if (maxpos <= minpos)
        return 0;

    for (k = 1; k < n; k++) {
        if (pos[k] - pos[k-1] > maxstep)
            maxstep = pos[k] - pos[k-1];
    }
    gwy_debug("maxstep %g", maxstep);
    res = (gint)ceil((*maxpos - *minpos)/maxstep) + 1;
    gwy_debug("estimated res %d", res);

    if (n % res != 0)
        return 0;

    return res;
}

/**
 * gwy_check_regular_2d_grid:
 * @coords: Array of @n coordinate pairs in plane.  You can also typecast #GwyXY or #GwyXYZ to doubles.
 * @stride: Actual number of double values in one block.  It must be at least 2 if @coords contains just alternating
 *          @x and @y.  If you pass an typecast #GwyXYZ array give stride as 3, etc.
 * @n: Number of items in @coords.
 * @tolerance: Relative distance from pixel center which is still considered OK.  Pass a negative value for some
 *             reasonable default. The maximum meaningful value is 0.5, beyond that the point would end up in
 *             a different pixel.
 * @xres: Location where to store the number of columns.
 * @yres: Location where to store the number of rows.
 * @xymin: Location where to store the minimum coordinates (top left corner).
 * @xystep: Location where to store the pixel size.
 *
 * Detects if points in plane form a regular rectangular grid oriented along the Cartesian axes.
 *
 * Points lying in one straight line are not considered to form a rectangle.
 *
 * When the function fails, i.e. the points do not form a regular grid, the values of output arguments are undefined.
 *
 * Returns: On success, a newly allocated array mapping grid indices (@i*@xres+@j) to indices in @coords.  %NULL is
 *          returned on failure.
 *
 * Since: 2.48
 **/
guint*
gwy_check_regular_2d_grid(const gdouble *coords, guint stride, guint n,
                          gdouble tolerance,
                          guint *pxres, guint *pyres,
                          GwyXY *xymin, GwyXY *xystep)
{
    gdouble xmin, xmax, ymin, ymax, dx, dy;
    gint xres, yres;
    guint k;
    gdouble *pos;
    guint *map;
    gboolean *encountered;
    gdouble matx[3], rhsx[2], maty[3], rhsy[3];

    g_return_val_if_fail(stride >= 2, NULL);
    g_return_val_if_fail(coords || !n, NULL);
    g_return_val_if_fail(pxres && pyres && xymin && xystep, NULL);

    if (n < 4)
        return NULL;

    if (tolerance < 0.0)
        tolerance = 0.05;

    pos = g_new(gdouble, n);
    gwy_debug("estimating yres from rows");
    for (k = 0; k < n; k++)
        pos[k] = coords[k*stride + 1];
    yres = estimate_regular_res(pos, n, &ymin, &ymax);

    gwy_debug("estimating xres from columns");
    for (k = 0; k < n; k++)
        pos[k] = coords[k*stride];
    xres = estimate_regular_res(pos, n, &xmin, &xmax);

    g_free(pos);

    if (yres) {
        xres = n/yres;
        gwy_debug("from rows xres %u, yres %u", xres, yres);
    }
    else if (xres) {
        yres = n/xres;
        gwy_debug("from columns xres %u, yres %u", xres, yres);
    }
    else
        return NULL;

    if (xres < 2 || yres < 2)
        return NULL;

    /* XXX: We could remove this condition but callers would need some means to tell if map[i] is an actual index or
     * something else, probably by putting G_MAXUINT there.   But this is an API change and current callers simply use
     * the values as indices.  A new function is needed. */
    if (xres*yres != n)
        return NULL;

    /* Widen the stripe by at most 1/2 at each side.  For large tolerance assume there is already a spread and widen
     * it less accordingly.  For exact coordinates this means differences from pixel centres are within the interval
     * [-tolerance/2, +tolerance/2], i.e. always safely smaller than tolerance in absolute value, regardless of the
     * tolerance.  So exact grids should always pass. */
    dx = (xmax - xmin)/(xres - 1 + tolerance);
    dy = (ymax - ymin)/(yres - 1 + tolerance);
    xmin -= 0.5*(1.0 - tolerance)*dx;
    xmax += 0.5*(1.0 - tolerance)*dx;
    ymin -= 0.5*(1.0 - tolerance)*dy;
    ymax += 0.5*(1.0 - tolerance)*dy;

    gwy_debug("x: [%g..%g] step %g", xmin, xmax, dx);
    gwy_debug("y: [%g..%g] step %g", ymin, ymax, dy);
    map = g_new(guint, n);
    encountered = g_new0(gboolean, n);
    gwy_clear(matx, 3);
    gwy_clear(maty, 3);
    gwy_clear(rhsx, 2);
    gwy_clear(rhsy, 2);
    for (k = 0; k < n; k++) {
        gdouble rawx = coords[k*stride + 0];
        gdouble rawy = coords[k*stride + 1];
        gdouble y = (rawy - ymin)/dy;
        gdouble x = (rawx - xmin)/dx;
        gint i = (gint)floor(y);
        gint j = (gint)floor(x);
        gdouble t;

        gwy_debug("(%g,%g) -> (%d,%d)", x, y, j, i);
        if (i < 0 || i >= yres || j < 0 || j >= xres) {
            g_critical("Points not inside estimated region?!");
            goto fail;
        }
        if (fabs(x - j - 0.5) > tolerance || fabs(y - i - 0.5) > tolerance) {
            gwy_debug("(%g,%g) too far from (%g,%g)", x, y, j+0.5, i+0.5);
            goto fail;
        }
        if (encountered[i*xres + j])
            goto fail;

        encountered[i*xres + j] = TRUE;
        map[i*xres + j] = k;

        t = j;
        matx[1] += t;
        matx[2] += t*t;
        rhsx[0] += rawx;
        rhsx[1] += t*rawx;

        t = i;
        maty[1] += t;
        maty[2] += t*t;
        rhsy[0] += rawy;
        rhsy[1] += t*rawy;
    }
    matx[0] = maty[0] = n;
    g_free(encountered);

    xymin->x = xmin;
    xymin->y = ymin;
    xystep->x = dx;
    xystep->y = dy;
    *pxres = xres;
    *pyres = yres;

    if (gwy_math_choleski_decompose(2, matx)) {
        gwy_math_choleski_solve(2, matx, rhsx);
        xystep->x = rhsx[1];
        xymin->x = rhsx[0] - 0.5*rhsx[1];
        gwy_debug("least-squares x-grid improvement to xoff=%g, xstep=%g", xymin->x, xystep->x);
    }
    if (gwy_math_choleski_decompose(2, maty)) {
        gwy_math_choleski_solve(2, maty, rhsy);
        xystep->y = rhsy[1];
        xymin->y = rhsy[0] - 0.5*rhsy[1];
        gwy_debug("least-squares y-grid improvement to yoff=%g, ystep=%g", xymin->y, xystep->y);
    }

    return map;

fail:
    g_free(map);
    g_free(encountered);
    return NULL;
}

/**
 * gwy_math_histogram:
 * @values: Values to make histogram from.
 * @n: Number of values in @values.
 * @min: Minimum value to consider (left edge of histogram).
 * @max: Maximum value to consider (right edge of histogram).
 * @nbins: Number of histogram bins (number of @counts items), a positive number.
 * @counts: Array where to store the counts.
 *
 * Counts the numbers of values falling into equal-sized bins.
 *
 * The value of @min must not be larger than @max.  The values may lie outside [@min,@max].  They are not counted in
 * the histogram, nor the returned total.
 *
 * Rounding rules for values exactly at the edge of two bins are arbitrary and must not be relied upon.
 *
 * Returns: The number of values inside the entire histogram, i.e. at most @n but possibly a reduced count.
 *
 * Since: 2.49
 **/
guint
gwy_math_histogram(const gdouble *values,
                   guint n,
                   gdouble min,
                   gdouble max,
                   guint nbins,
                   guint *counts)
{
    guint i, total = 0;
    gint nb;   /* Just a signed value */
    gdouble d;

    g_return_val_if_fail(nbins > 0, 0);
    g_return_val_if_fail(counts, 0);
    g_return_val_if_fail(values || !n, 0);
    g_return_val_if_fail(min <= max, 0);

    gwy_clear(counts, nbins);
    d = max - min;
    if (G_UNLIKELY(!(d > 0.0))) {
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:total) \
            private(i) \
            shared(values,n,min)
#endif
        for (i = 0; i < n; i++) {
            if (values[i] == min)
                total++;
        }
        counts[0] = total;
        return total;
    }

    gwy_clear(counts, nbins);
    d = nbins/d;
    nb = nbins;
#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            reduction(+:total) \
            private(i) \
            shared(counts,values,n,min,max,d,nb,nbins)
#endif
    {
        guint ifrom = gwy_omp_chunk_start(n), ito = gwy_omp_chunk_end(n);
        guint *tcounts = gwy_omp_if_threads_new0(counts, nbins);

        for (i = ifrom; i < ito; i++) {
            gdouble v = values[i];
            gint bi;

            if (v < min || v > max)
                continue;

            bi = (gint)floor((v - min)*d);
            if (G_LIKELY(bi >= 0 && bi < nb)) {
                tcounts[bi]++;
                total++;
            }
            else if (v == max) {
                tcounts[nbins-1]++;
                total++;
            }
        }
        gwy_omp_if_threads_sum_uint(counts, tcounts, nbins);
    }

    return total;
}

/**
 * gwy_xlnx_int:
 * @x: Value to calculate @x*log(@x) of.
 *
 * Calculates natural logarithm multiplied by the argument for integers.
 *
 * The value for zero @x is taken as the limit, i.e. zero.
 *
 * This function is useful for entropy calculations where values of @n*log(@n) can be evaulated a lot for small @n.
 * Therefore, values for small arguments are tabulated.  For large arguments the function is evaluated using the
 * standard log() function which is of course slower.
 *
 * Returns: Value of @x*log(@x).
 *
 * Since: 2.44
 **/
gdouble
gwy_xlnx_int(guint x)
{
    static const gdouble xlnx_table[] = {
        0.0,
        0.0,
        1.38629436111989061882,
        3.29583686600432907417,
        5.54517744447956247532,
        8.04718956217050187300,
        10.75055681536833000486,
        13.62137104338719313570,
        16.63553233343868742600,
        19.77502119602597444511,
        23.02585092994045684010,
        26.37684800078207598466,
        29.81887979745600372264,
        33.34434164699997756865,
        36.94680261461362060328,
        40.62075301653315098985,
        44.36141955583649980256,
        48.16462684895567336408,
        52.02669164213096445960,
        55.94434060416236874000,
        59.91464547107981986860,
        63.93497119219188292650,
        68.00293397388294877634,
        72.11636696637044288840,
        76.27329192835069487136,
        80.47189562170501873000,
    };

    /* Take the fast path quickly.  The slow path is slow anyway. */
    if (G_LIKELY(x < G_N_ELEMENTS(xlnx_table)))
        return xlnx_table[x];

    return x*log(x);
}

/**
 * gwy_sinc:
 * @x: Value to calculate sinc (cardinal sine) of.
 *
 * Calculates the sinc function.
 *
 * The sinc function is equal to sin(@x)/@x for non-zero @x, and defined to the limit 1 for zero @x.
 *
 * Returns: Value of sinc(@x).
 *
 * Since: 2.51
 **/
gdouble
gwy_sinc(gdouble x)
{
    if (G_LIKELY(fabs(x) > 3e-4))
        return sin(x)/x;
    return 1.0 - x*x/6.0;
}

/**
 * gwy_canonicalize_angle:
 * @phi: Angle to canonicalize, in radians.
 * @positive: %TRUE if a positive angle is requested, %FALSE for outputs symmetrical around zero.
 * @oriented: %TRUE for direction of a vector, %FALSE for the direction of a line (i.e. with no distinction between
 *            forward and backward direction).
 *
 * Canonicalizes an angle to requested interval.
 *
 * For @positive=%FALSE, @oriented=%FALSE the output interval is [-π/2,π/2].
 *
 * For @positive=%FALSE, @oriented=%TRUE the output interval is [-π,π].
 *
 * For @positive=%TRUE, @oriented=%FALSE the output interval is [0,π).
 *
 * For @positive=%TRUE, @oriented=%TRUE the output interval is [0,2π).
 *
 * Returns: Canonicalized angle, equivalent (in given sense) to @phi.
 *
 * Since: 2.50
 **/
gdouble
gwy_canonicalize_angle(gdouble phi, gboolean positive, gboolean oriented)
{
    if (oriented) {
        /* This can give anything from -2π up to 2π because fmod is based on rounding to zero. */
        phi = fmod(phi, GWY_TWO_PI);
        if (positive)
            return phi < 0.0 ? fmax(phi + GWY_TWO_PI, 0.0) : phi;

        if (phi < -G_PI)
            return phi + GWY_TWO_PI;
        if (phi > G_PI)
            return phi - GWY_TWO_PI;
        return phi;
    }

    /* This can give anything from -π up to π because fmod is based on rounding to zero. */
    phi = fmod(phi, G_PI);
    if (positive)
        return phi < 0.0 ? fmax(phi + G_PI, 0.0) : phi;

    if (phi < -G_PI_2)
        return phi + G_PI;
    if (phi > G_PI_2)
        return phi - G_PI;
    return phi;
}

/* Copyright (C) 1991, 1992, 1996, 1997, 1999 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Written by Douglas C. Schmidt (schmidt@ics.uci.edu).

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA  */

/* If you consider tuning this algorithm, you should consult first:
   Engineering a sort function; Jon Bentley and M. Douglas McIlroy;
   Software - Practice and Experience; Vol. 23 (11), 1249-1265, 1993.  */

/* The next 4 #defines implement a very fast in-line stack abstraction. */
#define PUSH(low, high) ((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define POP(low, high)  ((void) (--top, (low = top->lo), (high = top->hi)))
#define STACK_NOT_EMPTY (stack < top)

/* Order size using quicksort.  This implementation incorporates four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the next array partition to sort.  To save time,
   this maximum amount of space required to store an array of SIZE_MAX is allocated on the stack.  Assuming a 32-bit
   (64 bit) integer for size_t, this needs only 32 * sizeof(stack_node) == 256 bytes (for 64 bit: 1024 bytes). Pretty
   cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree. This reduces the probability of selecting a bad
   pivot value and eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving insertion sort to order the MAX_THRESH items within
   each partition. This is a big win, since insertion sort is faster for small, mostly sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the stack first, with the algorithm then
   concentrating on the smaller partition.  This *guarantees* no more than log(n) stack size is needed (actually O(1)
   in this case)!  */

/**
 * gwy_math_sort:
 * @n: Number of items in @array.
 * @array: Array of doubles to sort in place.
 *
 * Sorts an array of doubles using a quicksort algorithm.
 *
 * This is usually about twice as fast as the generic quicksort function thanks to specialization for doubles.
 **/
void
gwy_math_sort(gsize n, gdouble *array)
{
    /* Discontinue quicksort algorithm when partition gets below this size. This particular magic number was chosen to
     * work best on a Sun 4/260. */
    /* Specialization makes the insertion sort part relatively more efficient, after some benchmarking this seems be
     * about the best value on Athlon 64. */
    /* The stack needs log (total_elements) entries (we can even subtract log2(MAX_THRESH)).  Since total_elements has
     * type size_t, we get as upper bound for log (total_elements): bits per byte (CHAR_BIT) * sizeof(size_t).  */
    enum {
        MAX_THRESH = 12,
        LOG2_MAX_TRESH = 3,
        STACK_SIZE = CHAR_BIT*sizeof(gsize) - LOG2_MAX_TRESH,
    };

    /* Stack node declarations used to store unfulfilled partition obligations. */
    typedef struct {
        gdouble *lo;
        gdouble *hi;
    } stack_node;

    if (n < 2)
        /* Avoid lossage with unsigned arithmetic below.  */
        return;

    if (n > MAX_THRESH) {
        gdouble *lo = array;
        gdouble *hi = lo + (n - 1);
        stack_node stack[STACK_SIZE];
        stack_node *top = stack + 1;

        while (STACK_NOT_EMPTY) {
            gdouble *left_ptr;
            gdouble *right_ptr;

            /* Select median value from among LO, MID, and HI. Rearrange LO and HI so the three values are sorted.
             * This lowers the probability of picking a pathological pivot value and skips a comparison for both the
             * LEFT_PTR and RIGHT_PTR in the while loops. */

            gdouble *mid = lo + ((hi - lo) >> 1);

            if (*mid < *lo)
                DSWAP(*mid, *lo);

            if (*hi < *mid) {
                DSWAP(*mid, *hi);
                if (*mid < *lo)
                    DSWAP(*mid, *lo);
            }

            left_ptr  = lo + 1;
            right_ptr = hi - 1;

            /* Here's the famous ``collapse the walls'' section of quicksort. Gotta like those tight inner loops!
             * They are the main reason that this algorithm runs much faster than others. */
            do {
                while (*left_ptr < *mid)
                    left_ptr++;

                while (*mid < *right_ptr)
                    right_ptr--;

                if (left_ptr < right_ptr) {
                    DSWAP(*left_ptr, *right_ptr);
                    if (mid == left_ptr)
                        mid = right_ptr;
                    else if (mid == right_ptr)
                        mid = left_ptr;
                    left_ptr++;
                    right_ptr--;
                }
                else if (left_ptr == right_ptr) {
                    left_ptr++;
                    right_ptr--;
                    break;
                }
            }
            while (left_ptr <= right_ptr);

            /* Set up pointers for next iteration.  First determine whether left and right partitions are below the
             * threshold size.  If so, ignore one or both.  Otherwise, push the larger partition's bounds on the stack
             * and continue sorting the smaller one. */

            if ((gsize)(right_ptr - lo) <= MAX_THRESH) {
                if ((gsize)(hi - left_ptr) <= MAX_THRESH) {
                    /* Ignore both small partitions. */
                    POP(lo, hi);
                }
                else {
                    /* Ignore small left partition. */
                    lo = left_ptr;
                }
            }
            else if ((gsize)(hi - left_ptr) <= MAX_THRESH) {
                /* Ignore small right partition. */
                hi = right_ptr;
            }
            else if ((right_ptr - lo) > (hi - left_ptr)) {
                /* Push larger left partition indices. */
                PUSH(lo, right_ptr);
                lo = left_ptr;
            }
            else {
                /* Push larger right partition indices. */
                PUSH(left_ptr, hi);
                hi = right_ptr;
            }
        }
    }

    /* Once the BASE_PTR array is partially sorted by quicksort the rest is completely sorted using insertion sort,
     * since this is efficient for partitions below MAX_THRESH size. BASE_PTR points to the beginning of the array to
     * sort, and END_PTR points at the very last element in the array (*not* one beyond it!). */

    {
        gdouble *const end_ptr = array + (n - 1);
        gdouble *tmp_ptr = array;
        gdouble *thresh = MIN(end_ptr, array + MAX_THRESH);
        register gdouble *run_ptr;

        /* Find smallest element in first threshold and place it at the array's beginning.  This is the smallest array
         * element, and the operation speeds up insertion sort's inner loop. */
        for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++) {
            if (*run_ptr < *tmp_ptr)
                tmp_ptr = run_ptr;
        }

        if (tmp_ptr != array)
            DSWAP(*tmp_ptr, *array);

        /* Insertion sort, running from left-hand-side up to right-hand-side. */
        run_ptr = array + 1;
        while (++run_ptr <= end_ptr) {
            tmp_ptr = run_ptr - 1;
            while (*run_ptr < *tmp_ptr)
                tmp_ptr--;

            tmp_ptr++;
            if (tmp_ptr != run_ptr) {
                gdouble *hi, *lo;
                gdouble d;

                d = *run_ptr;
                for (hi = lo = run_ptr; --lo >= tmp_ptr; hi = lo)
                    *hi = *lo;
                *hi = d;
            }
        }
    }
}

/**
 * gwy_guint_sort:
 * @n: Number of items in @array.
 * @array: Array of #guint values to sort in place.
 *
 * Sorts an array of unsigned integers using a quicksort algorithm.
 *
 * This is usually about twice as fast as the generic quicksort function thanks to specialization for integers.
 *
 * Since: 2.50
 **/
void
gwy_guint_sort(gsize n, guint *array)
{
    /* Discontinue quicksort algorithm when partition gets below this size.
     * This particular magic number was chosen to work best on a Sun 4/260. */
    /* Specialization makes the insertion sort part relatively more
     * efficient, after some benchmarking this seems be about the best value
     * on Athlon 64. */
    /* The stack needs log (total_elements) entries (we can even subtract
     * log2(MAX_THRESH)).  Since total_elements has type size_t, we get as
     * upper bound for log (total_elements):
     * bits per byte (CHAR_BIT) * sizeof(size_t).  */
    enum {
        MAX_THRESH = 12,
        LOG2_MAX_TRESH = 3,
        STACK_SIZE = CHAR_BIT*sizeof(gsize) - LOG2_MAX_TRESH,
    };

    /* Stack node declarations used to store unfulfilled partition obligations.
     */
    typedef struct {
        guint *lo;
        guint *hi;
    } stack_node;

    if (n < 2)
        /* Avoid lossage with unsigned arithmetic below.  */
        return;

    if (n > MAX_THRESH) {
        guint *lo = array;
        guint *hi = lo + (n - 1);
        stack_node stack[STACK_SIZE];
        stack_node *top = stack + 1;

        while (STACK_NOT_EMPTY) {
            guint *left_ptr;
            guint *right_ptr;

            /* Select median value from among LO, MID, and HI. Rearrange
               LO and HI so the three values are sorted. This lowers the
               probability of picking a pathological pivot value and
               skips a comparison for both the LEFT_PTR and RIGHT_PTR in
               the while loops. */

            guint *mid = lo + ((hi - lo) >> 1);

            if (*mid < *lo)
                DSWAP(*mid, *lo);

            if (*hi < *mid) {
                DSWAP(*mid, *hi);
                if (*mid < *lo)
                    DSWAP(*mid, *lo);
            }

            left_ptr  = lo + 1;
            right_ptr = hi - 1;

            /* Here's the famous ``collapse the walls'' section of quicksort.
               Gotta like those tight inner loops!  They are the main reason
               that this algorithm runs much faster than others. */
            do {
                while (*left_ptr < *mid)
                    left_ptr++;

                while (*mid < *right_ptr)
                    right_ptr--;

                if (left_ptr < right_ptr) {
                    DSWAP(*left_ptr, *right_ptr);
                    if (mid == left_ptr)
                        mid = right_ptr;
                    else if (mid == right_ptr)
                        mid = left_ptr;
                    left_ptr++;
                    right_ptr--;
                }
                else if (left_ptr == right_ptr) {
                    left_ptr++;
                    right_ptr--;
                    break;
                }
            }
            while (left_ptr <= right_ptr);

          /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size.  If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

          if ((gsize)(right_ptr - lo) <= MAX_THRESH) {
              if ((gsize)(hi - left_ptr) <= MAX_THRESH)
                  /* Ignore both small partitions. */
                  POP(lo, hi);
              else
                  /* Ignore small left partition. */
                  lo = left_ptr;
          }
          else if ((gsize)(hi - left_ptr) <= MAX_THRESH)
              /* Ignore small right partition. */
              hi = right_ptr;
          else if ((right_ptr - lo) > (hi - left_ptr)) {
              /* Push larger left partition indices. */
              PUSH(lo, right_ptr);
              lo = left_ptr;
          }
          else {
              /* Push larger right partition indices. */
              PUSH(left_ptr, hi);
              hi = right_ptr;
          }
        }
    }

    /* Once the BASE_PTR array is partially sorted by quicksort the rest
       is completely sorted using insertion sort, since this is efficient
       for partitions below MAX_THRESH size. BASE_PTR points to the beginning
       of the array to sort, and END_PTR points at the very last element in
       the array (*not* one beyond it!). */

    {
        guint *const end_ptr = array + (n - 1);
        guint *tmp_ptr = array;
        guint *thresh = MIN(end_ptr, array + MAX_THRESH);
        register guint *run_ptr;

        /* Find smallest element in first threshold and place it at the
           array's beginning.  This is the smallest array element,
           and the operation speeds up insertion sort's inner loop. */

        for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++) {
            if (*run_ptr < *tmp_ptr)
                tmp_ptr = run_ptr;
        }

        if (tmp_ptr != array)
            DSWAP(*tmp_ptr, *array);

        /* Insertion sort, running from left-hand-side up to right-hand-side.
         */

        run_ptr = array + 1;
        while (++run_ptr <= end_ptr) {
            tmp_ptr = run_ptr - 1;
            while (*run_ptr < *tmp_ptr)
                tmp_ptr--;

            tmp_ptr++;
            if (tmp_ptr != run_ptr) {
                guint *hi, *lo;
                guint d;

                d = *run_ptr;
                for (hi = lo = run_ptr; --lo >= tmp_ptr; hi = lo)
                    *hi = *lo;
                *hi = d;
            }
        }
    }
}

/**
 * gwy_math_sort_with_index:
 * @n: Number of items in @array.
 * @array: Array of doubles to sort in place.
 * @index_array: Array of integer identifiers of the items that are permuted simultaneously with @array.
 *
 * Sorts an array of doubles using a quicksort algorithm, remembering the permutation.
 *
 * The simplest and probably most common use of @index_array is to fill it with numbers 0 to @n-1 before calling
 * gwy_math_sort().  After sorting, @index_array[@i] then contains the original position of the @i-th item of the
 * sorted array.
 *
 * Since: 2.50
 **/
/* FIXME: It is questionable whether it is still more efficient to use pointers instead of array indices when it
 * effectively doubles the number of variables.  This might force some variables from registers to memory... */
void
gwy_math_sort_with_index(gsize n, gdouble *array, guint *index_array)
{
    enum {
        MAX_THRESH = 12,
        LOG2_MAX_TRESH = 3,
        STACK_SIZE = CHAR_BIT*sizeof(gsize) - LOG2_MAX_TRESH,
    };

    /* Stack node declarations used to store unfulfilled partition obligations. */
    typedef struct {
        gdouble *lo;
        gdouble *hi;
        guint *loi;
        guint *hii;
    } stack_node;

    if (n < 2)
        /* Avoid lossage with unsigned arithmetic below.  */
        return;

    if (n > MAX_THRESH) {
        gdouble *lo = array;
        gdouble *hi = lo + (n - 1);
        guint *loi = index_array;
        guint *hii = loi + (n - 1);
        stack_node stack[STACK_SIZE];
        stack_node *top = stack + 1;

        while (STACK_NOT_EMPTY) {
            gdouble *left_ptr;
            gdouble *right_ptr;
            guint *left_ptri;
            guint *right_ptri;

            /* Select median value from among LO, MID, and HI. Rearrange
               LO and HI so the three values are sorted. This lowers the
               probability of picking a pathological pivot value and
               skips a comparison for both the LEFT_PTR and RIGHT_PTR in
               the while loops. */

            gdouble *mid = lo + ((hi - lo) >> 1);
            guint *midi = loi + ((hii - loi) >> 1);

            if (*mid < *lo) {
                DSWAP(*mid, *lo);
                ISWAP(*midi, *loi);
            }
            if (*hi < *mid) {
                DSWAP(*mid, *hi);
                ISWAP(*midi, *hii);

                if (*mid < *lo) {
                    DSWAP(*mid, *lo);
                    ISWAP(*midi, *loi);
                }
            }

          left_ptr  = lo + 1;
          right_ptr = hi - 1;
          left_ptri  = loi + 1;
          right_ptri = hii - 1;

          /* Here's the famous ``collapse the walls'' section of quicksort.
             Gotta like those tight inner loops!  They are the main reason
             that this algorithm runs much faster than others. */
          do {
              while (*left_ptr < *mid) {
                  left_ptr++;
                  left_ptri++;
              }

              while (*mid < *right_ptr) {
                  right_ptr--;
                  right_ptri--;
              }

              if (left_ptr < right_ptr) {
                  DSWAP(*left_ptr, *right_ptr);
                  ISWAP(*left_ptri, *right_ptri);
                  if (mid == left_ptr) {
                      mid = right_ptr;
                      midi = right_ptri;
                  }
                  else if (mid == right_ptr) {
                      mid = left_ptr;
                      midi = left_ptri;
                  }
                  left_ptr++;
                  left_ptri++;
                  right_ptr--;
                  right_ptri--;
              }
              else if (left_ptr == right_ptr) {
                  left_ptr++;
                  left_ptri++;
                  right_ptr--;
                  right_ptri--;
                  break;
              }
          }
          while (left_ptr <= right_ptr);

          /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size.  If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

          if ((gsize)(right_ptr - lo) <= MAX_THRESH) {
              if ((gsize)(hi - left_ptr) <= MAX_THRESH) {
                  /* Ignore both small partitions. */
                  --top;
                  lo = top->lo;
                  hi = top->hi;
                  loi = top->loi;
                  hii = top->hii;
              }
              else {
                  /* Ignore small left partition. */
                  lo = left_ptr;
                  loi = left_ptri;
              }
          }
          else if ((gsize)(hi - left_ptr) <= MAX_THRESH) {
              /* Ignore small right partition. */
              hi = right_ptr;
              hii = right_ptri;
          }
          else if ((right_ptr - lo) > (hi - left_ptr)) {
              /* Push larger left partition indices. */
              top->lo = lo;
              top->loi = loi;
              top->hi = right_ptr;
              top->hii = right_ptri;
              ++top;
              lo = left_ptr;
              loi = left_ptri;
          }
          else {
              /* Push larger right partition indices. */
              top->lo = left_ptr;
              top->loi = left_ptri;
              top->hi = hi;
              top->hii = hii;
              ++top;
              hi = right_ptr;
              hii = right_ptri;
          }
        }
    }

    /* Once the BASE_PTR array is partially sorted by quicksort the rest
       is completely sorted using insertion sort, since this is efficient
       for partitions below MAX_THRESH size. BASE_PTR points to the beginning
       of the array to sort, and END_PTR points at the very last element in
       the array (*not* one beyond it!). */

    {
        gdouble *const end_ptr = array + (n - 1);
        gdouble *tmp_ptr = array;
        guint *tmp_ptri = index_array;
        gdouble *thresh = MIN(end_ptr, array + MAX_THRESH);
        gdouble *run_ptr;
        guint *run_ptri;

        /* Find smallest element in first threshold and place it at the
           array's beginning.  This is the smallest array element,
           and the operation speeds up insertion sort's inner loop. */

        for (run_ptr = tmp_ptr + 1, run_ptri = tmp_ptri + 1;
             run_ptr <= thresh;
             run_ptr++, run_ptri++) {
            if (*run_ptr < *tmp_ptr) {
                tmp_ptr = run_ptr;
                tmp_ptri = run_ptri;
            }
        }

        if (tmp_ptr != array) {
            DSWAP(*tmp_ptr, *array);
            ISWAP(*tmp_ptri, *index_array);
        }

        /* Insertion sort, running from left-hand-side up to right-hand-side.
         */

        run_ptr = array + 1;
        run_ptri = index_array + 1;
        while (++run_ptr <= end_ptr) {
            tmp_ptr = run_ptr - 1;
            tmp_ptri = run_ptri;
            ++run_ptri;
            while (*run_ptr < *tmp_ptr) {
                tmp_ptr--;
                tmp_ptri--;
            }

            tmp_ptr++;
            tmp_ptri++;
            if (tmp_ptr != run_ptr) {
                gdouble *hi, *lo;
                guint *hii, *loi;
                gdouble d;
                guint i;

                d = *run_ptr;
                for (hi = lo = run_ptr; --lo >= tmp_ptr; hi = lo)
                    *hi = *lo;
                *hi = d;

                i = *run_ptri;
                for (hii = loi = run_ptri; --loi >= tmp_ptri; hii = loi)
                    *hii = *loi;
                *hii = i;
            }
        }
    }
}

/**
 * gwy_math_median_uncertainty:
 * @n: Number of items in @array.
 * @array: Array of doubles.  It is modified by this function.  All values are kept, but their positions in the array
 *         change.
 * @uarray: Array of value unvertainries.  It is modified by this function. All values are kept, but their positions
 *          in the array change.
 *
 * Find the uncertainty value corresponding to data median.
 *
 * Note that this is not the uncertainty arising from the calculation of the median.  It is just the uncertainty of
 * the single value that happens to be the data median.  As such, the function is not very useful.
 *
 * Since: 2.23
 *
 * Returns: The uncertainty of the median value.
 **/
gdouble
gwy_math_median_uncertainty(gsize n, gdouble *array, gdouble *uarray)
{
    gsize lo, hi;
    gsize median;
    gsize middle, ll, hh;

    lo = 0;
    hi = n - 1;
    median = n/2;
    while (TRUE) {
        if (hi <= lo)        /* One element only */
            return uarray[median];

        if (hi == lo + 1) {  /* Two elements only */
            if (array[lo] > array[hi]){
                DSWAP(array[lo], array[hi]);
                DSWAP(uarray[lo], uarray[hi]);
            }
            return uarray[median];
        }

        /* Find median of lo, middle and hi items; swap into position lo */
        middle = (lo + hi)/2;
        if (array[middle] > array[hi]){
            DSWAP(array[middle], array[hi]);
            DSWAP(uarray[middle], uarray[hi]);
        }
        if (array[lo] > array[hi]){
            DSWAP(array[lo], array[hi]);
            DSWAP(uarray[lo], uarray[hi]);
        }
        if (array[middle] > array[lo]){
            DSWAP(array[middle], array[lo]);
            DSWAP(uarray[middle], uarray[lo]);
        }

        /* Swap low item (now in position middle) into position (lo+1) */
        DSWAP(array[middle], array[lo + 1]);
        DSWAP(uarray[middle], uarray[lo + 1]);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = lo + 1;
        hh = hi;
        while (TRUE) {
            do {
                ll++;
            } while (array[lo] > array[ll]);
            do {
                hh--;
            } while (array[hh] > array[lo]);

            if (hh < ll)
                break;

            DSWAP(array[ll], array[hh]);
            DSWAP(uarray[ll], uarray[hh]);
        }
        /* Swap middle item (in position lo) back into correct position */
        DSWAP(array[lo], array[hh]);
        DSWAP(uarray[lo], uarray[hh]);

        /* Re-set active partition */
        if (hh <= median)
            lo = ll;
        if (hh >= median)
            hi = hh - 1;
    }
}

/************************** Documentation ****************************/

/**
 * SECTION:gwymath
 * @title: Math
 * @short_description: Mathematical utility functions
 * @see_also: #GwyNLFitter, non-linear least square fitter;
 *            <link linkend="libgwyddion-Math-Fallback">Math Fallback</link>,
 *            fallback mathematical functions
 *
 * Function gwy_math_humanize_numbers() deals with number representation.
 *
 * Nearest object finding functions gwy_math_find_nearest_line() and gwy_math_find_nearest_point() can be useful in
 * widget and vector layer implementation.
 *
 * And gwy_math_lin_solve(), gwy_math_lin_solve_rewrite(), and gwy_math_fit_polynom() are general purpose numeric
 * methods.
 **/

/**
 * ROUND:
 * @x: A double value.
 *
 * Rounds a number to nearest integer.  Use %GWY_ROUND instead.
 **/

/**
 * GWY_ROUND:
 * @x: A double value.
 *
 * Rounds a number to nearest integer.
 *
 * Since: 2.5
 **/

/**
 * GWY_SQRT3:
 *
 * The square root of 3.
 **/

/**
 * GWY_SQRT_PI:
 *
 * The square root of pi.
 **/

/**
 * GwyXY:
 * @x: X-coordinate.
 * @y: Y-coordinate.
 *
 * Representation of Cartesian coordinates in plane.
 *
 * Since: 2.45
 **/

/**
 * GwyXYZ:
 * @x: X-coordinate.
 * @y: Y-coordinate.
 * @z: Z-coordinate.
 *
 * Representation of Cartesian coordinates in space.
 *
 * Since: 2.45
 **/

/**
 * SECTION:gwymathfallback
 * @title: Math Fallback
 * @short_description: Fallback implementations of standard mathematical functions
 * @include: libgwyddion/gwymathfallback.h
 *
 * Fallback functions <function>gwy_math_fallback_<replaceable>foo</replaceable></function> are defined for
 * mathematical functions <function><replaceable>foo</replaceable></function> that might not be implemented on all
 * platforms and are commonly used in Gwyddion.  These functions are always defined (as <literal>static
 * inline</literal>), however, you should not use them as they can be less efficient or precise than the standard
 * functions.
 *
 * For each unavailable function (and only for those), this header file defines a replacement macro expanding to the
 * name of the fallback function. Therefore after including it, you can use for instance <function>cbrt</function>
 * regardless if the platform provides it or not. Note this header has to be included explicitly to avoid possible
 * inadvertent clashes with other definitions of <function>cbrt</function>.
 *
 * Since all replacement macros expand to names of functions, it is possible to take the address of any of them.
 **/

/**
 * gwy_math_fallback_cbrt:
 * @x: Floating point number.
 *
 * Fallback for the standard mathematical function <function>cbrt</function>.
 *
 * Returns: Cubic root of @x.
 *
 * Since: 2.9
 **/

/**
 * cbrt:
 *
 * Macro defined to gwy_math_fallback_cbrt() if the platform does not provide <function>cbrt</function>.
 **/

/**
 * gwy_math_fallback_pow10:
 * @x: Floating point number.
 *
 * Fallback for the standard mathematical function <function>pow10</function>.
 *
 * Returns: 10 raised to @x.
 *
 * Since: 2.9
 **/

/**
 * pow10:
 *
 * Macro defined to gwy_math_fallback_pow10() if the platform does not provide <function>pow10</function>.
 **/

/**
 * gwy_math_fallback_hypot:
 * @x: Floating point number.
 * @y: Floating point number.
 *
 * Fallback for the standard mathematical function <function>hypot</function>.
 *
 * Returns: Length of hypotenuse of a right-angle triangle with sides of lengths @x and @y.
 *
 * Since: 2.9
 **/

/**
 * hypot:
 *
 * Macro defined to gwy_math_fallback_hypot() if the platform does not provide <function>hypot</function>.
 **/

/**
 * gwy_math_fallback_acosh:
 * @x: Floating point number greater or equal to 1.0.
 *
 * Fallback for the standard mathematical function <function>acosh</function>.
 *
 * Returns: Inverse hyperbolic cosine of @x.
 *
 * Since: 2.9
 **/

/**
 * acosh:
 *
 * Macro defined to gwy_math_fallback_acosh() if the platform does not provide <function>acosh</function>.
 **/

/**
 * gwy_math_fallback_asinh:
 * @x: Floating point number.
 *
 * Fallback for the standard mathematical function <function>asinh</function>.
 *
 * Returns: Inverse hyperbolic sine of @x.
 *
 * Since: 2.9
 **/

/**
 * asinh:
 *
 * Macro defined to gwy_math_fallback_asinh() if the platform does not provide <function>asinh</function>.
 **/

/**
 * gwy_math_fallback_atanh:
 * @x: Floating point number in the range [-1, 1].
 *
 * Fallback for the standard mathematical function <function>atanh</function>.
 *
 * Returns: Inverse hyperbolic tangent of @x.
 *
 * Since: 2.9
 **/

/**
 * atanh:
 *
 * Macro defined to gwy_math_fallback_atanh() if the platform does not provide <function>atanh</function>.
 **/

/**
 * gwy_math_fallback_isnan:
 * @x: Floating point number.
 *
 * Fallback for the standard mathematical function <function>isnan</function>.
 *
 * Returns: %TRUE if @x is infinity, %FALSE otherwise.
 *
 * Since: 2.22
 **/

/**
 * gwy_isnan:
 *
 * Macro defined to working isnan() implementation, either a system one or gwy_math_fallback_isnan().
 *
 * Since: 2.22
 **/

/**
 * gwy_math_fallback_isinf:
 * @x: Floating point number.
 *
 * Fallback for the standard mathematical function <function>isinf</function>.
 *
 * Returns: %TRUE if @x is infinity, %FALSE otherwise.
 *
 * Since: 2.22
 **/

/**
 * gwy_isinf:
 *
 * Macro defined to working isinf() implementation, either a system one or gwy_math_fallback_isinf().
 *
 * Since: 2.22
 **/

/**
 * gwy_math_fallback_powi:
 * @x: Floating point number.
 * @i: Integer power.
 *
 * Fallback for the integer power function.
 *
 * It provides the same functionality as GCC's __builtin_powi() for finite fast math, without any precision guarantee.
 *
 * Returns: Value of @x raised to @i-th power.  If @i is zero, the return values is 1, even when @x is zero.
 *
 * Since: 2.53
 **/

/**
 * gwy_powi:
 *
 * Macro defined to working integer power implementation, either a compiler provided one or gwy_math_fallback_powi().
 *
 * Since: 2.53
 **/

/**
 * GwyPercentileInterpolationType:
 * @GWY_PERCENTILE_INTERPOLATION_LINEAR: Linear interpolation of the two nearest values.
 * @GWY_PERCENTILE_INTERPOLATION_LOWER: Round the rank down to an integer.
 * @GWY_PERCENTILE_INTERPOLATION_HIGHER: Round the rank up to an integer.
 * @GWY_PERCENTILE_INTERPOLATION_NEAREST: Round the rank to nearest integer.
 * @GWY_PERCENTILE_INTERPOLATION_MIDPOINT: Average of the two nearest values.
 *
 * Type of interpolation for percentile calculation.
 *
 * The interpolations are used when the percentile does not correspond exactly to a rank.
 *
 * Since: 2.50
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
