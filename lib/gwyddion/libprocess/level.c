/*
 *  $Id: level.c 23189 2021-02-23 17:32:52Z yeti-dn $
 *  Copyright (C) 2003-2020 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
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

#include <string.h>
#include <libgwyddion/gwymacros.h>
#include <libprocess/datafield.h>
#include <libprocess/level.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"

/**
 * gwy_data_field_fit_plane:
 * @data_field: A data field.
 * @pa: Where constant coefficient should be stored (or %NULL).
 * @pbx: Where x plane coefficient should be stored (or %NULL).
 * @pby: Where y plane coefficient should be stored (or %NULL).
 *
 * Fits a plane through a data field.
 *
 * The coefficients can be used for plane leveling using relation
 * data[i] := data[i] - (pa + pby*i + pbx*j);
 **/
void
gwy_data_field_fit_plane(GwyDataField *data_field,
                         gdouble *pa, gdouble *pbx, gdouble *pby)
{
    gdouble sumxi, sumxixi, sumyi, sumyiyi;
    gdouble sumsi, sumsixi, sumsiyi;
    gdouble bx, by;
    const gdouble *d, *drow;
    gint i, j, xres, yres;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    xres = data_field->xres;
    yres = data_field->yres;
    d = data_field->data;

    sumxi = (xres - 1.0)/2.0;
    sumxixi = (2.0*xres - 1.0)*(xres - 1.0)/6.0;
    sumyi = (yres - 1.0)/2.0;
    sumyiyi = (2.0*yres - 1.0)*(yres - 1.0)/6.0;
    sumsi = sumsixi = sumsiyi = 0.0;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
        reduction(+:sumsi,sumsixi,sumsiyi) \
        private(drow,i,j) \
        shared(d,xres,yres)
#endif
    for (i = 0; i < yres; i++) {
        drow = d + i*xres;
        for (j = 0; j < xres; j++) {
            sumsi += drow[j];
            sumsixi += drow[j] * j;
            sumsiyi += drow[j] * i;
        }
    }
    sumsi /= xres*yres;
    sumsixi /= xres*yres;
    sumsiyi /= xres*yres;

    bx = (sumsixi - sumsi*sumxi)/(sumxixi - sumxi*sumxi);
    by = (sumsiyi - sumsi*sumyi)/(sumyiyi - sumyi*sumyi);
    if (pbx)
        *pbx = bx;
    if (pby)
        *pby = by;
    if (pa)
        *pa = sumsi - bx*sumxi - by*sumyi;
}

/**
 * gwy_data_field_area_fit_plane:
 * @data_field: A data field
 * @mask: Mask of values to take values into account, or %NULL for full
 *        @data_field.  Values equal to 0.0 and below cause corresponding
 *        @data_field samples to be ignored, values equal to 1.0 and above
 *        cause inclusion of corresponding @data_field samples.  The behaviour
 *        for values inside (0.0, 1.0) is undefined (it may be specified
 *        in the future).
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @pa: Where constant coefficient should be stored (or %NULL).
 * @pbx: Where x plane coefficient should be stored (or %NULL).
 * @pby: Where y plane coefficient should be stored (or %NULL).
 *
 * Fits a plane through a rectangular part of a data field.
 *
 * The coefficients can be used for plane leveling using the same relation
 * as in gwy_data_field_fit_plane(), counting indices from area top left
 * corner.
 **/
void
gwy_data_field_area_fit_plane(GwyDataField *data_field,
                              GwyDataField *mask,
                              gint col, gint row, gint width, gint height,
                              gdouble *pa, gdouble *pbx, gdouble *pby)
{
    gwy_data_field_area_fit_plane_mask(data_field, mask, GWY_MASK_INCLUDE,
                                       col, row, width, height,
                                       pa, pbx, pby);
}

/**
 * gwy_data_field_area_fit_plane_mask:
 * @data_field: A data field
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @pa: Where constant coefficient should be stored (or %NULL).
 * @pbx: Where x plane coefficient should be stored (or %NULL).
 * @pby: Where y plane coefficient should be stored (or %NULL).
 *
 * Fits a plane through a rectangular part of a data field with masking.
 *
 * The coefficients can be used for plane leveling using the same relation
 * as in gwy_data_field_fit_plane(), counting indices from area top left
 * corner.
 *
 * Since: 2.56
 **/
void
gwy_data_field_area_fit_plane_mask(GwyDataField *data_field,
                                   GwyDataField *mask,
                                   GwyMaskingType masking,
                                   gint col, gint row,
                                   gint width, gint height,
                                   gdouble *pa, gdouble *pbx, gdouble *pby)
{
    const gdouble *datapos = NULL, *maskpos = NULL;
    gdouble sumx, sumy, sumz, sumxx, sumyy, sumxy, sumxz, sumyz;
    gdouble alpha1, alpha2, alpha3, beta1, beta2, gamma1;
    gdouble det;
    gdouble a, b, c;
    gdouble n;
    gint i, xres;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height)
        || !_gwy_data_field_check_mask(data_field, &mask, &masking))
        return;
    xres = data_field->xres;
    n = 0;

    /* try to return something reasonable even in degenerate cases */
    if (!width || !height)
        a = b = c = 0.0;
    else if (height == 1 && width == 1) {
        c = data_field->data[row*xres + col];
        a = b = 0.0;
    }
    else {
        sumx = sumy = sumz = sumxx = sumyy = sumxy = sumxz = sumyz = 0;
        datapos = data_field->data + row*xres + col;
        if (mask)
            maskpos = mask->data + row*xres + col;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
        reduction(+:n,sumx,sumy,sumz,sumxx,sumyy,sumxy,sumxz,sumyz) \
        private(i) \
        shared(datapos,maskpos,xres,width,height,masking)
#endif
        for (i = 0; i < height; i++) {
            const gdouble *drow = datapos + i*xres;
            guint j;

            if (masking == GWY_MASK_INCLUDE) {
                const gdouble *mrow = maskpos + i*xres;
                for (j = 0; j < width; j++) {
                    guint cc = (*mrow > 0.0);
                    gdouble x = cc*j;
                    gdouble y = cc*i;
                    gdouble z = cc*drow[j];

                    sumx += x;
                    sumy += y;
                    sumz += z;
                    sumxx += x*x;
                    sumyy += y*y;
                    sumxy += x*y;
                    sumxz += x*z;
                    sumyz += y*z;
                    n += cc;
                    mrow++;
                }
            }
            else if (masking == GWY_MASK_EXCLUDE) {
                const gdouble *mrow = maskpos + i*xres;
                for (j = 0; j < width; j++) {
                    guint cc = (*mrow <= 0.0);
                    gdouble x = cc*j;
                    gdouble y = cc*i;
                    gdouble z = cc*drow[j];

                    sumx += x;
                    sumy += y;
                    sumz += z;
                    sumxx += x*x;
                    sumyy += y*y;
                    sumxy += x*y;
                    sumxz += x*z;
                    sumyz += y*z;
                    n += cc;
                    mrow++;
                }
            }
            else {
                for (j = 0; j < width; j++) {
                    gdouble x = j;
                    gdouble y = i;
                    gdouble z = drow[j];

                    sumx += x;
                    sumy += y;
                    sumz += z;
                    sumxx += x*x;
                    sumyy += y*y;
                    sumxy += x*y;
                    sumxz += x*z;
                    sumyz += y*z;
                }
                n = width;
            }
        }

        det = (n*sumxx*sumyy) + (2*sumx*sumxy*sumy) - (sumx*sumx*sumyy)
                -(sumy*sumy*sumxx) - (n*sumxy*sumxy);

        /* try to return something reasonable in case of singularity */
        if (det == 0.0)
            a = b = c = 0.0;
        else {
            det = 1.0/det;

            alpha1 = (n*sumyy) - (sumy*sumy);
            alpha2 = (n*sumxx) - (sumx*sumx);
            alpha3 = (sumxx*sumyy) - (sumxy*sumxy);
            beta1 = (sumx*sumy) - (n*sumxy);
            beta2 = (sumx*sumxy) - (sumxx*sumy);
            gamma1 = (sumxy*sumy) - (sumx*sumyy);

            a = det*(alpha1*sumxz + beta1*sumyz + gamma1*sumz);
            b = det*(beta1*sumxz + alpha2*sumyz + beta2*sumz);
            c = det*(gamma1*sumxz + beta2*sumyz + alpha3*sumz);
        }
    }

    if (pbx)
        *pbx = a;
    if (pby)
        *pby = b;
    if (pa)
        *pa = c;
}

/**
 * gwy_data_field_fit_facet_plane:
 * @data_field: A data field.
 * @mfield: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use.
 * @pa: Where constant coefficient should be stored (or %NULL).
 * @pbx: Where x plane coefficient should be stored.
 * @pby: Where y plane coefficient should be stored.
 *
 * Calculates the inclination of a plane close to the dominant plane in a data
 * field.
 *
 * The dominant plane is determined by taking into account larger local slopes
 * with exponentially smaller weight.
 *
 * This is the basis of so-called facet levelling algorithm.  Usually, the
 * plane found by this method is subtracted using gwy_data_field_plane_level()
 * and the entire process is repeated until it converges.  A convergence
 * criterion may be sufficiently small values of the x and y plane
 * coefficients.  Note that since gwy_data_field_plane_level() uses pixel-based
 * lateral coordinates, the coefficients must be divided by
 * gwy_data_field_get_dx(data_field) and
 * gwy_data_field_get_dy(data_field) to obtain physical plane
 * coefficients.
 *
 * Returns: %TRUE if any plane was actually fitted; %FALSE if there was an
 *          insufficient number of unmasked pixels.
 *
 * Since: 2.37
 **/
gboolean
gwy_data_field_fit_facet_plane(GwyDataField *data_field,
                               GwyDataField *mask,
                               GwyMaskingType masking,
                               gdouble *pa, gdouble *pbx, gdouble *pby)
{
    const gdouble c = 1.0/20.0;

    gdouble *data, *row, *newrow;
    const gdouble *mdata, *mrow, *newmrow;
    gdouble vx, vy, q, sumvx, sumvy, sumvz, xr, yr, sigma2;
    gint xres, yres, n, i, j;

    if (!_gwy_data_field_check_mask(data_field, &mask, &masking))
        return FALSE;

    *pbx = *pby = 0.0;
    if (pa)
        *pa = 0.0;

    xres = data_field->xres;
    yres = data_field->yres;
    xr = data_field->xreal/xres;
    yr = data_field->yreal/yres;

    data = data_field->data;
    mdata = mask ? mask->data : NULL;

    sigma2 = 0.0;
    n = 0;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
        reduction(+:sigma2,n) \
        private(i,j,row,newrow,mrow,newmrow,vx,vy) \
        shared(data,mdata,xres,yres,xr,yr,masking)
#endif
    for (i = 1; i < yres; i++) {
        row = data + (i-1)*xres;
        newrow = data + i*xres;
        mrow = mdata + (i-1)*xres;
        newmrow = mdata + i*xres;
        for (j = 1; j < xres; j++) {
            if (masking == GWY_MASK_IGNORE
                || (masking == GWY_MASK_INCLUDE
                    && newmrow[j] >= 1.0 && mrow[j] >= 1.0
                    && newmrow[j-1] >= 1.0 && mrow[j-1] >= 1.0)
                || (masking == GWY_MASK_EXCLUDE
                    && newmrow[j] <= 0.0 && mrow[j] <= 0.0
                    && newmrow[j-1] <= 0.0 && mrow[j-1] <= 0.0)) {
                n++;
                vx = 0.5*(newrow[j] + row[j] - newrow[j-1] - row[j-1])/xr;
                vy = 0.5*(newrow[j-1] + newrow[j] - row[j-1] - row[j])/yr;
                sigma2 += vx*vx + vy*vy;
            }
        }
    }
    /* Do not try to level from some random pixel */
    gwy_debug("n=%d", n);
    if (n < 4)
        return FALSE;

    sigma2 = c*sigma2/n;

    sumvx = sumvy = sumvz = 0.0;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
        reduction(+:sumvx,sumvy,sumvz) \
        private(i,j,row,newrow,mrow,newmrow,vx,vy,q) \
        shared(data,mdata,xres,yres,xr,yr,masking,sigma2)
#endif
    for (i = 1; i < yres; i++) {
        row = data + (i-1)*xres;
        newrow = data + i*xres;
        mrow = mdata + (i-1)*xres;
        newmrow = mdata + i*xres;
        for (j = 1; j < xres; j++) {
            if (masking == GWY_MASK_IGNORE
                || (masking == GWY_MASK_INCLUDE
                    && newmrow[j] >= 1.0 && mrow[j] >= 1.0
                    && newmrow[j-1] >= 1.0 && mrow[j-1] >= 1.0)
                || (masking == GWY_MASK_EXCLUDE
                    && newmrow[j] <= 0.0 && mrow[j] <= 0.0
                    && newmrow[j-1] <= 0.0 && mrow[j-1] <= 0.0)) {
                vx = 0.5*(newrow[j] + row[j] - newrow[j-1] - row[j-1])/xr;
                vy = 0.5*(newrow[j-1] + newrow[j] - row[j-1] - row[j])/yr;
                /* XXX: I thought q alone (i.e., normal normalization) would
                 * give nice facet leveling, but the higher norm values has to
                 * be suppressed much more -- it seems */
                q = exp((vx*vx + vy*vy)/sigma2);
                sumvx += vx/q;
                sumvy += vy/q;
                sumvz += 1.0/q;
            }
        }
    }
    q = sumvz;
    *pbx = sumvx/q * xr;
    *pby = sumvy/q * yr;
    gwy_debug("beta=%g, sigma=%g sum=(%g, %g) q=%g b=(%g, %g)",
              sqrt(sigma2/c), sqrt(sigma2), sumvx, sumvy, q, *pbx, *pby);

    if (pa)
        *pa = -0.5*((*pbx)*xres + (*pby)*yres);

    return TRUE;
}

/**
 * gwy_data_field_plane_level:
 * @data_field: A data field.
 * @a: Constant coefficient.
 * @bx: X plane coefficient.
 * @by: Y plane coefficient.
 *
 * Subtracts plane from a data field.
 *
 * See gwy_data_field_fit_plane() for details.
 **/
void
gwy_data_field_plane_level(GwyDataField *data_field,
                           gdouble a, gdouble bx, gdouble by)
{
    gint i, j;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    for (i = 0; i < data_field->yres; i++) {
        gdouble *row = data_field->data + i*data_field->xres;
        gdouble rb = a + by*i;

        for (j = 0; j < data_field->xres; j++, row++)
            *row -= rb + bx*j;
    }

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_plane_rotate:
 * @data_field: A data field.
 * @xangle: Rotation angle in x direction (rotation along y axis, in radians).
 * @yangle: Rotation angle in y direction (rotation along x axis, in radians).
 * @interpolation: Interpolation type (can be only of two-point type).
 *
 * Performs rotation of plane along x and y axis.
 **/
void
gwy_data_field_plane_rotate(GwyDataField *data_field,
                            gdouble xangle,
                            gdouble yangle,
                            GwyInterpolationType interpolation)
{
    int k;
    GwyDataLine *l;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    if (xangle != 0) {
        l = gwy_data_line_new(data_field->xres, data_field->xreal, FALSE);
        for (k = 0; k < data_field->yres; k++) {
            gwy_data_field_get_row(data_field, l, k);
            gwy_data_line_rotate(l, -xangle, interpolation);
            gwy_data_field_set_row(data_field, l, k);
        }
        g_object_unref(l);
    }

    if (yangle != 0) {
        l = gwy_data_line_new(data_field->yres, data_field->yreal, FALSE);
        for (k = 0; k < data_field->xres; k++) {
            gwy_data_field_get_column(data_field, l, k);
            gwy_data_line_rotate(l, -yangle, interpolation);
            gwy_data_field_set_column(data_field, l, k);
        }
        g_object_unref(l);
    }

    gwy_data_field_invalidate(data_field);
}

#if 0
void
gwy_data_field_plane_true_rotate(GwyDataField *data_field,
                                 gdouble xangle,
                                 gdouble yangle,
                                 GwyInterpolationType interpolation)
{
    gdouble diag, dx, dy, phi, phi0, theta, tx, ty;
    gint xres, yres, txres, tyres, xbw, ybw, i;
    gdouble *data, *tdata;
    GwyDataField *tmp;

    if (xangle == 0 || yangle == 0) {
        gwy_data_field_plane_rotate(data_field, xangle, yangle, interpolation);
        return;
    }

    xres = data_field->xres;
    yres = data_field->yres;
    data = data_field->data;

    dx = tan(xangle);
    dy = tan(yangle);
    phi = atan2(dy, dx);
    theta = atan(hypot(dx, dy));
    phi0 = atan2(yres, xres);
    diag = hypot(xres, yres);
    tx = MAX(fabs(cos(-phi + phi0)), fabs(cos(-phi - phi0)));
    ty = MAX(fabs(sin(-phi + phi0)), fabs(sin(-phi - phi0)));
    txres = ((guint)GWY_ROUND(diag*tx + 2));
    tyres = ((guint)GWY_ROUND(diag*ty + 2));
    /* Keep parity to make the rotation less fuzzy */
    xbw = (txres - xres + 1)/2;
    if (xres + 2*xbw != txres)
        txres++;
    ybw = (tyres - yres + 1)/2;
    if (yres + 2*ybw != tyres)
        tyres++;

    /* Rotate to a temporary data field extended with border pixels */
    tmp = gwy_data_field_new(txres, tyres, 1.0, 1.0, FALSE);
    tdata = tmp->data;
    /* Copy */
    gwy_data_field_area_copy(data_field, tmp, 0, 0, xres, yres, xbw, ybw);
    /* Corners */
    gwy_data_field_area_fill(tmp, 0, 0, xbw, ybw,
                             data[0]);
    gwy_data_field_area_fill(tmp, xres + xbw, 0, xbw, ybw,
                             data[xres-1]);
    gwy_data_field_area_fill(tmp, 0, yres + ybw, xbw, ybw,
                             data[xres*(yres - 1)]);
    gwy_data_field_area_fill(tmp, xres + xbw, yres + ybw, xbw, ybw,
                             data[xres*yres - 1]);
    /* Sides */
    for (i = 0; i < ybw; i++)
        memcpy(tdata + i*txres + xbw, data,
               xres*sizeof(gdouble));
    for (i = 0; i < ybw; i++)
        memcpy(tdata + (yres + ybw + i)*txres + xbw, data + xres*(yres - 1),
               xres*sizeof(gdouble));
    for (i = 0; i < yres; i++) {
        gwy_data_field_area_fill(tmp, 0, ybw + i, xbw, 1,
                                 data[i*xres]);
        gwy_data_field_area_fill(tmp, xres + xbw, ybw + i, xbw, 1,
                                 data[i*xres + xres - 1]);
    }

    /* Rotate in xy to make the space rotation along y axis */
    gwy_data_field_rotate(tmp, -phi, interpolation);
    /* XXX: Still, individual gwy_data_line_rotate() can resample differently,
     * causing incompatible rows in the image.  And we cannot get the
     * resampling information from gwy_data_line_rotate(). */
    gwy_data_field_plane_rotate(tmp, theta, 0, GWY_INTERPOLATION_LINEAR);
    /* TODO:
     * recalculate xres
     * make samples square again
     */
    gwy_data_field_rotate(tmp, phi, interpolation);
    /* XXX: xbw is no longer correct border */
    gwy_data_field_area_copy(tmp, data_field, xbw, ybw, xres, yres, 0, 0);

    g_object_unref(tmp);
}
#endif

/**
 * gwy_data_field_fit_lines:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @degree: Fitted polynomial degree.
 * @exclude: If %TRUE, outside of area selected by @ulcol, @ulrow, @brcol,
 *           @brrow will be used for polynomial coefficients computation,
 *           instead of inside.
 * @orientation: Line orientation.
 *
 * Independently levels profiles on each row/column in a data field.
 *
 * Lines that have no intersection with area selected by @ulcol, @ulrow,
 * @brcol, @brrow are always leveled as a whole.  Lines that have intersection
 * with selected area, are leveled using polynomial coefficients computed only
 * from data inside (or outside for @exclude = %TRUE) the area.
 **/
void
gwy_data_field_fit_lines(GwyDataField *data_field,
                         gint col, gint row,
                         gint width, gint height,
                         gint degree,
                         gboolean exclude,
                         GwyOrientation orientation)
{

    gint i, j, xres, yres, res;
    gdouble real, coefs[4];
    GwyDataLine *hlp, *xdata = NULL, *ydata = NULL;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;

    xres = data_field->xres;
    yres = data_field->yres;
    res = (orientation == GWY_ORIENTATION_HORIZONTAL) ? xres : yres;
    real = (orientation == GWY_ORIENTATION_HORIZONTAL)
           ? data_field->xreal : data_field->yreal;
    hlp = gwy_data_line_new(res, real, FALSE);
    if (exclude) {
        xdata = gwy_data_line_new(res, real, FALSE);
        ydata = gwy_data_line_new(res, real, FALSE);
    }

    if (orientation == GWY_ORIENTATION_HORIZONTAL) {
        if (exclude) {
            for (i = j = 0; i < xres; i++) {
                if (i < col || i >= col + width)
                    xdata->data[j++] = i;
            }
        }

        for (i = 0; i < yres; i++) {
            gwy_data_field_get_row(data_field, hlp, i);
            if (i >= row && i < row + height) {
                if (exclude) {
                    gwy_assign(ydata->data, hlp->data, col);
                    gwy_assign(ydata->data + col, hlp->data + col + width,
                               xres - col - width);
                    gwy_math_fit_polynom(xres - width,
                                         xdata->data, ydata->data, degree,
                                         coefs);
                }
                else
                    gwy_data_line_part_fit_polynom(hlp, degree, coefs,
                                                   col, col + width);
            }
            else
                gwy_data_line_fit_polynom(hlp, degree, coefs);
            gwy_data_line_subtract_polynom(hlp, degree, coefs);
            gwy_data_field_set_row(data_field, hlp, i);
        }
    }
    else if (orientation == GWY_ORIENTATION_VERTICAL) {
        if (exclude) {
            for (i = j = 0; i < yres; i++) {
                if (i < row || i >= row + height)
                    xdata->data[j++] = i;
            }
        }

        for (i = 0; i < xres; i++) {
            gwy_data_field_get_column(data_field, hlp, i);
            if (i >= col && i < col + width) {
                if (exclude) {
                    gwy_assign(ydata->data, hlp->data, row);
                    gwy_assign(ydata->data + row, hlp->data + row + height,
                               yres - row - height);
                    gwy_math_fit_polynom(yres - height,
                                         xdata->data, ydata->data, degree,
                                         coefs);
                }
                else
                    gwy_data_line_part_fit_polynom(hlp, degree, coefs,
                                                   row, row + height);
            }
            else
                gwy_data_line_fit_polynom(hlp, degree, coefs);
            gwy_data_line_subtract_polynom(hlp, degree, coefs);
            gwy_data_field_set_column(data_field, hlp, i);
        }
    }
    g_object_unref(hlp);
    GWY_OBJECT_UNREF(xdata);
    GWY_OBJECT_UNREF(ydata);
}

/**
 * gwy_data_field_area_fit_polynom:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @col_degree: Degree of polynomial to fit column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to fit row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) to store the
 *          coefficients to, or %NULL (a fresh array is allocated then).
 *
 * Fits a two-dimensional polynomial to a rectangular part of a data field.
 *
 * The coefficients are stored by row into @coeffs, like data in a datafield.
 * Row index is y-degree, column index is x-degree.
 *
 * Note naive x^n y^m polynomial fitting is numerically unstable, therefore
 * this method works only up to @col_degree = @row_degree = 6.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 **/
gdouble*
gwy_data_field_area_fit_polynom(GwyDataField *data_field,
                                gint col, gint row,
                                gint width, gint height,
                                gint col_degree, gint row_degree,
                                gdouble *coeffs)
{
    gint i, j, size, xres, ns;
    gdouble *data, *sums, *m;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return NULL;
    g_return_val_if_fail(row_degree >= 0 && col_degree >= 0, NULL);
    xres = data_field->xres;

    data = data_field->data;
    size = (row_degree+1)*(col_degree+1);
    if (!coeffs)
        coeffs = g_new0(gdouble, size);
    else
        gwy_clear(coeffs, size);

    ns = (2*row_degree+1)*(2*col_degree+1);
    sums = g_new0(gdouble, ns);

#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(data,coeffs,sums,xres,row_degree,col_degree,row,col,width,height,size,ns)
#endif
    {
        gdouble *tsums = gwy_omp_if_threads_new0(sums, ns);
        gdouble *tcoeffs = gwy_omp_if_threads_new0(coeffs, size);
        gint ifrom = gwy_omp_chunk_start(height);
        gint ito = gwy_omp_chunk_end(height);

        for (i = row + ifrom; i < row + ito; i++) {
            for (j = col; j < col + width; j++) {
                gdouble ry = 1.0;
                gdouble z = data[i*xres + j];
                gint ix, iy;

                for (iy = 0; iy <= 2*row_degree; iy++) {
                    gdouble cx = 1.0;

                    for (ix = 0; ix <= 2*col_degree; ix++) {
                        sums[iy*(2*col_degree+1) + ix] += cx*ry;
                        cx *= j;
                    }
                    ry *= i;
                }

                ry = 1.0;
                for (iy = 0; iy <= row_degree; iy++) {
                    gdouble cx = 1.0;

                    for (ix = 0; ix <= col_degree; ix++) {
                        coeffs[iy*(col_degree+1) + ix] += cx*ry*z;
                        cx *= j;
                    }
                    ry *= i;
                }
            }
        }
        gwy_omp_if_threads_sum_double(coeffs, tcoeffs, size);
        gwy_omp_if_threads_sum_double(sums, tsums, ns);
    }

    m = g_new(gdouble, size*(size+1)/2);
    for (i = 0; i < size; i++) {
        gdouble *mrow = m + i*(i+1)/2;

        for (j = 0; j <= i; j++) {
            gint pow_x, pow_y;

            pow_x = i % (col_degree+1) + j % (col_degree+1);
            pow_y = i/(col_degree+1) + j/(col_degree+1);
            mrow[j] = sums[pow_y*(2*col_degree+1) + pow_x];
        }
    }

    if (!gwy_math_choleski_decompose(size, m))
        gwy_clear(coeffs, size);
    else
        gwy_math_choleski_solve(size, m, coeffs);

    g_free(m);
    g_free(sums);

    return coeffs;
}

/**
 * gwy_data_field_fit_polynom:
 * @data_field: A data field.
 * @col_degree: Degree of polynomial to fit column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to fit row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) to store the
 *          coefficients to, or %NULL (a fresh array is allocated then),
 *          see gwy_data_field_area_fit_polynom() for details.
 *
 * Fits a two-dimensional polynomial to a data field.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 **/
gdouble*
gwy_data_field_fit_polynom(GwyDataField *data_field,
                           gint col_degree, gint row_degree,
                           gdouble *coeffs)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    return gwy_data_field_area_fit_polynom(data_field, 0, 0,
                                           data_field->xres, data_field->yres,
                                           col_degree, row_degree, coeffs);
}

/**
 * gwy_data_field_area_subtract_polynom:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @col_degree: Degree of polynomial to subtract column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to subtract row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) with coefficients,
 *          see gwy_data_field_area_fit_polynom() for details.
 *
 * Subtracts a two-dimensional polynomial from a rectangular part of a data
 * field.
 **/
void
gwy_data_field_area_subtract_polynom(GwyDataField *data_field,
                                     gint col, gint row,
                                     gint width, gint height,
                                     gint col_degree, gint row_degree,
                                     const gdouble *coeffs)
{
    gint i, j, xres;
    gdouble *data;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(coeffs);
    g_return_if_fail(row_degree >= 0 && col_degree >= 0);
    xres = data_field->xres;
    data = data_field->data;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(data,coeffs,xres,row_degree,col_degree,row,col,width,height)
#endif
    for (i = row; i < row + height; i++) {
        for (j = col; j < col + width; j++) {
            gdouble ry = 1.0;
            gdouble z = data[i*xres + j];
            gint ix, iy;

            for (iy = 0; iy <= row_degree; iy++) {
                gdouble cx = 1.0;

                for (ix = 0; ix <= col_degree; ix++) {
                    /* FIXME: this is wrong, use Horner schema */
                    z -= coeffs[iy*(col_degree+1) + ix]*cx*ry;
                    cx *= j;
                }
                ry *= i;
            }

            data[i*xres + j] = z;
        }
    }

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_subtract_polynom:
 * @data_field: A data field.
 * @col_degree: Degree of polynomial to subtract column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to subtract row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) with coefficients,
 *          see gwy_data_field_area_fit_polynom() for details.
 *
 * Subtracts a two-dimensional polynomial from a data field.
 **/
void
gwy_data_field_subtract_polynom(GwyDataField *data_field,
                                gint col_degree, gint row_degree,
                                const gdouble *coeffs)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_subtract_polynom(data_field,
                                         0, 0,
                                         data_field->xres, data_field->yres,
                                         col_degree, row_degree, coeffs);
}

/* Calculate values of Legendre polynomials from 0 to @n in @x. */
static void
legendre_all(gdouble x,
             guint n,
             gdouble *p)
{
    guint m;

    p[0] = 1.0;
    if (n == 0)
        return;
    p[1] = x;
    if (n == 1)
        return;

    for (m = 2; m <= n; m++)
        p[m] = (x*(2*m - 1)*p[m-1] - (m - 1)*p[m-2])/m;
}

/**
 * gwy_data_field_area_fit_legendre:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @col_degree: Degree of polynomial to fit column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to fit row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) to store the
 *          coefficients to, or %NULL (a fresh array is allocated then).
 *
 * Fits two-dimensional Legendre polynomial to a rectangular part of a data
 * field.
 *
 * The @col_degree and @row_degree parameters limit the maximum powers of x and
 * y exactly as if simple powers were fitted, therefore if you do not intend to
 * interpret contents of @coeffs youself, the only difference is that this
 * method is much more numerically stable.
 *
 * The coefficients are organized exactly like in
 * gwy_data_field_area_fit_polynom(), but they are not coefficients of
 * x^n y^m, instead they are coefficients of P_n(x) P_m(x), where P are
 * Legendre polynomials.  The polynomials are evaluated in coordinates where
 * first row (column) corresponds to -1.0, and the last row (column) to 1.0.
 *
 * Note the polynomials are normal Legendre polynomials that are not exactly
 * orthogonal on a discrete point set (if their degrees are equal mod 2).
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 **/
gdouble*
gwy_data_field_area_fit_legendre(GwyDataField *data_field,
                                 gint col, gint row,
                                 gint width, gint height,
                                 gint col_degree, gint row_degree,
                                 gdouble *coeffs)
{
    gint r, c, i, j, size, maxsize, xres, col_n, row_n;
    gint isize, jsize, thissize;
    gdouble *data, *m, *pmx, *pmy, *sumsx, *sumsy, *rhs;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return NULL;
    g_return_val_if_fail(row_degree >= 0 && col_degree >= 0, NULL);
    xres = data_field->xres;
    data = data_field->data;
    col_n = col_degree + 1;
    row_n = row_degree + 1;

    size = col_n*row_n;
    /* The maximum necessary matrix size (order), it is approximately four
     * times smaller thanks to separation of even and odd polynomials */
    maxsize = ((col_n + 1)/2)*((row_n + 1)/2);
    if (!coeffs)
        coeffs = g_new0(gdouble, size);
    else
        gwy_clear(coeffs, size);

    sumsx = g_new0(gdouble, col_n*col_n);
    sumsy = g_new0(gdouble, row_n*row_n);
    rhs = g_new(gdouble, maxsize);
    m = g_new(gdouble, MAX(maxsize*(maxsize + 1)/2, col_n + row_n));
    /* pmx, pmy and m are not needed at the same time, reuse it */
    pmx = m;
    pmy = m + col_n;

    /* Calculate <P_m(x) P_n(y) z(x,y)> (normalized to complete area) */
    for (r = 0; r < height; r++) {
        legendre_all(2*r/(height - 1.0) - 1.0, row_degree, pmy);
        for (c = 0; c < width; c++) {
            gdouble z = data[(row + r)*xres + (col + c)];

            legendre_all(2*c/(width - 1.0) - 1.0, col_degree, pmx);
            for (i = 0; i < row_n; i++) {
                for (j = 0; j < col_n; j++)
                    coeffs[i*col_n + j] += z*pmx[j]*pmy[i];
            }
        }
    }

    /* Calculate <P_m(x) P_a(x)> (normalized to single row).
     * 3/4 of these values are zeroes, but it only takes O(width) time. */
    for (c = 0; c < width; c++) {
        legendre_all(2*c/(width - 1.0) - 1.0, col_degree, pmx);
        for (i = 0; i < col_n; i++) {
            for (j = 0; j < col_n; j++)
                sumsx[i*col_n + j] += pmx[i]*pmx[j];
        }
    }

    /* Calculate <P_n(y) P_b(y)> (normalized to single column)
     * 3/4 of these values are zeroes, but it only takes O(height) time. */
    for (r = 0; r < height; r++) {
        legendre_all(2*r/(height - 1.0) - 1.0, row_degree, pmy);
        for (i = 0; i < row_n; i++) {
            for (j = 0; j < row_n; j++)
                sumsy[i*row_n + j] += pmy[i]*pmy[j];
        }
    }

    /* (Even, Even) */
    isize = (row_n + 1)/2;
    jsize = (col_n + 1)/2;
    thissize = jsize*isize;
    /* This is always true */
    if (thissize) {
        /* Construct the submatrix */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize);
            gint iy = 2*(i/jsize);
            gdouble *mrow = m + i*(i + 1)/2;

            for (j = 0; j <= i; j++) {
                gint jx = 2*(j % jsize);
                gint jy = 2*(j/jsize);

                mrow[j] = sumsx[ix*col_n + jx]*sumsy[iy*row_n + jy];
            }
        }
        /* Construct the subrhs */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize);
            gint iy = 2*(i/jsize);

            rhs[i] = coeffs[iy*col_n + ix];
        }
        /* Solve */
        if (!gwy_math_choleski_decompose(thissize, m)) {
            gwy_clear(coeffs, size);
            goto fail;
        }
        gwy_math_choleski_solve(thissize, m, rhs);
        /* Copy back */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize);
            gint iy = 2*(i/jsize);

            coeffs[iy*col_n + ix] = rhs[i];
        }
    }

    /* (Even, Odd) */
    isize = (row_n + 1)/2;
    jsize = col_n/2;
    thissize = jsize*isize;
    if (thissize) {
        /* Construct the submatrix */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize) + 1;
            gint iy = 2*(i/jsize);
            gdouble *mrow = m + i*(i + 1)/2;

            for (j = 0; j <= i; j++) {
                gint jx = 2*(j % jsize) + 1;
                gint jy = 2*(j/jsize);

                mrow[j] = sumsx[ix*col_n + jx]*sumsy[iy*row_n + jy];
            }
        }
        /* Construct the subrhs */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize) + 1;
            gint iy = 2*(i/jsize);

            rhs[i] = coeffs[iy*col_n + ix];
        }
        /* Solve */
        if (!gwy_math_choleski_decompose(thissize, m)) {
            gwy_clear(coeffs, size);
            goto fail;
        }
        gwy_math_choleski_solve(thissize, m, rhs);
        /* Copy back */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize) + 1;
            gint iy = 2*(i/jsize);

            coeffs[iy*col_n + ix] = rhs[i];
        }
    }

    /* (Odd, Even) */
    isize = row_n/2;
    jsize = (col_n + 1)/2;
    thissize = jsize*isize;
    if (thissize) {
        /* Construct the submatrix */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize);
            gint iy = 2*(i/jsize) + 1;
            gdouble *mrow = m + i*(i + 1)/2;

            for (j = 0; j <= i; j++) {
                gint jx = 2*(j % jsize);
                gint jy = 2*(j/jsize) + 1;

                mrow[j] = sumsx[ix*col_n + jx]*sumsy[iy*row_n + jy];
            }
        }
        /* Construct the subrhs */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize);
            gint iy = 2*(i/jsize) + 1;

            rhs[i] = coeffs[iy*col_n + ix];
        }
        /* Solve */
        if (!gwy_math_choleski_decompose(thissize, m)) {
            gwy_clear(coeffs, size);
            goto fail;
        }
        gwy_math_choleski_solve(thissize, m, rhs);
        /* Copy back */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize);
            gint iy = 2*(i/jsize) + 1;

            coeffs[iy*col_n + ix] = rhs[i];
        }
    }

    /* (Odd, Odd) */
    isize = row_n/2;
    jsize = col_n/2;
    thissize = jsize*isize;
    if (thissize) {
        /* Construct the submatrix */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize) + 1;
            gint iy = 2*(i/jsize) + 1;
            gdouble *mrow = m + i*(i + 1)/2;

            for (j = 0; j <= i; j++) {
                gint jx = 2*(j % jsize) + 1;
                gint jy = 2*(j/jsize) + 1;

                mrow[j] = sumsx[ix*col_n + jx]*sumsy[iy*row_n + jy];
            }
        }
        /* Construct the subrhs */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize) + 1;
            gint iy = 2*(i/jsize) + 1;

            rhs[i] = coeffs[iy*col_n + ix];
        }
        /* Solve */
        if (!gwy_math_choleski_decompose(thissize, m)) {
            gwy_clear(coeffs, size);
            goto fail;
        }
        gwy_math_choleski_solve(thissize, m, rhs);
        /* Copy back */
        for (i = 0; i < thissize; i++) {
            gint ix = 2*(i % jsize) + 1;
            gint iy = 2*(i/jsize) + 1;

            coeffs[iy*col_n + ix] = rhs[i];
        }
    }

fail:
    g_free(m);
    g_free(rhs);
    g_free(sumsx);
    g_free(sumsy);

    return coeffs;
}

/**
 * gwy_data_field_fit_legendre:
 * @data_field: A data field.
 * @col_degree: Degree of polynomial to fit column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to fit row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) to store the
 *          coefficients to, or %NULL (a fresh array is allocated then).
 *
 * Fits two-dimensional Legendre polynomial to a data field.
 *
 * See gwy_data_field_area_fit_legendre() for details.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 **/
gdouble*
gwy_data_field_fit_legendre(GwyDataField *data_field,
                            gint col_degree, gint row_degree,
                            gdouble *coeffs)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    return gwy_data_field_area_fit_legendre(data_field, 0, 0,
                                            data_field->xres, data_field->yres,
                                            col_degree, row_degree, coeffs);
}

/**
 * gwy_data_field_area_subtract_legendre:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @col_degree: Degree of polynomial to subtract column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to subtract row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) with coefficients,
 *          see gwy_data_field_area_fit_legendre() for details.
 *
 * Subtracts a two-dimensional Legendre polynomial fit from a rectangular part
 * of a data field.
 *
 * Due to the transform of coordinates to [-1,1] x [-1,1], this method can be
 * used on an area of dimensions different than the area the coefficients were
 * calculated for.
 **/
void
gwy_data_field_area_subtract_legendre(GwyDataField *data_field,
                                      gint col, gint row,
                                      gint width, gint height,
                                      gint col_degree, gint row_degree,
                                      const gdouble *coeffs)
{
    gint r, c, i, j, xres, col_n, row_n;
    gdouble *data, *pmx, *pmy;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(coeffs);
    g_return_if_fail(row_degree >= 0 && col_degree >= 0);
    xres = data_field->xres;
    data = data_field->data;
    col_n = col_degree + 1;
    row_n = row_degree + 1;

    pmx = g_new0(gdouble, col_n + row_n);
    pmy = pmx + col_n;

    for (r = 0; r < height; r++) {
        legendre_all(2*r/(height - 1.0) - 1.0, row_degree, pmy);
        for (c = 0; c < width; c++) {
            gdouble z = data[(row + r)*xres + (col + c)];

            legendre_all(2*c/(width - 1.0) - 1.0, col_degree, pmx);
            for (i = 0; i < row_n; i++) {
                for (j = 0; j < col_n; j++)
                    z -= coeffs[i*col_n + j]*pmx[j]*pmy[i];
            }

            data[(row + r)*xres + (col + c)] = z;
        }
    }

    g_free(pmx);

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_subtract_legendre:
 * @data_field: A data field.
 * @col_degree: Degree of polynomial to subtract column-wise (x-coordinate).
 * @row_degree: Degree of polynomial to subtract row-wise (y-coordinate).
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+1) with coefficients,
 *          see gwy_data_field_area_fit_legendre() for details.
 *
 * Subtracts a two-dimensional Legendre polynomial fit from a data field.
 **/
void
gwy_data_field_subtract_legendre(GwyDataField *data_field,
                                 gint col_degree, gint row_degree,
                                 const gdouble *coeffs)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_subtract_legendre(data_field,
                                          0, 0,
                                          data_field->xres, data_field->yres,
                                          col_degree, row_degree, coeffs);
}

/**
 * gwy_data_field_area_fit_poly_max:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @max_degree: Maximum total polynomial degree, that is the maximum of m+n
 *              in x^n y^m terms.
 * @coeffs: An array of size (@max_degree+1)*(@max_degree+2)/2 to store the
 *          coefficients to, or %NULL (a fresh array is allocated then).
 *
 * Fits two-dimensional polynomial with limited total degree to a rectangular
 * part of a data field.
 *
 * See gwy_data_field_area_fit_legendre() for description.  This function
 * differs by limiting the total maximum degree, while
 * gwy_data_field_area_fit_legendre() limits the maximum degrees in horizontal
 * and vertical directions independently.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 **/
gdouble*
gwy_data_field_area_fit_poly_max(GwyDataField *data_field,
                                 gint col, gint row,
                                 gint width, gint height,
                                 gint max_degree,
                                 gdouble *coeffs)
{
    gint r, c, i, j, size, xres, degree_n;
    gint ix, jx, iy, jy;
    gdouble *data, *m, *pmx, *pmy, *sumsx, *sumsy;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return NULL;
    g_return_val_if_fail(max_degree >= 0, NULL);
    xres = data_field->xres;
    data = data_field->data;
    degree_n = max_degree + 1;
    size = degree_n*(degree_n + 1)/2;
    g_return_val_if_fail(width*height > size, NULL);
    /* The maximum necessary matrix size (order), it is approximately four
     * times smaller thanks to separation of even and odd polynomials */
    if (!coeffs)
        coeffs = g_new0(gdouble, size);
    else
        gwy_clear(coeffs, size);

    sumsx = g_new0(gdouble, degree_n*degree_n);
    sumsy = g_new0(gdouble, degree_n*degree_n);
    m = g_new(gdouble, MAX(size*(size + 1)/2, 2*degree_n));
    /* pmx, pmy and m are not needed at the same time, reuse it */
    pmx = m;
    pmy = m + degree_n;

    /* Calculate <P_m(x) P_n(y) z(x,y)> (normalized to complete area) */
    for (r = 0; r < height; r++) {
        legendre_all(2*r/(height - 1.0) - 1.0, max_degree, pmy);
        for (c = 0; c < width; c++) {
            gdouble z = data[(row + r)*xres + (col + c)];

            legendre_all(2*c/(width - 1.0) - 1.0, max_degree, pmx);
            for (i = 0; i < degree_n; i++) {
                for (j = 0; j < degree_n - i; j++)
                    coeffs[i*(2*degree_n + 1 - i)/2 + j] += z*pmx[j]*pmy[i];
            }
        }
    }

    /* Calculate <P_m(x) P_a(x)> (normalized to single row).
     * 3/4 of these values are zeroes, but it only takes O(width) time. */
    for (c = 0; c < width; c++) {
        legendre_all(2*c/(width - 1.0) - 1.0, max_degree, pmx);
        for (i = 0; i < degree_n; i++) {
            for (j = 0; j < degree_n; j++)
                sumsx[i*degree_n + j] += pmx[i]*pmx[j];
        }
    }

    /* Calculate <P_n(y) P_b(y)> (normalized to single column)
     * 3/4 of these values are zeroes, but it only takes O(height) time. */
    for (r = 0; r < height; r++) {
        legendre_all(2*r/(height - 1.0) - 1.0, max_degree, pmy);
        for (i = 0; i < degree_n; i++) {
            for (j = 0; j < degree_n; j++)
                sumsy[i*degree_n + j] += pmy[i]*pmy[j];
        }
    }

    /* Construct the matrix */
    for (iy = 0; iy < degree_n; iy++) {
        for (jy = 0; jy < degree_n - iy; jy++) {
            gdouble *mrow;

            i = iy*(2*degree_n + 1 - iy)/2 + jy;
            mrow = m + i*(i + 1)/2;
            for (ix = 0; ix < degree_n; ix++) {
                for (jx = 0; jx < degree_n - ix; jx++) {
                    j = ix*(2*degree_n + 1 - ix)/2 + jx;
                    /* It is easier to go through all the coeffs and ignore
                     * the upper right triangle than to construct conditions
                     * directly for jy, jy, etc. */
                    if (j > i)
                        continue;
                    mrow[j] = sumsx[jy*degree_n + jx]*sumsy[iy*degree_n + ix];
                }
            }
        }
    }
    /* Solve */
    if (!gwy_math_choleski_decompose(size, m)) {
        gwy_clear(coeffs, size);
        goto fail;
    }
    gwy_math_choleski_solve(size, m, coeffs);

fail:
    g_free(m);
    g_free(sumsx);
    g_free(sumsy);

    return coeffs;
}

/**
 * gwy_data_field_fit_poly_max:
 * @data_field: A data field.
 * @max_degree: Maximum total polynomial degree, that is the maximum of m+n
 *              in x^n y^m terms.
 * @coeffs: An array of size (@max_degree+1)*(@max_degree+2)/2 to store the
 *          coefficients to, or %NULL (a fresh array is allocated then).
 *
 * Fits two-dimensional polynomial with limited total degree to a data field.
 *
 * See gwy_data_field_area_fit_poly_max() for details.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 **/
gdouble*
gwy_data_field_fit_poly_max(GwyDataField *data_field,
                            gint max_degree,
                            gdouble *coeffs)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    return gwy_data_field_area_fit_poly_max(data_field, 0, 0,
                                            data_field->xres, data_field->yres,
                                            max_degree, coeffs);
}

/**
 * gwy_data_field_area_subtract_poly_max:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @max_degree: Maximum total polynomial degree, that is the maximum of m+n
 *              in x^n y^m terms.
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+2)/2 with
 *          coefficients, see gwy_data_field_area_fit_poly_max() for details.
 *
 * Subtracts a two-dimensional polynomial with limited total degree from a
 * rectangular part of a data field.
 *
 * Due to the transform of coordinates to [-1,1] x [-1,1], this method can be
 * used on an area of dimensions different than the area the coefficients were
 * calculated for.
 **/
void
gwy_data_field_area_subtract_poly_max(GwyDataField *data_field,
                                      gint col, gint row,
                                      gint width, gint height,
                                      gint max_degree,
                                      const gdouble *coeffs)
{
    gint r, c, i, j, xres, degree_n;
    gdouble *data, *pmx, *pmy;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(coeffs);
    g_return_if_fail(max_degree >= 0);
    xres = data_field->xres;
    data = data_field->data;
    degree_n = max_degree + 1;

    pmx = g_new0(gdouble, 2*degree_n);
    pmy = pmx + degree_n;

    for (r = 0; r < height; r++) {
        legendre_all(2*r/(height - 1.0) - 1.0, max_degree, pmy);
        for (c = 0; c < width; c++) {
            gdouble z = data[(row + r)*xres + (col + c)];

            legendre_all(2*c/(width - 1.0) - 1.0, max_degree, pmx);
            for (i = 0; i < degree_n; i++) {
                for (j = 0; j < degree_n - i; j++)
                    z -= coeffs[i*(2*degree_n + 1 - i)/2 + j]*pmx[j]*pmy[i];
            }

            data[(row + r)*xres + (col + c)] = z;
        }
    }

    g_free(pmx);

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_subtract_poly_max:
 * @data_field: A data field.
 * @max_degree: Maximum total polynomial degree, that is the maximum of m+n
 *              in x^n y^m terms.
 * @coeffs: An array of size (@row_degree+1)*(@col_degree+2)/2 with
 *          coefficients, see gwy_data_field_area_fit_poly_max() for details.
 *
 * Subtracts a two-dimensional polynomial with limited total degree from
 * a data field.
 **/
void
gwy_data_field_subtract_poly_max(GwyDataField *data_field,
                                 gint max_degree,
                                 const gdouble *coeffs)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_subtract_poly_max(data_field,
                                          0, 0,
                                          data_field->xres, data_field->yres,
                                          max_degree, coeffs);
}

/**
 * gwy_data_field_area_fit_poly:
 * @data_field: A data field.
 * @mask_field: Mask of values to take values into account, or %NULL for full
 *        @data_field.  Values equal to 0.0 and below cause corresponding
 *        @data_field samples to be ignored, values equal to 1.0 and above
 *        cause inclusion of corresponding @data_field samples.  The behaviour
 *        for values inside (0.0, 1.0) is undefined (it may be specified
 *        in the future).
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @nterms: The number of polynomial terms to take into account (half the
 *          number of items in @term_powers).
 * @term_powers: Array of size 2*@nterms describing the terms to fit.  Each
 *               terms is described by a couple of powers (powerx, powery).
 * @exclude: Interpret values @w in the mask as 1.0-@w.
 * @coeffs: Array of size @nterms to store the coefficients to, or %NULL to
 *          allocate a new array.
 *
 * Fit a given set of polynomial terms to a rectangular part of a data field.
 *
 * The polynomial coefficients correspond to normalized coordinates that
 * are always from the interval [-1,1] where -1 corresponds to the left/topmost
 * pixel and 1 corresponds to the bottom/rightmost pixel of the area.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 *
 * Since: 2.11
 **/
gdouble*
gwy_data_field_area_fit_poly(GwyDataField *data_field,
                             GwyDataField *mask_field,
                             gint col, gint row,
                             gint width, gint height,
                             gint nterms,
                             const gint *term_powers,
                             gboolean exclude,
                             gdouble *coeffs)
{
    const gdouble *data, *mask;
    gint xres;
    guint np;
    gdouble *m;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height)
        || !_gwy_data_field_check_mask(data_field, &mask_field, NULL))
        return NULL;
    g_return_val_if_fail(nterms >= 0, NULL);
    xres = data_field->xres;

    if (!nterms)
        return coeffs;

    data = data_field->data;
    mask = mask_field ? mask_field->data : NULL;

    if (!coeffs)
        coeffs = g_new0(gdouble, nterms);
    else
        gwy_clear(coeffs, nterms);

    np = nterms*(nterms + 1)/2;
    m = g_new0(gdouble, np);

#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
        shared(data,mask,m,coeffs,xres,np,width,height,row,col,nterms,exclude,term_powers)
#endif
    {
        gint ifrom = gwy_omp_chunk_start(height);
        gint ito = gwy_omp_chunk_end(height);
        gdouble *tm = gwy_omp_if_threads_new0(m, np);
        gdouble *tcoeffs = gwy_omp_if_threads_new0(coeffs, nterms);
        gdouble *p, *py;
        gint i, j, k, l;

        p = g_new(gdouble, 2*nterms);
        py = p + nterms;
        for (i = ifrom; i < ito; i++) {
            gdouble y = 2*i/(height - 1.0) - 1.0;

            for (k = 0; k < nterms; k++)
                py[k] = gwy_powi(y, term_powers[2*k + 1]);

            for (j = 0; j < width; j++) {
                gdouble x = 2*j/(width - 1.0) - 1.0;
                gdouble z = data[(row + i)*xres + (col + j)];
                gdouble *mm = tm;

                if (mask) {
                    gdouble w = mask[(row + i)*xres + (col + j)];
                    if ((exclude && w > 0.0)
                        || (!exclude && w <= 0.0))
                        continue;
                }

                for (k = 0; k < nterms; k++)
                    p[k] = gwy_powi(x, term_powers[2*k]) * py[k];

                for (k = 0; k < nterms; k++) {
                    for (l = 0; l <= k; l++, mm++)
                        *mm += p[k]*p[l];
                    tcoeffs[k] += z*p[k];
                }
            }
        }
        g_free(p);

        gwy_omp_if_threads_sum_double(coeffs, tcoeffs, nterms);
        gwy_omp_if_threads_sum_double(m, tm, np);
    }

    if (!gwy_math_choleski_decompose(nterms, m))
        gwy_clear(coeffs, nterms);
    else
        gwy_math_choleski_solve(nterms, m, coeffs);

    g_free(m);

    return coeffs;
}

/**
 * gwy_data_field_fit_poly:
 * @data_field: A data field.
 * @mask_field: Mask of values to take values into account, or %NULL for full
 *        @data_field.  Values equal to 0.0 and below cause corresponding
 *        @data_field samples to be ignored, values equal to 1.0 and above
 *        cause inclusion of corresponding @data_field samples.  The behaviour
 *        for values inside (0.0, 1.0) is undefined (it may be specified
 *        in the future).
 * @nterms: The number of polynomial terms to take into account (half the
 *          number of items in @term_powers).
 * @term_powers: Array of size 2*@nterms describing the terms to fit.  Each
 *               terms is described by a couple of powers (powerx, powery).
 * @exclude: Interpret values @w in the mask as 1.0-@w.
 * @coeffs: Array of size @nterms to store the coefficients to, or %NULL to
 *          allocate a new array.
 *
 * Fit a given set of polynomial terms to a data field.
 *
 * Returns: Either @coeffs if it was not %NULL, or a newly allocated array
 *          with coefficients.
 *
 * Since: 2.11
 **/
gdouble*
gwy_data_field_fit_poly(GwyDataField *data_field,
                        GwyDataField *mask_field,
                        gint nterms,
                        const gint *term_powers,
                        gboolean exclude,
                        gdouble *coeffs)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    return gwy_data_field_area_fit_poly(data_field, mask_field,
                                        0, 0,
                                        data_field->xres, data_field->yres,
                                        nterms, term_powers, exclude,
                                        coeffs);
}

/**
 * gwy_data_field_area_subtract_poly:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @nterms: The number of polynomial terms to take into account (half the
 *          number of items in @term_powers).
 * @term_powers: Array of size 2*@nterms describing the fitted terms.  Each
 *               terms is described by a couple of powers (powerx, powery).
 * @coeffs: Array of size @nterms to store with the coefficients.
 *
 * Subtract a given set of polynomial terms from a rectangular part of a data
 * field.
 *
 * Since: 2.11
 **/
void
gwy_data_field_area_subtract_poly(GwyDataField *data_field,
                                  gint col, gint row,
                                  gint width, gint height,
                                  gint nterms,
                                  const gint *term_powers,
                                  const gdouble *coeffs)
{
    gdouble *data;
    gint xres;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(nterms >= 0);
    g_return_if_fail(coeffs);

    if (!nterms)
        return;

    xres = data_field->xres;
    data = data_field->data;

#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            shared(data,coeffs,xres,width,height,row,col,nterms,term_powers)
#endif
    {
        gdouble *py = g_new(gdouble, nterms);
        gint ifrom = gwy_omp_chunk_start(height);
        gint ito = gwy_omp_chunk_end(height);
        gint i, j, k;

        for (i = ifrom; i < ito; i++) {
            gdouble y = 2*i/(height - 1.0) - 1.0;

            for (k = 0; k < nterms; k++)
                py[k] = gwy_powi(y, term_powers[2*k + 1]);

            for (j = 0; j < width; j++) {
                gdouble x = 2*j/(width - 1.0) - 1.0;
                gdouble z = data[(row + i)*xres + (col + j)];

                for (k = 0; k < nterms; k++)
                    z -= coeffs[k] * gwy_powi(x, term_powers[2*k]) * py[k];

                data[(row + i)*xres + (col + j)] = z;
            }
        }
        g_free(py);
    }

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_subtract_poly:
 * @data_field: A data field.
 * @nterms: The number of polynomial terms to take into account (half the
 *          number of items in @term_powers).
 * @term_powers: Array of size 2*@nterms describing the fitter terms.  Each
 *               terms is described by a couple of powers (powerx, powery).
 * @coeffs: Array of size @nterms to store with the coefficients.
 *
 * Subtract a given set of polynomial terms from a data field.
 *
 * Since: 2.11
 **/
void
gwy_data_field_subtract_poly(GwyDataField *data_field,
                             gint nterms,
                             const gint *term_powers,
                             const gdouble *coeffs)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_subtract_poly(data_field,
                                      0, 0,
                                      data_field->xres, data_field->yres,
                                      nterms, term_powers, coeffs);
}

/**
 * gwy_data_field_area_fit_local_planes:
 * @data_field: A data field.
 * @size: Neighbourhood size (must be at least 2).  It is centered around
 *        each pixel, unless @size is even when it sticks to the right.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @nresults: The number of requested quantities.
 * @types: The types of requested quantities.
 * @results: An array to store quantities to, may be %NULL to allocate a new
 *           one which must be freed by caller then.  If any item is %NULL,
 *           a new data field is allocated for it, existing data fields
 *           are resized to @width x @height.
 *
 * Fits a plane through neighbourhood of each sample in a rectangular part
 * of a data field.
 *
 * The sample is always in the origin of its local (x,y) coordinate system,
 * even if the neighbourhood is not centered about it (e.g. because sample
 * is on the edge of data field).  Z-coordinate is however not centered,
 * that is @GWY_PLANE_FIT_A is normal mean value.
 *
 * Returns: An array of data fields with requested quantities, that is
 *          @results unless it was %NULL and a new array was allocated.
 **/
GwyDataField**
gwy_data_field_area_fit_local_planes(GwyDataField *data_field,
                                     gint size,
                                     gint col, gint row,
                                     gint width, gint height,
                                     gint nresults,
                                     const GwyPlaneFitQuantity *types,
                                     GwyDataField **results)
{
    gdouble xreal, yreal, qx, qy, asymshfit;
    gint xres, yres, ri, i, j;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return NULL;
    g_return_val_if_fail(size > 1, NULL);
    if (!nresults)
        return NULL;
    g_return_val_if_fail(nresults > 0, NULL);
    g_return_val_if_fail(types, NULL);
    for (ri = 0; ri < nresults; ri++) {
        g_return_val_if_fail(types[ri] >= GWY_PLANE_FIT_A
                             && types[ri] <= GWY_PLANE_FIT_S0_REDUCED,
                             NULL);
        g_return_val_if_fail(!results
                             || !results[ri]
                             || GWY_IS_DATA_FIELD(results[ri]),
                             NULL);
    }
    if (!results)
        results = g_new0(GwyDataField*, nresults);

    /* Allocate output data fields or fix their dimensions */
    xres = data_field->xres;
    yres = data_field->yres;
    qx = data_field->xreal/xres;
    qy = data_field->yreal/yres;
    xreal = qx*width;
    yreal = qy*height;
    for (ri = 0; ri < nresults; ri++) {
        if (!results[ri])
            results[ri] = gwy_data_field_new(width, height, xreal, yreal,
                                             FALSE);
        else {
            gwy_data_field_resample(results[ri], width, height,
                                    GWY_INTERPOLATION_NONE);
            gwy_data_field_set_xreal(results[ri], xreal);
            gwy_data_field_set_yreal(results[ri], yreal);
        }
    }

    /* Fit local planes */
    asymshfit = (1 - size % 2)/2.0;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j,ri) \
            shared(data_field,results,xres,yres,col,row,width,height,size,asymshfit,qx,qy,types,nresults)
#endif
    for (i = 0; i < height; i++) {
        gint ifrom = MAX(0, i + row - (size-1)/2);
        gint ito = MIN(yres-1, i + row + size/2);

        /* Prevent fitting plane through just one pixel on bottom edge when
         * size == 2 */
        if (G_UNLIKELY(ifrom == ito) && ifrom)
            ifrom--;

        for (j = 0; j < width; j++) {
            gint jfrom = MAX(0, j + col - (size-1)/2);
            gint jto = MIN(xres-1, j + col + size/2);
            gdouble *drect;
            gdouble sumz, sumzx, sumzy, sumzz, sumx, sumy, sumxx, sumxy, sumyy;
            gdouble n, bx, by, s0, s0r, det, shift;
            gdouble coeffs[GWY_PLANE_FIT_S0_REDUCED + 1];
            gint ii, jj;

            /* Prevent fitting plane through just one pixel on right edge when
             * size == 2 */
            if (G_UNLIKELY(jfrom == jto) && jfrom)
                jfrom--;

            drect = data_field->data + ifrom*xres + jfrom;
            /* Compute sums with origin in top left corner */
            sumz = sumzx = sumzy = sumzz = 0.0;
            for (ii = 0; ii <= ito - ifrom; ii++) {
                gdouble *drow = drect + xres*ii;

                for (jj = 0; jj <= jto - jfrom; jj++) {
                    sumz += drow[jj];
                    sumzx += drow[jj]*jj;
                    sumzy += drow[jj]*ii;
                    sumzz += drow[jj]*drow[jj];
                }
            }
            n = (ito - ifrom + 1)*(jto - jfrom + 1);
            sumx = n*(jto - jfrom)/2.0;
            sumy = n*(ito - ifrom)/2.0;
            sumxx = sumx*(2*(jto - jfrom) + 1)/3.0;
            sumyy = sumy*(2*(ito - ifrom) + 1)/3.0;
            sumxy = sumx*sumy/n;

            /* Move origin to pixel, including in z coordinate, remembering
             * average z value in shift */
            shift = ifrom - (i + row + asymshfit);
            sumxy += shift*sumx;
            sumyy += shift*(2*sumy + n*shift);
            sumzy += shift*sumz;
            sumy += n*shift;

            shift = jfrom - (j + col + asymshfit);
            sumxx += shift*(2*sumx + n*shift);
            sumxy += shift*sumy;
            sumzx += shift*sumz;
            sumx += n*shift;

            shift = -sumz/n;
            sumzx += shift*sumx;
            sumzy += shift*sumy;
            sumzz += shift*(2*sumz + n*shift);
            /* sumz = 0.0;  unused */

            /* Compute coefficients */
            det = sumxx*sumyy - sumxy*sumxy;
            bx = (sumzx*sumyy - sumxy*sumzy)/det;
            by = (sumzy*sumxx - sumxy*sumzx)/det;
            s0 = sumzz - bx*sumzx - by*sumzy;
            s0r = s0/(1.0 + bx*bx/qx/qx + by*by/qy/qy);

            coeffs[GWY_PLANE_FIT_A] = -shift;
            coeffs[GWY_PLANE_FIT_BX] = bx;
            coeffs[GWY_PLANE_FIT_BY] = by;
            coeffs[GWY_PLANE_FIT_ANGLE] = atan2(by, bx);
            coeffs[GWY_PLANE_FIT_SLOPE] = sqrt(bx*bx + by*by);
            coeffs[GWY_PLANE_FIT_S0] = s0;
            coeffs[GWY_PLANE_FIT_S0_REDUCED] = s0r;

            for (ri = 0; ri < nresults; ri++)
                results[ri]->data[width*i + j] = coeffs[types[ri]];
        }
    }

    for (ri = 0; ri < nresults; ri++)
        gwy_data_field_invalidate(results[ri]);

    return results;
}

/**
 * gwy_data_field_fit_local_planes:
 * @data_field: A data field.
 * @size: Neighbourhood size.
 * @nresults: The number of requested quantities.
 * @types: The types of requested quantities.
 * @results: An array to store quantities to.
 *
 * Fits a plane through neighbourhood of each sample in a data field.
 *
 * See gwy_data_field_area_fit_local_planes() for details.
 *
 * Returns: An array of data fields with requested quantities.
 **/
GwyDataField**
gwy_data_field_fit_local_planes(GwyDataField *data_field,
                                gint size,
                                gint nresults,
                                const GwyPlaneFitQuantity *types,
                                GwyDataField **results)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    return gwy_data_field_area_fit_local_planes(data_field, size,
                                                0, 0,
                                                data_field->xres,
                                                data_field->yres,
                                                nresults, types, results);
}

/**
 * gwy_data_field_area_local_plane_quantity:
 * @data_field: A data field.
 * @size: Neighbourhood size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @type: The type of requested quantity.
 * @result: A data field to store result to, or %NULL to allocate a new one.
 *
 * Convenience function to get just one quantity from
 * gwy_data_field_area_fit_local_planes().
 *
 * Returns: @result if it isn't %NULL, otherwise a newly allocated data field.
 **/
GwyDataField*
gwy_data_field_area_local_plane_quantity(GwyDataField *data_field,
                                         gint size,
                                         gint col, gint row,
                                         gint width, gint height,
                                         GwyPlaneFitQuantity type,
                                         GwyDataField *result)
{
    gwy_data_field_area_fit_local_planes(data_field, size,
                                         col, row, width, height,
                                         1, &type, &result);
    return result;
}

/**
 * gwy_data_field_local_plane_quantity:
 * @data_field: A data field.
 * @size: Neighbourhood size.
 * @type: The type of requested quantity.
 * @result: A data field to store result to, or %NULL to allocate a new one.
 *
 * Convenience function to get just one quantity from
 * gwy_data_field_fit_local_planes().
 *
 * Returns: @result if it isn't %NULL, otherwise a newly allocated data field.
 **/
GwyDataField*
gwy_data_field_local_plane_quantity(GwyDataField *data_field,
                                    gint size,
                                    GwyPlaneFitQuantity type,
                                    GwyDataField *result)
{
    gwy_data_field_fit_local_planes(data_field, size, 1, &type, &result);

    return result;
}

/************************** Documentation ****************************/

/**
 * SECTION:level
 * @title: level
 * @short_description: Leveling and background removal
 **/

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
