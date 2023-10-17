/*
 *  $Id: filters.c 25308 2023-04-21 13:16:28Z yeti-dn $
 *  Copyright (C) 2003-2021 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
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

#include <string.h>
#include <libgwyddion/gwymacros.h>
#include <libprocess/filters.h>
#include <libprocess/arithmetic.h>
#include <libprocess/stats.h>
#include <libprocess/inttrans.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"

static void thin_data_field(GwyDataField *data_field);

/**
 * gwy_data_field_normalize:
 * @data_field: A data field.
 *
 * Normalizes data in a data field to range 0.0 to 1.0.
 *
 * It is equivalent to gwy_data_field_renormalize(@data_field, 1.0, 0.0);
 *
 * If @data_field is filled with only one value, it is changed to 0.0.
 **/
void
gwy_data_field_normalize(GwyDataField *data_field)
{
    gdouble min, max;
    gdouble *p;
    gint xres, yres, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    gwy_data_field_get_min_max(data_field, &min, &max);
    if (min == max) {
        gwy_data_field_clear(data_field);
        return;
    }
    if (!min) {
        if (max != 1.0)
            gwy_data_field_multiply(data_field, 1.0/max);
        return;
    }

    /* The general case */
    max -= min;
    xres = data_field->xres;
    yres = data_field->yres;
    for (i = xres*yres, p = data_field->data; i; i--, p++)
        *p = (*p - min)/max;

    /* We can transform some stats */
    data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED) | CBIT(ARF) | CBIT(ART);
    CVAL(data_field, MIN) = 0.0;
    CVAL(data_field, MAX) = 1.0;
    CVAL(data_field, SUM) /= (CVAL(data_field, SUM) - xres*yres*min)/max;
    CVAL(data_field, RMS) /= max;
    CVAL(data_field, MED) = (CVAL(data_field, MED) - min)/max;
    CVAL(data_field, ART) = (CVAL(data_field, ART) - min)/max;
    CVAL(data_field, ARF) = (CVAL(data_field, ARF) - min)/max;
    /* There is a transformation formula for MSQ, but it can be prone to ugly cancellation errors. */
}

/**
 * gwy_data_field_renormalize:
 * @data_field: A data field.
 * @range: New data interval size.
 * @offset: New data interval offset.
 *
 * Transforms data in a data field with linear function to given range.
 *
 * When @range is positive, the new data range is (@offset, @offset+@range); when @range is negative, the new data
 * range is (@offset-@range, @offset). In neither case the data are flipped, negative range only means different
 * selection of boundaries.
 *
 * When @range is zero, this method is equivalent to gwy_data_field_fill(@data_field, @offset).
 **/
void
gwy_data_field_renormalize(GwyDataField *data_field,
                           gdouble range,
                           gdouble offset)
{
    gdouble min, max, v;
    gdouble *d;
    gint xres, yres, n, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    if (!range) {
        gwy_data_field_fill(data_field, offset);
        return;
    }

    gwy_data_field_get_min_max(data_field, &min, &max);
    if (min == max) {
        gwy_data_field_fill(data_field, offset);
        return;
    }

    if ((range > 0 && min == offset && min + range == max) || (range < 0 && max == offset && min - range == max))
        return;

    /* The general case */
    xres = data_field->xres;
    yres = data_field->yres;
    n = xres*yres;
    d = data_field->data;

    if (range > 0) {
        max -= min;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i) \
            shared(d,n,min,max,range,offset)
#endif
        for (i = 0; i < n; i++)
            d[i] = (d[i] - min)/max*range + offset;

        /* We can transform stats */
        data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED);
        CVAL(data_field, MIN) = offset;
        CVAL(data_field, MAX) = offset + range;
        v = CVAL(data_field, SUM);
        CVAL(data_field, SUM) = (v - n*min)/max*range + offset*n;
        CVAL(data_field, RMS) = CVAL(data_field, RMS)/max*range;
        CVAL(data_field, MED) = (CVAL(data_field, MED) - min)/max*range + offset;
        /* FIXME: we can recompute ARF and ART too */
        /* There is a transformation formula for MSQ, but it can be prone to ugly cancellation errors. */
    }
    else {
        min = max - min;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i) \
            shared(d,n,min,max,range,offset)
#endif
        for (i = 0; i < n; i++)
            d[i] = (max - d[i])/min*range + offset;

        /* We can transform stats */
        data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED);
        CVAL(data_field, MIN) = offset + range;
        CVAL(data_field, MAX) = offset;
        v = CVAL(data_field, SUM);
        CVAL(data_field, SUM) = (xres*yres*max - v)/min*range + offset*xres*yres;
        CVAL(data_field, RMS) = CVAL(data_field, RMS)/min*(-range);
        CVAL(data_field, MED) = (max - CVAL(data_field, MED))/min*range + offset;
        /* FIXME: we can recompute ARF and ART too */
        /* There is a transformation formula for MSQ, but it can be prone to ugly cancellation errors. */
    }
}

/**
 * gwy_data_field_area_renormalize:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @range: New data interval size.
 * @offset: New data interval offset.
 *
 * Transforms data in a part of a data field with linear function to given range.
 *
 * When @range is positive, the new data range is (@offset, @offset+@range); when @range is negative, the new data
 * range is (@offset-@range, @offset). In neither case the data are flipped, negative range only means different
 * selection of boundaries.
 *
 * When @range is zero, this method is equivalent to gwy_data_field_fill(@data_field, @offset).
 *
 * Since: 2.45
 **/
void
gwy_data_field_area_renormalize(GwyDataField *dfield,
                                gint col, gint row,
                                gint width, gint height,
                                gdouble range,
                                gdouble offset)
{
    gdouble min, max;
    gdouble *d, *r;
    gint xres, yres, i, j;

    if (!_gwy_data_field_check_area(dfield, col, row, width, height))
        return;

    xres = dfield->xres;
    yres = dfield->yres;
    if (col == 0 && row == 0 && width == xres && height == yres) {
        gwy_data_field_renormalize(dfield, range, offset);
        return;
    }

    if (!range) {
        gwy_data_field_area_fill(dfield, col, row, width, height, offset);
        return;
    }

    gwy_data_field_area_get_min_max_mask(dfield, NULL, GWY_MASK_IGNORE, col, row, width, height, &min, &max);
    if (min == max) {
        gwy_data_field_area_fill(dfield, col, row, width, height, offset);
        return;
    }

    if ((range > 0 && min == offset && min + range == max) || (range < 0 && max == offset && min - range == max))
        return;

    /* The general case */
    d = dfield->data;
    if (range > 0) {
        max -= min;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j,r) \
            shared(d,width,height,col,row,xres,min,max,range,offset)
#endif
        for (i = 0; i < height; i++) {
            r = d + (i + row)*xres + col;
            for (j = 0; j < width; j++)
                r[j] = (r[j] - min)/max*range + offset;
        }
    }
    else {
        min = max - min;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j,r) \
            shared(d,width,height,col,row,xres,min,max,range,offset)
#endif
        for (i = 0; i < height; i++) {
            r = d + (i + row)*xres + col;
            for (j = 0; j < width; j++)
                r[j] = (max - r[j])/min*range + offset;
        }
    }

    gwy_data_field_invalidate(dfield);
}

/**
 * gwy_data_field_threshold:
 * @data_field: A data field.
 * @threshval: Threshold value.
 * @bottom: Lower replacement value.
 * @top: Upper replacement value.
 *
 * Tresholds values of a data field.
 *
 * Values smaller than @threshold are set to value @bottom, values higher than @threshold or equal to it are set to
 * value @top
 *
 * Returns: The total number of values above threshold.
 **/
gint
gwy_data_field_threshold(GwyDataField *data_field,
                         gdouble threshval, gdouble bottom, gdouble top)
{
    gint i, n, tot = 0;
    gdouble *d;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0);

    n = data_field->xres * data_field->yres;
    d = data_field->data;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:tot) \
            private(i) \
            shared(d,n,threshval,bottom,top)
#endif
    for (i = 0; i < n; i++) {
        if (d[i] < threshval)
            d[i] = bottom;
        else {
            d[i] = top;
            tot++;
        }
    }

    /* We can precompute stats */
    data_field->cached = CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED);
    if (tot == n)
        CVAL(data_field, MIN) = CVAL(data_field, MAX) = top;
    else if (tot == 0)
        CVAL(data_field, MIN) = CVAL(data_field, MAX) = bottom;
    else {
        CVAL(data_field, MIN) = MIN(top, bottom);
        CVAL(data_field, MAX) = MAX(top, bottom);
    }
    CVAL(data_field, SUM) = tot*top + (n - tot)*bottom;
    CVAL(data_field, RMS) = (top - bottom)*(top - bottom) * tot/(gdouble)n * (n - tot)/(gdouble)n;
    /* FIXME: may be incorrect for tot == n/2(?) */
    CVAL(data_field, MED) = tot > n/2 ? top : bottom;

    return tot;
}


/**
 * gwy_data_field_area_threshold:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @threshval: Threshold value.
 * @bottom: Lower replacement value.
 * @top: Upper replacement value.
 *
 * Tresholds values of a rectangular part of a data field.
 *
 * Values smaller than @threshold are set to value @bottom, values higher than @threshold or equal to it are set to
 * value @top
 *
 * Returns: The total number of values above threshold.
 **/
gint
gwy_data_field_area_threshold(GwyDataField *data_field,
                              gint col, gint row, gint width, gint height,
                              gdouble threshval, gdouble bottom, gdouble top)
{
    gint xres, i, j, tot = 0;
    gdouble *d;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return 0;

    d = data_field->data;
    xres = data_field->xres;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:tot) \
            private(i,j) \
            shared(d,width,height,xres,col,row,threshval,bottom,top)
#endif
    for (i = 0; i < height; i++) {
        gdouble *drow = d + (row + i)*xres + col;

        for (j = 0; j < width; j++) {
            if (drow[j] < threshval)
                drow[j] = bottom;
            else {
                drow[j] = top;
                tot++;
            }
        }
    }
    gwy_data_field_invalidate(data_field);

    return tot;
}

/**
 * gwy_data_field_clamp:
 * @data_field: A data field.
 * @bottom: Lower limit value.
 * @top: Upper limit value.
 *
 * Limits data field values to a range.
 *
 * Returns: The number of changed values, i.e., values that were outside [@bottom, @top].
 **/
gint
gwy_data_field_clamp(GwyDataField *data_field,
                     gdouble bottom, gdouble top)
{
    gint i, n, tot = 0;
    gdouble *d = data_field->data;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0);
    g_return_val_if_fail(bottom <= top, 0);

    n = data_field->xres * data_field->yres;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:tot) \
            private(i) \
            shared(d,n,bottom,top)
#endif
    for (i = 0; i < n; i++) {
        if (d[i] < bottom) {
            d[i] = bottom;
            tot++;
        }
        else if (d[i] > top) {
            d[i] = top;
            tot++;
        }
    }
    if (tot) {
        /* We can precompute stats */
        data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(MED);
        CVAL(data_field, MIN) = MAX(bottom, CVAL(data_field, MIN));
        CVAL(data_field, MAX) = MIN(top, CVAL(data_field, MAX));
        if (CTEST(data_field, MED) && (CVAL(data_field, MED) < bottom || CVAL(data_field, MED) > top))
            data_field->cached &= ~CBIT(MED);
    }

    return tot;
}

/**
 * gwy_data_field_area_clamp:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @bottom: Lower limit value.
 * @top: Upper limit value.
 *
 * Limits values in a rectangular part of a data field to a range.
 *
 * Returns: The number of changed values, i.e., values that were outside [@bottom, @top].
 **/
gint
gwy_data_field_area_clamp(GwyDataField *data_field,
                          gint col, gint row,
                          gint width, gint height,
                          gdouble bottom, gdouble top)
{
    gint xres, i, j, tot = 0;
    gdouble *d;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return 0;

    d = data_field->data;
    xres = data_field->xres;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:tot) \
            private(i,j) \
            shared(d,xres,width,height,col,row,bottom,top)
#endif
    for (i = 0; i < height; i++) {
        gdouble *drow = d + (row + i)*xres + col;

        for (j = 0; j < width; j++) {
            if (drow[j] < bottom) {
                drow[j] = bottom;
                tot++;
            }
            else if (drow[j] > top) {
                drow[j] = top;
                tot++;
            }
        }
    }
    if (tot)
        gwy_data_field_invalidate(data_field);

    return tot;
}

/**
 * gwy_data_field_area_gather:
 * @data_field: A data field.
 * @result: A data field to put the result to, it may be @data_field itself.
 * @buffer: A data field to use as a scratch area, its size must be at least @width*@height.  May be %NULL to allocate
 *          a private temporary buffer.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @hsize: Horizontal size of gathered area.  The area is centered around each sample if @hsize is odd, it extends one
 *         pixel more to the right if @hsize is even.
 * @vsize: Vertical size of gathered area.  The area is centered around each sample if @vsize is odd, it extends one
 *         pixel more down if @vsize is even.
 * @average: %TRUE to divide resulting sums by the number of involved samples to get averages instead of sums.
 *
 * Sums or averages values in reactangular areas around each sample in a data field.
 *
 * When the gathered area extends out of calculation area, only samples from their intersection are taken into the
 * local sum (or average).
 *
 * There are no restrictions on values of @hsize and @vsize with regard to @width and @height, but they have to be
 * positive.
 *
 * The result is calculated by the means of two-dimensional rolling sums. One one hand it means the calculation time
 * depends linearly on (@width + @hsize)*(@height + @vsize) instead of @width*@hsize*@height*@vsize.  On the other
 * hand it means absolute rounding errors of all output values are given by the largest input values, that is relative
 * precision of results small in absolute value may be poor.
 **/
void
gwy_data_field_area_gather(GwyDataField *data_field,
                           GwyDataField *result,
                           GwyDataField *buffer,
                           gint hsize,
                           gint vsize,
                           gboolean average,
                           gint col, gint row,
                           gint width, gint height)
{
    const gdouble *srow, *trow;
    gdouble *drow, *r;
    gint xres, yres, i, j, m;
    gint hs2p, hs2m, vs2p, vs2m;
    gdouble v;

    g_return_if_fail(hsize > 0 && vsize > 0);
    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    xres = data_field->xres;
    yres = data_field->yres;
    g_return_if_fail(GWY_IS_DATA_FIELD(result));
    g_return_if_fail(result->xres == xres && result->yres == yres);
    if (buffer) {
        g_return_if_fail(GWY_IS_DATA_FIELD(buffer));
        g_return_if_fail(buffer->xres*buffer->yres >= width*height);
        g_object_ref(buffer);
    }
    else
        buffer = gwy_data_field_new(width, height, 1.0, 1.0, FALSE);

    /* Extension to the left and to the right (for asymmetric sizes extend to the right more) */
    hs2m = (hsize - 1)/2;
    hs2p = hsize/2;
    vs2m = (vsize - 1)/2;
    vs2p = vsize/2;

    /* Row-wise sums */
    /* FIXME: This is inefficient, split the inner loops to explicitly according to the conditions inside */
    for (i = 0; i < height; i++) {
        srow = data_field->data + (i + row)*xres + col;
        drow = buffer->data + i*width;

        /* Left half */
        drow[0] = 0.0;
        m = MIN(hs2p, width-1);
        for (j = 0; j <= m; j++)
            drow[0] += srow[j];
        for (j = 1; j < width/2; j++) {
            v = ((j + hs2p < width ? srow[j + hs2p] : 0.0) - (j-1 - hs2m >= 0 ? srow[j-1 - hs2m] : 0.0));
            drow[j] = drow[j-1] + v;
        }

        /* Right half */
        drow[width-1] = 0.0;
        m = width-1 - MIN(hs2m, width-1);
        for (j = width-1; j >= m; j--)
            drow[width-1] += srow[j];
        for (j = width-2; j >= width/2; j--) {
            v = ((j - hs2m >= 0 ? srow[j - hs2m] : 0.0) - (j+1 + hs2p < width ? srow[j+1 + hs2p] : 0.0));
            drow[j] = drow[j+1] + v;
        }
    }

    /* Column-wise sums (but iterate row-wise to access memory linearly) */
    /* Top half */
    r = result->data;
    drow = r + row*xres + col;
    for (j = 0; j < width; j++)
        drow[j] = 0.0;
    m = MIN(vs2p, height-1);
    for (i = 0; i <= m; i++) {
        srow = buffer->data + i*width;
        for (j = 0; j < width; j++)
            drow[j] += srow[j];
    }
    for (i = 1; i < height/2; i++) {
        drow = r + (i + row)*xres + col;
        if (i + vs2p < height) {
            srow = buffer->data + (i + vs2p)*width;
            if (i-1 - vs2m >= 0) {
                trow = buffer->data + (i-1 - vs2m)*width;
                for (j = 0; j < width; j++)
                    drow[j] = *(drow + j - xres) + (srow[j] - trow[j]);
            }
            else {
                for (j = 0; j < width; j++)
                    drow[j] = *(drow + j - xres) + srow[j];
            }
        }
        else {
            if (G_UNLIKELY(i-1 - vs2m >= 0)) {
                g_warning("Me thinks pure subtraction cannot occur.");
                trow = buffer->data + (i-1 - vs2m)*width;
                for (j = 0; j < width; j++)
                    drow[j] = *(drow + j - xres) - trow[j];
            }
            else {
                for (j = 0; j < width; j++)
                    drow[j] = *(drow + j - xres);
            }
        }
    }

    /* Bottom half */
    drow = r + (height-1 + row)*xres + col;
    for (j = 0; j < width; j++)
        drow[j] = 0.0;
    m = height-1 - MIN(vs2m, height-1);
    for (i = height-1; i >= m; i--) {
        srow = buffer->data + i*width;
        for (j = 0; j < width; j++)
            drow[j] += srow[j];
    }
    for (i = height-2; i >= height/2; i--) {
        drow = r + (i + row)*xres + col;
        if (i+1 + vs2p < height) {
            srow = buffer->data + (i+1 + vs2p)*width;
            if (G_LIKELY(i - vs2m >= 0)) {
                trow = buffer->data + (i - vs2m)*width;
                for (j = 0; j < width; j++)
                    drow[j] = drow[j + xres] + (trow[j] - srow[j]);
            }
            else {
                g_warning("Me thinks pure subtraction cannot occur.");
                for (j = 0; j < width; j++)
                    drow[j] = drow[j + xres] - srow[j];
            }
        }
        else {
            if (i - vs2m >= 0) {
                trow = buffer->data + (i - vs2m)*width;
                for (j = 0; j < width; j++)
                    drow[j] = drow[j + xres] + trow[j];
            }
            else {
                for (j = 0; j < width; j++)
                    drow[j] = drow[j + xres];
            }
        }
    }

    gwy_data_field_invalidate(result);
    g_object_unref(buffer);

    if (!average)
        return;

    /* Divide sums by the numbers of pixels that entered them */
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(r,width,height,xres,row,col,vsize,hsize,vs2m,vs2p,hs2m,hs2p)
#endif
    for (i = 0; i < height; i++) {
        gint iw;

        if (i <= vs2m)
            iw = vs2p + 1 + i;
        else if (i >= height-1 - vs2p)
            iw = vs2m + height - i;
        else
            iw = vsize;
        iw = MIN(iw, height);

        for (j = 0; j < width; j++) {
            gint jw;

            if (j <= hs2m)
                jw = hs2p + 1 + j;
            else if (j >= width-1 - hs2p)
                jw = hs2m + width - j;
            else
                jw = hsize;
            jw = MIN(jw, width);

            r[(i + row)*xres + j + col] /= iw*jw;
        }
    }
}

/**
 * gwy_data_field_area_filter_mean:
 * @data_field: A data field to apply the filter to.
 * @size: Averaged area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with mean filter of size @size.
 *
 * This method is a simple gwy_data_field_area_gather() wrapper, so the kernel is square.  Use convolution
 * gwy_data_field_area_ext_convolve() to perform a mean filter with different, for instance circular, kernel.
 **/
void
gwy_data_field_area_filter_mean(GwyDataField *data_field,
                                gint size,
                                gint col, gint row,
                                gint width, gint height)
{
    gwy_data_field_area_gather(data_field, data_field, NULL, size, size, TRUE, col, row, width, height);
}

/**
 * gwy_data_field_filter_mean:
 * @data_field: A data field to apply the filter to.
 * @size: Averaged area size.
 *
 * Filters a data field with mean filter of size @size.
 *
 * This method is a simple gwy_data_field_area_gather() wrapper, so the kernel is square.  Use convolution
 * gwy_data_field_area_ext_convolve() to perform a mean filter with different, for instance circular, kernel.
 **/
void
gwy_data_field_filter_mean(GwyDataField *data_field,
                           gint size)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_mean(data_field, size, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_area_filter_rms:
 * @data_field: A data field to apply RMS filter to.
 * @size: Area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with RMS filter of size @size.
 *
 * RMS filter computes root mean square in given area.
 **/
void
gwy_data_field_area_filter_rms(GwyDataField *data_field,
                               gint size,
                               gint col, gint row,
                               gint width, gint height)
{
    GwyDataField *avg2, *buffer;
    gint i, j;
    const gdouble *arow;
    gdouble *drow;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(size > 0);

    if (size == 1) {
        gwy_data_field_clear(data_field);
        return;
    }

    avg2 = gwy_data_field_area_extract(data_field, col, row, width, height);
    for (i = 0; i < width*height; i++)
        avg2->data[i] *= avg2->data[i];

    buffer = gwy_data_field_new_alike(avg2, FALSE);
    gwy_data_field_area_gather(avg2, avg2, buffer, size, size, TRUE, 0, 0, width, height);
    gwy_data_field_area_gather(data_field, data_field, buffer, size, size, TRUE, col, row, width, height);
    g_object_unref(buffer);

    for (i = 0; i < height; i++) {
        arow = avg2->data + i*width;
        drow = data_field->data + (i + row)*data_field->xres + col;
        for (j = 0; j < width; j++) {
            drow[j] = arow[j] - drow[j]*drow[j];
            drow[j] = sqrt(MAX(drow[j], 0.0));
        }
    }
    g_object_unref(avg2);

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_filter_rms:
 * @data_field: A data field to apply RMS filter to.
 * @size: Area size.
 *
 * Filters a data field with RMS filter.
 **/
void
gwy_data_field_filter_rms(GwyDataField *data_field,
                          gint size)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_rms(data_field, size, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_filter_canny:
 * @data_field: A data field to apply the filter to.
 * @threshold: Slope detection threshold (range 0..1).
 *
 * Filters a rectangular part of a data field with canny edge detector filter.
 **/
void
gwy_data_field_filter_canny(GwyDataField *data_field,
                            gdouble threshold)
{
    GwyDataField *sobel_horizontal;
    GwyDataField *sobel_vertical;
    gint i, j, k;
    gint xres, yres, n;
    gdouble angle, min, max;
    gint pass;
    gdouble *data, *svdata, *shdata;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    sobel_horizontal = gwy_data_field_duplicate(data_field);
    sobel_vertical = gwy_data_field_duplicate(data_field);

    gwy_data_field_filter_sobel(sobel_horizontal, GWY_ORIENTATION_HORIZONTAL);
    gwy_data_field_filter_sobel(sobel_vertical, GWY_ORIENTATION_VERTICAL);

    data = data_field->data;
    xres = data_field->xres;
    yres = data_field->yres;
    n = xres*yres;
    svdata = sobel_vertical->data;
    shdata = sobel_horizontal->data;
    for (k = 0; k < n; k++)
        data[k] = fabs(shdata[k]) + fabs(svdata[k]);

    gwy_data_field_invalidate(data_field);

    gwy_data_field_get_min_max(data_field, &min, &max);
    threshold = min + (max - min)*threshold;

    /* We do not need sobel array more, so use sobel_horizontal to store data
     * results. */
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j,pass,angle) \
            shared(data,shdata,svdata,xres,yres,threshold)
#endif
    for (i = 1; i < yres-1; i++) {
        for (j = 1; j < xres-1; j++) {
            pass = 0;
            if (data[i*xres + j] > threshold) {
                angle = atan2(svdata[i*xres + j], shdata[i*xres + j]);

                if (angle < 0.3925 || angle > 5.8875 || (angle > 2.7475 && angle < 3.5325)) {
                    if (data[j + 1 + xres*i] > threshold)
                        pass = 1;
                }
                else if ((angle > 1.178 && angle < 1.9632) || (angle > 4.318 && angle < 5.1049)) {
                    if (data[j + 1 + xres*(i + 1)] > threshold)
                        pass = 1;
                }
                else {
                    if (data[j + xres*(i + 1)] > threshold)
                        pass = 1;
                }
            }
            shdata[i*xres + j] = pass;
        }
    }
    /*result is now in sobel_horizontal field*/
    gwy_data_field_copy(sobel_horizontal, data_field, FALSE);

    g_object_unref(sobel_horizontal);
    g_object_unref(sobel_vertical);

    /*thin the lines*/
    thin_data_field(data_field);
}

/**
 * gwy_data_field_filter_slope:
 * @data_field: A data field to apply the filter to.
 * @xder: Data field where the x-derivative is to be stored, or %NULL if you are only interested in the y-derivative.
 *        It will be resized to match @data_field.
 * @yder: Data field where the y-derivative is to be stored, or %NULL if you are only interested in the x-derivative.
 *        It will be resized to match @data_field.
 *
 * Calculates x and y derivaties for an entire field.
 *
 * The derivatives are in physical units (not pixel-wise) and calculated from simple symmetrical differences, except at
 * the edges where the differences are one-sided.
 *
 * Since: 2.37
 **/
void
gwy_data_field_filter_slope(GwyDataField *data_field,
                            GwyDataField *xder,
                            GwyDataField *yder)
{
    guint xres, yres, i, j;
    gdouble dx, dy;
    const gdouble *d;
    gdouble *bx, *by;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(!xder || GWY_IS_DATA_FIELD(xder));
    g_return_if_fail(!yder || GWY_IS_DATA_FIELD(yder));
    if (!xder && !yder)
        return;

    xres = data_field->xres;
    yres = data_field->yres;
    if (xder) {
        gwy_data_field_resample(xder, xres, yres, GWY_INTERPOLATION_NONE);
        xder->xreal = data_field->xreal;
        xder->yreal = data_field->yreal;
        xder->xoff = data_field->xoff;
        xder->yoff = data_field->yoff;
        gwy_si_unit_assign(gwy_data_field_get_si_unit_xy(xder), gwy_data_field_get_si_unit_xy(data_field));
        gwy_si_unit_divide(gwy_data_field_get_si_unit_z(data_field), gwy_data_field_get_si_unit_xy(data_field),
                           gwy_data_field_get_si_unit_z(xder));
    }
    if (yder) {
        gwy_data_field_resample(yder, xres, yres, GWY_INTERPOLATION_NONE);
        yder->xreal = data_field->xreal;
        yder->yreal = data_field->yreal;
        yder->xoff = data_field->xoff;
        yder->yoff = data_field->yoff;
        gwy_si_unit_assign(gwy_data_field_get_si_unit_xy(xder), gwy_data_field_get_si_unit_xy(data_field));
        gwy_si_unit_divide(gwy_data_field_get_si_unit_z(data_field), gwy_data_field_get_si_unit_xy(data_field),
                           gwy_data_field_get_si_unit_z(xder));
    }
    dx = gwy_data_field_get_dx(data_field);
    dy = gwy_data_field_get_dy(data_field);
    d = data_field->data;
    bx = xder ? xder->data : NULL;
    by = yder ? yder->data : NULL;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(d,bx,by,dx,dy,xres,yres)
#endif
    for (i = 0; i < yres; i++) {
        const gdouble *row = d + i*xres, *prev = row - xres, *next = row + xres;
        gdouble *bxrow = (bx && xres > 1) ? bx + i*xres : NULL;
        gdouble *byrow = (by && yres > 1) ? by + i*xres : NULL;

        for (j = 0; j < xres; j++) {
            gdouble xd, yd;

            if (bxrow) {
                if (!j)
                    xd = row[j + 1] - row[j];
                else if (j == xres-1)
                    xd = row[j] - row[j - 1];
                else
                    xd = (row[j + 1] - row[j - 1])/2;

                bxrow[j] = xd/dx;
            }

            if (byrow) {
                if (!i)
                    yd = next[j] - row[j];
                else if (i == yres-1)
                    yd = row[j] - prev[j];
                else
                    yd = (next[j] - prev[j])/2;

                byrow[j] = yd/dy;
            }
        }
    }

    if (xder) {
        if (xres == 1)
            gwy_data_field_clear(xder);
        else
            gwy_data_field_invalidate(xder);
    }
    if (yder) {
        if (yres == 1)
            gwy_data_field_clear(yder);
        else
            gwy_data_field_invalidate(yder);
    }
}

/**
 * gwy_data_field_filter_gauss_step:
 * @data_field: A data field to apply the filter to.
 * @sigma: Gaussian filter width (in pixels).
 *
 * Processes a data field with Gaussian step detection filter.
 *
 * The filter is a multi-directional combination of convolutions with Gaussian multiplied by a signed step function.
 *
 * The resulting values correspond roughly to the step height around the pixel.
 *
 * Since: 2.54
 **/
void
gwy_data_field_filter_gauss_step(GwyDataField *data_field, gdouble sigma)
{
    GwyDataField *extfield, *kernel, *freal, *fimag, *kreal, *kimag, *gauss;
    guint ext, extx, exty, extxres, extyres, extn, xres, yres, i, n, idir, ndirs = 6;
    gdouble *fre, *fim, *g, *k, *kre, *kim, *f, *d;
    gdouble gs;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(sigma >= 0.0);
    if (sigma < 1e-9) {
        gwy_data_field_clear(data_field);
        return;
    }

    ext = GWY_ROUND(3.0*sigma) + 1;
    xres = data_field->xres;
    yres = data_field->yres;
    n = xres*yres;
    extx = gwy_fft_find_nice_size(xres + ext) - xres;
    exty = gwy_fft_find_nice_size(yres + ext) - yres;
    extfield = gwy_data_field_extend(data_field, extx/2, extx - extx/2, exty/2, exty - exty/2,
                                     GWY_EXTERIOR_BORDER_EXTEND, 0.0, FALSE);

    extxres = extfield->xres;
    extyres = extfield->yres;

    extn = extxres*extyres;
    freal = gwy_data_field_new_alike(extfield, FALSE);
    fimag = gwy_data_field_new_alike(extfield, FALSE);

    kernel = gwy_data_field_new_alike(extfield, FALSE);
    k = kernel->data;
    kreal = gwy_data_field_new_alike(extfield, FALSE);
    kimag = gwy_data_field_new_alike(extfield, FALSE);
    kre = kreal->data;
    kim = kimag->data;

    /* Subtract background on scale much larger than sigma. */
    gauss = gwy_data_field_new_alike(extfield, FALSE);
    g = gwy_data_field_get_data(gauss);
    gs = 0.0;
#ifdef _OPENMP
#pragma omp parallel for if (gwy_threads_are_enabled()) default(none) \
            reduction(+:gs) \
            shared(extxres,extyres,g,sigma) \
            private(i)
#endif
    for (i = 0; i < extyres; i++) {
        guint j, m = i*extxres;
        gdouble y = (i < (extyres + 1)/2 ? i : extyres-i)/(10.0*sigma);
        for (j = 0; j < extxres; j++, m++) {
            gdouble x = (j < (extxres + 1)/2 ? j : extxres-j)/(10.0*sigma);

            g[m] = exp(-x*x-y*y);
            gs += g[m];
        }
    }
    gwy_data_field_multiply(gauss, 1.0/gs);

    gwy_data_field_copy(extfield, kernel, FALSE);
    /* FIXME: We could do this in frequency space and directly modify freal
     * and fimag. */
    gwy_data_field_fft_convolve(kernel, gauss);
    gwy_data_field_subtract_fields(extfield, extfield, kernel);

    gwy_data_field_2dfft_raw(extfield, NULL, freal, fimag, GWY_TRANSFORM_DIRECTION_FORWARD);
    fre = freal->data;
    fim = fimag->data;

    /* Prepare a Gaussian. */
#ifdef _OPENMP
#pragma omp parallel for if (gwy_threads_are_enabled()) default(none) \
            shared(extxres,extyres,g,sigma) \
            private(i)
#endif
    for (i = 0; i < extyres; i++) {
        guint j, m = i*extxres;
        gdouble y = (i < (extyres + 1)/2 ? i : extyres-i)/sigma;
        for (j = 0; j < extxres; j++, m++) {
            gdouble x = (j < (extxres + 1)/2 ? j : extxres-j)/sigma;

            g[m] = exp(-x*x-y*y);
        }
    }
    g[0] = 0.0;
    gs = gwy_data_field_get_sum(gauss);

    /* Convolve with the Gaussian multiplied by rotated signum function. */
    gwy_data_field_clear(extfield);
    f = extfield->data;
    for (idir = 0; idir < ndirs; idir++) {
        gdouble phi = (idir + 0.5)/ndirs*G_PI;
        gdouble cphi = cos(phi), sphi = sin(phi);

#ifdef _OPENMP
#pragma omp parallel for if (gwy_threads_are_enabled()) default(none) \
            shared(extxres,extyres,g,k,sigma,cphi,sphi) \
            private(i)
#endif
        for (i = 0; i < extyres; i++) {
            gdouble y = (i < (extyres + 1)/2 ? i/sigma : (extyres-i)/(-sigma));
            guint j, m = i*extxres;
            for (j = 0; j < extxres; j++, m++) {
                gdouble x = (j < (extxres + 1)/2 ? j/sigma : (extxres-j)/(-sigma));
                gdouble v = x*cphi + y*sphi;

                if (G_UNLIKELY(fabs(v) < 1e-9))
                    k[m] = 0.0;
                else
                    k[m] = g[m]*((v > 0.0) - 0.5);
            }
        }

        gwy_data_field_2dfft_raw(kernel, NULL, kreal, kimag, GWY_TRANSFORM_DIRECTION_FORWARD);

#ifdef _OPENMP
#pragma omp parallel for if (gwy_threads_are_enabled()) default(none) \
            shared(extn,kre,kim,fre,fim) \
            private(i)
#endif
        for (i = 0; i < extn; i++) {
            gdouble re = kre[i]*fre[i] - kim[i]*fim[i];
            gdouble im = kre[i]*fim[i] + kim[i]*fre[i];
            kre[i] = re;
            kim[i] = im;
        }
        gwy_data_field_2dfft_raw(kreal, kimag, kernel, NULL, GWY_TRANSFORM_DIRECTION_BACKWARD);
        for (i = 0; i < extn; i++)
            f[i] += k[i]*k[i];
    }

    gwy_data_field_area_copy(extfield, data_field, extx/2, exty/2, xres, yres, 0, 0);

    d = data_field->data;
#ifdef _OPENMP
#pragma omp parallel for if (gwy_threads_are_enabled()) default(none) \
            shared(n,d,gs) \
            private(i)
#endif
    for (i = 0; i < n; i++)
        d[i] = sqrt(d[i]*n*9.0)/gs;
    gwy_data_field_invalidate(data_field);

    g_object_unref(extfield);
    g_object_unref(kreal);
    g_object_unref(kimag);
    g_object_unref(freal);
    g_object_unref(fimag);
    g_object_unref(kernel);
    g_object_unref(gauss);
}

/**
 * gwy_data_field_area_filter_conservative:
 * @data_field: A data field to apply the filter to.
 * @size: Filtered area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with conservative denoise filter.
 **/
void
gwy_data_field_area_filter_conservative(GwyDataField *data_field,
                                        gint size,
                                        gint col, gint row,
                                        gint width, gint height)
{
    gint xres, yres, i, j, ii, jj;
    gdouble *data, *hlpdata;
    GwyDataField *hlp_df;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(size > 0);
    xres = data_field->xres;
    yres = data_field->yres;
    if (size == 1)
        return;
    if (size > width || size > height) {
        g_warning("Kernel size larger than field area size.");
        return;
    }

    hlp_df = gwy_data_field_new(width, height, 1.0, 1.0, FALSE);

    data = data_field->data;
    hlpdata = hlp_df->data;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j,ii,jj) \
            shared(data,hlpdata,xres,yres,width,height,col,row,size)
#endif
    for (i = 0; i < height; i++) {
        gint ifrom = MAX(0, i + row - (size-1)/2);
        gint ito = MIN(yres-1, i + row + size/2);

        for (j = 0; j < width; j++) {
            gint jfrom = MAX(0, j + col - (size-1)/2);
            gint jto = MIN(xres-1, j + col + size/2);
            gdouble maxval = -G_MAXDOUBLE, minval = G_MAXDOUBLE;

            for (ii = 0; ii <= ito - ifrom; ii++) {
                gdouble *drow = data + (ifrom + ii)*xres + jfrom;

                for (jj = 0; jj <= jto - jfrom; jj++) {
                    if (i + row == ii + ifrom && j + col == jj + jfrom)
                        continue;

                    if (drow[jj] < minval)
                        minval = drow[jj];
                    if (drow[jj] > maxval)
                        maxval = drow[jj];
                }
            }

            hlpdata[i*width + j] = CLAMP(data[(i + row)*xres + j + col], minval, maxval);
        }
    }
    /* fix bottom right corner for size == 2 */
    if (size == 2)
        hlpdata[height*width - 1] = data[(row + height-1)*xres + col + width-1];

    gwy_data_field_area_copy(hlp_df, data_field, 0, 0, width, height, col, row);
    g_object_unref(hlp_df);
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_filter_conservative:
 * @data_field: A data field to apply the filter to.
 * @size: Filtered area size.
 *
 * Filters a data field with conservative denoise filter.
 **/
void
gwy_data_field_filter_conservative(GwyDataField *data_field,
                                   gint size)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_conservative(data_field, size, 0, 0, data_field->xres, data_field->yres);
}

/* Return 8*number-of-neighbours + number-of-neighbour-segments.
 * Segments must be separated by unset segments, so it is zero for all set
 * neighbourhood. */
static guint
neighbour_segments(const gdouble *d, gint xres, gint k)
{
    static const guint neightable[0x100] = {
        0,  9,  9,  17, 9,  18, 17, 25, 9,  17, 18, 25, 18, 26, 26, 33,
        9,  18, 18, 26, 17, 26, 25, 33, 18, 26, 27, 34, 26, 34, 34, 41,
        9,  18, 18, 26, 18, 27, 26, 34, 17, 25, 26, 33, 26, 34, 34, 41,
        18, 27, 27, 35, 26, 35, 34, 42, 26, 34, 35, 42, 34, 42, 42, 49,
        9,  18, 18, 26, 18, 27, 26, 34, 18, 26, 27, 34, 27, 35, 35, 42,
        18, 27, 27, 35, 26, 35, 34, 42, 27, 35, 36, 43, 35, 43, 43, 50,
        17, 26, 26, 34, 26, 35, 34, 42, 25, 33, 34, 41, 34, 42, 42, 49,
        26, 35, 35, 43, 34, 43, 42, 50, 34, 42, 43, 50, 42, 50, 50, 57,
        9,  18, 18, 26, 18, 27, 26, 34, 18, 26, 27, 34, 27, 35, 35, 42,
        17, 26, 26, 34, 25, 34, 33, 41, 26, 34, 35, 42, 34, 42, 42, 49,
        18, 27, 27, 35, 27, 36, 35, 43, 26, 34, 35, 42, 35, 43, 43, 50,
        26, 35, 35, 43, 34, 43, 42, 50, 34, 42, 43, 50, 42, 50, 50, 57,
        17, 26, 26, 34, 26, 35, 34, 42, 26, 34, 35, 42, 35, 43, 43, 50,
        25, 34, 34, 42, 33, 42, 41, 49, 34, 42, 43, 50, 42, 50, 50, 57,
        25, 34, 34, 42, 34, 43, 42, 50, 33, 41, 42, 49, 42, 50, 50, 57,
        33, 42, 42, 50, 41, 50, 49, 57, 41, 49, 50, 57, 49, 57, 57, 64,
    };
    guint b = 0;

    if (d[k-xres-1] > 0.0)
        b |= 1;
    if (d[k-xres] > 0.0)
        b |= 2;
    if (d[k-xres+1] > 0.0)
        b |= 4;
    if (d[k-1] > 0.0)
        b |= 8;
    if (d[k+1] > 0.0)
        b |= 16;
    if (d[k+xres-1] > 0.0)
        b |= 32;
    if (d[k+xres] > 0.0)
        b |= 64;
    if (d[k+xres+1] > 0.0)
        b |= 128;

    return neightable[b];
}

static gboolean
pixel_thinnable(const gdouble *data, gint xres, gint k)
{
    guint neighval;

    neighval = neighbour_segments(data, xres, k);
    /* One contiguous neighbour segment. */
    if (neighval % 8 != 1)
        return FALSE;

    /* Two to six neighbours. */
    neighval /= 8;
    if (neighval < 2 || neighval > 6)
        return FALSE;

    /* We could pull the first parts into the table in neighbour_segments(), but I am currently too lazy. */
    if (data[k + 1] > 0.0 && data[k - 1] > 0.0 && data[k + xres] > 0.0
        && neighbour_segments(data, xres, k+xres) % 8 == 1)
        return FALSE;

    if (data[k + 1] > 0.0 && data[k - xres] > 0.0 && data[k + xres] > 0.0
        && neighbour_segments(data, xres, k+1) % 8 == 1)
        return FALSE;

    return TRUE;
}

static gint
thinstep(GwyDataField *data_field,
         GwyDataField *buffer)
{
    gint i, j, ch;
    gint xres = data_field->xres, yres = data_field->yres;
    gdouble *data = data_field->data;
    gdouble *bdata = buffer->data;

    gwy_data_field_clear(buffer);
    ch = 0;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:ch) \
            private(i,j) \
            shared(data,bdata,xres,yres)
#endif
    for (i = 2; i < yres-2; i++) {
        for (j = 2; j < xres-2; j++) {
            gint k = i*xres + j;
            if (data[k] > 0.0 && pixel_thinnable(data, xres, k)) {
                ch++;
                bdata[k] = 1.0;
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(data,bdata,xres,yres)
#endif
    for (i = 2; i < yres-1; i++) {
        for (j = 2; j < xres-1; j++) {
            gint k = i*xres + j;
            if (bdata[k] > 0.0)
                data[k] = 0.0;
        }
    }
    gwy_data_field_invalidate(data_field);

    return ch;
}

/* XXX: Could be an alternative public function to gwy_data_field_thin().  It is more aggressive and loses some
 * branches, but this can be a good thing.  But for that it would be nice to implement it using a pixel queue so that
 * we do not have to repeatedly scan the entire data field. */
static void
thin_data_field(GwyDataField *data_field)
{
    GwyDataField *buffer;
    gint xres, yres;

    xres = data_field->xres;
    yres = data_field->yres;
    gwy_data_field_area_clear(data_field, 0, 0, xres, 1);
    gwy_data_field_area_clear(data_field, 0, 0, 1, yres);
    gwy_data_field_area_clear(data_field, xres-1, 0, 1, yres);
    gwy_data_field_area_clear(data_field, 0, yres-1, xres, 1);

    buffer = gwy_data_field_new_alike(data_field, FALSE);
    while (thinstep(data_field, buffer))
        ;
    g_object_unref(buffer);

    gwy_data_field_invalidate(data_field);
}

/**
 * kuwahara_block:
 * @a: points to a 5x5 matrix (array of 25 doubles)
 *
 * Computes a new value of the center pixel according to the Kuwahara filter.
 *
 * Return: Filtered value.
 */
static gdouble
kuwahara_block(const gdouble *a)
{
   static const gint r1[] = { 0, 1, 2, 5, 6, 7, 10, 11, 12 };
   static const gint r2[] = { 2, 3, 4, 7, 8, 9, 12, 13, 14 };
   static const gint r3[] = { 12, 13, 14, 17, 18, 19, 22, 23, 24 };
   static const gint r4[] = { 10, 11, 12, 15, 16, 17, 20, 21, 22 };
   gdouble mean1 = 0.0, mean2 = 0.0, mean3 = 0.0, mean4 = 0.0;
   gdouble var1 = 0.0, var2 = 0.0, var3 = 0.0, var4 = 0.0;
   gint i;

   for (i = 0; i < 9; i++) {
       mean1 += a[r1[i]]/9.0;
       mean2 += a[r2[i]]/9.0;
       mean3 += a[r3[i]]/9.0;
       mean4 += a[r4[i]]/9.0;
       var1 += a[r1[i]]*a[r1[i]]/9.0;
       var2 += a[r2[i]]*a[r2[i]]/9.0;
       var3 += a[r3[i]]*a[r3[i]]/9.0;
       var4 += a[r4[i]]*a[r4[i]]/9.0;
   }

   var1 -= mean1 * mean1;
   var2 -= mean2 * mean2;
   var3 -= mean3 * mean3;
   var4 -= mean4 * mean4;

   if (var1 <= var2 && var1 <= var3 && var1 <= var4)
       return mean1;
   if (var2 <= var3 && var2 <= var4 && var2 <= var1)
       return mean2;
   if (var3 <= var4 && var3 <= var1 && var3 <= var2)
       return mean3;
   if (var4 <= var1 && var4 <= var2 && var4 <= var3)
       return mean4;
   return 0.0;
}

/**
 * gwy_data_field_area_filter_kuwahara:
 * @data_field: A data filed to apply Kuwahara filter to.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with a Kuwahara (edge-preserving smoothing) filter.
 **/
void
gwy_data_field_area_filter_kuwahara(GwyDataField *data_field,
                                    gint col, gint row,
                                    gint width, gint height)
{
    gint i, j, xres, yres;
    gdouble *buffer, *d;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;

    buffer = g_new(gdouble, width*height);
    xres = data_field->xres;
    yres = data_field->yres;
    d = data_field->data;

#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(d,buffer,xres,yres,width,height,row,col)
#endif
    {
        gdouble kernel[25];
        gint ifrom = gwy_omp_chunk_start(height);
        gint ito = gwy_omp_chunk_end(height);

        for (i = ifrom; i < ito; i++) {
            for (j = 0; j < width; j++) {
                gint x, y, ii, jj, ctr = 0;

                for (y = -2; y <= 2; y++) {
                    ii = CLAMP(row + i + y, 0, yres-1);
                    for (x = -2; x <= 2; x++) {
                        jj = CLAMP(col + j + x, 0, xres-1);
                        kernel[ctr++] = d[ii*xres + jj];
                    }
                }
                buffer[i*width + j] = kuwahara_block(kernel);
            }
        }
    }

    for (i = 0; i < height; i++)
        gwy_assign(d + xres*(row + i) + col, buffer + i*width, width);

    g_free(buffer);
}

/**
 * gwy_data_field_filter_kuwahara:
 * @data_field: A data field to apply Kuwahara filter to.
 *
 * Filters a data field with Kuwahara filter.
 **/
void
gwy_data_field_filter_kuwahara(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_kuwahara(data_field, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_shade:
 * @data_field: A data field.
 * @target_field: A data field to put the shade to.  It will be resized to
 *                match @data_field.
 * @theta: Shading angle (in radians, from north pole).
 * @phi: Shade orientation in xy plane (in radians, counterclockwise).
 *
 * Shades a data field.
 **/
void
gwy_data_field_shade(GwyDataField *data_field,
                     GwyDataField *target_field,
                     gdouble theta, gdouble phi)
{
    gint i, j, xres, yres;
    gdouble max, maxval, v, cphi, sphi;
    gdouble *tdata;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(target_field));

    xres = data_field->xres;
    yres = data_field->yres;
    gwy_data_field_resample(target_field, xres, yres, GWY_INTERPOLATION_NONE);
    tdata = target_field->data;

    max = -G_MAXDOUBLE;
    cphi = cos(phi);
    sphi = sin(phi);
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(max:max) \
            private(i,j,v) \
            shared(data_field,tdata,xres,yres,cphi,sphi)
#endif
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++) {
            v = _gwy_data_field_xder(data_field, j, i)*cphi + _gwy_data_field_yder(data_field, j, i)*sphi;
            tdata[j + xres*i] = -v;
            if (max < v)
                max = v;
        }
    }

    maxval = theta/max;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i) \
            shared(tdata,xres,yres,max,maxval)
#endif
    for (i = 0; i < xres*yres; i++)
        tdata[i] = max - fabs(maxval - tdata[i]);

    gwy_data_field_invalidate(target_field);
}

/**
 * gwy_data_field_filter_harris:
 * @x_gradient: Data field with pre-calculated horizontal derivative.
 * @y_gradient: Data field with pre-calculated vertical derivative.
 * @result: Data field for the result.
 * @neighbourhood: Neighbourhood size.
 * @alpha: Sensitivity paramter (the squared trace is multiplied by it).
 *
 * Applies Harris corner detection filter to a pair of gradient data fields.
 *
 * All passed data field must have the same size.
 **/
void
gwy_data_field_filter_harris(GwyDataField *x_gradient,
                             GwyDataField *y_gradient,
                             GwyDataField *result,
                             gint neighbourhood,
                             gdouble alpha)
{
    gdouble mult, sigma;
    gint yres, xres, i, j;
    GwyDataField *xx, *xy, *yy;
    gdouble *xxdata, *xydata, *yydata, *xg, *yg, *r;

    g_return_if_fail(GWY_IS_DATA_FIELD(x_gradient));
    g_return_if_fail(GWY_IS_DATA_FIELD(y_gradient));
    g_return_if_fail(GWY_IS_DATA_FIELD(result));
    yres = result->yres;
    xres = result->xres;
    g_return_if_fail(x_gradient->xres == xres);
    g_return_if_fail(x_gradient->yres == yres);
    g_return_if_fail(y_gradient->xres == xres);
    g_return_if_fail(y_gradient->yres == yres);
    g_return_if_fail(neighbourhood > 0);

    if (2*neighbourhood >= MIN(xres, yres)) {
        gwy_data_field_clear(result);
        return;
    }

    mult = (fabs(gwy_data_field_get_max(x_gradient) - gwy_data_field_get_min(x_gradient))
            + fabs(gwy_data_field_get_max(y_gradient) - gwy_data_field_get_min(y_gradient)));
    mult = 1.0/(mult*mult);

    xx = result;
    xy = gwy_data_field_new_alike(result, TRUE);
    yy = gwy_data_field_new_alike(result, TRUE);
    xxdata = xx->data;
    xydata = xy->data;
    yydata = yy->data;
    xg = x_gradient->data;
    yg = y_gradient->data;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(xg,yg,xxdata,xydata,yydata,xres,yres,neighbourhood,mult)
#endif
    for (i = neighbourhood; i < yres - neighbourhood; i++) {
         for (j = neighbourhood; j < xres - neighbourhood; j++) {
             gint k = i*xres + j;
             gdouble vx = xg[k], vy = yg[k];
             xxdata[k] = vx*vx*mult;
             xydata[k] = vx*vy*mult;
             yydata[k] = vy*vy*mult;
         }
    }

    sigma = neighbourhood/5.0;
    gwy_data_field_filter_gaussian(xx, sigma);
    gwy_data_field_filter_gaussian(xy, sigma);
    gwy_data_field_filter_gaussian(yy, sigma);

    r = result->data;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(r,xxdata,xydata,yydata,xres,yres,neighbourhood,alpha)
#endif
    for (i = neighbourhood; i < yres - neighbourhood; i++) {
         for (j = neighbourhood; j < xres - neighbourhood; j++) {
             gint k = i*xres + j;
             gdouble pxx = xxdata[k], pxy = xydata[k], pyy = yydata[k];
             gdouble det = pxx*pyy - pxy*pxy, trace = pxx + pyy;
             r[k] = det - alpha*trace*trace;
          }
    }

    gwy_data_field_area_clear(result, 0, 0, xres, neighbourhood);
    gwy_data_field_area_clear(result, 0, neighbourhood, neighbourhood, yres-2*neighbourhood);
    gwy_data_field_area_clear(result, xres-neighbourhood, neighbourhood, neighbourhood, yres-2*neighbourhood);
    gwy_data_field_area_clear(result, 0, yres-neighbourhood, xres, neighbourhood);

    /* xx = result */
    g_object_unref(xy);
    g_object_unref(yy);
}

/************************** Documentation ****************************/

/**
 * SECTION:filters
 * @title: filters
 * @short_description: Convolution and other 2D data filters
 *
 * Filters are point-wise operations, such as thresholding, or generally local operations producing a value based on
 * the data in the vicinity of each point: gradients, step detectors and convolutions.  Some simple common point-wise
 * operations, e.g. value inversion, are also found in base #GwyDataField methods.
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
