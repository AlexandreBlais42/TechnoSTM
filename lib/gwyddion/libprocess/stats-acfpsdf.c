/*
 *  $Id: stats-acfpsdf.c 25319 2023-05-09 12:00:36Z yeti-dn $
 *  Copyright (C) 2003-2022 David Necas (Yeti), Petr Klapetek.
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
#include <libprocess/simplefft.h>
#include <libprocess/datafield.h>
#include <libprocess/filters.h>
#include <libprocess/correct.h>
#include <libprocess/arithmetic.h>
#include <libprocess/level.h>
#include <libprocess/stats.h>
#include <libprocess/linestats.h>
#include <libprocess/inttrans.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"
#include "gwyfftw.h"

static void
gwy_data_field_area_func_fft(GwyDataField *data_field,
                             GwyDataLine *target_line,
                             GwyFFTAreaFunc func,
                             gint col, gint row,
                             gint width, gint height,
                             GwyOrientation orientation,
                             GwyInterpolationType interpolation,
                             gint nstats)
{
    GwyDataLine *din, *dout;
    fftw_plan plan;
    gdouble *in, *out, *drow, *dcol;
    gint i, j, xres, yres, res = 0;
    gdouble avg;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));
    xres = data_field->xres;
    yres = data_field->yres;
    g_return_if_fail(orientation == GWY_ORIENTATION_HORIZONTAL || orientation == GWY_ORIENTATION_VERTICAL);

    if (orientation == GWY_ORIENTATION_VERTICAL) {
        res = gwy_fft_find_nice_size(2*yres);
        gwy_data_line_resample(target_line, height, GWY_INTERPOLATION_NONE);
    }
    else {
        res = gwy_fft_find_nice_size(2*xres);
        gwy_data_line_resample(target_line, width, GWY_INTERPOLATION_NONE);
    }
    gwy_data_line_clear(target_line);
    gwy_data_line_set_offset(target_line, 0.0);

    din = gwy_data_line_new(res, 1.0, FALSE);
    dout = gwy_data_line_new(res, 1.0, FALSE);
    in = gwy_data_line_get_data(din);
    out = gwy_data_line_get_data(dout);
    plan = gwy_fftw_plan_r2r_1d(res, in, out, FFTW_R2HC, FFTW_ESTIMATE);
    g_return_if_fail(plan);

    if (orientation == GWY_ORIENTATION_VERTICAL) {
        for (i = 0; i < width; i++) {
            dcol = data_field->data + row*xres + (i + col);
            avg = gwy_data_field_area_get_avg(data_field, NULL, col+i, row, 1, height);
            for (j = 0; j < height; j++)
                in[j] = dcol[j*xres] - avg;
            func(plan, din, dout, target_line);
        }
        gwy_data_line_set_real(target_line, gwy_data_field_itor(data_field, height));
        gwy_data_line_multiply(target_line, 1.0/(res*width));
    }
    else {
        for (i = 0; i < height; i++) {
            drow = data_field->data + (i + row)*xres + col;
            avg = gwy_data_field_area_get_avg(data_field, NULL, col, row+i, width, 1);
            for (j = 0; j < width; j++)
                in[j] = drow[j] - avg;
            func(plan, din, dout, target_line);
        }
        gwy_data_line_set_real(target_line, gwy_data_field_jtor(data_field, width));
        gwy_data_line_multiply(target_line, 1.0/(res*height));
    }

    fftw_destroy_plan(plan);
    g_object_unref(din);
    g_object_unref(dout);

    if (nstats > 1)
        gwy_data_line_resample(target_line, nstats, interpolation);
}

/**
 * gwy_data_field_area_acf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @orientation: Orientation of lines (ACF is simply averaged over the other orientation).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, @width (@height) is used.
 *
 * Calculates one-dimensional autocorrelation function of a rectangular part of a data field.
 **/
void
gwy_data_field_area_acf(GwyDataField *data_field,
                        GwyDataLine *target_line,
                        gint col, gint row,
                        gint width, gint height,
                        GwyOrientation orientation,
                        GwyInterpolationType interpolation,
                        gint nstats)
{
    gwy_data_field_area_func_fft(data_field, target_line, &do_fft_acf, col, row, width, height,
                                 orientation, interpolation, nstats);
    /* Set proper units */
    _gwy_copy_si_unit(data_field->si_unit_xy, &target_line->si_unit_x);
    gwy_si_unit_power(gwy_data_field_get_si_unit_z(data_field), 2, gwy_data_line_get_si_unit_y(target_line));
}

/**
 * gwy_data_field_acf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @orientation: Orientation of lines (ACF is simply averaged over the other orientation).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, data field width (height) is
 * used.
 *
 * Calculates one-dimensional autocorrelation function of a data field.
 **/
void
gwy_data_field_acf(GwyDataField *data_field,
                   GwyDataLine *target_line,
                   GwyOrientation orientation,
                   GwyInterpolationType interpolation,
                   gint nstats)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_acf(data_field, target_line, 0, 0, data_field->xres, data_field->yres,
                            orientation, interpolation, nstats);
}

/**
 * gwy_data_field_area_hhcf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @orientation: Orientation of lines (HHCF is simply averaged over the other orientation).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, @width (@height) is used.
 *
 * Calculates one-dimensional autocorrelation function of a rectangular part of a data field.
 **/
void
gwy_data_field_area_hhcf(GwyDataField *data_field,
                         GwyDataLine *target_line,
                         gint col, gint row,
                         gint width, gint height,
                         GwyOrientation orientation,
                         GwyInterpolationType interpolation,
                         gint nstats)
{
    gwy_data_field_area_func_fft(data_field, target_line, &do_fft_hhcf, col, row, width, height,
                                 orientation, interpolation, nstats);

    /* Set proper units */
    _gwy_copy_si_unit(data_field->si_unit_xy, &target_line->si_unit_x);
    gwy_si_unit_power(gwy_data_field_get_si_unit_z(data_field), 2, gwy_data_line_get_si_unit_y(target_line));
}

/**
 * gwy_data_field_hhcf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @orientation: Orientation of lines (HHCF is simply averaged over the other orientation).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, data field width (height) is
 *          used.
 *
 * Calculates one-dimensional autocorrelation function of a data field.
 **/
void
gwy_data_field_hhcf(GwyDataField *data_field,
                    GwyDataLine *target_line,
                    GwyOrientation orientation,
                    GwyInterpolationType interpolation,
                    gint nstats)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_hhcf(data_field, target_line, 0, 0, data_field->xres, data_field->yres,
                             orientation, interpolation, nstats);
}

/**
 * gwy_data_field_area_psdf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @orientation: Orientation of lines (PSDF is simply averaged over the other orientation).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @windowing: Windowing type to use.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, data field width (height) is
 *          used.
 *
 * Calculates one-dimensional power spectrum density function of a rectangular part of a data field.
 **/
void
gwy_data_field_area_psdf(GwyDataField *data_field,
                         GwyDataLine *target_line,
                         gint col, gint row,
                         gint width, gint height,
                         GwyOrientation orientation,
                         GwyInterpolationType interpolation,
                         GwyWindowingType windowing,
                         gint nstats)
{
    GwyDataField *re_field, *im_field;
    GwySIUnit *xyunit, *zunit, *lineunit;
    gdouble *re, *im, *target;
    gint i, j, xres, yres, size;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));
    xres = data_field->xres;
    yres = data_field->yres;
    size = (orientation == GWY_ORIENTATION_HORIZONTAL) ? width : height;
    g_return_if_fail(size >= 4);
    g_return_if_fail(orientation == GWY_ORIENTATION_HORIZONTAL || orientation == GWY_ORIENTATION_VERTICAL);

    if (nstats < 1)
        nstats = size/2 - 1;
    gwy_data_line_resample(target_line, size/2, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(target_line);
    gwy_data_line_set_offset(target_line, 0.0);

    re_field = gwy_data_field_new(width, height, 1.0, 1.0, FALSE);
    im_field = gwy_data_field_new(width, height, 1.0, 1.0, FALSE);
    target = target_line->data;
    if (orientation == GWY_ORIENTATION_VERTICAL) {
        gwy_data_field_area_1dfft(data_field, NULL, re_field, im_field, col, row, width, height,
                                  orientation, windowing, GWY_TRANSFORM_DIRECTION_FORWARD, interpolation, TRUE, 2);
        re = re_field->data;
        im = im_field->data;
        for (i = 0; i < width; i++) {
            for (j = 0; j < size/2; j++)
                target[j] += re[j*width + i]*re[j*width + i] + im[j*width + i]*im[j*width + i];
        }
        gwy_data_line_multiply(target_line, data_field->yreal/yres/(2*G_PI*width));
        gwy_data_line_set_real(target_line, G_PI*yres/data_field->yreal);
    }
    else {
        gwy_data_field_area_1dfft(data_field, NULL, re_field, im_field, col, row, width, height,
                                  orientation, windowing, GWY_TRANSFORM_DIRECTION_FORWARD, interpolation, TRUE, 2);
        re = re_field->data;
        im = im_field->data;
        for (i = 0; i < height; i++) {
            for (j = 0; j < size/2; j++)
                target[j] += re[i*width + j]*re[i*width + j] + im[i*width + j]*im[i*width + j];
        }
        gwy_data_line_multiply(target_line, data_field->xreal/xres/(2*G_PI*height));
        gwy_data_line_set_real(target_line, G_PI*xres/data_field->xreal);
    }

    gwy_data_line_set_offset(target_line, target_line->real/target_line->res);
    gwy_data_line_resize(target_line, 1, target_line->res);
    gwy_data_line_resample(target_line, nstats, interpolation);

    g_object_unref(re_field);
    g_object_unref(im_field);

    /* Set proper units */
    xyunit = gwy_data_field_get_si_unit_xy(data_field);
    zunit = gwy_data_field_get_si_unit_z(data_field);
    lineunit = gwy_data_line_get_si_unit_x(target_line);
    gwy_si_unit_power(xyunit, -1, lineunit);
    lineunit = gwy_data_line_get_si_unit_y(target_line);
    gwy_si_unit_power(zunit, 2, lineunit);
    gwy_si_unit_multiply(lineunit, xyunit, lineunit);
}

/**
 * gwy_data_field_psdf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @orientation: Orientation of lines (PSDF is simply averaged over the other orientation).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @windowing: Windowing type to use.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, data field width (height) is
 *          used.
 *
 * Calculates one-dimensional power spectrum density function of a data field.
 **/
void
gwy_data_field_psdf(GwyDataField *data_field,
                    GwyDataLine *target_line,
                    GwyOrientation orientation,
                    GwyInterpolationType interpolation,
                    GwyWindowingType windowing,
                    gint nstats)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_psdf(data_field, target_line, 0, 0, data_field->xres, data_field->yres,
                             orientation, interpolation, windowing, nstats);
}

/**
 * gwy_data_field_area_rpsdf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @windowing: Windowing type to use.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, data field width (height) is
 *          used.
 *
 * Calculates radial power spectrum density function of a rectangular part of a data field.
 *
 * Since: 2.7
 **/
void
gwy_data_field_area_rpsdf(GwyDataField *data_field,
                          GwyDataLine *target_line,
                          gint col, gint row,
                          gint width, gint height,
                          GwyInterpolationType interpolation,
                          GwyWindowingType windowing,
                          gint nstats)
{
    GwyDataField *re_field, *im_field;
    GwySIUnit *xyunit, *zunit, *lineunit;
    gdouble *re, *im;
    gint i, j, k, xres, yres;
    gdouble xreal, yreal, r;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));
    g_return_if_fail(width >= 4 && height >= 4);
    xres = data_field->xres;
    yres = data_field->yres;
    xreal = data_field->xreal;
    yreal = data_field->yreal;

    re_field = gwy_data_field_new(width, height, width*xreal/xres, height*yreal/yres, FALSE);
    im_field = gwy_data_field_new_alike(re_field, FALSE);
    gwy_data_field_area_2dfft(data_field, NULL, re_field, im_field, col, row, width, height,
                              windowing, GWY_TRANSFORM_DIRECTION_FORWARD, interpolation, TRUE, 2);
    re = re_field->data;
    im = im_field->data;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            k = i*width + j;
            re[k] = re[k]*re[k] + im[k]*im[k];
        }
    }
    g_object_unref(im_field);

    gwy_data_field_fft_postprocess(re_field, TRUE);
    r = 0.5*MAX(re_field->xreal, re_field->yreal);
    gwy_data_field_angular_average(re_field, target_line, NULL, GWY_MASK_IGNORE, 0.0, 0.0, r, nstats ? nstats+1 : 0);
    g_object_unref(re_field);
    /* Get rid of the zero first element which is bad for logscale. */
    nstats = target_line->res-1;
    gwy_data_line_resize(target_line, 1, nstats+1);
    target_line->off += target_line->real/nstats;

    /* Postprocess does not use angular coordinates, fix that. */
    target_line->real *= 2.0*G_PI;
    target_line->off *= 2.0*G_PI;
    r = xreal*yreal/(2.0*G_PI*width*height) * target_line->real/nstats;
    for (k = 0; k < nstats; k++)
        target_line->data[k] *= r*(k + 1);

    /* Set proper value units */
    xyunit = gwy_data_field_get_si_unit_xy(data_field);
    zunit = gwy_data_field_get_si_unit_z(data_field);
    lineunit = gwy_data_line_get_si_unit_x(target_line);
    gwy_si_unit_power(xyunit, -1, lineunit);
    lineunit = gwy_data_line_get_si_unit_y(target_line);
    gwy_si_unit_power(zunit, 2, lineunit);
    gwy_si_unit_multiply(lineunit, xyunit, lineunit);
}

/**
 * gwy_data_field_rpsdf:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to requested width.
 * @interpolation: Interpolation to use when @nstats is given and requires resampling.
 * @windowing: Windowing type to use.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, data field width (height) is
 *          used.
 *
 * Calculates radial power spectrum density function of a data field.
 *
 * Since: 2.7
 **/
void
gwy_data_field_rpsdf(GwyDataField *data_field,
                     GwyDataLine *target_line,
                     GwyInterpolationType interpolation,
                     GwyWindowingType windowing,
                     gint nstats)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_rpsdf(data_field, target_line, 0, 0, data_field->xres, data_field->yres,
                              interpolation, windowing, nstats);
}

/**
 * gwy_data_field_area_racf:
 * @data_field: A data field.
 * @target_line: A data line to store the autocorrelation function to.  It will be resampled to requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @nstats: The number of samples to take on the autocorrelation function.  If nonpositive, a suitable resolution is
 *          chosen automatically.
 *
 * Calculates radially averaged autocorrelation function of a rectangular part of a data field.
 *
 * Since: 2.22
 **/
void
gwy_data_field_area_racf(GwyDataField *data_field,
                         GwyDataLine *target_line,
                         gint col, gint row,
                         gint width, gint height,
                         gint nstats)
{
    GwyDataField *acf_field;
    gint size;
    gdouble r;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));
    g_return_if_fail(width >= 4 && height >= 4);

    size = MIN(width, height)/2;
    if (nstats < 1)
        nstats = size;

    acf_field = gwy_data_field_new(2*size - 1, 2*size - 1, 1.0, 1.0, FALSE);
    gwy_data_field_area_2dacf(data_field, acf_field, col, row, width, height, size, size);
    r = 0.5*MAX(acf_field->xreal, acf_field->yreal);
    gwy_data_field_angular_average(acf_field, target_line, NULL, GWY_MASK_IGNORE, 0.0, 0.0, r, nstats);
    g_object_unref(acf_field);
}

/**
 * gwy_data_field_racf:
 * @data_field: A data field.
 * @target_line: A data line to store the autocorrelation function to.  It will be resampled to requested width.
 * @nstats: The number of samples to take on the autocorrelation function.  If nonpositive, a suitable resolution is
 *          chosen automatically.
 *
 * Calculates radially averaged autocorrelation function of a data field.
 *
 * Since: 2.22
 **/
void
gwy_data_field_racf(GwyDataField *data_field,
                    GwyDataLine *target_line,
                    gint nstats)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_racf(data_field, target_line, 0, 0, data_field->xres, data_field->yres, nstats);
}

static void
fix_levelling_degree(gint *level)
{
    if (*level > 2) {
        g_warning("Levelling degree %u is not supported, changing to 2.", *level);
        *level = 2;
    }
}

static void
execute_2d_acf(GwyDataField *extfield,
               fftw_plan plan, fftw_complex *cbuf, guint cstride,
               guint width, guint height)
{
    guint i, j, xsize = extfield->xres, ysize = extfield->yres;
    gdouble *row1, *row2, *extdata = extfield->data;
    fftw_complex *c;
    gdouble re;

    gwy_data_field_area_clear(extfield, width, 0, xsize - width, height);
    gwy_data_field_area_clear(extfield, 0, height, xsize, ysize - height);
    gwy_fftw_execute(plan);

    c = cbuf;
    for (i = 0; i < ysize; i++) {
        row1 = extdata + i*xsize;
        row2 = extdata + ((ysize - i) % ysize)*xsize + xsize-1;
        re = gwycreal(*c)*gwycreal(*c) + gwycimag(*c)*gwycimag(*c);
        *(row1++) = re;
        c++;
        for (j = 1; j < cstride; j++) {
            re = gwycreal(*c)*gwycreal(*c) + gwycimag(*c)*gwycimag(*c);
            *(row1++) = re;
            *(row2--) = re;
            c++;
        }
    }
    gwy_fftw_execute(plan);
}

static void
extract_2d_acf_real(GwyDataField *field,
                    fftw_complex *cbuf, guint cstride, guint ysize)
{
    guint i, j, yrange, xrange, txres = field->xres, tyres = field->yres;
    gdouble *row1, *row2;
    fftw_complex *c;

    xrange = txres/2 + 1;
    yrange = tyres/2 + 1;
    for (i = 0; i < tyres; i++) {
        c = cbuf + ((i + ysize - (yrange-1)) % ysize)*cstride;
        row1 = field->data + i*txres;
        for (j = xrange-1; j < txres; j++) {
            row1[j] = gwycreal(*c);
            c++;
        }

        row1 = field->data + (tyres-1 - i)*txres;
        row2 = field->data + i*txres + txres-1;
        for (j = 0; j < xrange-1; j++)
            *(row1++) = *(row2--);
    }
}

/**
 * gwy_data_field_area_2dacf_mask:
 * @data_field: A data field.
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use (has any effect only with non-%NULL @mask).
 * @target_field: A data field to store the result to.  It will be resampled to (2@xrange-1)×(2@yrange-1).
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @xrange: Horizontal correlation range.  Non-positive value means the default range of half of @data_field width
 *          will be used.
 * @yrange: Vertical correlation range.  Non-positive value means the default range of half of @data_field height will
 *          be used.
 * @weights: Field to store the denominators to (or %NULL).  It will be resized like @target_field.  The denominators
 *           are integers equal to the number of terms that contributed to each value.  They are suitable as fitting
 *           weights if the ACF is fitted.
 *
 * Calculates two-dimensional autocorrelation function of a data field area.
 *
 * The resulting data field has the correlation corresponding to (0,0) in the centre.
 *
 * The maximum possible values of @xrange and @yrange are @data_field width and height, respectively.  However, as the
 * values for longer distances are calculated from smaller number of data points they become increasingly bogus,
 * therefore the default range is half of the size.
 *
 * Since: 2.50
 **/
void
gwy_data_field_area_2dacf_mask(GwyDataField *field,
                               GwyDataField *target,
                               GwyDataField *mask,
                               GwyMaskingType masking,
                               gint col, gint row,
                               gint width, gint height,
                               gint xrange, gint yrange,
                               GwyDataField *weights)
{
    GwyDataField *extfield;
    fftw_plan plan;
    gint i, j, xres, yres, xsize, ysize, cstride, txres, tyres, qi;
    gdouble xreal, yreal, thresh;
    fftw_complex *cbuf;
    gdouble *trow, *qrow, *mrow, *drow;

    if (!_gwy_data_field_check_area(field, col, row, width, height)
        || !_gwy_data_field_check_mask(field, &mask, &masking))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(target));
    g_return_if_fail(!weights || GWY_IS_DATA_FIELD(weights));
    xres = field->xres;
    yres = field->yres;
    if (xrange <= 0)
        xrange = width/2;
    if (yrange <= 0)
        yrange = height/2;
    g_return_if_fail(xrange <= width && yrange <= height);
    xreal = field->xreal;
    yreal = field->yreal;

    xsize = gwy_fft_find_nice_size(width + xrange);
    ysize = gwy_fft_find_nice_size(height + yrange);
    cstride = xsize/2 + 1;

    txres = 2*xrange - 1;
    tyres = 2*yrange - 1;
    gwy_data_field_resample(target, txres, tyres, GWY_INTERPOLATION_NONE);

    if (weights) {
        gwy_data_field_resample(weights, txres, tyres, GWY_INTERPOLATION_NONE);
        g_object_ref(weights);
    }
    else
        weights = gwy_data_field_new_alike(target, FALSE);

    cbuf = gwy_fftw_new_complex(cstride*ysize);
    extfield = gwy_data_field_new(xsize, ysize, 1.0, 1.0, FALSE);
    plan = gwy_fftw_plan_dft_r2c_2d(ysize, xsize, extfield->data, cbuf, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    if (mask) {
        /* Calculate unnormalised ACF of the mask, i.e. the denominators. */
        for (i = 0; i < height; i++) {
            mrow = mask->data + (i + row)*xres + col;
            trow = extfield->data + i*xsize;
            if (masking == GWY_MASK_INCLUDE) {
                for (j = 0; j < width; j++)
                    trow[j] = (mrow[j] > 0.0);
            }
            else {
                for (j = 0; j < width; j++)
                    trow[j] = (mrow[j] <= 0.0);
            }
        }
        execute_2d_acf(extfield, plan, cbuf, cstride, width, height);
        extract_2d_acf_real(weights, cbuf, cstride, ysize);
        gwy_data_field_multiply(weights, 1.0/(xsize*ysize));

        /* Calculate unnormalised ACF of the premultiplied image, i.e. the
         * numerators. */
        for (i = 0; i < height; i++) {
            mrow = mask->data + (i + row)*xres + col;
            drow = field->data + (i + row)*xres + col;
            trow = extfield->data + i*xsize;
            if (masking == GWY_MASK_INCLUDE) {
                for (j = 0; j < width; j++)
                    trow[j] = (mrow[j] > 0.0) * drow[j];
            }
            else {
                for (j = 0; j < width; j++)
                    trow[j] = (mrow[j] <= 0.0) * drow[j];
            }
        }
    }
    else {
        qrow = weights->data;
        for (j = 0; j < txres; j++)
            qrow[j] = width - ABS(j - (xrange-1));
        for (i = 1; i < tyres; i++) {
            qi = height - ABS(i - (yrange-1));
            trow = weights->data + i*txres;
            for (j = 0; j < txres; j++)
                trow[j] = qrow[j] * qi;
        }
        for (j = 0; j < txres; j++)
            qrow[j] *= height - (yrange - 1);

        gwy_data_field_area_copy(field, extfield, col, row, width, height, 0, 0);
    }
    execute_2d_acf(extfield, plan, cbuf, cstride, width, height);
    extract_2d_acf_real(target, cbuf, cstride, ysize);
    gwy_data_field_multiply(target, 1.0/(xsize*ysize));

    fftw_destroy_plan(plan);
    fftw_free(cbuf);

    if (mask) {
        gwy_data_field_resample(extfield, txres, tyres, GWY_INTERPOLATION_NONE);
        thresh = 1.000001*GWY_ROUND(log(width*height));
        for (i = 0; i < tyres; i++) {
            qrow = weights->data + i*txres;
            trow = target->data + i*txres;
            mrow = extfield->data + i*txres;
            for (j = 0; j < txres; j++) {
                if (qrow[j] > thresh) {
                    mrow[j] = 0.0;
                    trow[j] /= qrow[j];
                }
                else {
                    mrow[j] = 1.0;
                    trow[j] = 0.0;
                }
            }
        }
        gwy_data_field_laplace_solve(target, extfield, -1, 1.0);
    }
    else {
        gwy_data_field_divide_fields(target, target, weights);
    }
    g_object_unref(extfield);

    target->xreal = xreal*txres/xres;
    target->yreal = yreal*tyres/yres;
    target->xoff = -0.5*target->xreal;
    target->yoff = -0.5*target->yreal;

    _gwy_copy_si_unit(field->si_unit_xy, &target->si_unit_xy);
    gwy_si_unit_power(gwy_data_field_get_si_unit_z(field), 2, gwy_data_field_get_si_unit_z(target));

    weights->xreal = xreal*txres/xres;
    weights->yreal = yreal*tyres/yres;
    weights->xoff = -0.5*weights->xreal;
    weights->yoff = -0.5*weights->yreal;

    _gwy_copy_si_unit(field->si_unit_xy, &weights->si_unit_xy);
    _gwy_copy_si_unit(NULL, &weights->si_unit_z);

    gwy_data_field_invalidate(target);
    gwy_data_field_invalidate(weights);
    g_object_unref(weights);
}

/**
 * gwy_data_field_area_2dacf:
 * @data_field: A data field.
 * @target_field: A data field to store the result to.  It will be resampled to (2@xrange-1)×(2@yrange-1).
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @xrange: Horizontal correlation range.  Non-positive value means the default range of half of @data_field width
 *          will be used.
 * @yrange: Vertical correlation range.  Non-positive value means the default range of half of @data_field height will
 *          be used.
 *
 * Calculates two-dimensional autocorrelation function of a data field area.
 *
 * The resulting data field has the correlation corresponding to (0,0) in the centre.
 *
 * The maximum possible values of @xrange and @yrange are @data_field width and height, respectively.  However, as the
 * values for longer distances are calculated from smaller number of data points they become increasingly bogus,
 * therefore the default range is half of the size.
 *
 * Since: 2.7
 **/
void
gwy_data_field_area_2dacf(GwyDataField *field,
                          GwyDataField *target,
                          gint col, gint row,
                          gint width, gint height,
                          gint xrange, gint yrange)
{
    gwy_data_field_area_2dacf_mask(field, target, NULL, GWY_MASK_IGNORE, col, row, width, height, xrange, yrange,
                                   NULL);
}

/**
 * gwy_data_field_2dacf:
 * @data_field: A data field.
 * @target_field: A data field to store the result to.
 *
 * Calculates two-dimensional autocorrelation function of a data field.
 *
 * See gwy_data_field_area_2dacf() for details.  Parameters missing (not adjustable) in this function are set to their
 * default values.
 *
 * Since: 2.7
 **/
void
gwy_data_field_2dacf(GwyDataField *data_field,
                     GwyDataField *target_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    gwy_data_field_area_2dacf(data_field, target_field, 0, 0, data_field->xres, data_field->yres, 0, 0);
}

/* Assumes @plan acts on buf->data. */
static void
execute_2d_cacf(GwyDataField *buf, GwyDataField *target,
                fftw_plan plan, fftw_complex *cbuf, guint cstride)
{
    guint i, j, xres = buf->xres, yres = buf->yres;
    gdouble *row1, *row2, *data;
    fftw_complex *c;
    gdouble re;

    gwy_fftw_execute(plan);
    c = cbuf;
    data = buf->data;
    for (i = 0; i < yres; i++) {
        row1 = data + i*xres;
        row2 = data + ((yres - i) % yres)*xres + xres-1;
        re = gwycreal(*c)*gwycreal(*c) + gwycimag(*c)*gwycimag(*c);
        *(row1++) = re;
        c++;
        for (j = 1; j < cstride; j++) {
            re = gwycreal(*c)*gwycreal(*c) + gwycimag(*c)*gwycimag(*c);
            *(row1++) = re;
            *(row2--) = re;
            c++;
        }
    }

    gwy_fftw_execute(plan);
    c = cbuf;
    data = target->data;
    for (i = 0; i < yres; i++) {
        row1 = data + i*xres;
        row2 = data + ((yres - i) % yres)*xres + xres-1;
        re = gwycreal(*c);
        *(row1++) = re;
        c++;
        for (j = 1; j < cstride; j++) {
            re = gwycreal(*c);
            *(row1++) = re;
            *(row2--) = re;
            c++;
        }
    }
}

/* Extract real parts.  Callers might not like negative PSDF much but it is the unbiased estimate. */
static void
extract_2d_fft_real(GwyDataField *target,
                    fftw_plan plan, fftw_complex *cbuf, guint cstride)
{
    guint i, j, xres = target->xres, yres = target->yres;
    gdouble *row1, *row2, *data;
    fftw_complex *c;
    gdouble re;

    gwy_fftw_execute(plan);
    c = cbuf;
    data = target->data;
    for (i = 0; i < yres; i++) {
        row1 = data + i*xres;
        row2 = data + ((yres - i) % yres)*xres + xres-1;
        re = gwycreal(*c);
        *(row1++) = re;
        c++;
        for (j = 1; j < cstride; j++) {
            re = gwycreal(*c);
            *(row1++) = re;
            *(row2--) = re;
            c++;
        }
    }
}

/* Ensure we produce non-negative output even in the presence of rounding errors.  Callers might not like negative
 * PSDF much. */
static void
extract_2d_fft_cnorm(GwyDataField *target,
                     fftw_plan plan, fftw_complex *cbuf, guint cstride)
{
    guint i, j, xres = target->xres, yres = target->yres;
    gdouble *row1, *row2, *data;
    fftw_complex *c;
    gdouble re;

    gwy_fftw_execute(plan);
    c = cbuf;
    data = target->data;
    for (i = 0; i < yres; i++) {
        row1 = data + i*xres;
        row2 = data + ((yres - i) % yres)*xres + xres-1;
        re = gwycreal(*c)*gwycreal(*c) + gwycimag(*c)*gwycimag(*c);
        *(row1++) = re;
        c++;
        for (j = 1; j < cstride; j++) {
            re = gwycreal(*c)*gwycreal(*c) + gwycimag(*c)*gwycimag(*c);
            *(row1++) = re;
            *(row2--) = re;
            c++;
        }
    }
}

/**
 * gwy_data_field_area_2dpsdf_mask:
 * @field: A data field.
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use (has any effect only with non-%NULL @mask).
 * @windowing: Windowing type to use.
 * @level: The first polynomial degree to keep in the area; lower degrees than @level are subtracted.
 * @target_field: A data field to store the result to.  It will be resampled to @width×@height.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Calculates two-dimensional power spectrum density function of a data field area.
 *
 * The resulting data field has the spectrum density corresponding zero frequency (0,0) in the centre.
 *
 * Only @level values 0 (no levelling) and 1 (subtract the mean value) used to be available. For SPM data, you usually
 * wish to pass 1.  Since 2.56 you can also pass 2 for mean plane subtraction.
 *
 * The reduction of the total energy by windowing is compensated by multiplying the PSDF to make its sum of squares
 * equal to the input data sum of squares.
 *
 * Do not assume the PSDF values are all positive, when masking is in effect. The PSDF should still have the correct
 * integral, but it will be contaminated with noise, both positive and negative.
 *
 * Since: 2.51
 **/
void
gwy_data_field_area_2dpsdf_mask(GwyDataField *field,
                                GwyDataField *target,
                                GwyDataField *mask,
                                GwyMaskingType masking,
                                gint col, gint row,
                                gint width, gint height,
                                GwyWindowingType windowing,
                                gint level)
{
    GwyDataField *buf, *weights;
    fftw_plan plan;
    gint i, cstride, n;
    gdouble rms, newrms, r, thresh, avg, bx, by;
    fftw_complex *cbuf;
    gdouble *b, *m, *w, *t;

    if (!_gwy_data_field_check_area(field, col, row, width, height)
        || !_gwy_data_field_check_mask(field, &mask, &masking))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(target));
    fix_levelling_degree(&level);

    /* We cannot pad to a nice-for-FFT size because we want to calculate cyclic ACF. */
    n = width*height;
    cstride = width/2 + 1;
    gwy_data_field_resample(target, width, height, GWY_INTERPOLATION_NONE);
    cbuf = gwy_fftw_new_complex(cstride*height);

    if (mask) {
        buf = gwy_data_field_new_alike(target, FALSE);
        b = buf->data;
        plan = gwy_fftw_plan_dft_r2c_2d(height, width, b, cbuf, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

        /* Level and window the area. */
        gwy_data_field_area_copy(field, target, col, row, width, height, 0, 0);
        weights = gwy_data_field_new_alike(target, FALSE);
        gwy_data_field_area_copy(mask, weights, col, row, width, height, 0, 0);
        if (level > 1) {
            gwy_data_field_area_fit_plane_mask(target, weights, masking, 0, 0, width, height, &avg, &bx, &by);
            gwy_data_field_plane_level(target, avg, bx, by);
        }
        else if (level) {
            avg = gwy_data_field_area_get_avg_mask(target, weights, masking, 0, 0, width, height);
            gwy_data_field_add(target, -avg);
        }
        rms = gwy_data_field_area_get_rms_mask(target, weights, masking, 0, 0, width, height);
        gwy_data_field_fft_window(target, windowing);

        newrms = gwy_data_field_area_get_rms_mask(target, weights, masking, 0, 0, width, height);
        if (newrms > 0.0)
            gwy_data_field_multiply(target, rms/newrms);

        /* Calculate unnormalised CACF of the mask, i.e. the denominators. */
        m = mask->data;
        if (masking == GWY_MASK_INCLUDE) {
            for (i = 0; i < n; i++)
                b[i] = (m[i] > 0.0);
        }
        else {
            for (i = 0; i < n; i++)
                b[i] = (m[i] <= 0.0);
        }
        gwy_data_field_multiply_fields(target, buf, target);

        execute_2d_cacf(buf, weights, plan, cbuf, cstride);

        /* Calculate unnormalised CACF of the premultiplied image, i.e. the numerators. */
        gwy_data_field_copy(target, buf, FALSE);
        execute_2d_cacf(buf, target, plan, cbuf, cstride);

        /* Divide, with interpolation of missing values. */
        thresh = 1.000001*GWY_ROUND(log(n));
        w = weights->data;
        t = target->data;
        for (i = 0; i < n; i++) {
            if (w[i] > thresh) {
                b[i] = t[i]/w[i];
                w[i] = 0.0;
            }
            else {
                b[i] = 0.0;
                w[i] = 1.0;
            }
        }
        gwy_data_field_laplace_solve(buf, weights, -1, 1.0);

        extract_2d_fft_real(target, plan, cbuf, cstride);
        g_object_unref(weights);
        g_object_unref(buf);
    }
    else {
        t = target->data;
        plan = gwy_fftw_plan_dft_r2c_2d(height, width, t, cbuf, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

        /* Level and window the area. */
        gwy_data_field_area_copy(field, target, col, row, width, height, 0, 0);
        if (level > 1) {
            gwy_data_field_fit_plane(target, &avg, &bx, &by);
            gwy_data_field_plane_level(target, avg, bx, by);
        }
        else if (level) {
            avg = gwy_data_field_get_avg(target);
            gwy_data_field_add(target, -avg);
        }
        rms = gwy_data_field_get_rms(target);
        gwy_data_field_fft_window(target, windowing);
        newrms = gwy_data_field_get_rms(target);
        if (newrms > 0.0)
            gwy_data_field_multiply(target, rms/newrms);

        /* Do the FFT and gather squared Fourier coeffs. */
        extract_2d_fft_cnorm(target, plan, cbuf, cstride);
        gwy_data_field_multiply(target, 1.0/n);
    }
    gwy_data_field_2dfft_humanize(target);

    fftw_destroy_plan(plan);
    fftw_free(cbuf);

    gwy_si_unit_power(gwy_data_field_get_si_unit_xy(field), -1, gwy_data_field_get_si_unit_xy(target));
    gwy_si_unit_power_multiply(gwy_data_field_get_si_unit_xy(field), 2,
                               gwy_data_field_get_si_unit_z(field), 2,
                               gwy_data_field_get_si_unit_z(target));

    gwy_data_field_set_xreal(target, 2.0*G_PI/gwy_data_field_get_dx(field));
    gwy_data_field_set_yreal(target, 2.0*G_PI/gwy_data_field_get_dy(field));

    r = (width + 1 - width % 2)/2.0;
    gwy_data_field_set_xoffset(target, -gwy_data_field_jtor(target, r));

    r = (height + 1 - height % 2)/2.0;
    gwy_data_field_set_yoffset(target, -gwy_data_field_itor(target, r));

    gwy_data_field_multiply(target, 1.0/(target->xreal*target->yreal));
    gwy_data_field_invalidate(target);
}

/**
 * gwy_data_field_2dpsdf:
 * @field: A data field.
 * @windowing: Windowing type to use.
 * @level: The first polynomial degree to keep in the area; lower degrees than @level are subtracted.  Note only
 *         values 0, 1, and 2 are available at present. For SPM data, you usually wish to pass 1.
 * @target_field: A data field to store the result to.  It will be resampled to the same size as @data_field.
 *
 * Calculates two-dimensional power spectrum density function of a data field.
 *
 * See gwy_data_field_area_2dpsdf_mask() for details and discussion.
 *
 * Since: 2.51
 **/
void
gwy_data_field_2dpsdf(GwyDataField *field,
                      GwyDataField *target,
                      GwyWindowingType windowing,
                      gint level)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(field));
    gwy_data_field_area_2dpsdf_mask(field, target, NULL, GWY_MASK_IGNORE, 0, 0, field->xres, field->yres,
                                    windowing, level);
}

/**
 * gwy_data_field_psdf_to_angular_spectrum:
 * @psdf: A data field containing 2D spectral density, humanized (i.e. with zero frequency in the centre).
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, a suitable number is chosen
 *          automatically.
 *
 * Transforms 2D power spectral density to an angular spectrum.
 *
 * Returns: A new one-dimensional data line with the angular spectrum.
 *
 * Since: 2.56
 **/
GwyDataLine*
gwy_data_field_psdf_to_angular_spectrum(GwyDataField *psdf,
                                        gint nstats)
{
    gint xres, yres, nn, k;
    GwyDataLine *angspec;
    const gdouble *d;
    gdouble xreal, yreal, asum, wsum;
    gdouble *a;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(psdf), NULL);

    xres = psdf->xres;
    yres = psdf->yres;
    xreal = psdf->xreal;
    yreal = psdf->yreal;
    d = psdf->data;
    nn = xres*yres;

    if (nstats < 1) {
        nstats = floor(5.49*cbrt(nn) + 0.5);
        nstats = MAX(nstats, 2);
        nstats += (nstats & 1);
    }

    angspec = gwy_data_line_new(nstats, 2.0*G_PI, TRUE);
    gwy_data_line_set_offset(angspec, -G_PI/nstats);
    a = angspec->data;
#ifdef _OPENMP
#pragma omp parallel for if (gwy_threads_are_enabled()) default(none) \
            shared(xres,yres,xreal,yreal,nstats,d,a) private(k)
#endif
    for (k = 0; k < nstats; k++) {
        const gdouble step = 0.14589803375031551;
        gdouble alpha0 = k*2.0*G_PI/nstats, s = 0.0;
        gint kk = 0, mm = 0;

        for (kk = -5; kk <= 5; kk++) {
            gdouble alpha = alpha0 + kk/10.0*2.0*G_PI/nstats;
            /* Now we have to transform the real-space angle to an image in the reciprocal space.  For non-square
             * pixels it changes. */
            gdouble Kxsa = xreal/xres * sin(alpha);
            gdouble Kyca = yreal/yres * cos(alpha);
            gdouble h = sqrt(Kxsa*Kxsa + Kyca*Kyca);
            gdouble cb = Kyca/h, sb = Kxsa/h;
            gint m = 0;

            while (TRUE) {
                /* Zero frequency is always at res/2 in humanized data. */
                gdouble x = xres/2 + 0.5 + m*step*cb;
                gdouble y = yres/2 + 0.5 - m*step*sb;
                gint i = floor(y), j = floor(x);
                if (i < 0 || i > yres-1 || j < 0 || j > xres-1)
                    break;

                s += d[i*xres + j];
                m++;
            }
            mm += m;
        }
        if (mm)
            a[k] = s/mm;
    }

    /* The result now has some weird normalisation depending on how we hit and missed the pixels.  Preserve the
     * integral which must be sigma^2 (but we only make half of the curve). */
    wsum = gwy_data_field_get_sum(psdf) * xreal/xres * yreal/yres;
    asum = gwy_data_line_get_sum(angspec) * 2.0*G_PI/nstats;
    if (asum > 0.0)
        gwy_data_line_multiply(angspec, wsum/asum);

    gwy_si_unit_power_multiply(gwy_data_field_get_si_unit_z(psdf), 1,
                               gwy_data_field_get_si_unit_xy(psdf), 2,
                               gwy_data_line_get_si_unit_y(angspec));

    return angspec;
}

/* Does not really belong here, but is is used only by functions from this
 * source file, so... */

/**
 * gwy_data_field_angular_average:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to @nstats size.
 * @mask: Mask of pixels to include from/exclude in the averaging, or %NULL for full @data_field.
 * @masking: Masking mode to use.  See the introduction for description of masking modes.
 * @x: X-coordinate of the averaging disc origin, in real coordinates including offsets.
 * @y: Y-coordinate of the averaging disc origin, in real coordinates including offsets.
 * @r: Radius, in real coordinates.  It determines the real length of the resulting line.
 * @nstats: The number of samples the resulting line should have.  A non-positive value means the sampling will be
 *          determined automatically.
 *
 * Performs angular averaging of a part of a data field.
 *
 * The result of such averaging is an radial profile, starting from the disc centre.
 *
 * The function does not guarantee that @target_line will have exactly @nstats samples upon return.  A smaller number
 * of samples than requested may be calculated for instance if either central or outer part of the disc is excluded by
 * masking.
 *
 * Since: 2.42
 **/
void
gwy_data_field_angular_average(GwyDataField *data_field,
                               GwyDataLine *target_line,
                               GwyDataField *mask,
                               GwyMaskingType masking,
                               gdouble x,
                               gdouble y,
                               gdouble r,
                               gint nstats)
{
    gint ifrom, ito, jfrom, jto, k, kfrom, kto, xres, yres;
    gdouble xreal, yreal, dx, dy, xoff, yoff, h;
    const gdouble *d, *m;
    gdouble *target, *weight;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));
    g_return_if_fail(r >= 0.0);
    xres = data_field->xres;
    yres = data_field->yres;
    if (masking == GWY_MASK_IGNORE)
        mask = NULL;
    else if (!mask)
        masking = GWY_MASK_IGNORE;

    if (mask) {
        g_return_if_fail(GWY_IS_DATA_FIELD(mask));
        g_return_if_fail(mask->xres == xres);
        g_return_if_fail(mask->yres == yres);
    }

    xreal = data_field->xreal;
    yreal = data_field->yreal;
    xoff = data_field->xoff;
    yoff = data_field->yoff;
    g_return_if_fail(x >= xoff && x <= xoff + xreal);
    g_return_if_fail(y >= yoff && y <= yoff + yreal);
    /* Just for integer overflow; we limit i and j ranges explicitly later. */
    r = MIN(r, hypot(xreal, yreal));
    x -= xoff;
    y -= yoff;

    dx = xreal/xres;
    dy = yreal/yres;

    /* Prefer sampling close to the shorter step. */
    if (nstats < 1) {
        h = 2.0*dx*dy/(dx + dy);
        nstats = GWY_ROUND(r/h);
        nstats = MAX(nstats, 1);
    }
    h = r/nstats;

    d = data_field->data;
    m = mask ? mask->data : NULL;

    gwy_data_line_resample(target_line, nstats, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(target_line);
    gwy_data_field_copy_units_to_data_line(data_field, target_line);
    target_line->real = h*nstats;
    target_line->off = 0.0;
    target = target_line->data;
    /* Just return something for single-point lines. */
    if (nstats < 2 || r == 0.0) {
        /* NB: gwy_data_field_get_dval_real() does not use offsets. */
        target[0] = gwy_data_field_get_dval_real(data_field, x, y, GWY_INTERPOLATION_ROUND);
        return;
    }

    ifrom = (gint)floor(gwy_data_field_rtoi(data_field, y - r));
    ifrom = MAX(ifrom, 0);
    ito = (gint)ceil(gwy_data_field_rtoi(data_field, y + r));
    ito = MIN(ito, yres-1);

    jfrom = (gint)floor(gwy_data_field_rtoj(data_field, x - r));
    jfrom = MAX(jfrom, 0);
    jto = (gint)ceil(gwy_data_field_rtoj(data_field, x + r));
    jto = MIN(jto, xres-1);

    weight = g_new0(gdouble, nstats);
#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            shared(d,m,target,weight,ifrom,jfrom,ito,jto,masking,xres,yres,nstats,h,x,y,dx,dy)
#endif
    {
        gint tifrom = gwy_omp_chunk_start(ito+1 - ifrom) + ifrom;
        gint tito = gwy_omp_chunk_end(ito+1 - ifrom) + ifrom;
        gdouble *ttarget = gwy_omp_if_threads_new0(target, nstats);
        gdouble *tweight = gwy_omp_if_threads_new0(weight, nstats);
        gint i, j;

        for (i = tifrom; i < tito; i++) {
            gdouble yy = (i + 0.5)*dy - y;
            for (j = jfrom; j <= jto; j++) {
                gdouble xx = (j + 0.5)*dx - x;
                gdouble v = d[i*xres + j];
                gdouble rr;
                gint kk;

                if ((masking == GWY_MASK_INCLUDE && m[i*xres + j] <= 0.0)
                    || (masking == GWY_MASK_EXCLUDE && m[i*xres + j] >= 1.0))
                    continue;

                rr = sqrt(xx*xx + yy*yy)/h;
                kk = floor(rr);
                if (kk+1 >= nstats) {
                    if (kk+1 == nstats) {
                        ttarget[kk] += v;
                        tweight[kk] += 1.0;
                    }
                    continue;
                }

                rr -= kk;
                if (rr <= 0.5)
                    rr = 2.0*rr*rr;
                else
                    rr = 1.0 - 2.0*(1.0 - rr)*(1.0 - rr);

                ttarget[kk] += (1.0 - rr)*v;
                ttarget[kk+1] += rr*v;
                tweight[kk] += 1.0 - rr;
                tweight[kk+1] += rr;
            }
        }

        gwy_omp_if_threads_sum_double(weight, tweight, nstats);
        gwy_omp_if_threads_sum_double(target, ttarget, nstats);
    }

    /* Get rid of initial and trailing no-data segment. */
    for (kfrom = 0; kfrom < nstats; kfrom++) {
        if (weight[kfrom])
            break;
    }
    for (kto = nstats-1; kto > kfrom; kto--) {
        if (weight[kto])
            break;
    }
    if (kto - kfrom < 2) {
        /* XXX: This is not correct.  We do not care. */
        target_line->real = h;
        target[0] = gwy_data_field_get_dval_real(data_field, x, y, GWY_INTERPOLATION_ROUND);
        return;
    }

    if (kfrom != 0 || kto != nstats-1) {
        nstats = kto+1 - kfrom;
        gwy_data_line_resize(target_line, kfrom, kto+1);
        target = target_line->data;
        target_line->off = kfrom*h;
        memmove(weight, weight + kfrom, nstats*sizeof(gdouble));
    }
    g_assert(weight[0]);
    g_assert(weight[nstats-1]);

    /* Fill holes where we have no weight, this can occur near the start if large nstats is requested. */
    kfrom = -1;
    for (k = 0; k < nstats; k++) {
        if (weight[k]) {
            target[k] /= weight[k];
            if (kfrom+1 != k) {
                gdouble first = target[kfrom];
                gdouble last = target[k];
                gint j;
                for (j = kfrom+1; j < k; j++) {
                    gdouble w = (j - kfrom)/(gdouble)(k - kfrom);
                    target[j] = w*last + (1.0 - w)*first;
                }
            }
            kfrom = k;
        }
    }

    g_free(weight);
}

/**************************************************************************
 *
 * Masked ACF, HHCF, PSDF
 * DOI 10.1016/j.ultramic.2012.08.002
 *
 **************************************************************************/

static void
row_assign_mask(GwyDataField *mask,
                guint col,
                guint row,
                guint width,
                GwyMaskingType masking,
                gdouble *out)
{
    const gdouble *m = mask->data + row*mask->xres + col;
    guint j;

    if (masking == GWY_MASK_INCLUDE) {
        for (j = width; j; j--, out++, m++)
            *out = (*m > 0.0);
    }
    else {
        for (j = width; j; j--, out++, m++)
            *out = (*m <= 0.0);
    }
}

static void
row_accumulate(gdouble *accum,
               const gdouble *data,
               guint size)
{
    guint j;

    for (j = size; j; j--, accum++, data++)
        *accum += *data;
}

/* FFTW calculates unnormalised DFT so we divide the result of the first transformation with (1/√size)² = 1/size and
 * keep the second transfrom as-is to obtain exactly g_k.
 *
 * Here we deviate from the paper and try to smoothly interpolate the missing values to reduce spurious high-frequency
 * content.  It helps sometimes... */
static void
row_divide_nonzero_with_laplace(const gdouble *numerator,
                                const gdouble *denominator,
                                gdouble *out,
                                guint size, guint thresh)
{
    GwyDataLine *line, *mask;
    guint j;
    gboolean have_zero = FALSE;

    for (j = 0; j < size; j++) {
        if (denominator[j] > thresh)
            out[j] = numerator[j]/denominator[j];
        else {
            out[j] = 0.0;
            have_zero = TRUE;
        }
    }

    if (!have_zero)
        return;

    line = gwy_data_line_new(size, size, FALSE);
    gwy_assign(line->data, out, size);
    mask = gwy_data_line_new(size, size, FALSE);
    for (j = 0; j < size; j++)
        mask->data[j] = (denominator[j] == 0);

    gwy_data_line_correct_laplace(line, mask);
    gwy_assign(out, line->data, size);

    g_object_unref(line);
    g_object_unref(mask);
}

static void
row_accum_cnorm(gdouble *accum,
                const fftw_complex *fftc,
                guint size,
                gdouble q)
{
    gdouble *out = accum, *out2 = accum + (size-1);
    gdouble re = gwycreal(*fftc), im = gwycimag(*fftc);
    gdouble v;
    guint j;

    q /= size;
    v = q*(re*re + im*im);
    *out += v;
    out++, fftc++;
    for (j = (size + 1)/2 - 1; j; j--, fftc++, out++, out2--) {
        re = gwycreal(*fftc);
        im = gwycimag(*fftc);
        v = q*(re*re + im*im);
        *out += v;
        *out2 += v;
    }
    if (size % 2 == 0) {
        re = gwycreal(*fftc);
        im = gwycimag(*fftc);
        v = q*(re*re + im*im);
        *out += v;
    }
}

static void
row_extfft_accum_cnorm(fftw_plan plan,
                       gdouble *fftr,
                       gdouble *accum,
                       fftw_complex *fftc,
                       guint size,
                       guint width,
                       gdouble q)
{
    gwy_clear(fftr + width, size - width);
    gwy_fftw_execute(plan);
    row_accum_cnorm(accum, fftc, size, q);
}

/* Calculate the product A*B+AB*, equal to 2*(Re A Re B + Im A Im B), of two R2HC outputs (the result is added to @out
 * including the redundant even terms). */
static void
row_accum_cprod(const fftw_complex *fftca,
                const fftw_complex *fftcb,
                gdouble *out,
                guint size,
                gdouble q)
{
    gdouble *out2 = out + size-1;
    gdouble rea = gwycreal(*fftca), ima = gwycimag(*fftca),
            reb = gwycreal(*fftcb), imb = gwycimag(*fftcb), v;
    guint j;

    q *= 2.0/size;
    v = q*(rea*reb + ima*imb);
    *out += v;
    out++, fftca++, fftcb++;
    for (j = (size + 1)/2 - 1; j; j--, out++, fftca++, fftcb++, out2--) {
        rea = gwycreal(*fftca);
        ima = gwycimag(*fftca);
        reb = gwycreal(*fftcb);
        imb = gwycimag(*fftcb);
        v = q*(rea*reb + ima*imb);
        *out += v;
        *out2 += v;
    }
    if (size % 2 == 0) {
        rea = gwycreal(*fftca);
        ima = gwycimag(*fftca);
        reb = gwycreal(*fftcb);
        imb = gwycimag(*fftcb);
        v = q*(rea*reb + ima*imb);
        *out += v;
    }
}

/* Used in cases when we expect the imaginary part to be zero but do not want to bother with specialised DCT. */
static void
row_extfft_extract_re(fftw_plan plan,
                      gdouble *fftr,
                      gdouble *out,
                      fftw_complex *fftc,
                      guint size,
                      guint width)
{
    guint j;

    gwy_assign(fftr, out, size);
    gwy_fftw_execute(plan);
    for (j = 0; j < width; j++)
        out[j] = gwycreal(fftc[j]);
}

static void
row_extfft_symmetrise_re(fftw_plan plan,
                         gdouble *fftr,
                         gdouble *out,
                         fftw_complex *fftc,
                         guint size)
{
    gdouble *out2 = out + size-1;
    guint j;

    gwy_assign(fftr, out, size);
    gwy_fftw_execute(plan);

    *out = gwycreal(*fftc);
    out++, fftc++;
    for (j = (size + 1)/2 - 1; j; j--, fftc++, out++, out2--)
        *out = *out2 = gwycreal(*fftc);
    if (size % 2 == 0)
        *out = gwycreal(*fftc);
}

static void
row_accumulate_vk(const gdouble *data,
                  gdouble *v,
                  guint size)
{
    const gdouble *data2 = data + (size-1);
    gdouble sum = 0.0;
    guint j;

    v += size-1;
    for (j = size; j; j--, data++, data2--, v--) {
        sum += (*data)*(*data) + (*data2)*(*data2);
        *v += sum;
    }
}

/* Level a row of data by subtracting the mean value. */
static void
row_level1(const gdouble *in,
           gdouble *out,
           guint n)
{
    gdouble sumz = 0.0;
    const gdouble *pdata = in;
    guint i;

    for (i = n; i; i--, pdata++)
        sumz += *pdata;

    sumz /= n;
    for (i = n; i; i--, in++, out++)
        *out = *in - sumz;
}

static void
row_level2(const gdouble *in,
           gdouble *out,
           guint n)
{
    gdouble sumx = 0.0, sumxx = 0.0, sumz = 0.0, sumxz = 0.0, a, b;
    const gdouble *pdata = in;
    guint i;

    if (n < 2) {
        gwy_clear(out, n);
        return;
    }

    for (i = n; i; i--, pdata++) {
        gdouble z = *pdata;
        gdouble x = i - 0.5*n;
        sumz += z;
        sumxz += x*z;
        sumx += x;
        sumxx += x*x;
    }

    {
        gdouble matrix[3], rhs[2];

        matrix[0] = n;
        matrix[1] = sumx;
        matrix[2] = sumxx;
        rhs[0] = sumz;
        rhs[1] = sumxz;
        gwy_math_choleski_decompose(2, matrix);
        gwy_math_choleski_solve(2, matrix, rhs);
        a = rhs[0];
        b = rhs[1];
    }

    pdata = in;
    for (i = n; i; i--, pdata++, out++) {
        gdouble z = *pdata;
        gdouble x = i - 0.5*n;
        *out = z - a - b*x;
    }
}

/* Level a row of data by subtracting the mean value of data under mask and clear (set to zero) all data not under
 * mask.  Note how the zeroes nicely ensure that the subsequent functions Just Work(TM) and don't need to know we use
 * masking at all. */
static guint
row_level1_mask(const gdouble *in,
                gdouble *out,
                guint n,
                const gdouble *m,
                GwyMaskingType masking)
{
    gdouble sumz = 0.0, a;
    const gdouble *pdata = in, *mdata = m;
    guint i, nd = 0;

    if (masking == GWY_MASK_INCLUDE) {
        for (i = n; i; i--, pdata++, mdata++) {
            guint c = (*mdata > 0.0);
            gdouble z = *pdata;
            sumz += c*z;
            nd += c;
        }
    }
    else {
        for (i = n; i; i--, pdata++, mdata++) {
            guint c = (*mdata <= 0.0);
            gdouble z = *pdata;
            sumz += c*z;
            nd += c;
        }
    }

    if (!nd) {
        gwy_clear(out, n);
        return nd;
    }

    a = sumz/nd;
    pdata = in;
    mdata = m;
    if (masking == GWY_MASK_INCLUDE) {
        for (i = n; i; i--, pdata++, mdata++, out++) {
            guint c = (*mdata > 0.0);
            gdouble z = *pdata;
            *out = c*(z - a);
        }
    }
    else {
        for (i = n; i; i--, pdata++, mdata++, out++) {
            guint c = (*mdata <= 0.0);
            gdouble z = *pdata;
            *out = c*(z - a);
        }
    }

    return nd;
}

static guint
row_level2_mask(const gdouble *in,
                gdouble *out,
                guint n,
                const gdouble *m,
                GwyMaskingType masking)
{
    gdouble sumx = 0.0, sumxx = 0.0, sumz = 0.0, sumxz = 0.0, a, b;
    const gdouble *pdata = in, *mdata = m;
    guint i, nd = 0;

    if (masking == GWY_MASK_INCLUDE) {
        for (i = n; i; i--, pdata++, mdata++) {
            guint c = (*mdata > 0.0);
            gdouble z = c*(*pdata);
            gdouble x = c*(i - 0.5*n);
            sumz += z;
            sumxz += x*z;
            sumx += x;
            sumxx += x*x;
            nd += c;
        }
    }
    else {
        for (i = n; i; i--, pdata++, mdata++) {
            guint c = (*mdata <= 0.0);
            gdouble z = c*(*pdata);
            gdouble x = c*(i - 0.5*n);
            sumz += z;
            sumxz += x*z;
            sumx += x;
            sumxx += x*x;
            nd += c;
        }
    }

    if (nd < 2) {
        gwy_clear(out, n);
        return nd;
    }

    {
        gdouble matrix[3], rhs[2];

        matrix[0] = nd;
        matrix[1] = sumx;
        matrix[2] = sumxx;
        rhs[0] = sumz;
        rhs[1] = sumxz;
        gwy_math_choleski_decompose(2, matrix);
        gwy_math_choleski_solve(2, matrix, rhs);
        a = rhs[0];
        b = rhs[1];
    }

    pdata = in;
    mdata = m;
    if (masking == GWY_MASK_INCLUDE) {
        for (i = n; i; i--, pdata++, mdata++, out++) {
            guint c = (*mdata > 0.0);
            gdouble z = *pdata;
            gdouble x = i - 0.5*n;
            *out = (z - a - b*x)*c;
        }
    }
    else {
        for (i = n; i; i--, pdata++, mdata++, out++) {
            guint c = (*mdata <= 0.0);
            gdouble z = *pdata;
            gdouble x = i - 0.5*n;
            *out = (z - a - b*x)*c;
        }
    }

    return nd;
}

/* Window a row using a sampled windowing function. */
static void
row_window(gdouble *data, const gdouble *window, guint n)
{
    guint i;

    for (i = n; i; i--, data++, window++)
        *data *= *window;
}

/* Level and count the number of valid data in a row */
static guint
row_level_and_count(const gdouble *in,
                    gdouble *out,
                    guint width,
                    GwyDataField *mask,
                    GwyMaskingType masking,
                    guint maskcol,
                    guint maskrow,
                    guint level)
{
    guint i, count;
    const gdouble *m;

    if (!mask || masking == GWY_MASK_IGNORE) {
        if (level > 1)
            row_level2(in, out, width);
        else if (level)
            row_level1(in, out, width);
        else
            gwy_assign(out, in, width);
        return width;
    }

    m = mask->data + mask->xres*maskrow + maskcol;
    if (level > 1)
        return row_level2_mask(in, out, width, m, masking);
    else
        return row_level1_mask(in, out, width, m, masking);

    count = 0;
    if (masking == GWY_MASK_INCLUDE) {
        for (i = width; i; i--, in++, out++, m++) {
            guint c = (*m > 0.0);
            *out = c*(*in);
            count += c;
        }
    }
    else {
        for (i = width; i; i--, in++, out++, m++) {
            guint c = (*m <= 0.0);
            *out = c*(*in);
            count += c;
        }
    }
    return count;
}

static void
set_cf_units(GwyDataField *field,
             GwyDataLine *line,
             GwyDataLine *weights)
{
    gwy_data_field_copy_units_to_data_line(field, line);
    if (line->si_unit_y)
        gwy_si_unit_power(line->si_unit_y, 2, line->si_unit_y);

    if (weights) {
        gwy_data_field_copy_units_to_data_line(field, weights);
        gwy_si_unit_set_from_string(gwy_data_line_get_si_unit_y(weights), NULL);
    }
}

/**
 * gwy_data_field_area_row_acf:
 * @field: A two-dimensional data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use (has any effect only with non-%NULL @mask).
 * @level: The first polynomial degree to keep in the rows, lower degrees than @level are subtracted.
 * @weights: Line to store the denominators to (or %NULL).  It will be resized to match the returned line.  The
 *           denominators are integers equal to the number of terms that contributed to each value.  They are suitable
 *           as fitting weights if the ACF is fitted.
 *
 * Calculates the row-wise autocorrelation function (ACF) of a field.
 *
 * The calculated ACF has the natural number of points, i.e. @width.
 *
 * Masking is performed by omitting all terms that contain excluded pixels. Since different rows contain different
 * numbers of pixels, the resulting ACF values are calculated as a weighted sums where weight of each row's
 * contribution is proportional to the number of contributing terms.  In other words, the weighting is fair: each
 * contributing pixel has the same influence on the result.
 *
 * Only @level values 0 (no levelling) and 1 (subtract the mean value) used to be available. For SPM data, you usually
 * wish to pass 1.  Since 2.56 you can also pass 2 for mean line subtraction.
 *
 * Returns: A new one-dimensional data line with the ACF.
 **/
GwyDataLine*
gwy_data_field_area_row_acf(GwyDataField *field,
                            GwyDataField *mask,
                            GwyMaskingType masking,
                            guint col, guint row,
                            guint width, guint height,
                            guint level,
                            GwyDataLine *weights)
{
    GwyDataLine *line = NULL;
    fftw_complex *fftc;
    const gdouble *base;
    gdouble *fftr, *accum_data, *accum_mask;
    guint nfullrows = 0, nemptyrows = 0;
    guint size, cstride, i, j;
    fftw_plan plan;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(field), NULL);
    g_return_val_if_fail(!mask || GWY_IS_DATA_FIELD(mask), NULL);
    if (masking == GWY_MASK_IGNORE)
        mask = NULL;
    if (mask) {
        g_return_val_if_fail(mask->xres == field->xres
                             && mask->yres == field->yres, NULL);
    }
    g_return_val_if_fail(!weights || GWY_IS_DATA_LINE(weights), NULL);

    /* Transform size must be at least twice the data size for zero padding. An even size is necessary due to
     * alignment constraints in FFTW. Using this size for all buffers is a bit excessive but safe. */
    line = gwy_data_line_new(width, 1.0, TRUE);
    size = gwy_fft_find_nice_size((width + 1)/2*4);
    /* The innermost (contiguous) dimension of R2C the complex output is slightly larger than the real input.  Note
     * @cstride is measured in fftw_complex, multiply it by 2 for doubles. */
    cstride = size/2 + 1;
    base = field->data + row*field->xres + col;
    fftr = gwy_fftw_new_real(size);
    accum_data = g_new(gdouble, 2*size);
    accum_mask = accum_data + size;
    fftc = gwy_fftw_new_complex(cstride);
    plan = gwy_fftw_plan_dft_r2c_1d(size, fftr, fftc,
                                    FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    gwy_clear(accum_data, size);
    gwy_clear(accum_mask, size);

    /* Gather squared Fourier coefficients for all rows. */
    for (i = 0; i < height; i++) {
        guint count = row_level_and_count(base + i*field->xres, fftr, width, mask, masking, col, row + i, level);
        if (!count) {
            nemptyrows++;
            continue;
        }

        /* Calculate and gather squared Fourier coefficients of the data. */
        row_extfft_accum_cnorm(plan, fftr, accum_data, fftc, size, width, 1.0);

        if (count == width) {
            nfullrows++;
            continue;
        }

        /* Calculate and gather squared Fourier coefficients of the mask. */
        row_assign_mask(mask, col, row + i, width, masking, fftr);
        row_extfft_accum_cnorm(plan, fftr, accum_mask, fftc, size, width, 1.0);
    }

    /* Numerator of G_k, i.e. FFT of squared data Fourier coefficients. */
    row_extfft_extract_re(plan, fftr, accum_data, fftc, size, width);

    /* Denominator of G_k, i.e. FFT of squared mask Fourier coefficients. Don't perform the FFT if there were no
     * partial rows. */
    if (nfullrows + nemptyrows < height)
        row_extfft_extract_re(plan, fftr, accum_mask, fftc, size, width);

    for (j = 0; j < width; j++) {
        /* Denominators must be rounded to integers because they are integers and this permits to detect zeroes in the
         * denominator. */
        accum_mask[j] = GWY_ROUND(accum_mask[j]) + nfullrows*(width - j);
    }
    row_divide_nonzero_with_laplace(accum_data, accum_mask, line->data, line->res, GWY_ROUND(log(width*height)));

    line->real = gwy_data_field_get_dx(field)*line->res;
    /* line->off = -0.5*line->real/line->res; */

    if (weights) {
        gwy_data_line_resample(weights, line->res, GWY_INTERPOLATION_NONE);
        gwy_data_line_set_real(weights, line->real);
        gwy_data_line_set_offset(weights, line->off);
        gwy_assign(weights->data, accum_mask, weights->res);
    }

    fftw_destroy_plan(plan);
    fftw_free(fftc);
    g_free(accum_data);
    fftw_free(fftr);

    set_cf_units(field, line, weights);
    return line;
}

static void
normalise_window_square(gdouble *window, guint n)
{
    gdouble s = 0.0;
    guint i;

    for (i = 0; i < n; i++)
        s += window[i]*window[i];

    if (!s)
        return;

    s = sqrt(n/s);
    for (i = 0; i < n; i++)
        window[i] *= s;
}

/**
 * gwy_data_field_area_row_psdf:
 * @field: A two-dimensional data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use (has any effect only with non-%NULL @mask).
 * @windowing: Windowing type to use.
 * @level: The first polynomial degree to keep in the rows; lower degrees than @level are subtracted.
 *
 * Calculates the row-wise power spectrum density function (PSDF) of a rectangular part of a field.
 *
 * The calculated PSDF has the natural number of points that follows from DFT, i.e. @width/2+1.
 *
 * The reduction of the total energy by windowing is compensated by multiplying the PSDF to make its sum of squares
 * equal to the input data sum of squares.
 *
 * Masking is performed by omitting all terms that contain excluded pixels. Since different rows contain different
 * numbers of pixels, the resulting PSDF is calculated as a weighted sum where each row's weight is proportional to
 * the number of contributing pixels.  In other words, the weighting is fair: each contributing pixel has the same
 * influence on the result.
 *
 * Only @level values 0 (no levelling) and 1 (subtract the mean value) used to be available. For SPM data, you usually
 * wish to pass 1.  Since 2.56 you can also pass 2 for mean line subtraction.
 *
 * Do not assume the PSDF values are all positive, when masking is in effect. The PSDF should still have the correct
 * integral, but it will be contaminated with noise, both positive and negative.
 *
 * Returns: A new one-dimensional data line with the PSDF.
 **/
GwyDataLine*
gwy_data_field_area_row_psdf(GwyDataField *field,
                             GwyDataField *mask,
                             GwyMaskingType masking,
                             guint col, guint row,
                             guint width, guint height,
                             GwyWindowingType windowing,
                             guint level)
{
    GwyDataLine *line = NULL;
    fftw_complex *fftc;
    const gdouble *base;
    gdouble *fftr, *accum_data, *accum_mask, *window;
    guint nfullrows = 0, nemptyrows = 0;
    guint size, cstride, i, j;
    fftw_plan plan;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(field), NULL);
    g_return_val_if_fail(!mask || GWY_IS_DATA_FIELD(mask), NULL);
    if (masking == GWY_MASK_IGNORE)
        mask = NULL;
    if (mask) {
        g_return_val_if_fail(mask->xres == field->xres && mask->yres == field->yres, NULL);
    }
    fix_levelling_degree(&level);

    /* The innermost (contiguous) dimension of R2C the complex output is
     * slightly larger than the real input.  Note @cstride is measured in
     * fftw_complex, multiply it by 2 for doubles. */
    cstride = width/2 + 1;
    /* An even size is necessary due to alignment constraints in FFTW.
     * Using this size for all buffers is a bit excessive but safe. */
    line = gwy_data_line_new(cstride, 1.0, TRUE);
    size = (width + 3)/4*4;
    base = field->data + row*field->xres + col;
    fftr = gwy_fftw_new_real(size);
    accum_data = g_new(gdouble, 2*size + width);
    accum_mask = accum_data + size;
    window = accum_data + 2*size;
    fftc = gwy_fftw_new_complex(cstride);

    gwy_clear(accum_data, size);
    gwy_clear(accum_mask, size);

    for (j = 0; j < width; j++)
        window[j] = 1.0;
    gwy_fft_window(width, window, windowing);
    normalise_window_square(window, width);
    plan = gwy_fftw_plan_dft_r2c_1d(width, fftr, fftc, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    for (i = 0; i < height; i++) {
        guint count = row_level_and_count(base + i*field->xres, fftr, width, mask, masking, col, row + i, level);
        if (!count) {
            nemptyrows++;
            continue;
        }

        /* Calculate and gather squared Fourier coefficients of the data. */
        row_window(fftr, window, width);
        row_extfft_accum_cnorm(plan, fftr, accum_data, fftc, width, width, 1.0);

        if (count == width) {
            nfullrows++;
            continue;
        }

        /* Calculate and gather squared Fourier coefficients of the mask. */
        row_assign_mask(mask, col, row + i, width, masking, fftr);
        row_extfft_accum_cnorm(plan, fftr, accum_mask, fftc, width, width, 1.0);
    }

    /* Numerator of A_k, i.e. FFT of squared data Fourier coefficients. */
    row_extfft_symmetrise_re(plan, fftr, accum_data, fftc, width);

    /* Denominator of A_k, i.e. FFT of squared mask Fourier coefficients. Don't perform the FFT if there were no
     * partial rows. */
    if (nfullrows + nemptyrows < height)
        row_extfft_symmetrise_re(plan, fftr, accum_mask, fftc, width);

    for (j = 0; j < width; j++) {
        /* Denominators must be rounded to integers because they are integers and this permits to detect zeroes in the
         * denominator. */
        accum_mask[j] = GWY_ROUND(accum_mask[j]) + nfullrows*width;
    }
    row_divide_nonzero_with_laplace(accum_data, accum_mask, fftr, width, GWY_ROUND(log(width*height)));

    /* The transform is the other way round – for complex numbers.  Since it is in fact a DCT here we don't care and
     * run it as a forward transform. */
    gwy_fftw_execute(plan);
    for (j = 0; j < line->res; j++)
        line->data[j] = gwycreal(fftc[j]);

    fftw_destroy_plan(plan);
    fftw_free(fftc);
    fftw_free(fftr);
    g_free(accum_data);

    gwy_data_line_multiply(line, gwy_data_field_get_dx(field)/(2*G_PI));
    line->real = G_PI/gwy_data_field_get_dx(field);
    /* line->off = -0.5*line->real/line->res; */

    gwy_si_unit_power(gwy_data_field_get_si_unit_xy(field), -1, gwy_data_line_get_si_unit_x(line));
    gwy_si_unit_power_multiply(gwy_data_field_get_si_unit_xy(field), 1,
                               gwy_data_field_get_si_unit_z(field), 2,
                               gwy_data_line_get_si_unit_y(line));
    return line;
}

/**
 * gwy_data_field_area_row_hhcf:
 * @field: A two-dimensional data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use (has any effect only with non-%NULL @mask).
 * @level: The first polynomial degree to keep in the rows, lower degrees than @level are subtracted.
 * @weights: Line to store the denominators to (or %NULL).  It will be resized to match the returned line.  The
 *           denominators are integers equal to the number of terms that contributed to each value.  They are suitable
 *           as fitting weights if the HHCF is fitted.
 *
 * Calculates the row-wise height-height correlation function (HHCF) of a rectangular part of a field.
 *
 * The calculated HHCF has the natural number of points, i.e. @width.
 *
 * Masking is performed by omitting all terms that contain excluded pixels. Since different rows contain different
 * numbers of pixels, the resulting HHCF values are calculated as a weighted sums where weight of each row's
 * contribution is proportional to the number of contributing terms.  In other words, the weighting is fair: each
 * contributing pixel has the same influence on the result.
 *
 * Only @level values 0 (no levelling) and 1 (subtract the mean value) used to be available.  There is no difference
 * between them for HHCF. Since 2.56 you can also pass 2 for mean line subtraction.
 *
 * Returns: A new one-dimensional data line with the HHCF.
 **/
GwyDataLine*
gwy_data_field_area_row_hhcf(GwyDataField *field,
                             GwyDataField *mask,
                             GwyMaskingType masking,
                             guint col, guint row,
                             guint width, guint height,
                             guint level,
                             GwyDataLine *weights)
{
    GwyDataLine *line = NULL;
    fftw_complex *fftc, *tmp;
    const gdouble *base, *q;
    gdouble *fftr, *accum_data, *accum_mask, *accum_v, *p;
    guint nfullrows = 0, nemptyrows = 0;
    guint size, cstride, i, j;
    fftw_plan plan;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(field), NULL);
    g_return_val_if_fail(!mask || GWY_IS_DATA_FIELD(mask), NULL);
    if (masking == GWY_MASK_IGNORE)
        mask = NULL;
    if (mask) {
        g_return_val_if_fail(mask->xres == field->xres && mask->yres == field->yres, NULL);
    }
    g_return_val_if_fail(!weights || GWY_IS_DATA_LINE(weights), NULL);
    fix_levelling_degree(&level);

    /* Transform size must be at least twice the data size for zero padding.
     * An even size is necessary due to alignment constraints in FFTW.
     * Using this size for all buffers is a bit excessive but safe. */
    line = gwy_data_line_new(width, 1.0, TRUE);
    size = gwy_fft_find_nice_size((width + 1)/2*4);
    /* The innermost (contiguous) dimension of R2C the complex output is
     * slightly larger than the real input.  Note @cstride is measured in
     * fftw_complex, multiply it by 2 for doubles. */
    cstride = size/2 + 1;
    base = field->data + row*field->xres + col;
    fftr = gwy_fftw_new_real(size);
    accum_data = g_new(gdouble, 3*size);
    accum_mask = accum_data + size;
    accum_v = accum_data + 2*size;
    fftc = gwy_fftw_new_complex(cstride);
    tmp = g_new(fftw_complex, cstride);
    plan = gwy_fftw_plan_dft_r2c_1d(size, fftr, fftc, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    gwy_clear(accum_data, size);
    gwy_clear(accum_mask, size);
    gwy_clear(accum_v, size);

    // Gather V_ν-2|Z_ν|² for all rows, except that for full rows we actually gather just -2|Z_ν|² because v_k can be
    // calculated without DFT.
    for (i = 0; i < height; i++) {
        guint count = row_level_and_count(base + i*field->xres, fftr, width, mask, masking, col, row + i, level);
        if (!count) {
            nemptyrows++;
            continue;
        }

        /* Calculate v_k before FFT destroys the input levelled/filtered data. */
        if (count == width)
            row_accumulate_vk(fftr, accum_v, width);
        else {
            // For partial rows, we will need the data later to calculate FFT of their squares.  Save them to the line
            // that conveniently has the right size.
            gwy_assign(line->data, fftr, width);
        }

        /* Calculate and gather -2 times squared Fourier coefficients. */
        row_extfft_accum_cnorm(plan, fftr, accum_data, fftc, size, width, -2.0);

        if (count == width) {
            nfullrows++;
            continue;
        }

        /* First calculate U_ν (Fourier cofficients of squared data).  Save them to tmp. */
        q = line->data;
        p = fftr;
        for (j = width; j; j--, p++, q++)
            *p = (*q)*(*q);
        gwy_clear(fftr + width, size - width);
        gwy_fftw_execute(plan);
        gwy_assign(tmp, fftc, cstride);

        /* Mask.  We need the intermediate result C_ν to combine it with U_ν. */
        row_assign_mask(mask, col, row + i, width, masking, fftr);
        gwy_clear(fftr + width, size - width);
        gwy_fftw_execute(plan);

        /* Accumulate V_ν (calculated from C_ν and U_ν) to accum_data. */
        row_accum_cprod(tmp, fftc, accum_data, size, 1.0);

        /* And accumulate squared mask Fourier coeffs |C_ν|². */
        row_accum_cnorm(accum_mask, fftc, size, 1.0);
    }

    /* Numerator of H_k, excluding non-DFT data in v_k. */
    row_extfft_extract_re(plan, fftr, accum_data, fftc, size, width);
    /* Combine it with v_k to get the full numerator in accum_data. */
    row_accumulate(accum_data, accum_v, width);

    // Denominator of H_k, i.e. FFT of squared mask Fourier coefficients. Don't perform the FFT if there were no
    // partial rows.
    if (nfullrows + nemptyrows < height)
        row_extfft_extract_re(plan, fftr, accum_mask, fftc, size, width);

    for (j = 0; j < width; j++) {
        // Denominators must be rounded to integers because they are integers and this permits to detect zeroes in the
        // denominator.
        accum_mask[j] = GWY_ROUND(accum_mask[j]) + nfullrows*(width - j);
    }
    row_divide_nonzero_with_laplace(accum_data, accum_mask, line->data, line->res, GWY_ROUND(log(width*height)));

    line->real = gwy_data_field_get_dx(field)*line->res;
    /* line->off = -0.5*line->real/line->res; */

    if (weights) {
        gwy_data_line_resample(weights, line->res, GWY_INTERPOLATION_NONE);
        gwy_data_line_set_real(weights, line->real);
        gwy_data_line_set_offset(weights, line->off);
        gwy_assign(weights->data, accum_mask, weights->res);
    }

    fftw_destroy_plan(plan);
    fftw_free(fftc);
    fftw_free(fftr);
    g_free(accum_data);
    g_free(tmp);

    set_cf_units(field, line, weights);
    return line;
}

/* Recalculate area excess based on second-order expansion to the true one, assuming the distribution is exponential.
 * */
static inline gdouble
asg_correction(gdouble ex)
{
    if (ex < 1e-3)
        return ex*(1.0 - ex*(1.0 - 3.0*ex*(1.0 - 5.0*ex*(1.0 - 7.0*ex*(1.0 - 9.0*ex*(1.0 - 11.0*ex))))));

    return sqrt(0.5*G_PI*ex) * exp(0.5/ex) * erfc(sqrt(0.5/ex));
}

/**
 * gwy_data_field_area_row_asg:
 * @field: A two-dimensional data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use (has any effect only with non-%NULL @mask).
 * @level: The first polynomial degree to keep in the rows, lower degrees than @level are subtracted.
 *
 * Calculates the row-wise area scale graph (ASG) of a rectangular part of a field.
 *
 * The calculated ASG has the natural number of points, i.e. @width-1.
 *
 * The ASG represents the apparent area excess (ratio of surface and projected area minus one) observed at given
 * length scale.  The quantity calculated by this function serves a similar purpose as ASME B46.1 area scale graph but
 * is defined differently, based on the HHCF.  See gwy_data_field_area_row_hhcf() for details of its calculation.
 *
 * Only @level values 0 (no levelling) and 1 (subtract the mean value) used to be available.  There is no difference
 * between them for HHCF. Since 2.56 you can also pass 2 for mean line subtraction.
 *
 * Returns: A new one-dimensional data line with the ASG.
 **/
GwyDataLine*
gwy_data_field_area_row_asg(GwyDataField *field,
                            GwyDataField *mask,
                            GwyMaskingType masking,
                            guint col, guint row,
                            guint width, guint height,
                            guint level)
{
    GwyDataLine *hhcf;
    GwyDataLine *line = NULL;
    gdouble dx;
    guint i;

    hhcf = gwy_data_field_area_row_hhcf(field, mask, masking, col, row, width, height, level, NULL);

    g_return_val_if_fail(hhcf, NULL);
    dx = hhcf->real/hhcf->res;
    if (hhcf->res < 2) {
        line = gwy_data_line_new(1, dx, TRUE);
        g_object_unref(hhcf);
        return line;
    }

    line = gwy_data_line_new(hhcf->res - 1, dx*(hhcf->res - 1), FALSE);
    line->off = 0.5*dx;

    for (i = 0; i < line->res; i++) {
        gdouble t = (i + 0.5)*dx + line->off;
        line->data[i] = asg_correction(hhcf->data[i+1]/(t*t));
    }

    _gwy_copy_si_unit(field->si_unit_xy, &line->si_unit_x);
    g_object_unref(hhcf);

    return line;
}

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
