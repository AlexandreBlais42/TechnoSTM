/*
 *  $Id: inttrans.c 25190 2023-01-06 12:32:34Z yeti-dn $
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

#include "config.h"
#include <string.h>
#include <libgwyddion/gwymacros.h>
#include <libgwyddion/gwymath.h>
#include <libprocess/arithmetic.h>
#include <libprocess/inttrans.h>
#include <libprocess/linestats.h>
#include <libprocess/simplefft.h>
#include <libprocess/level.h>
#include <libprocess/stats.h>
#include "gwyprocessinternal.h"
#include "gwyfftw.h"

#ifdef HAVE_SINCOS
#define _gwy_sincos sincos
#else
static inline void
_gwy_sincos(gdouble x, gdouble *s, gdouble *c)
{
    *s = sin(x);
    *c = cos(x);
}
#endif

static void gwy_data_line_fft_do          (GwyDataLine *rsrc,
                                           GwyDataLine *isrc,
                                           GwyDataLine *rdest,
                                           GwyDataLine *idest,
                                           GwyTransformDirection direction);
static void gwy_data_line_fft_real_do     (GwyDataLine *rsrc,
                                           GwyDataLine *ibuf,
                                           GwyDataLine *rdest,
                                           GwyDataLine *idest,
                                           GwyTransformDirection direction);
static void zoom_fft_1d_do                (const gdouble *rein,
                                           const gdouble *imin,
                                           gint n,
                                           gdouble *reout,
                                           gdouble *imout,
                                           gint m,
                                           gdouble f0,
                                           gdouble f1);
static void zoom_fft_2d_do                (const gdouble *rein,
                                           const gdouble *imin,
                                           gint nx,
                                           gint ny,
                                           gdouble *reout,
                                           gdouble *imout,
                                           gint mx,
                                           gint my,
                                           gdouble fx0,
                                           gdouble fy0,
                                           gdouble fx1,
                                           gdouble fy1);
static void gwy_data_field_area_2dfft_real(GwyDataField *ra,
                                           GwyDataField *rb,
                                           GwyDataField *ib,
                                           gint col,
                                           gint row,
                                           gint width,
                                           gint height,
                                           GwyWindowingType windowing,
                                           GwyTransformDirection direction,
                                           gboolean preserverms,
                                           gint level);
static void gwy_data_field_area_xfft      (GwyDataField *ra,
                                           GwyDataField *ia,
                                           GwyDataField *rb,
                                           GwyDataField *ib,
                                           gint col,
                                           gint row,
                                           gint width,
                                           gint height,
                                           GwyWindowingType windowing,
                                           GwyTransformDirection direction,
                                           gboolean preserverms,
                                           gint level);
static void gwy_data_field_xfft_do        (GwyDataField *rin,
                                           GwyDataField *iin,
                                           GwyDataField *rout,
                                           GwyDataField *iout,
                                           GwyTransformDirection direction);
static void gwy_data_field_area_xfft_real (GwyDataField *ra,
                                           GwyDataField *rb,
                                           GwyDataField *ib,
                                           gint col,
                                           gint row,
                                           gint width,
                                           gint height,
                                           GwyWindowingType windowing,
                                           GwyTransformDirection direction,
                                           gboolean preserverms,
                                           gint level);
static void gwy_data_field_xfft_real_do   (GwyDataField *rin,
                                           GwyDataField *ibuf,
                                           GwyDataField *rout,
                                           GwyDataField *iout,
                                           GwyTransformDirection direction);
static void gwy_data_field_area_yfft      (GwyDataField *ra,
                                           GwyDataField *ia,
                                           GwyDataField *rb,
                                           GwyDataField *ib,
                                           gint col,
                                           gint row,
                                           gint width,
                                           gint height,
                                           GwyWindowingType windowing,
                                           GwyTransformDirection direction,
                                           gboolean preserverms,
                                           gint level);
static void gwy_data_field_yfft_do        (GwyDataField *rin,
                                           GwyDataField *iin,
                                           GwyDataField *rout,
                                           GwyDataField *iout,
                                           GwyTransformDirection direction);
static void gwy_data_field_area_yfft_real (GwyDataField *ra,
                                           GwyDataField *rb,
                                           GwyDataField *ib,
                                           gint col,
                                           gint row,
                                           gint width,
                                           gint height,
                                           GwyWindowingType windowing,
                                           GwyTransformDirection direction,
                                           gboolean preserverms,
                                           gint level);
static void gwy_data_field_yfft_real_do   (GwyDataField *rin,
                                           GwyDataField *ibuf,
                                           GwyDataField *rout,
                                           GwyDataField *iout,
                                           GwyTransformDirection direction);
static void gwy_level_simple              (gint n,
                                           gint stride,
                                           gdouble *data,
                                           gint level);
static void gwy_preserve_rms_simple       (gint nsrc,
                                           gint stridesrc,
                                           const gdouble *src1,
                                           const gdouble *src2,
                                           gint ndata,
                                           gint stridedata,
                                           gdouble *data1,
                                           gdouble *data2);

/**
 * gwy_data_line_fft:
 * @rsrc: Real input data line.
 * @isrc: Imaginary input data line.
 * @rdest: Real output data line.  It will be resized to the size of the input data line.
 * @idest: Imaginary output data line.  It will be resized to the size of the input data line.
 * @windowing: Windowing type to use.
 * @direction: FFT direction.
 * @interpolation: Interpolation type. Ignored since 2.8 as no resampling is performed.
 * @preserverms: %TRUE to preserve RMS value while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract line (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Calculates Fast Fourier Transform of a data line.
 *
 * A windowing or data leveling can be applied if requested.
 **/
void
gwy_data_line_fft(GwyDataLine *rsrc, GwyDataLine *isrc,
                  GwyDataLine *rdest, GwyDataLine *idest,
                  GwyWindowingType windowing,
                  GwyTransformDirection direction,
                  GwyInterpolationType interpolation,
                  gboolean preserverms,
                  gint level)
{
    g_return_if_fail(GWY_IS_DATA_LINE(rsrc));

    gwy_data_line_part_fft(rsrc, isrc, rdest, idest, 0, rsrc->res,
                           windowing, direction, interpolation, preserverms, level);
}

/**
 * gwy_data_line_part_fft:
 * @rsrc: Real input data line.
 * @isrc: Imaginary input data line. Since 2.7 it can be %NULL for real-to-complex transforms.
 * @rdest: Real output data line, it will be resized to @len.
 * @idest: Imaginary output data line, it will be resized to @len.
 * @from: The index in input lines to start from (inclusive).
 * @len: Lenght of data line part, it must be at least 2.
 * @windowing: Windowing type to use.
 * @direction: FFT direction.
 * @interpolation: Interpolation type. Ignored since 2.8 as no resampling is performed.
 * @preserverms: %TRUE to preserve RMS value while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract line (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Calculates Fast Fourier Transform of a part of a data line.
 *
 * A windowing or data leveling can be applied if requested.
 **/
void
gwy_data_line_part_fft(GwyDataLine *rsrc, GwyDataLine *isrc,
                       GwyDataLine *rdest, GwyDataLine *idest,
                       gint from, gint len,
                       GwyWindowingType windowing,
                       GwyTransformDirection direction,
                       G_GNUC_UNUSED GwyInterpolationType interpolation,
                       gboolean preserverms,
                       gint level)
{
    GwyDataLine *rbuf, *ibuf;

    g_return_if_fail(GWY_IS_DATA_LINE(rsrc));
    g_return_if_fail(!isrc || GWY_IS_DATA_LINE(isrc));
    if (isrc) {
        g_return_if_fail(!gwy_data_line_check_compatibility(rsrc, isrc, GWY_DATA_COMPATIBILITY_RES));
    }
    g_return_if_fail(GWY_IS_DATA_LINE(rdest));
    g_return_if_fail(GWY_IS_DATA_LINE(idest));
    g_return_if_fail(level >= 0 && level <= 2);
    g_return_if_fail(from >= 0 && len >= 2 && from + len <= rsrc->res);

    gwy_data_line_resample(rdest, len, GWY_INTERPOLATION_NONE);
    gwy_data_line_resample(idest, len, GWY_INTERPOLATION_NONE);

    rbuf = gwy_data_line_part_extract(rsrc, from, len);
    gwy_level_simple(len, 1, rbuf->data, level);
    gwy_data_line_fft_window(rbuf, windowing);

    if (isrc) {
        ibuf = gwy_data_line_part_extract(isrc, from, len);
        gwy_level_simple(len, 1, ibuf->data, level);
        gwy_data_line_fft_window(ibuf, windowing);
        gwy_data_line_fft_do(rbuf, ibuf, rdest, idest, direction);
        if (preserverms) {
            gwy_preserve_rms_simple(len, 1, rsrc->data + from, isrc->data + from,
                                    len, 1, rdest->data, idest->data);
        }
    }
    else {
        ibuf = gwy_data_line_new_alike(rbuf, FALSE);
        gwy_data_line_fft_real_do(rbuf, ibuf, rdest, idest, direction);
        if (preserverms) {
            gwy_preserve_rms_simple(len, 1, rsrc->data + from, NULL,
                                    len, 1, rdest->data, idest->data);
        }
    }

    g_object_unref(rbuf);
    g_object_unref(ibuf);
}

/**
 * gwy_data_line_fft_raw:
 * @rsrc: Real input data line.
 * @isrc: Imaginary input data line.  Since 2.7 it can be %NULL for real-to-complex transform.
 * @rdest: Real output data line.  It will be resized to the size of the input data line.
 * @idest: Imaginary output data line.  It will be resized to the size of the input data line.
 * @direction: FFT direction.
 *
 * Calculates Fast Fourier Transform of a data line.
 *
 * No leveling, windowing nor scaling is performed.
 *
 * The normalisation of FFT is symmetrical, so transformations in both directions are unitary.
 *
 * Since 2.8 the dimensions need not to be from the set of sizes returned by gwy_fft_find_nice_size().
 *
 * Since: 2.1
 **/
void
gwy_data_line_fft_raw(GwyDataLine *rsrc,
                      GwyDataLine *isrc,
                      GwyDataLine *rdest,
                      GwyDataLine *idest,
                      GwyTransformDirection direction)
{
    g_return_if_fail(GWY_IS_DATA_LINE(rsrc));
    g_return_if_fail(!isrc || GWY_IS_DATA_LINE(isrc));
    if (isrc) {
        g_return_if_fail(!gwy_data_line_check_compatibility(rsrc, isrc, GWY_DATA_COMPATIBILITY_RES));
    }
    g_return_if_fail(GWY_IS_DATA_LINE(rdest));
    g_return_if_fail(GWY_IS_DATA_LINE(idest));

    gwy_data_line_resample(rdest, rsrc->res, GWY_INTERPOLATION_NONE);
    gwy_data_line_resample(idest, rsrc->res, GWY_INTERPOLATION_NONE);

    if (isrc)
        g_object_ref(isrc);
    else
        isrc = gwy_data_line_new_alike(rsrc, TRUE);

    gwy_data_line_fft_do(rsrc, isrc, rdest, idest, direction);
    g_object_unref(isrc);
}

static void
gwy_data_line_fft_do(GwyDataLine *rsrc,
                     GwyDataLine *isrc,
                     GwyDataLine *rdest,
                     GwyDataLine *idest,
                     GwyTransformDirection direction)
{
    fftw_iodim dims[1], howmany_dims[1];
    gdouble *rbuf = gwy_fftw_new_real(rsrc->res);
    gdouble *ibuf = gwy_fftw_new_real(isrc->res);
    fftw_plan plan;

    dims[0].n = rsrc->res;
    dims[0].is = 1;
    dims[0].os = 1;
    howmany_dims[0].n = 1;
    howmany_dims[0].is = rsrc->res;
    howmany_dims[0].os = rsrc->res;
    /* Backward direction is equivalent to switching real and imaginary parts */
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD) {
        plan = gwy_fftw_plan_guru_split_dft(1, dims, 1, howmany_dims, rbuf, ibuf, rdest->data, idest->data,
                                            FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    }
    else {
        plan = gwy_fftw_plan_guru_split_dft(1, dims, 1, howmany_dims, ibuf, rbuf, idest->data, rdest->data,
                                            FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    }
    g_return_if_fail(plan);
    gwy_assign(rbuf, rsrc->data, rsrc->res);
    gwy_assign(ibuf, isrc->data, isrc->res);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(ibuf);
    fftw_free(rbuf);

    gwy_data_line_multiply(rdest, 1.0/sqrt(rsrc->res));
    gwy_data_line_multiply(idest, 1.0/sqrt(rsrc->res));
}

static void
gwy_data_line_fft_real_do(GwyDataLine *rsrc,
                          GwyDataLine *ibuf,
                          GwyDataLine *rdest,
                          GwyDataLine *idest,
                          GwyTransformDirection direction)
{
    fftw_iodim dims[1], howmany_dims[1];
    fftw_plan plan;
    gint j;

    dims[0].n = rsrc->res;
    dims[0].is = 1;
    dims[0].os = 1;
    howmany_dims[0].n = 1;
    howmany_dims[0].is = rsrc->res;
    howmany_dims[0].os = rsrc->res;
    /* Backward direction is equivalent to switching real and imaginary parts */
    plan = gwy_fftw_plan_guru_split_dft_r2c(1, dims, 1, howmany_dims, ibuf->data, rdest->data, idest->data,
                                            FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    g_return_if_fail(plan);
    /* R2C destroys input, and especially, the planner destroys input too */
    gwy_data_line_copy(rsrc, ibuf);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* Complete the missing half of transform.  */
    for (j = rsrc->res/2 + 1; j < rsrc->res; j++) {
        rdest->data[j] = rdest->data[rsrc->res - j];
        idest->data[j] = -idest->data[rsrc->res - j];
    }

    gwy_data_line_multiply(rdest, 1.0/sqrt(rsrc->res));
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD)
        gwy_data_line_multiply(idest, 1.0/sqrt(rsrc->res));
    else
        gwy_data_line_multiply(idest, -1.0/sqrt(rsrc->res));
}

/**
 * gwy_data_line_zoom_fft:
 * @rsrc: Real input data line.
 * @isrc: Imaginary input data line.  It can be %NULL for real-to-complex transform.
 * @rdest: Real output data line.  It will be resized to @m samples.
 * @idest: Imaginary output data line.  It will be resized to @m samples.
 * @f0: The first spatial frequency, measured in DFT frequency steps.
 * @f1: The last spatial frequency, measured in DFT frequency steps.
 * @m: The number of frequencies to compute. It must be at least 2.
 *
 * Computes Zoom FFT of a data line.
 *
 * The output is DFTs, but computed for an arbitrary linear sequence of frequencies. The frequencies do not have to
 * be in any relation to the data sampling step.
 *
 * The first item of output corresponds exactly to @f0 and the last exactly to @f1.  So the frequency sampling step
 * will be (@f1 − @f0)/(@m − 1), instead of the more usual division by @m. To follow the usual Gwyddion conventions,
 * the output data line real size will be (@f1 − @f0)/(@m − 1)*@m. If it seems confusing, just take the output as
 * indexed by integers and work with that.
 *
 * Frequency step of one corresponds to the normal DFT frequency step. Therefore, passing @f0=0, @f1=@n-1 (where @rsrc
 * has @n points) and @m=@n reproduces the usual DFT, except more slowly. The result is normalised as raw FFT and the
 * units of the output data lines are unchanged.
 *
 * The transform direction is always forward.  Windowing or other preprocessing need to be done separately beforehand.
 * They would be usually once, but followed by any number of (Zoom) FFTs.
 *
 * Since: 2.61
 **/
void
gwy_data_line_zoom_fft(GwyDataLine *rsrc,
                       GwyDataLine *isrc,
                       GwyDataLine *rdest,
                       GwyDataLine *idest,
                       gint m,
                       gdouble f0,
                       gdouble f1)
{
    g_return_if_fail(GWY_IS_DATA_LINE(rsrc));
    g_return_if_fail(!isrc || GWY_IS_DATA_LINE(isrc));
    if (isrc) {
        g_return_if_fail(!gwy_data_line_check_compatibility(rsrc, isrc, GWY_DATA_COMPATIBILITY_RES));
    }
    g_return_if_fail(GWY_IS_DATA_LINE(rdest));
    g_return_if_fail(GWY_IS_DATA_LINE(idest));
    g_return_if_fail(m > 1);

    gwy_data_line_resample(rdest, m, GWY_INTERPOLATION_NONE);
    gwy_data_line_resample(idest, m, GWY_INTERPOLATION_NONE);

    zoom_fft_1d_do(rsrc->data, isrc ? isrc->data : NULL, rsrc->res,
                   rdest->data, idest->data, m, f0, f1);

    rdest->real = idest->real = (f1 - f0)*m/(m - 1.0);
    rdest->off = idest->off = f0;
}

/**
 * gwy_data_field_zoom_fft:
 * @rsrc: Real input data field.
 * @isrc: Imaginary input data field.  It can be %NULL for real-to-complex transform.
 * @rdest: Real output data field.  It will be resized to @mx × @my samples.
 * @idest: Imaginary output data field.  It will be resized to @mx × @my samples.
 * @fx0: The first horizontal spatial frequency, measured in DFT frequency steps.
 * @fy0: The first vertical spatial frequency, measured in DFT frequency steps.
 * @fx1: The last horizontal spatial frequency, measured in DFT frequency steps.
 * @fy1: The last vetical spatial frequency, measured in DFT frequency steps.
 * @mx: The number of horizontal frequencies to compute. It must be at least 2.
 * @my: The number of vertical frequencies to compute. It must be at least 2.
 *
 * Computes Zoom FFT of a data field.
 *
 * The output is DFTs, but computed for an arbitrary 2D Cartesian grid of frequencies along @x and @y. The frequencies
 * do not have to be in any relation to the data sampling step.
 *
 * The top-left pixel of output corresponds exactly to (@fx0,@fy0) and the bottom right exactly to (@fx1,@fy1).  So
 * the frequency sampling steps will be (@fx1 − @fx0)/(@mx − 1) and (@fy1 − @fy0)/(@my − 1), instead of the more usual
 * division by @mx and @my. To follow the usual Gwyddion conventions, the output data field real size will be (@fx1
 * − @fx0)/(@mx − 1)*@mx along @x, and similarly along @y. If it seems confusing, just take the output as indexed by
 * integers and work with that.
 *
 * Frequency step of one corresponds to the normal DFT frequency step. Therefore, passing @fx0=0, @fx1=@xres–1,
 * @fy0=0, @fy1=@yres–1 (where @rsrc has @xres × @yres points), @mx=@xres and @my=@yres reproduces the usual DFT,
 * except more slowly. The result is normalised as raw FFT and the units of the output data fields are unchanged.
 *
 * The transform direction is always forward.  Windowing or other preprocessing need to be done separately beforehand.
 * They would be usually once, but followed by any number of (Zoom) FFTs.
 *
 * Since: 2.62
 **/
void
gwy_data_field_zoom_fft(GwyDataField *rsrc,
                        GwyDataField *isrc,
                        GwyDataField *rdest,
                        GwyDataField *idest,
                        gint mx,
                        gint my,
                        gdouble fx0,
                        gdouble fy0,
                        gdouble fx1,
                        gdouble fy1)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(rsrc));
    g_return_if_fail(!isrc || GWY_IS_DATA_FIELD(isrc));
    if (isrc) {
        g_return_if_fail(!gwy_data_field_check_compatibility(rsrc, isrc, GWY_DATA_COMPATIBILITY_RES));
    }
    g_return_if_fail(GWY_IS_DATA_FIELD(rdest));
    g_return_if_fail(GWY_IS_DATA_FIELD(idest));
    g_return_if_fail(mx > 1);
    g_return_if_fail(my > 1);

    gwy_data_field_resample(rdest, mx, my, GWY_INTERPOLATION_NONE);
    gwy_data_field_resample(idest, mx, my, GWY_INTERPOLATION_NONE);

    zoom_fft_2d_do(rsrc->data, isrc ? isrc->data : NULL, rsrc->xres, rsrc->yres,
                   rdest->data, idest->data, mx, my, fx0, fy0, fx1, fy1);

    rdest->xreal = idest->xreal = (fx1 - fx0)*mx/(mx - 1.0);
    rdest->yreal = idest->yreal = (fy1 - fy0)*my/(my - 1.0);
    rdest->xoff = idest->xoff = fx0;
    rdest->yoff = idest->yoff = fy0;
}

/* Precompute the factors w_k = exp(-2πik²/(ND)) */
static void
make_chirp(fftw_complex *w, gint m, gint n, gint size, gdouble D)
{
    gint minsize = m + n - 1;
    gint k, mm = MIN(m, n);
    gdouble q = -G_PI/n*D;

    /* NB: Gwydion has swapped forward and backward signs! This basically means q has the opposite sign than it would
     * normally have. */
    gwycreal(w[0]) = 1.0;
    gwycimag(w[0]) = 0.0;
    for (k = 1; k < mm; k++) {
        gdouble s, c;
        _gwy_sincos(q*k*k, &s, &c);
        gwycreal(w[k]) = gwycreal(w[size - k]) = c;
        gwycimag(w[k]) = gwycimag(w[size - k]) = s;
        //w[k] = w[size - k] = cexp(q*k*k*I);
    }
    /* Only one of the two following actually does something, depending on which of m and n is larger. */
    for (k = n; k < m; k++) {
        gdouble s, c;
        _gwy_sincos(q*k*k, &s, &c);
        gwycreal(w[k]) = c;
        gwycimag(w[k]) = s;
        //w[k] = cexp(q*k*k*I);
    }
    for (k = m; k < n; k++) {
        gdouble s, c;
        _gwy_sincos(q*k*k, &s, &c);
        gwycreal(w[size - k]) = c;
        gwycimag(w[size - k]) = s;
        //w[size - k] = cexp(q*k*k*I);
    }
    gwy_clear(w + m, size - minsize);
}

static void
make_cexp(fftw_complex *cs, gint n, gdouble q)
{
    gint k;

    gwycreal(cs[0]) = 1.0;
    gwycimag(cs[0]) = 0.0;
    for (k = 1; k < n; k++) {
        gdouble s, c;
        _gwy_sincos(q*k, &s, &c);
        gwycreal(cs[k]) = c;
        gwycimag(cs[k]) = s;
    }
}

#define cmulre(a,b) (gwycreal(a)*gwycreal(b) - gwycimag(a)*gwycimag(b))
#define cmulim(a,b) (gwycimag(a)*gwycreal(b) + gwycreal(a)*gwycimag(b))
#define cjmulre(a,b) (gwycreal(a)*gwycreal(b) + gwycimag(a)*gwycimag(b))
#define cjmulim(a,b) (gwycimag(a)*gwycreal(b) - gwycreal(a)*gwycimag(b))

static void
zoom_fft_1d_do(const gdouble *rein, const gdouble *imin, gint n,
               gdouble *reout, gdouble *imout, gint m, gdouble f0, gdouble f1)
{
    /* The range of chirp coefficient w indices is -n+1,-n+2,...,m-2,m-1.
     * The range of data d indices is 0,1,...,n-1.
     * We only need convolution result to be correct for indices 0,1,...,m-1.
     * This means we do not have to pad, i.e. the minimum transform length is m+n-1. */
    gint minsize = m + n - 1;
    gint size = gwy_fft_find_nice_size(minsize);
    fftw_complex *cs = g_new(fftw_complex, n);
    fftw_complex *x = gwy_fftw_new_complex(size);
    fftw_complex *w = gwy_fftw_new_complex(size);
    fftw_complex *fx = gwy_fftw_new_complex(size);
    fftw_complex *fw = gwy_fftw_new_complex(size);
    fftw_plan fplan, bplan;
    gdouble q, D = (f1 - f0)/(m - 1);
    gint k;

    fplan = gwy_fftw_plan_dft_1d(size, x, fw, FFTW_FORWARD, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    bplan = gwy_fftw_plan_dft_1d(size, fx, x, FFTW_BACKWARD, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    /* Precompute the factors w_k = exp(-2πik²/(ND)) */
    make_chirp(w, m, n, size, D);

    /* Transform premultiplied data. */
    make_cexp(cs, n, 2.0*G_PI*f0/n);
    gwycreal(x[0]) = rein[0];
    gwycimag(x[0]) = imin ? imin[0] : 0.0;
    for (k = 1; k < n; k++) {
        fftw_complex t, d;
        gwycreal(t) = cjmulre(cs[k], w[size - k]);
        gwycimag(t) = cjmulim(cs[k], w[size - k]);
        gwycreal(d) = rein[k];
        gwycimag(d) = imin ? imin[k] : 0.0;
        gwycreal(x[k]) = cmulre(t, d);
        gwycimag(x[k]) = cmulim(t, d);
        //x[k] = cexp(q*k*I)*conj(w[size - k])*d[k];
    }
    g_free(cs);
    gwy_clear(x + n, size - n);
    fftw_execute(fplan);
    gwy_assign(fx, fw, size);

    /* Transform chirp w. */
    gwy_assign(x, w, size);
    fftw_execute(fplan);
    fftw_destroy_plan(fplan);

    /* Multiply and transform back. */
    for (k = 0; k < size; k++) {
        gdouble tre = cmulre(fx[k], fw[k]), tim = cmulim(fx[k], fw[k]);
        gwycreal(fx[k]) = tre;
        gwycimag(fx[k]) = tim;
        //fx[k] *= fw[k];
    }
    fftw_execute(bplan);
    fftw_destroy_plan(bplan);
    fftw_free(fx);
    fftw_free(fw);

    /* And finally post-multiply by conj(w). */
    q = 1.0/size/sqrt(n);
    for (k = 0; k < m; k++) {
        reout[k] = q*cjmulre(x[k], w[k]);
        imout[k] = q*cjmulim(x[k], w[k]);
        //f[k] = q*x[k]*conj(w[k]);
    }
    fftw_free(x);
    fftw_free(w);
}

static void
zoom_fft_2d_do(const gdouble *rein, const gdouble *imin, gint nx, gint ny,
               gdouble *reout, gdouble *imout, gint mx, gint my,
               gdouble fx0, gdouble fy0, gdouble fx1, gdouble fy1)
{
    /* The range of chirp coefficient w indices is -n+1,-n+2,...,m-2,m-1.
     * The range of data d indices is 0,1,...,n-1.
     * We only need convolution result to be correct for indices 0,1,...,m-1.
     * This means we do not have to pad, i.e. the minimum transform length is m+n-1. */
    gint minxsize = mx + nx - 1;
    gint minysize = my + ny - 1;
    gint xsize = gwy_fft_find_nice_size(minxsize);
    gint ysize = gwy_fft_find_nice_size(minysize);
    gint size = xsize*ysize;
    fftw_complex *w_x = g_new(fftw_complex, xsize);
    fftw_complex *w_y = g_new(fftw_complex, ysize);
    fftw_complex *cs_x = g_new(fftw_complex, nx);
    fftw_complex *cs_y = g_new(fftw_complex, ny);
    fftw_complex *x = gwy_fftw_new_complex(size);
    fftw_complex *fx = gwy_fftw_new_complex(size);
    fftw_complex *fw = gwy_fftw_new_complex(size);
    fftw_plan fplan, bplan;
    gdouble q, Dx = (fx1 - fx0)/(mx - 1), Dy = (fy1 - fy0)/(my - 1);
    gint k, kx, ky;

    fplan = gwy_fftw_plan_dft_2d(xsize, ysize, x, fw, FFTW_FORWARD, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    bplan = gwy_fftw_plan_dft_2d(xsize, ysize, fx, x, FFTW_BACKWARD, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    /* Precompute the factors w_k = exp(-2πik²/(ND))
     * This part is entirely factorable. So we only compute 1D arrays along x and y and then basically make outer
     * products. */
    make_chirp(w_x, mx, nx, xsize, Dx);
    make_chirp(w_y, my, ny, ysize, Dy);

    /* Premultiply data.  We first combine cexp(qx*kx*I)*conj(wx[xsize - kx]) to a single x-factor, and similarly for
     * y. We put it to cs[] because it is needed only here, whereas w[] will be still needed later. */
    make_cexp(cs_x, nx, 2.0*G_PI*fx0/nx);
    make_cexp(cs_y, ny, 2.0*G_PI*fy0/ny);
    for (kx = 1; kx < nx; kx++) {
        fftw_complex t;
        gwycreal(t) = cjmulre(cs_x[kx], w_x[xsize - kx]);
        gwycimag(t) = cjmulim(cs_x[kx], w_x[xsize - kx]);
        gwycreal(cs_x[kx]) = gwycreal(t);
        gwycimag(cs_x[kx]) = gwycimag(t);
    }
    for (ky = 1; ky < ny; ky++) {
        fftw_complex t;
        gwycreal(t) = cjmulre(cs_y[ky], w_y[ysize - ky]);
        gwycimag(t) = cjmulim(cs_y[ky], w_y[ysize - ky]);
        gwycreal(cs_y[ky]) = gwycreal(t);
        gwycimag(cs_y[ky]) = gwycimag(t);
    }
    for (ky = 0; ky < ny; ky++) {
        for (kx = 0; kx < nx; kx++) {
            /* The arrays have different sizes! */
            gint kkin = nx*ky + kx;
            gint kkx = xsize*ky + kx;
            fftw_complex t, d;
            gwycreal(t) = cmulre(cs_x[kx], cs_y[ky]);
            gwycimag(t) = cmulim(cs_x[kx], cs_y[ky]);
            gwycreal(d) = rein[kkin];
            gwycimag(d) = imin ? imin[kkin] : 0.0;
            gwycreal(x[kkx]) = cmulre(t, d);
            gwycimag(x[kkx]) = cmulim(t, d);
            //x[ky,kx] = cexp(qx*kx*I)*cexp(qy*ky*I)*conj(wx[xsize - kx])*conj(wy[ysize - ky])*d[ky,kx];
        }
        gwy_clear(x + xsize*ky + nx, xsize - nx);
    }
    g_free(cs_x);
    g_free(cs_y);
    gwy_clear(x + xsize*ny, xsize*(ysize - ny));
    fftw_execute(fplan);
    gwy_assign(fx, fw, size);

    /* Transform chirp w (must form it as an outer product of w_x and w_y). */
    for (ky = 0; ky < ysize; ky++) {
        for (kx = 0; kx < xsize; kx++) {
            gint kk = xsize*ky + kx;
            gwycreal(x[kk]) = cmulre(w_x[kx], w_y[ky]);
            gwycimag(x[kk]) = cmulim(w_x[kx], w_y[ky]);
            //x[ky,kx] = w[ky,kx]
        }
    }
    fftw_execute(fplan);
    fftw_destroy_plan(fplan);

    /* Multiply and transform back. */
    for (k = 0; k < size; k++) {
        gdouble tre = cmulre(fx[k], fw[k]), tim = cmulim(fx[k], fw[k]);
        gwycreal(fx[k]) = tre;
        gwycimag(fx[k]) = tim;
        //fx[k] *= fw[k];
    }
    fftw_execute(bplan);
    fftw_destroy_plan(bplan);
    fftw_free(fx);
    fftw_free(fw);

    /* And finally post-multiply by conj(w). The arrays have different sizes; we extract only a part of the result. */
    q = 1.0/size/sqrt(nx*ny);
    for (ky = 0; ky < my; ky++) {
        for (kx = 0; kx < mx; kx++) {
            gint kkx = xsize*ky + kx;
            gint kkout = mx*ky + kx;
            fftw_complex t;
            gwycreal(t) = cmulre(w_x[kx], w_y[ky]);
            gwycimag(t) = cmulim(w_x[kx], w_y[ky]);
            reout[kkout] = q*cjmulre(x[kkx], t);
            imout[kkout] = q*cjmulim(x[kkx], t);
            //f[ky,kx] = q*x[ky,kx]*conj(w[ky,kx]);
        }
    }
    fftw_free(x);
    fftw_free(w_x);
    fftw_free(w_y);
}

static void
gwy_data_field_2dfft_prepare(GwyDataField *dfield,
                             gint level,
                             GwyWindowingType windowing,
                             gboolean preserverms,
                             gdouble *rms)
{
    gdouble a, bx, by;

    if (level == 2) {
        gwy_data_field_fit_plane(dfield, &a, &bx, &by);
        gwy_data_field_plane_level(dfield, a, bx, by);
    }
    else if (level == 1)
        gwy_data_field_add(dfield, -gwy_data_field_get_avg(dfield));
    if (preserverms) {
        a = gwy_data_field_get_rms(dfield);
        *rms = hypot(*rms, a);
    }
    gwy_data_field_fft_window(dfield, windowing);
}

/**
 * gwy_data_field_area_2dfft:
 * @rin: Real input data field.
 * @iin: Imaginary input data field.  It can be %NULL for real-to-complex transform which can be somewhat faster than
 *       complex-to-complex transform.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns), must be at least 2.
 * @height: Area height (number of rows), must be at least 2.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @interpolation: Interpolation type. Ignored since 2.8 as no resampling is performed.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract plane (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Calculates 2D Fast Fourier Transform of a rectangular area of a data field.
 *
 * If requested a windowing and/or leveling is applied to preprocess data to
 * obtain reasonable results.
 **/
void
gwy_data_field_area_2dfft(GwyDataField *rin, GwyDataField *iin,
                          GwyDataField *rout, GwyDataField *iout,
                          gint col, gint row,
                          gint width, gint height,
                          GwyWindowingType windowing,
                          GwyTransformDirection direction,
                          G_GNUC_UNUSED GwyInterpolationType interpolation,
                          gboolean preserverms, gint level)
{
    gint j, k;
    GwyDataField *rbuf, *ibuf;
    gdouble *out_rdata, *out_idata;
    gdouble rmsa = 0.0, rmsb;

    if (!iin) {
        gwy_data_field_area_2dfft_real(rin, rout, iout, col, row, width, height,
                                       windowing, direction, preserverms, level);
        return;
    }

    if (!_gwy_data_field_check_area(rin, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iin));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));
    g_return_if_fail(rin->xres == iin->xres && rin->yres == iin->yres);
    g_return_if_fail(level >= 0 && level <= 2);
    g_return_if_fail(width >= 2 && height >= 2);

    rbuf = gwy_data_field_area_extract(rin, col, row, width, height);
    gwy_data_field_2dfft_prepare(rbuf, level, windowing, preserverms, &rmsa);

    ibuf = gwy_data_field_area_extract(iin, col, row, width, height);
    gwy_data_field_2dfft_prepare(ibuf, level, windowing, preserverms, &rmsa);

    gwy_data_field_2dfft_raw(rbuf, ibuf, rout, iout, direction);

    if (preserverms) {
        out_idata = iout->data;
        out_rdata = rout->data;

        /* Ignore coefficient [0,0] */
        rmsb = -(out_rdata[0]*out_rdata[0] + out_idata[0]*out_idata[0]);
        for (j = 0; j < height; j++) {
            for (k = 0; k < width; k++)
                rmsb += (out_rdata[j*width + k]*out_rdata[j*width + k] + out_idata[j*width + k]*out_idata[j*width + k]);
        }
        rmsb = sqrt(rmsb)/(width*height);
        if (rmsb > 0.0) {
            gwy_data_field_multiply(rout, rmsa/rmsb);
            gwy_data_field_multiply(iout, rmsa/rmsb);
        }
    }

    g_object_unref(rbuf);
    g_object_unref(ibuf);

    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

static void
field_fft_2d_c2c(GwyDataField *rin, GwyDataField *iin,
                 GwyDataField *rout, GwyDataField *iout,
                 GwyTransformDirection direction)
{
    gint xres = rin->xres, yres = rin->yres;
    fftw_complex *in, *out;
    const gdouble *rindata = rin->data, *iindata = iin->data;
    gdouble *routdata = rout->data, *ioutdata = iout->data;
    fftw_plan plan;
    guint flags = FFTW_DESTROY_INPUT | FFTW_ESTIMATE;
    /* XXX: We have the sign reversed with respect to FFTW. */
    gint sign = (direction == GWY_TRANSFORM_DIRECTION_FORWARD ? FFTW_BACKWARD : FFTW_FORWARD);
    gdouble q;
    gint i, n;

    n = xres*yres;
    in = gwy_fftw_new_complex(n);
    out = gwy_fftw_new_complex(n);
    plan = gwy_fftw_plan_dft_2d(yres, xres, in, out, sign, flags);
    for (i = 0; i < n; i++) {
        in[i][0] = rindata[i];
        in[i][1] = iindata[i];
    }
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(in);
    q = 1.0/sqrt(n);
    for (i = 0; i < n; i++) {
        routdata[i] = q*out[i][0];
        ioutdata[i] = q*out[i][1];
    }
    fftw_free(out);

    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

static void
field_fft_2d_r2c(GwyDataField *rin,
                 GwyDataField *rout, GwyDataField *iout,
                 GwyTransformDirection direction)
{
    gint xres = rin->xres, yres = rin->yres, xres2 = xres/2;
    gint cstride = xres2 + 1;
    fftw_complex *out;
    gdouble *in, *routdata = rout->data, *ioutdata = iout->data;
    gdouble *rrow, *irow, *rrow2, *irow2;
    const fftw_complex *crow;
    fftw_plan plan;
    gdouble q;
    gint i, j;

    /* The planner may destroy input.  Use rout as a temporary input buffer. */
    in = routdata;
    out = gwy_fftw_new_complex(cstride*yres);
    plan = gwy_fftw_plan_dft_r2c_2d(yres, xres, in, out, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    gwy_assign(in, rin->data, xres*yres);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* Expand the R2C data to full-sized fields using the Hermitean symmetry. The zeroth row and column are not
     * mirrored; the central row and column might be (sort of), depending on parity.
     * XXX: We also have the sign reversed with respect to FFTW. */
    q = 1.0/sqrt(xres*yres);
    /* The zeroth row. */
    crow = out;
    rrow = routdata;
    irow = ioutdata;
    rrow[0] = q*crow[0][0];
    irow[0] = -q*crow[0][1];   /* Should be actually zero */
    if (xres % 2) {
        for (j = 1; j <= xres2; j++) {
            rrow[j] = rrow[xres-j] = q*crow[j][0];
            irow[j] = -q*crow[j][1];
            irow[xres-j] = q*crow[j][1];
        }
    }
    else {
        for (j = 1; j < xres2; j++) {
            rrow[j] = rrow[xres-j] = q*crow[j][0];
            irow[j] = -q*crow[j][1];
            irow[xres-j] = q*crow[j][1];
        }
        rrow[xres2] = q*crow[xres2][0];
        irow[xres2] = -q*crow[xres2][1];   /* Should be actually zero */
    }
    /* Remaining yres-1 rows. */
    for (i = 1; i < yres; i++) {
        rrow = routdata + i*xres;
        irow = ioutdata + i*xres;
        rrow2 = routdata + (yres-i)*xres;
        irow2 = ioutdata + (yres-i)*xres;
        crow = out + i*cstride;

        rrow[0] = q*crow[0][0];
        irow[0] = -q*crow[0][1];   /* Should be actually zero */
        if (xres % 2) {
            for (j = 1; j <= xres2; j++) {
                rrow[j] = rrow2[xres-j] = q*crow[j][0];
                irow[j] = -q*crow[j][1];
                irow2[xres-j] = q*crow[j][1];
            }
        }
        else {
            for (j = 1; j < xres2; j++) {
                rrow[j] = rrow2[xres-j] = q*crow[j][0];
                irow[j] = -q*crow[j][1];
                irow2[xres-j] = q*crow[j][1];
            }
            rrow[xres2] = q*crow[xres2][0];
            irow[xres2] = -q*crow[xres2][1];   /* Should be actually zero */
        }
    }
    fftw_free(out);

    /* Backward R2C is a silly case, but implement it, mainly because the API has always accepted both transform
     * directions.  Inverse transform is the complex conjugation of forward transform of conjugated input.  Since the
     * input is real, we just conjugate the output. */
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD) {
        for (i = 0; i < xres*yres; i++)
            ioutdata[i] = -ioutdata[i];
    }

    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

static void
field_fft_2d_c2r(GwyDataField *rin, GwyDataField *iin,
                 GwyDataField *rout,
                 GwyTransformDirection direction)
{
    gint xres = rin->xres, yres = rin->yres, xres2 = xres/2;
    gint cstride = xres2 + 1;
    fftw_complex *in;
    const gdouble *rindata = rin->data, *iindata = iin->data;
    gdouble *routdata = rout->data;
    fftw_plan plan;
    gdouble qr, qi;
    gint i, j;

    /* Forward C2R is a silly case, but implement it, mainly because the API has always accepted both transform
     * directions.  Inverse transform is the complex conjugation of forward transform of conjugated input.  Since the
     * output is real, we just conjugate the input. */
    qr = 1.0/sqrt(xres*yres);
    /* XXX: We have the sign reversed with respect to FFTW. */
    if (direction == GWY_TRANSFORM_DIRECTION_FORWARD)
        qi = -qr;
    else
        qi = qr;

    in = gwy_fftw_new_complex(cstride*yres);
    plan = gwy_fftw_plan_dft_c2r_2d(yres, xres, in, routdata, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    g_assert(plan);
    /* Use half of input fields, assuming the Hermitean symmetry.  Do not
     * attempt to enforce zeros in imaginary parts either.  */
    for (i = 0; i < yres; i++) {
        const gdouble *rrow = rindata + i*xres, *irow = iindata + i*xres;
        fftw_complex *crow = in + i*cstride;

        for (j = 0; j < cstride; j++) {
            crow[j][0] = qr*rrow[j];
            crow[j][1] = -qi*irow[j];
        }
    }
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(in);

    gwy_data_field_invalidate(rout);
}

/**
 * gwy_data_field_2dfft_raw:
 * @rin: Real input data field.
 * @iin: Imaginary input data field.  It can be %NULL for real-to-complex transform.
 * @rout: Real output data field, it will be resized to @rin size.
 * @iout: Imaginary output data field, it will be resized to @rin size.
 * @direction: FFT direction.  It should be %GWY_TRANSFORM_DIRECTION_FORWARD for real-to-complex transforms and
 *             %GWY_TRANSFORM_DIRECTION_BACKWARD for complex-to-real transforms.
 *
 * Calculates 2D Fast Fourier Transform of a data field.
 *
 * No leveling, windowing nor scaling is performed.
 *
 * The normalisation of FFT is symmetrical, so transformations in both directions are unitary.
 *
 * Since 2.8 the dimensions need not to be from the set of sizes returned by gwy_fft_find_nice_size().
 *
 * Lateral dimensions, offsets and units are unchanged.  See gwy_data_field_fft_postprocess() for that.
 *
 * Since 2.53 @iout can be %NULL for complex-to-real transforms.  Note that this means Hermitean symmetry of the input
 * data is assumed, i.e. about half of the input is ignored.  If you want to extract the real part of a complex
 * transform, you must pass a non-%NULL @iout.
 *
 * Since: 2.1
 **/
void
gwy_data_field_2dfft_raw(GwyDataField *rin,
                         GwyDataField *iin,
                         GwyDataField *rout,
                         GwyDataField *iout,
                         GwyTransformDirection direction)
{
    gint xres, yres;

    g_return_if_fail(GWY_IS_DATA_FIELD(rin));
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    xres = rin->xres;
    yres = rin->yres;

    g_return_if_fail(!iin || GWY_IS_DATA_FIELD(iin));
    if (iin) {
        g_return_if_fail(iin->xres == rin->xres);
        g_return_if_fail(iin->yres == rin->yres);
    }
    g_return_if_fail(!iout || GWY_IS_DATA_FIELD(iout));
    /* We could also special-case R2R transforms, but they are not commonly
     * needed in Gwyddion code. */
    g_return_if_fail(iin || iout);

    gwy_data_field_resample(rout, xres, yres, GWY_INTERPOLATION_NONE);
    if (iout)
        gwy_data_field_resample(iout, xres, yres, GWY_INTERPOLATION_NONE);

    if (iin && iout)
        field_fft_2d_c2c(rin, iin, rout, iout, direction);
    else if (iout)
        field_fft_2d_r2c(rin, rout, iout, direction);
    else
        field_fft_2d_c2r(rin, iin, rout, direction);
}

/**
 * gwy_data_field_area_2dfft_real:
 * @rin: Real input data field.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns), must be at least 2.
 * @height: Area height (number of rows), must be at least 2.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract plane (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Calculates 2D Fast Fourier Transform of a rectangular area of a data field.
 *
 * As the input is only real, the computation can be a somewhat faster than gwy_data_field_2dfft().
 **/
static void
gwy_data_field_area_2dfft_real(GwyDataField *rin,
                               GwyDataField *rout, GwyDataField *iout,
                               gint col, gint row,
                               gint width, gint height,
                               GwyWindowingType windowing,
                               GwyTransformDirection direction,
                               gboolean preserverms, gint level)
{
    GwyDataField *rbuf;
    gdouble *out_rdata, *out_idata;
    gdouble rmsa = 0.0, rmsb;

    if (!_gwy_data_field_check_area(rin, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));
    g_return_if_fail(width >= 2 && height >= 2);

    rbuf = gwy_data_field_area_extract(rin, col, row, width, height);
    gwy_data_field_2dfft_prepare(rbuf, level, windowing, preserverms, &rmsa);

    gwy_data_field_2dfft_raw(rbuf, NULL, rout, iout, direction);

    if (preserverms) {
        gint k;

        /* Ignore coefficient [0,0] */
        out_rdata = rout->data;
        out_idata = iout->data;
        rmsb = -(out_rdata[0]*out_rdata[0] + out_idata[0]*out_idata[0]);
        for (k = height*width; k; k--, out_rdata++, out_idata++)
            rmsb += (*out_rdata)*(*out_rdata) + (*out_idata)*(*out_idata);
        rmsb = sqrt(rmsb/(width*height));
        if (rmsb > 0.0) {
            gwy_data_field_multiply(rout, rmsa/rmsb);
            gwy_data_field_multiply(iout, rmsa/rmsb);
        }
    }

    g_object_unref(rbuf);
}

/**
 * gwy_data_field_2dfft:
 * @rin: Real input data field.
 * @iin: Imaginary input data field.  It can be %NULL for real-to-complex transform which can be somewhat faster than
 *       complex-to-complex transform.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @interpolation: Interpolation type.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract plane (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Calculates 2D Fast Fourier Transform of a rectangular a data field.
 *
 * If requested a windowing and/or leveling is applied to preprocess data to obtain reasonable results.
 *
 * Lateral dimensions, offsets and units are unchanged.  See gwy_data_field_fft_postprocess() for that.
 **/
void
gwy_data_field_2dfft(GwyDataField *rin, GwyDataField *iin,
                     GwyDataField *rout, GwyDataField *iout,
                     GwyWindowingType windowing,
                     GwyTransformDirection direction,
                     GwyInterpolationType interpolation,
                     gboolean preserverms, gint level)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(rin));

    if (!iin) {
        gwy_data_field_area_2dfft_real(rin, rout, iout, 0, 0, rin->xres, rin->yres,
                                       windowing, direction, preserverms, level);
    }
    else {
        gwy_data_field_area_2dfft(rin, iin, rout, iout, 0, 0, rin->xres, rin->yres,
                                  windowing, direction, interpolation, preserverms, level);
    }
}

/**
 * gwy_data_field_2dfft_humanize_in_place:
 * @data_field: A data field to (de)humanize.
 *
 * (De)humanizes a data field with Fourier coefficients in-place.
 *
 * This method can be only used for even-sized data fields and then it is an involutory operation.
 **/
static void
gwy_data_field_2dfft_humanize_in_place(GwyDataField *data_field)
{
    gint i, j, im, jm, xres, yres;
    gdouble *data;

    data = data_field->data;
    xres = data_field->xres;
    yres = data_field->yres;
    im = yres/2;
    jm = xres/2;

    for (i = 0; i < im; i++) {
        for (j = 0; j < jm; j++) {
            GWY_SWAP(gdouble, data[j + i*xres], data[(j + jm) + (i + im)*xres]);
            GWY_SWAP(gdouble, data[j + (i + im)*xres], data[(j + jm) + i*xres]);
        }
    }

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_2dfft_humanize:
 * @data_field: A data field.
 *
 * Rearranges 2D FFT output to a human-friendly form.
 *
 * Top-left, top-right, bottom-left and bottom-right sub-rectangles are swapped to obtain a humanized 2D FFT output
 * with (0,0) in the centre.
 *
 * More precisely, for even field dimensions the equally-sized blocks starting with the Nyquist frequency and with the
 * zero frequency (constant component) will exchange places.  For odd field dimensions, the block containing the zero
 * frequency is one item larger and the constant component will actually end up in the exact centre.
 *
 * Also note if both dimensions are even, this function is involutory and identical to
 * gwy_data_field_2dfft_dehumanize().  However, if any dimension is odd, gwy_data_field_2dfft_humanize() and
 * gwy_data_field_2dfft_dehumanize() are different, therefore they must be paired properly.
 **/
void
gwy_data_field_2dfft_humanize(GwyDataField *data_field)
{
    GwyDataField *tmp;
    gint im, jm, xres, yres;

    xres = data_field->xres;
    yres = data_field->yres;
    jm = data_field->xres/2;
    im = data_field->yres/2;

    if (xres == 2*jm && yres == 2*im) {
        gwy_data_field_2dfft_humanize_in_place(data_field);
        return;
    }

    tmp = gwy_data_field_new_alike(data_field, FALSE);
    gwy_data_field_area_copy(data_field, tmp, 0, 0, xres-jm, yres-im, jm, im);
    gwy_data_field_area_copy(data_field, tmp, xres-jm, 0, jm, yres-im, 0, im);
    gwy_data_field_area_copy(data_field, tmp, 0, yres-im, xres-jm, im, jm, 0);
    gwy_data_field_area_copy(data_field, tmp, xres-jm, yres-im, jm, im, 0, 0);
    gwy_data_field_copy(tmp, data_field, FALSE);
    g_object_unref(tmp);
}

/**
 * gwy_data_field_2dfft_dehumanize:
 * @data_field: A data field.
 *
 * Rearranges 2D FFT output back from the human-friendly form.
 *
 * Top-left, top-right, bottom-left and bottom-right sub-rectangles are swapped to reshuffle a humanized 2D FFT output
 * back into the natural positions.
 *
 * See gwy_data_field_2dfft_humanize() for discussion.
 *
 * Since: 2.8
 **/
void
gwy_data_field_2dfft_dehumanize(GwyDataField *data_field)
{
    GwyDataField *tmp;
    gint im, jm, xres, yres;

    xres = data_field->xres;
    yres = data_field->yres;
    jm = data_field->xres/2;
    im = data_field->yres/2;

    if (xres == 2*jm && yres == 2*im) {
        gwy_data_field_2dfft_humanize_in_place(data_field);
        return;
    }

    tmp = gwy_data_field_new_alike(data_field, FALSE);
    gwy_data_field_area_copy(data_field, tmp, 0, 0, jm, im, xres-jm, yres-im);
    gwy_data_field_area_copy(data_field, tmp, jm, 0, xres-jm, im, 0, yres-im);
    gwy_data_field_area_copy(data_field, tmp, 0, im, jm, yres-im, xres-jm, 0);
    gwy_data_field_area_copy(data_field, tmp, jm, im, xres-jm, yres-im, 0, 0);
    gwy_data_field_copy(tmp, data_field, FALSE);
    g_object_unref(tmp);
}

/**
 * gwy_data_field_fft_postprocess:
 * @data_field: A data field.
 * @humanize: %TRUE to rearrange data to have the frequency origin in the centre.
 *
 * Updates units, dimensions and offsets for a 2D FFT-processed field.
 *
 * The field is expected to have dimensions and units of the original direct-space data.  The lateral units and
 * resolutions are updated to correspond to its Fourier transform.
 *
 * The real dimensions are set for spatial frequencies, not wavevectors. For wavevector lateral coordinates, mutiply
 * all real dimensions and offsets by 2*%G_PI.
 *
 * If @humanize is %TRUE gwy_data_field_2dfft_humanize() is applied to the field data and the lateral offsets are set
 * accordingly.  Otherwise the offsets are cleared.
 *
 * Value units are kept intact.
 *
 * Since: 2.38
 **/
void
gwy_data_field_fft_postprocess(GwyDataField *dfield,
                               gboolean humanize)
{
    GwySIUnit *xyunit;
    gint res;
    gdouble r;

    xyunit = gwy_data_field_get_si_unit_xy(dfield);
    gwy_si_unit_power(xyunit, -1, xyunit);

    gwy_data_field_set_xreal(dfield, 1.0/gwy_data_field_get_dx(dfield));
    gwy_data_field_set_yreal(dfield, 1.0/gwy_data_field_get_dy(dfield));

    if (!humanize) {
        gwy_data_field_invalidate(dfield);
        gwy_data_field_set_xoffset(dfield, 0.0);
        gwy_data_field_set_yoffset(dfield, 0.0);
        return;
    }

    gwy_data_field_2dfft_humanize(dfield);

    res = dfield->xres;
    r = (res + 1 - res % 2)/2.0;
    gwy_data_field_set_xoffset(dfield, -gwy_data_field_jtor(dfield, r));

    res = dfield->yres;
    r = (res + 1 - res % 2)/2.0;
    gwy_data_field_set_yoffset(dfield, -gwy_data_field_itor(dfield, r));
}

/**
 * gwy_data_field_area_1dfft:
 * @rin: Real input data field.
 * @iin: Imaginary input data field.  It can be %NULL for real-to-complex transform which can be somewhat faster than
 *       complex-to-complex transform.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns), must be at least 2 for horizontal transforms.
 * @height: Area height (number of rows), must be at least 2 for vertical transforms.
 * @orientation: Orientation: pass %GWY_ORIENTATION_HORIZONTAL to transform rows, %GWY_ORIENTATION_VERTICAL to
 *               transform columns.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @interpolation: Interpolation type. Ignored since 2.8 as no resampling is performed.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract lines (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Transforms all rows or columns in a rectangular part of a data field with Fast Fourier Transform.
 *
 * If requested a windowing and/or leveling is applied to preprocess data to obtain reasonable results.
 **/
void
gwy_data_field_area_1dfft(GwyDataField *rin, GwyDataField *iin,
                          GwyDataField *rout, GwyDataField *iout,
                          gint col, gint row,
                          gint width, gint height,
                          GwyOrientation orientation,
                          GwyWindowingType windowing,
                          GwyTransformDirection direction,
                          G_GNUC_UNUSED GwyInterpolationType interpolation,
                          gboolean preserverms,
                          gint level)
{
    if (orientation == GWY_ORIENTATION_HORIZONTAL) {
        if (!iin) {
            gwy_data_field_area_xfft_real(rin, rout, iout, col, row, width, height,
                                          windowing, direction, preserverms, level);
        }
        else {
            gwy_data_field_area_xfft(rin, iin, rout, iout, col, row, width, height,
                                     windowing, direction, preserverms, level);
        }
    }
    else if (orientation == GWY_ORIENTATION_VERTICAL) {
        if (!iin) {
            gwy_data_field_area_yfft_real(rin, rout, iout, col, row, width, height,
                                          windowing, direction, preserverms, level);
        }
        else {
            gwy_data_field_area_yfft(rin, iin, rout, iout, col, row, width, height,
                                     windowing, direction, preserverms, level);
        }
    }
    else {
        g_assert_not_reached();
    }
}

/**
 * gwy_data_field_1dfft:
 * @rin: Real input data field.
 * @iin: Imaginary input data field.  It can be %NULL for real-to-complex transform which can be somewhat faster than
 *       complex-to-complex transform.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @orientation: Orientation: pass %GWY_ORIENTATION_HORIZONTAL to transform rows, %GWY_ORIENTATION_VERTICAL to
 *               transform columns.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @interpolation: Interpolation type. Ignored since 2.8 as no resampling is performed.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract line (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Transforms all rows or columns in a data field with Fast Fourier Transform.
 *
 * If requested a windowing and/or leveling is applied to preprocess data to obtain reasonable results.
 **/
void
gwy_data_field_1dfft(GwyDataField *rin, GwyDataField *iin,
                     GwyDataField *rout, GwyDataField *iout,
                     GwyOrientation orientation,
                     GwyWindowingType windowing,
                     GwyTransformDirection direction,
                     G_GNUC_UNUSED GwyInterpolationType interpolation,
                     gboolean preserverms,
                     gint level)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(rin));

    if (orientation == GWY_ORIENTATION_HORIZONTAL) {
        if (!iin) {
            gwy_data_field_area_xfft_real(rin, rout, iout, 0, 0, rin->xres, rin->yres,
                                          windowing, direction, preserverms, level);
        }
        else {
            gwy_data_field_area_xfft(rin, iin, rout, iout, 0, 0, rin->xres, rin->yres,
                                     windowing, direction, preserverms, level);
        }
    }
    else if (orientation == GWY_ORIENTATION_VERTICAL) {
        if (!iin) {
            gwy_data_field_area_yfft_real(rin, rout, iout, 0, 0, rin->xres, rin->yres,
                                          windowing, direction, preserverms, level);
        }
        else {
            gwy_data_field_area_yfft(rin, iin, rout, iout, 0, 0, rin->xres, rin->yres,
                                     windowing, direction, preserverms, level);
        }
    }
    else {
        g_assert_not_reached();
    }
}

/**
 * gwy_data_field_1dfft_raw:
 * @rin: Real input data field.
 * @iin: Imaginary input data field.  It can be %NULL for real-to-complex transform.
 * @rout: Real output data field, it will be resized to @rin size.
 * @iout: Imaginary output data field, it will be resized to @rin size.
 * @orientation: Orientation: pass %GWY_ORIENTATION_HORIZONTAL to transform rows, %GWY_ORIENTATION_VERTICAL to
 *               transform columns.
 * @direction: FFT direction.
 *
 * Transforms all rows or columns in a data field with Fast Fourier Transform.
 *
 * No leveling, windowing nor scaling is performed.
 *
 * The normalisation of FFT is symmetrical, so transformations in both directions are unitary.
 *
 * Since 2.8 the dimensions need not to be from the set of sizes returned by gwy_fft_find_nice_size().
 *
 * Since: 2.1
 **/
void
gwy_data_field_1dfft_raw(GwyDataField *rin,
                         GwyDataField *iin,
                         GwyDataField *rout,
                         GwyDataField *iout,
                         GwyOrientation orientation,
                         GwyTransformDirection direction)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(rin));
    g_return_if_fail(!iin || GWY_IS_DATA_FIELD(iin));
    if (iin) {
        g_return_if_fail(!gwy_data_field_check_compatibility(rin, iin, GWY_DATA_COMPATIBILITY_RES));
    }
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));

    gwy_data_field_resample(rout, rin->xres, rin->yres, GWY_INTERPOLATION_NONE);
    gwy_data_field_resample(iout, rin->xres, rin->yres, GWY_INTERPOLATION_NONE);
    if (orientation == GWY_ORIENTATION_HORIZONTAL) {
        if (iin)
            gwy_data_field_xfft_do(rin, iin, rout, iout, direction);
        else {
            iin = gwy_data_field_new_alike(rin, FALSE);
            gwy_data_field_xfft_real_do(rin, iin, rout, iout, direction);
            g_object_unref(iin);
        }
    }
    else if (orientation == GWY_ORIENTATION_VERTICAL) {
        if (iin)
            gwy_data_field_yfft_do(rin, iin, rout, iout, direction);
        else {
            iin = gwy_data_field_new_alike(rin, FALSE);
            gwy_data_field_yfft_real_do(rin, iin, rout, iout, direction);
            g_object_unref(iin);
        }
    }
    else {
        g_assert_not_reached();
    }
}

/**
 * gwy_data_field_area_xfft:
 * @ra: Real input data field.
 * @ia: Imaginary input data field.  It can be %NULL for real-to-complex transform which can be somewhat faster than
 *      complex-to-complex transform.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns), must be at least 2.
 * @height: Area height (number of rows).
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract
 *         lines (the number can be interpreted as the first polynomial degree
 *         to keep, but only the enumerated three values are available).
 *
 * Transforms all rows in a data field with Fast Fourier Transform.
 *
 * If requested a windowing and/or leveling is applied to preprocess data to
 * obtain reasonable results.
 **/
static void
gwy_data_field_area_xfft(GwyDataField *rin, GwyDataField *iin,
                         GwyDataField *rout, GwyDataField *iout,
                         gint col, gint row,
                         gint width, gint height,
                         GwyWindowingType windowing,
                         GwyTransformDirection direction,
                         gboolean preserverms, gint level)
{
    gint k;
    GwyDataField *rbuf, *ibuf;

    if (!_gwy_data_field_check_area(rin, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iin));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));
    g_return_if_fail(rin->xres == iin->xres && rin->yres == rout->yres);
    g_return_if_fail(level >= 0 && level <= 2);
    g_return_if_fail(width >= 2);

    gwy_data_field_resample(rout, width, height, GWY_INTERPOLATION_NONE);
    gwy_data_field_resample(iout, width, height, GWY_INTERPOLATION_NONE);

    rbuf = gwy_data_field_area_extract(rin, col, row, width, height);
    if (level) {
        for (k = 0; k < height; k++)
            gwy_level_simple(width, 1, rbuf->data + k*width, level);
    }
    gwy_data_field_fft_window_1d(rbuf, GWY_ORIENTATION_HORIZONTAL, windowing);

    ibuf = gwy_data_field_area_extract(iin, col, row, width, height);
    if (level) {
        for (k = 0; k < height; k++)
            gwy_level_simple(width, 1, ibuf->data + k*width, level);
    }
    gwy_data_field_fft_window_1d(ibuf, GWY_ORIENTATION_HORIZONTAL, windowing);

    gwy_data_field_xfft_do(rbuf, ibuf, rout, iout, direction);

    if (preserverms) {
        for (k = 0; k < height; k++)
            gwy_preserve_rms_simple(width, 1,
                                    rin->data + rin->xres*(row + k) + col,
                                    iin->data + iin->xres*(row + k) + col,
                                    width, 1,
                                    rout->data + width*k,
                                    iout->data + width*k);
    }

    g_object_unref(rbuf);
    g_object_unref(ibuf);
}

static void
gwy_data_field_xfft_do(GwyDataField *rin,
                       GwyDataField *iin,
                       GwyDataField *rout,
                       GwyDataField *iout,
                       GwyTransformDirection direction)
{
    fftw_iodim dims[1], howmany_dims[1];
    fftw_plan plan;

    dims[0].n = rin->xres;
    dims[0].is = 1;
    dims[0].os = 1;
    howmany_dims[0].n = rin->yres;
    howmany_dims[0].is = rin->xres;
    howmany_dims[0].os = rin->xres;
    /* Backward direction is equivalent to switching real and imaginary parts */
    /* XXX: Planner destroys input, we have to either allocate memory or use in-place transform.  In some cases caller
     * could provide us with already allocated buffers. */
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD) {
        plan = gwy_fftw_plan_guru_split_dft(1, dims, 1, howmany_dims, rout->data, iout->data, rout->data, iout->data,
                                            FFTW_ESTIMATE);
    }
    else {
        plan = gwy_fftw_plan_guru_split_dft(1, dims, 1, howmany_dims, iout->data, rout->data, iout->data, rout->data,
                                            FFTW_ESTIMATE);
    }
    g_return_if_fail(plan);
    gwy_data_field_copy(rin, rout, FALSE);
    gwy_data_field_copy(iin, iout, FALSE);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);

    gwy_data_field_multiply(rout, 1.0/sqrt(rin->xres));
    gwy_data_field_multiply(iout, 1.0/sqrt(rin->xres));
    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

/**
 * gwy_data_field_area_xfft_real:
 * @rin: Real input data field.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns), must be at least 2.
 * @height: Area height (number of rows).
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract lines (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Transforms all rows in a data real field with Fast Fourier Transform.
 *
 * As the input is only real, the computation can be a somewhat faster than gwy_data_field_xfft().
 **/
static void
gwy_data_field_area_xfft_real(GwyDataField *rin, GwyDataField *rout,
                              GwyDataField *iout,
                              gint col, gint row,
                              gint width, gint height,
                              GwyWindowingType windowing,
                              GwyTransformDirection direction,
                              gboolean preserverms, gint level)
{
    gint k;
    GwyDataField *rbuf, *ibuf;

    if (!_gwy_data_field_check_area(rin, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));
    g_return_if_fail(level >= 0 && level <= 2);
    g_return_if_fail(width >= 2);

    gwy_data_field_resample(rout, width, height, GWY_INTERPOLATION_NONE);
    gwy_data_field_resample(iout, width, height, GWY_INTERPOLATION_NONE);

    rbuf = gwy_data_field_area_extract(rin, col, row, width, height);
    if (level) {
        for (k = 0; k < height; k++)
            gwy_level_simple(width, 1, rbuf->data + k*width, level);
    }
    gwy_data_field_fft_window_1d(rbuf, GWY_ORIENTATION_HORIZONTAL, windowing);

    ibuf = gwy_data_field_new_alike(rbuf, FALSE);

    gwy_data_field_xfft_real_do(rbuf, ibuf, rout, iout, direction);

    if (preserverms) {
        for (k = 0; k < height; k++) {
            gwy_preserve_rms_simple(width, 1,
                                    rin->data + rin->xres*(row + k) + col, NULL,
                                    width, 1,
                                    rout->data + width*k,
                                    iout->data + width*k);
        }
    }

    g_object_unref(rbuf);
    g_object_unref(ibuf);
}

static void
gwy_data_field_xfft_real_do(GwyDataField *rin,
                            GwyDataField *ibuf,
                            GwyDataField *rout,
                            GwyDataField *iout,
                            GwyTransformDirection direction)
{
    fftw_iodim dims[1], howmany_dims[1];
    fftw_plan plan;
    gint j, k;

    dims[0].n = rin->xres;
    dims[0].is = 1;
    dims[0].os = 1;
    howmany_dims[0].n = rin->yres;
    howmany_dims[0].is = rin->xres;
    howmany_dims[0].os = rin->xres;
    plan = gwy_fftw_plan_guru_split_dft_r2c(1, dims, 1, howmany_dims, ibuf->data, rout->data, iout->data,
                                            FFTW_ESTIMATE);
    /* R2C destroys input, and especially, the planner destroys input too */
    gwy_data_field_copy(rin, ibuf, FALSE);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* Complete the missing half of transform.  */
    for (k = 0; k < rin->yres; k++) {
        gdouble *re, *im;

        re = rout->data + k*rin->xres;
        im = iout->data + k*rin->xres;
        for (j = rin->xres/2 + 1; j < rin->xres; j++) {
            re[j] = re[rin->xres - j];
            im[j] = -im[rin->xres - j];
        }
    }

    gwy_data_field_multiply(rout, 1.0/sqrt(rin->xres));
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD)
        gwy_data_field_multiply(iout, 1.0/sqrt(rin->xres));
    else
        gwy_data_field_multiply(iout, -1.0/sqrt(rin->xres));
    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

/**
 * gwy_data_field_area_yfft:
 * @ra: Real input data field.
 * @ia: Imaginary input data field.  It can be %NULL for real-to-complex transform which can be somewhat faster than
 *      complex-to-complex transform.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows), must be at least 2.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract lines (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Transforms all columns in a data field with Fast Fourier Transform.
 *
 * If requested a windowing and/or leveling is applied to preprocess data to
 * obtain reasonable results.
 **/
static void
gwy_data_field_area_yfft(GwyDataField *rin, GwyDataField *iin,
                         GwyDataField *rout, GwyDataField *iout,
                         gint col, gint row,
                         gint width, gint height,
                         GwyWindowingType windowing,
                         GwyTransformDirection direction,
                         gboolean preserverms, gint level)
{
    gint k;
    GwyDataField *rbuf, *ibuf;

    if (!_gwy_data_field_check_area(rin, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iin));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));
    g_return_if_fail(rin->xres == iin->xres && rin->yres == rout->yres);
    g_return_if_fail(level >= 0 && level <= 2);
    g_return_if_fail(height >= 2);

    gwy_data_field_resample(rout, width, height, GWY_INTERPOLATION_NONE);
    gwy_data_field_resample(iout, width, height, GWY_INTERPOLATION_NONE);

    rbuf = gwy_data_field_area_extract(rin, col, row, width, height);
    if (level) {
        for (k = 0; k < width; k++)
            gwy_level_simple(height, width, rbuf->data + k, level);
    }
    gwy_data_field_fft_window_1d(rbuf, GWY_ORIENTATION_VERTICAL, windowing);

    ibuf = gwy_data_field_area_extract(iin, col, row, width, height);
    if (level) {
        for (k = 0; k < width; k++)
            gwy_level_simple(height, width, ibuf->data + k, level);
    }
    gwy_data_field_fft_window_1d(ibuf, GWY_ORIENTATION_VERTICAL, windowing);

    gwy_data_field_yfft_do(rbuf, ibuf, rout, iout, direction);

    if (preserverms) {
        for (k = 0; k < width; k++) {
            gwy_preserve_rms_simple(height, rin->xres,
                                    rin->data + rin->xres*row + col + k,
                                    iin->data + iin->xres*row + col + k,
                                    height, width,
                                    rout->data + k,
                                    iout->data + k);
        }
    }

    g_object_unref(rbuf);
    g_object_unref(ibuf);
}

static void
gwy_data_field_yfft_do(GwyDataField *rin,
                       GwyDataField *iin,
                       GwyDataField *rout,
                       GwyDataField *iout,
                       GwyTransformDirection direction)
{
    fftw_iodim dims[1], howmany_dims[1];
    fftw_plan plan;

    dims[0].n = rin->yres;
    dims[0].is = rin->xres;
    dims[0].os = rin->xres;
    howmany_dims[0].n = rin->xres;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 1;
    /* Backward direction is equivalent to switching real and imaginary parts */
    /* XXX: Planner destroys input, we have to either allocate memory or use in-place transform.  In some cases caller
     * could provide us with already allocated buffers. */
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD) {
        plan = gwy_fftw_plan_guru_split_dft(1, dims, 1, howmany_dims, rout->data, iout->data, rout->data, iout->data,
                                            FFTW_ESTIMATE);
    }
    else {
        plan = gwy_fftw_plan_guru_split_dft(1, dims, 1, howmany_dims, iout->data, rout->data, iout->data, rout->data,
                                            FFTW_ESTIMATE);
    }
    g_return_if_fail(plan);
    gwy_data_field_copy(rin, rout, FALSE);
    gwy_data_field_copy(iin, iout, FALSE);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);

    gwy_data_field_multiply(rout, 1.0/sqrt(rin->yres));
    gwy_data_field_multiply(iout, 1.0/sqrt(rin->yres));
    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

/**
 * gwy_data_field_area_yfft_real:
 * @ra: Real input data field.
 * @rout: Real output data field, it will be resized to area size.
 * @iout: Imaginary output data field, it will be resized to area size.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows), must be at least 2.
 * @windowing: Windowing type.
 * @direction: FFT direction.
 * @preserverms: %TRUE to preserve RMS while windowing.
 * @level: 0 to perform no leveling, 1 to subtract mean value, 2 to subtract lines (the number can be interpreted as
 *         the first polynomial degree to keep, but only the enumerated three values are available).
 *
 * Transforms all columns in a data real field with Fast Fourier Transform.
 *
 * As the input is only real, the computation can be a somewhat faster than gwy_data_field_yfft().
 **/
static void
gwy_data_field_area_yfft_real(GwyDataField *rin, GwyDataField *rout,
                              GwyDataField *iout,
                              gint col, gint row,
                              gint width, gint height,
                              GwyWindowingType windowing,
                              GwyTransformDirection direction,
                              gboolean preserverms, gint level)
{
    gint k;
    GwyDataField *rbuf, *ibuf;

    if (!_gwy_data_field_check_area(rin, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(rout));
    g_return_if_fail(GWY_IS_DATA_FIELD(iout));
    g_return_if_fail(level >= 0 && level <= 2);
    g_return_if_fail(height >= 2);

    gwy_data_field_resample(rout, width, height, GWY_INTERPOLATION_NONE);
    gwy_data_field_resample(iout, width, height, GWY_INTERPOLATION_NONE);

    rbuf = gwy_data_field_area_extract(rin, col, row, width, height);
    if (level) {
        for (k = 0; k < width; k++)
            gwy_level_simple(height, width, rbuf->data + k, level);
    }
    gwy_data_field_fft_window_1d(rbuf, GWY_ORIENTATION_VERTICAL, windowing);

    ibuf = gwy_data_field_new_alike(rbuf, FALSE);

    gwy_data_field_yfft_real_do(rbuf, ibuf, rout, iout, direction);

    if (preserverms) {
        for (k = 0; k < width; k++) {
            gwy_preserve_rms_simple(height, rin->xres,
                                    rin->data + rin->xres*row + col + k, NULL,
                                    height, width,
                                    rout->data + k,
                                    iout->data + k);
        }
    }

    g_object_unref(rbuf);
    g_object_unref(ibuf);
}

static void
gwy_data_field_yfft_real_do(GwyDataField *rin,
                            GwyDataField *ibuf,
                            GwyDataField *rout,
                            GwyDataField *iout,
                            GwyTransformDirection direction)
{
    fftw_iodim dims[1], howmany_dims[1];
    fftw_plan plan;
    gint j, k;

    dims[0].n = rin->yres;
    dims[0].is = rin->xres;
    dims[0].os = rin->xres;
    howmany_dims[0].n = rin->xres;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 1;
    plan = gwy_fftw_plan_guru_split_dft_r2c(1, dims, 1, howmany_dims, ibuf->data, rout->data, iout->data,
                                            FFTW_ESTIMATE);
    /* R2C destroys input, and especially, the planner destroys input too */
    gwy_data_field_copy(rin, ibuf, FALSE);
    gwy_fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* Complete the missing half of transform.  */
    for (k = 0; k < rin->xres; k++) {
        gdouble *re, *im;

        re = rout->data + k;
        im = iout->data + k;
        for (j = rin->yres/2 + 1; j < rin->yres; j++) {
            re[rin->xres*j] = re[rin->xres*(rin->yres - j)];
            im[rin->xres*j] = -im[rin->xres*(rin->yres - j)];
        }
    }

    gwy_data_field_multiply(rout, 1.0/sqrt(rin->yres));
    if (direction == GWY_TRANSFORM_DIRECTION_BACKWARD)
        gwy_data_field_multiply(iout, 1.0/sqrt(rin->yres));
    else
        gwy_data_field_multiply(iout, -1.0/sqrt(rin->yres));
    gwy_data_field_invalidate(rout);
    gwy_data_field_invalidate(iout);
}

static void
gwy_level_simple(gint n,
                 gint stride,
                 gdouble *data,
                 gint level)
{
    gdouble sumxi, sumxixi, sumsi, sumsixi, a, b;
    gdouble *pdata;
    gint i;

    level = MIN(level, n);

    if (!level)
        return;

    if (level == 1) {
        sumsi = 0.0;
        pdata = data;
        for (i = n; i; i--, pdata += stride)
            sumsi += *pdata;

        a = sumsi/n;
        pdata = data;
        for (i = n; i; i--, pdata += stride)
            *pdata -= a;

        return;
    }

    g_return_if_fail(level == 2);

    /* These are already averages, not sums */
    sumxi = (n + 1.0)/2.0;
    sumxixi = (2.0*n + 1.0)*(n + 1.0)/6.0;

    sumsi = sumsixi = 0.0;

    pdata = data;
    for (i = n; i; i--, pdata += stride) {
        sumsi += *pdata;
        sumsixi += *pdata * i;
    }
    sumsi /= n;
    sumsixi /= n;

    b = (sumsixi - sumsi*sumxi)/(sumxixi - sumxi*sumxi);
    a = (sumsi*sumxixi - sumxi*sumsixi)/(sumxixi - sumxi*sumxi);

    pdata = data;
    sumsi = 0;
    for (i = n; i; i--, pdata += stride) {
        *pdata -= a + b*i;
        sumsi += *pdata;
    }
}

static void
gwy_preserve_rms_simple(gint nsrc,
                        gint stridesrc,
                        const gdouble *src1,
                        const gdouble *src2,
                        gint ndata,
                        gint stridedata,
                        gdouble *data1,
                        gdouble *data2)
{
    gdouble sum2, sum0, sum02, a, b, q;
    gdouble *pdata;
    gint i;

    /* Calculate original RMS */
    sum0 = sum02 = 0.0;
    for (i = nsrc; i; i--, src1 += stridesrc) {
        sum0 += *src1;
        sum02 += *src1 * *src1;
    }
    a = sum02 - sum0*sum0/nsrc;
    if (src2) {
        sum0 = sum02 = 0.0;
        for (i = nsrc; i; i--, src1 += stridesrc) {
            sum0 += *src2;
            sum02 += *src2 * *src2;
        }
        a += sum02 - sum0*sum0/nsrc;
    }
    if (a <= 0.0)
        return;
    a = sqrt(a/nsrc);

    /* Calculare new RMS ignoring 0th elements that correspond to constants */
    sum2 = 0.0;
    for (i = ndata-1, pdata = data1 + 1; i; i--, pdata += stridedata)
        sum2 += *pdata * *pdata;
    for (i = ndata-1, pdata = data2 + 1; i; i--, pdata += stridedata)
        sum2 += *pdata * *pdata;
    if (sum2 == 0.0)
        return;
    b = sqrt(sum2/ndata);

    /* Multiply output to get the same RMS */
    q = a/b;
    for (i = ndata, pdata = data1; i; i--, pdata += stridedata)
        *pdata *= q;
    for (i = ndata, pdata = data2; i; i--, pdata += stridedata)
        *pdata *= q;
}

static GwyDataLine*
resample_dline_for_1d_fft_filter(GwyDataLine *dline, gint res,
                                 GwyInterpolationType interpolation)
{
    GwyDataLine *half = gwy_data_line_new_resampled(dline, (res + 1)/2, interpolation);
    GwyDataLine *full = gwy_data_line_new(res, res, FALSE);
    gint i;

    /* Fill the full line symmetrically.  The central element may we written twice, but with the same value. */
    gwy_assign(full->data, half->data, half->res);
    for (i = 0; i < half->res; i++)
        full->data[res-1 - i] = half->data[i];

    g_object_unref(half);
    return full;
}

/**
 * gwy_data_field_fft_filter_1d:
 * @data_field: A data field to filter.
 * @result_field: A data field to store the result to.  It will be resampled to @data_field's size.
 * @weights: Filter weights for the lower half of the spectrum (the other half is symmetric).  Its size can be
 *           arbitrary, it will be interpolated.
 * @orientation: Filter direction.
 * @interpolation: The interpolation to use for resampling.
 *
 * Performs 1D FFT filtering of a data field.
 **/
void
gwy_data_field_fft_filter_1d(GwyDataField *data_field,
                             GwyDataField *result_field,
                             GwyDataLine *weights,
                             GwyOrientation orientation,
                             GwyInterpolationType interpolation)
{
    GwyDataField *iresult_field, *hlp_rdfield, *hlp_idfield;
    GwyDataLine *w;
    gint i, j, xres, yres;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(result_field));
    g_return_if_fail(GWY_IS_DATA_LINE(weights));

    yres = data_field->yres;
    xres = data_field->xres;
    gwy_data_field_resample(result_field, xres, yres, GWY_INTERPOLATION_NONE);

    hlp_rdfield = gwy_data_field_new_alike(data_field, TRUE);
    hlp_idfield = gwy_data_field_new_alike(data_field, TRUE);
    iresult_field = gwy_data_field_new_alike(data_field, TRUE);

    gwy_data_field_1dfft_raw(data_field, NULL, hlp_rdfield, hlp_idfield,
                             orientation, GWY_TRANSFORM_DIRECTION_FORWARD);

    if (orientation == GWY_ORIENTATION_VERTICAL)
        w = resample_dline_for_1d_fft_filter(weights, yres, interpolation);
    else
        w = resample_dline_for_1d_fft_filter(weights, xres, interpolation);

    for (i = 0; i < yres; i++) {
        gdouble *rrow = hlp_rdfield->data + i*xres;
        gdouble *irow = hlp_idfield->data + i*xres;

        if (orientation == GWY_ORIENTATION_VERTICAL) {
            gdouble wi = w->data[i];
            for (j = 0; j < xres; j++) {
                rrow[j] *= wi;
                irow[j] *= wi;
            }
        }
        else {
            gdouble *wrow = w->data;
            for (j = 0; j < xres; j++) {
                rrow[j] *= wrow[j];
                irow[j] *= wrow[j];
            }
        }
    }
    g_object_unref(w);

    gwy_data_field_1dfft_raw(hlp_rdfield, hlp_idfield, result_field, iresult_field,
                             orientation, GWY_TRANSFORM_DIRECTION_BACKWARD);
    g_object_unref(iresult_field);
    g_object_unref(hlp_rdfield);
    g_object_unref(hlp_idfield);
}

/**
 * gwy_data_line_fft_window:
 * @line: A one-dimensional data line.
 * @windowing: Windowing type to use.
 *
 * Performs windowing of a data line in preparation for FFT.
 *
 * Since: 2.62
 **/
void
gwy_data_line_fft_window(GwyDataLine *line,
                         GwyWindowingType windowing)
{
    g_return_if_fail(GWY_IS_DATA_LINE(line));
    g_return_if_fail(windowing <= GWY_WINDOWING_KAISER25);
    if (windowing <= GWY_WINDOWING_NONE)
        return;
    gwy_fft_window(line->res, line->data, windowing);
}

/**
 * gwy_data_field_fft_window:
 * @field: A two-dimensional data field.
 * @windowing: Windowing type to use.
 *
 * Performs two-dimensional windowing of a data field in preparation for 2D FFT.
 *
 * The same windowing function is used row-wise and column-wise.
 *
 * Since: 2.62
 **/
void
gwy_data_field_fft_window(GwyDataField *field,
                          GwyWindowingType windowing)
{
    gint xres, yres, j, i;
    gdouble *data, *wx = NULL, *wy = NULL;

    g_return_if_fail(GWY_IS_DATA_FIELD(field));
    g_return_if_fail(windowing <= GWY_WINDOWING_KAISER25);
    if (windowing <= GWY_WINDOWING_NONE)
        return;

    xres = field->xres;
    yres = field->yres;
    data = field->data;
    if (xres == yres) {
        wx = wy = g_new(gdouble, xres);
        for (i = 0; i < xres; i++)
            wx[i] = 1.0;
        gwy_fft_window(xres, wx, windowing);
    }
    else {
        wx = g_new(gdouble, xres + yres);
        wy = wx + xres;
        for (i = 0; i < xres + yres; i++)
            wx[i] = 1.0;
        gwy_fft_window(xres, wx, windowing);
        gwy_fft_window(yres, wy, windowing);
    }

    for (i = 0; i < yres; i++) {
        gdouble qy = wy[i];
        for (j = 0; j < xres; j++)
            data[i*xres + j] *= qy*wx[j];
    }

    g_free(wx);
    gwy_data_field_invalidate(field);
}

/**
 * gwy_data_field_fft_window_1d:
 * @field: A two-dimensional data field.
 * @orientation: Windowing orientation (the same as corresponding FFT orientation).
 * @windowing: Windowing type to use.
 *
 * Performs row-wise or column-wise windowing of a data field in preparation for 1D FFT.
 *
 * Since: 2.62
 **/
void
gwy_data_field_fft_window_1d(GwyDataField *field,
                             GwyOrientation orientation,
                             GwyWindowingType windowing)
{
    gint xres, yres, wres, j, i;
    gdouble *data, *w = NULL;

    g_return_if_fail(GWY_IS_DATA_FIELD(field));
    g_return_if_fail(windowing <= GWY_WINDOWING_KAISER25);
    g_return_if_fail(orientation == GWY_ORIENTATION_HORIZONTAL || orientation == GWY_ORIENTATION_VERTICAL);
    if (windowing <= GWY_WINDOWING_NONE)
        return;

    xres = field->xres;
    yres = field->yres;
    data = field->data;
    wres = (orientation == GWY_ORIENTATION_VERTICAL ? yres : xres);
    w = g_new(gdouble, wres);
    for (i = 0; i < wres; i++)
        w[i] = 1.0;
    gwy_fft_window(wres, w, windowing);

    if (orientation == GWY_ORIENTATION_VERTICAL) {
        for (i = 0; i < yres; i++) {
            gdouble q = w[i];
            for (j = 0; j < xres; j++)
                data[i*xres + j] *= q;
        }
    }
    else {
        for (i = 0; i < yres; i++) {
            for (j = 0; j < xres; j++)
                data[i*xres + j] *= w[j];
        }
    }

    g_free(w);
    gwy_data_field_invalidate(field);
}

/************************** Documentation ****************************/

/**
 * SECTION:inttrans
 * @title: inttrans
 * @short_description: FFT and other integral transforms
 *
 * There are two main groups of FFT functions.
 *
 * High-level functions such as gwy_data_field_2dfft(), gwy_data_line_fft() can perform windowing, leveling and other
 * pre- and postprocessing. This makes them suitable for calculation of spectral densities and other statistical
 * characteristics.
 *
 * Low-level functions have <literal>raw</literal> appended to their name: gwy_data_field_2dfft_raw(),
 * gwy_data_line_fft_raw().  They perform no other operations on the data beside the transform itself. This makes them
 * suitable for applications where both forward and inverse transform is performed.
 *
 * Both types of functions wrap <ulink url="http://fftw.org/">FFTW3</ulink> routines.
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
