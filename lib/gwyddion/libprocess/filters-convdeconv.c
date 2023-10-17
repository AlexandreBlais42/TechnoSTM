/*
 *  $Id: filters-convdeconv.c 24753 2022-03-29 13:41:56Z yeti-dn $
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
#include <libprocess/filters.h>
#include <libprocess/stats.h>
#include <libprocess/linestats.h>
#include <libprocess/grains.h>
#include <libprocess/inttrans.h>
#include <libprocess/simplefft.h>
#include <libprocess/arithmetic.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"
#include "gwyfftw.h"

typedef struct {
    /* Image-sized. */
    fftw_complex *fideal;      /* Precalculated, for P = ... */
    fftw_complex *fmeas;       /* Precalculated, for P = ... */
    fftw_complex *fpsf;        /* P = ..., in each step */
    fftw_complex *cbuf;        /* FFT scratch buffer*/
    GwyDataField *psf;         /* Real-space PSF */
    fftw_plan pplan;
    gdouble sigma_scale;
    guint cstride;
    /* Area for PSF width measurement. */
    gint row;
    gint col;
    gint width;
    gint height;
} PSFSigmaOptData;

typedef struct {
    GwyDataField *matrix;
    GwyDataField *rhs;
    GwyDataField *tf;
    GwyDataField *buf;
    GwyDataField *ffield;
    GwyDataField *vfield;
    GwyDataField *wfield;
    fftw_plan fplan;
    fftw_plan bplan;
    fftw_complex *fmat;
    fftw_complex *fvec;
    gdouble matrix00;
    gdouble mnorm;
    gint want_txres;
    gint want_tyres;
    gint txres;
    gint tyres;
    gint xsize;
    gint ysize;
    gint border;
    gint cstride;
} PSFConjGradData;

typedef gdouble (*DoubleArrayFunc)(const gdouble *results);

static void
gwy_data_field_area_convolve_3x3(GwyDataField *data_field,
                                 const gdouble *kernel,
                                 gint col, gint row,
                                 gint width, gint height)
{
    gdouble *rm, *rc, *rp;
    gdouble t, v;
    gint xres, i, j;

    xres = data_field->xres;
    rp = data_field->data + row*xres + col;

    /* Special-case width == 1 to avoid complications below.  It's silly but the API guarantees it. */
    if (width == 1) {
        t = rp[0];
        for (i = 0; i < height; i++) {
            rc = rp = data_field->data + (row + i)*xres + col;
            if (i < height-1)
                rp += xres;

            v = (kernel[0] + kernel[1] + kernel[2])*t
                + (kernel[3] + kernel[4] + kernel[5])*rc[0]
                + (kernel[6] + kernel[7] + kernel[8])*rp[0];
            t = rc[0];
            rc[0] = v;
        }
        gwy_data_field_invalidate(data_field);

        return;
    }

    /* NB: This is not trivial to parallelise because of the one-pass method where we just keep a single-row buffer.
     * With multiple threads we have to ensure they do not write what we still want to read. */
    rm = g_new(gdouble, width);
    gwy_assign(rm, rp, width);

    for (i = 0; i < height; i++) {
        rc = rp;
        if (i < height-1)
            rp += xres;
        v = (kernel[0] + kernel[1])*rm[0] + kernel[2]*rm[1] + (kernel[3] + kernel[4])*rc[0] + kernel[5]*rc[1]
            + (kernel[6] + kernel[7])*rp[0] + kernel[8]*rp[1];
        t = rc[0];
        rc[0] = v;
        if (i < height-1) {
            for (j = 1; j < width-1; j++) {
                v = kernel[0]*rm[j-1] + kernel[1]*rm[j] + kernel[2]*rm[j+1] + kernel[3]*t + kernel[4]*rc[j]
                    + kernel[5]*rc[j+1] + kernel[6]*rp[j-1] + kernel[7]*rp[j] + kernel[8]*rp[j+1];
                rm[j-1] = t;
                t = rc[j];
                rc[j] = v;
            }
            v = kernel[0]*rm[j-1] + (kernel[1] + kernel[2])*rm[j] + kernel[3]*t + (kernel[4] + kernel[5])*rc[j]
                + kernel[6]*rp[j-1] + (kernel[7] + kernel[8])*rp[j];
        }
        else {
            for (j = 1; j < width-1; j++) {
                v = kernel[0]*rm[j-1] + kernel[1]*rm[j] + kernel[2]*rm[j+1] + kernel[3]*t + kernel[4]*rc[j]
                    + kernel[5]*rc[j+1] + kernel[6]*t + kernel[7]*rc[j] + kernel[8]*rc[j+1];
                rm[j-1] = t;
                t = rc[j];
                rc[j] = v;
            }
            v = kernel[0]*rm[j-1] + (kernel[1] + kernel[2])*rm[j] + kernel[3]*t + (kernel[4] + kernel[5])*rc[j]
                + kernel[6]*t + (kernel[7] + kernel[8])*rc[j];
        }
        rm[j-1] = t;
        rm[j] = rc[j];
        rc[j] = v;
    }

    g_free(rm);
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_area_convolve:
 * @data_field: A data field to convolve.  It must be at least as large as 1/3 of @kernel_field in each dimension.
 * @kernel_field: Kenrel field to convolve @data_field with.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Convolves a rectangular part of a data field with given kernel.
 *
 * Note that the convolution is done by summation and can be slow for large kernels.
 **/
void
gwy_data_field_area_convolve(GwyDataField *data_field,
                             GwyDataField *kernel_field,
                             gint col, gint row,
                             gint width, gint height)
{
    gint xres, yres, kxres, kyres, i, j;
    GwyDataField *hlp_df;
    gdouble *d, *h, *k;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(kernel_field));

    xres = data_field->xres;
    yres = data_field->yres;
    kxres = kernel_field->xres;
    kyres = kernel_field->yres;

    if (kxres == 3 && kyres == 3) {
        gwy_data_field_area_convolve_3x3(data_field, kernel_field->data, col, row, width, height);
        return;
    }

    hlp_df = gwy_data_field_new(width, height, 1.0, 1.0, FALSE);
    d = data_field->data;
    k = kernel_field->data;
    h = hlp_df->data;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(d,k,h,xres,yres,kxres,kyres,col,row,width,height)
#endif
    for (i = row; i < row + height; i++) {
        for (j = col; j < col + width; j++) {
            gdouble v = 0.0;
            gint m, n, ii, jj;

            for (m = -kyres/2; m < kyres - kyres/2; m++) {
                ii = i + m;
                if (G_UNLIKELY(ii < 0))
                    ii = -ii-1;
                else if (G_UNLIKELY(ii >= yres))
                    ii = 2*yres-1 - ii;

                for (n = -kxres/2; n < kxres - kxres/2; n++) {
                    jj = j + n;
                    if (G_UNLIKELY(jj < 0))
                        jj = -jj-1;
                    else if (G_UNLIKELY(jj >= xres))
                        jj = 2*xres-1 - jj;

                    v += d[ii*xres + jj] * k[kxres*(m + kyres/2) + n + kxres/2];
                }
            }
            h[(i - row)*width + (j - col)] = v;
        }
    }
    gwy_data_field_area_copy(hlp_df, data_field, 0, 0, width, height, col, row);
    g_object_unref(hlp_df);

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_convolve:
 * @data_field: A data field to convolve.  It must be at least as large as 1/3 of @kernel_field in each dimension.
 * @kernel_field: Kenrel field to convolve @data_field with.
 *
 * Convolves a data field with given kernel.
 *
 * Note that the convolution is done by summation and can be slow for large kernels.
 **/
void
gwy_data_field_convolve(GwyDataField *data_field,
                        GwyDataField *kernel_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_convolve(data_field, kernel_field, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_fft_convolve:
 * @data_field: A data field to convolve.
 * @kernel_field: Kenrel field to convolve @data_field with.  It must have the same size as @data_field.
 *
 * Convolves a data field with given kernel of the same size using FFT.
 *
 * This is a simple FFT-based convolution done by multiplication in the frequency domain.
 *
 * This is a somewhat low-level function.  There is no padding or boundary treatment; images are considered periodic.
 * The result is normalised as if the convolution was done by summation and the physical units of @data_field are
 * unchanged.
 *
 * Also note that in order to obtain unshifted result, the kernel needs to be centered around the top left corner.
 * You can use gwy_data_field_2dfft_dehumanize() to transform a centered kernel.
 *
 * Since: 2.54
 **/
void
gwy_data_field_fft_convolve(GwyDataField *data_field,
                            GwyDataField *kernel_field)
{
    fftw_plan fplan, bplan;
    fftw_complex *fieldc, *kernelc;
    gdouble *rbuf;
    guint i, xres, yres, n, cstride;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(kernel_field));
    xres = data_field->xres;
    yres = data_field->yres;
    g_return_if_fail(kernel_field->xres == xres);
    g_return_if_fail(kernel_field->yres == yres);

    cstride = xres/2 + 1;
    n = xres*yres;
    rbuf = gwy_fftw_new_real(n);
    fieldc = gwy_fftw_new_complex(cstride*yres);
    kernelc = gwy_fftw_new_complex(cstride*yres);

    fplan = gwy_fftw_plan_dft_r2c_2d(yres, xres, rbuf, fieldc, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    bplan = gwy_fftw_plan_dft_c2r_2d(yres, xres, fieldc, rbuf, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    gwy_assign(rbuf, kernel_field->data, xres*yres);
    gwy_fftw_execute(fplan);
    gwy_assign(kernelc, fieldc, cstride*yres);

    gwy_assign(rbuf, data_field->data, xres*yres);
    gwy_fftw_execute(fplan);

    for (i = 0; i < cstride*yres; i++) {
        gdouble re = (gwycreal(fieldc[i])*gwycreal(kernelc[i]) - gwycimag(fieldc[i])*gwycimag(kernelc[i]));
        gdouble im = (gwycreal(fieldc[i])*gwycimag(kernelc[i]) + gwycimag(fieldc[i])*gwycreal(kernelc[i]));
        gwycreal(fieldc[i]) = re;
        gwycimag(fieldc[i]) = im;
    }
    gwy_fftw_execute(bplan);

    gwy_assign(data_field->data, rbuf, n);
    gwy_data_field_invalidate(data_field);
    gwy_data_field_multiply(data_field, 1.0/n);

    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);
    fftw_free(fieldc);
    fftw_free(kernelc);
    fftw_free(rbuf);
}

static void
gwy_data_field_area_hconvolve(GwyDataField *data_field,
                              GwyDataLine *kernel_line,
                              gint col, gint row,
                              gint width, gint height)
{
    gint xres, kres, mres, k0;
    const gdouble *kernel;
    gdouble *data;

    xres = data_field->xres;
    kres = kernel_line->res;
    data = data_field->data;
    kernel = kernel_line->data;
    mres = 2*width;
    k0 = (kres/2 + 1)*mres;

#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            shared(data,kernel,xres,kres,mres,k0,col,row,width,height)
#endif
    {
        gdouble *drow, *buf = g_new(gdouble, kres);
        gint ifrom = gwy_omp_chunk_start(height);
        gint ito = gwy_omp_chunk_end(height);
        gint i, j, k, pos;
        gdouble d;

        for (i = ifrom; i < ito; i++) {
            drow = data + (row + i)*xres + col;
            /* Initialize with a triangluar sums, mirror-extend */
            gwy_clear(buf, kres);
            for (j = 0; j < kres; j++) {
                k = (j - kres/2 + k0) % mres;
                d = drow[k < width ? k : mres-1 - k];
                for (k = 0; k <= j; k++)
                    buf[k] += kernel[j - k]*d;
            }
            pos = 0;
            /* Middle part and tail with mirror extension again, we do some O(1/2*k^2) of useless work here by not
             * separating the tail */
            for (j = 0; j < width; j++) {
                drow[j] = buf[pos];
                buf[pos] = 0.0;
                pos = (pos + 1) % kres;
                k = (j + kres - kres/2 + k0) % mres;
                d = drow[G_LIKELY(k < width) ? k : mres-1 - k];
                for (k = pos; k < kres; k++)
                    buf[k] += kernel[kres-1 - (k - pos)]*d;
                for (k = 0; k < pos; k++)
                    buf[k] += kernel[pos-1 - k]*d;
            }
        }

        g_free(buf);
    }
}

static void
gwy_data_field_area_vconvolve(GwyDataField *data_field,
                              GwyDataLine *kernel_line,
                              gint col, gint row,
                              gint width, gint height)
{
    gint kres, xres, mres, k0;
    const gdouble *kernel;
    gdouble *data;

    xres = data_field->xres;
    kres = kernel_line->res;
    data = data_field->data;
    kernel = kernel_line->data;
    mres = 2*height;
    k0 = (kres/2 + 1)*mres;

    /* This looks like a bad memory access pattern.  And for small kernels it indeed is (we should iterate row-wise
     * and directly calculate the sums). For large kernels this is mitigated by the maximum possible amount of work
     * done per a data field access. */
#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            shared(data,kernel,xres,kres,mres,k0,col,row,width,height)
#endif
    {
        gdouble *dcol, *buf = g_new(gdouble, kres);
        gint jfrom = gwy_omp_chunk_start(width);
        gint jto = gwy_omp_chunk_end(width);
        gint i, j, k, pos;
        gdouble d;

        for (j = jfrom; j < jto; j++) {
            dcol = data + row*xres + (col + j);
            /* Initialize with a triangluar sums, mirror-extend */
            gwy_clear(buf, kres);
            for (i = 0; i < kres; i++) {
                k = (i - kres/2 + k0) % mres;
                d = dcol[k < height ? k*xres : (mres-1 - k)*xres];
                for (k = 0; k <= i; k++)
                    buf[k] += kernel[i - k]*d;
            }
            pos = 0;
            /* Middle part and tail with mirror extension again, we do some O(1/2*k^2) of useless work here by not
             * separating the tail */
            for (i = 0; i < height; i++) {
                dcol[i*xres] = buf[pos];
                buf[pos] = 0.0;
                pos = (pos + 1) % kres;
                k = (i + kres - kres/2 + k0) % mres;
                d = dcol[G_LIKELY(k < height) ? k*xres : (mres-1 - k)*xres];
                for (k = pos; k < kres; k++)
                    buf[k] += kernel[kres-1 - (k - pos)]*d;
                for (k = 0; k < pos; k++)
                    buf[k] += kernel[pos-1 - k]*d;
            }
        }

        g_free(buf);
    }
}

/**
 * gwy_data_field_area_convolve_1d:
 * @data_field: A data field to convolve.  It must be at least as large as 1/3 of @kernel_field in the corresponding
 *              dimension.
 * @kernel_line: Kernel line to convolve @data_field with.
 * @orientation: Filter orientation (%GWY_ORIENTATION_HORIZONTAL for row-wise convolution, %GWY_ORIENTATION_VERTICAL
 *               for column-wise convolution).
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Convolves a rectangular part of a data field with given linear kernel.
 *
 * For large separable kernels it can be more efficient to use a sequence of horizontal and vertical convolutions
 * instead one 2D convolution.
 *
 * Since: 2.4
 **/
void
gwy_data_field_area_convolve_1d(GwyDataField *data_field,
                                GwyDataLine *kernel_line,
                                GwyOrientation orientation,
                                gint col, gint row,
                                gint width, gint height)
{
    gint kres;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_LINE(kernel_line));

    kres = kernel_line->res;
    if (kres == 1) {
        gwy_data_field_area_multiply(data_field, col, row, width, height, kernel_line->data[0]);
        return;
    }

    if (orientation == GWY_ORIENTATION_HORIZONTAL)
        gwy_data_field_area_hconvolve(data_field, kernel_line, col, row, width, height);
    else
        gwy_data_field_area_vconvolve(data_field, kernel_line, col, row, width, height);

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_convolve_1d:
 * @data_field: A data field to convolve.  It must be at least as large as 1/3 of @kernel_field in the corresponding
 *              dimension.
 * @kernel_line: Kenrel line to convolve @data_field with.
 * @orientation: Filter orientation (see gwy_data_field_area_convolve_1d()).
 *
 * Convolves a data field with given linear kernel.
 *
 * Since: 2.4
 **/
void
gwy_data_field_convolve_1d(GwyDataField *data_field,
                           GwyDataLine *kernel_line,
                           GwyOrientation orientation)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_convolve_1d(data_field, kernel_line, orientation, 0, 0, data_field->xres, data_field->yres);
}

static void
ensure_defined_exterior(GwyExteriorType *exterior, gdouble *fill_value)
{
    if (*exterior == GWY_EXTERIOR_UNDEFINED) {
        g_warning("Do not use GWY_EXTERIOR_UNDEFINED for convolutions and correlations.  "
                  "Fixing to zero-filled exterior.");
        *exterior = GWY_EXTERIOR_FIXED_VALUE;
        *fill_value = 0.0;
    }
}

static void
row_convolve_direct(const GwyDataField *field,
                    guint col, guint row,
                    guint width, guint height,
                    GwyDataField *target,
                    guint targetcol, guint targetrow,
                    const GwyDataLine *kernel,
                    RowExtendFunc extend_row,
                    gdouble fill_value)
{
    guint xres = field->xres;
    guint kres = kernel->res;
    const gdouble *kdata = kernel->data;
    guint size = width + kres - 1;
    guint extend_left, extend_right, i, j, k;
    gdouble *extdata = g_new(gdouble, size);

    make_symmetrical_extension(width, size, &extend_left, &extend_right);

    /* The direct method is used only if kres ≪ res.  Don't bother optimising the boundaries, just make the inner loop
     * tight. */
    for (i = 0; i < height; i++) {
        gdouble *trow = target->data + (targetrow + i)*target->xres + targetcol;
        extend_row(field->data + (row + i)*xres, extdata, col, width, xres, extend_left, extend_right, fill_value);
        for (j = 0; j < width; j++) {
            const gdouble *d = extdata + extend_left + kres/2 + j;
            gdouble v = 0.0;
            for (k = 0; k < kres; k++, d--)
                v += kdata[k] * *d;
            trow[j] = v;
        }
    }

    g_free(extdata);
}

static inline void
complex_multiply_with(fftw_complex *a, const fftw_complex *b)
{
    gdouble re = (*a)[0]*(*b)[0] - (*a)[1]*(*b)[1];

    (*a)[1] = (*a)[1]*(*b)[0] + (*a)[0]*(*b)[1];
    (*a)[0] = re;
}

static void
row_convolve_fft(GwyDataField *field,
                 guint col, guint row,
                 guint width, guint height,
                 GwyDataField *target,
                 guint targetcol, guint targetrow,
                 GwyDataLine *kernel,
                 RowExtendFunc extend_row,
                 gdouble fill_value)
{
    guint xres = field->xres, kres = kernel->res;
    guint size = gwy_fft_find_nice_size(width + kres - 1);
    // The innermost (contiguous) dimension of R2C the complex output is slightly larger than the real input.  Note
    // @cstride is measured in fftw_complex, multiply it by 2 for doubles.
    guint cstride = size/2 + 1;
    guint extend_left, extend_right, i, j;
    gdouble *extdata = gwy_fftw_new_real(size);
    fftw_complex *datac = gwy_fftw_new_complex(cstride);
    fftw_complex *kernelc = gwy_fftw_new_complex(cstride);
    fftw_plan dplan, cplan;
    gdouble q;

    /* The R2C plan for transforming the extended data row (or kernel). */
    dplan = gwy_fftw_plan_dft_r2c_1d(size, extdata, datac, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    /* The C2R plan the backward transform of the convolution of each row. */
    cplan = gwy_fftw_plan_dft_c2r_1d(size, datac, extdata, FFTW_ESTIMATE);

    // Transform the kernel.
    extend_kernel_row(kernel->data, kres, extdata, size);
    gwy_fftw_execute(dplan);
    gwy_assign(kernelc, datac, cstride);

    // Convolve rows
    make_symmetrical_extension(width, size, &extend_left, &extend_right);
    q = 1.0/size;
    for (i = 0; i < height; i++) {
        extend_row(field->data + (row + i)*xres, extdata, col, width, xres, extend_left, extend_right, fill_value);
        gwy_fftw_execute(dplan);
        for (j = 0; j < cstride; j++) {
            complex_multiply_with(datac + j, kernelc + j);
            datac[j][0] *= q;
            datac[j][1] *= q;
        }
        gwy_fftw_execute(cplan);
        gwy_assign(target->data + (targetrow + i)*target->xres + targetcol, extdata + extend_left, width);
    }

    fftw_destroy_plan(cplan);
    fftw_destroy_plan(dplan);
    fftw_free(extdata);
    fftw_free(datac);
    fftw_free(extdata);
}

/**
 * gwy_data_field_area_ext_row_convolve:
 * @field: A two-dimensional data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @target: A two-dimensional data field where the result will be placed. It may be @field for an in-place
 *          modification.
 * @kernel: Kernel to convolve @field with.
 * @exterior: Exterior pixels handling.
 * @fill_value: The value to use with %GWY_EXTERIOR_FIXED_VALUE exterior.
 * @as_integral: %TRUE for normalisation and units as a convolution integral, %FALSE as a sum.
 *
 * Convolve a field row-wise with a one-dimensional kernel.
 *
 * Pixel dimensions of @target may match either @field or just the rectangular area.  In the former case the result is
 * written in the same rectangular area; in the latter case the result fills the entire @target.
 *
 * The convolution is performed with the kernel centred on the respective field pixels.  For an odd-sized kernel this
 * holds precisely.  For an even-sized kernel this means the kernel centre is placed 0.5 pixel to the left (towards
 * lower column indices) from the respective field pixel.
 *
 * See gwy_data_field_extend() for what constitutes the exterior and how it is handled.
 *
 * If @as_integral is %FALSE the function performs a simple discrete convolution sum and the value units of @target
 * are set to product of @field and @kernel units.
 *
 * If @as_integral is %TRUE the function approximates a convolution integral. In this case @kernel should be a sampled
 * continuous transfer function. The units of value @target are set to product of @field and @kernel value units and
 * @field lateral units.  Furthermore, the discrete sum is multiplied by the pixel size (i.e. d@x in the integral).
 *
 * In either case, the lateral units and pixel size of @kernel are assumed to be the same as for a @field's row
 * (albeit not checked), because the convolution does not make sense otherwise.
 *
 * Since: 2.49
 **/
void
gwy_data_field_area_ext_row_convolve(GwyDataField *field,
                                     guint col, guint row,
                                     guint width, guint height,
                                     GwyDataField *target,
                                     GwyDataLine *kernel,
                                     GwyExteriorType exterior,
                                     gdouble fill_value,
                                     gboolean as_integral)
{
    guint xres, yres;
    guint targetcol, targetrow;
    GwySIUnit *funit, *kunit, *tunit;
    RowExtendFunc extend_row;
    gdouble dx, dy;

    if (!_gwy_data_field_check_area(field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(target));
    g_return_if_fail(GWY_IS_DATA_LINE(kernel));
    xres = field->xres;
    yres = field->yres;
    g_return_if_fail((target->xres == xres && target->yres == yres)
                     || (target->xres == width && target->yres == height));
    targetcol = (target->xres == xres) ? col : 0;
    targetrow = (target->yres == yres) ? row : 0;

    ensure_defined_exterior(&exterior, &fill_value);
    if (!(extend_row = _gwy_get_row_extend_func(exterior)))
        return;

    if (width <= 12 || kernel->res <= 3.0*(log(width) - 1.0)) {
        row_convolve_direct(field, col, row, width, height,
                            target, targetcol, targetrow,
                            kernel, extend_row, fill_value);
    }
    else {
        row_convolve_fft(field, col, row, width, height,
                         target, targetcol, targetrow,
                         kernel, extend_row, fill_value);
    }

    dx = field->xreal/field->xres;
    dy = field->yreal/field->yres;
    if (target != field) {
        _gwy_copy_si_unit(field->si_unit_xy, &target->si_unit_xy);
        target->xreal = dx*target->xres;
        target->yreal = dy*target->yres;
    }

    funit = gwy_data_field_get_si_unit_z(field);
    kunit = gwy_data_line_get_si_unit_y(kernel);
    tunit = gwy_data_field_get_si_unit_z(target);
    gwy_si_unit_multiply(funit, kunit, tunit);
    if (as_integral) {
        funit = gwy_data_field_get_si_unit_xy(field);
        gwy_si_unit_multiply(tunit, funit, tunit);
        gwy_data_field_multiply(target, dx);
    }

    gwy_data_field_invalidate(target);
}

/**
 * multiconvolve_direct:
 * @field: A two-dimensional data field.
 * @col: First ROI column.
 * @row: First RIO row.
 * @width: ROI width.
 * @height: RIO height.
 * @target: A two-dimensional data field where the result will be placed. It may be @field itself.
 * @targetcol: Column to place the result into @target.
 * @targetrow: Target to place the result into @target.
 * @kernel: Array of @nkernel equally-sized kernel.
 * @nkernel: Number of items in @kernel.
 * @combine_results: Function to combine results of individual convolutions to the final result put to @target.  May
 *                   be %NULL if @nkernel is 1.
 * @extend_rect: Rectangle extending method.
 * @fill_value: The value to use with fixed-value exterior.
 *
 * Performs convolution of a field with a number of equally-sized kenrels, combining the results of individual
 * convolutions into a single value.
 */
static void
multiconvolve_direct(GwyDataField *field,
                     guint col, guint row,
                     guint width, guint height,
                     GwyDataField *target,
                     guint targetcol, guint targetrow,
                     GwyDataField **kernel,
                     guint nkernel,
                     DoubleArrayFunc combine_results,
                     RectExtendFunc extend_rect,
                     gdouble fill_value)
{
    guint xres, yres, kxres, kyres, xsize, ysize;
    guint extend_left, extend_right, extend_up, extend_down;
    gdouble *extdata;
    guint kno, i, j, ik, jk;

    g_return_if_fail(nkernel);
    g_return_if_fail(kernel);
    g_return_if_fail(nkernel == 1 || combine_results);

    xres = field->xres;
    yres = field->yres;
    kxres = kernel[0]->xres;
    kyres = kernel[0]->yres;
    for (kno = 1; kno < nkernel; kno++) {
        g_return_if_fail(kernel[kno]->xres == kxres && kernel[kno]->yres == kyres);
    }

    xsize = width + kxres - 1;
    ysize = height + kyres - 1;
    extdata = g_new(gdouble, xsize*ysize);
    make_symmetrical_extension(width, xsize, &extend_left, &extend_right);
    make_symmetrical_extension(height, ysize, &extend_up, &extend_down);

    extend_rect(field->data, xres, extdata, xsize,
                col, row, width, height, xres, yres,
                extend_left, extend_right, extend_up, extend_down, fill_value);

    /* The direct method is used only if kres ≪ res.  Don't bother optimising
     * the boundaries, just make the inner loop tight. */
    if (nkernel == 1) {
        const gdouble *kdata = kernel[0]->data;
        for (i = 0; i < height; i++) {
            gdouble *trow = target->data + ((targetrow + i)*target->xres + targetcol);
            for (j = 0; j < width; j++) {
                const gdouble *id = extdata + (extend_up + kyres/2 + i)*xsize;
                gdouble v = 0.0;
                for (ik = 0; ik < kyres; ik++, id -= xsize) {
                    const gdouble *jd = id + extend_left + kxres/2 + j;
                    const gdouble *krow = kdata + ik*kxres;
                    for (jk = 0; jk < kxres; jk++, jd--)
                        v += krow[jk] * *jd;
                }
                trow[j] = v;
            }
        }
    }
    else {
        for (i = 0; i < height; i++) {
            gdouble *trow = target->data + ((targetrow + i)*target->xres + targetcol);
            for (j = 0; j < width; j++) {
                gdouble results[nkernel];
                for (kno = 0; kno < nkernel; kno++) {
                    const gdouble *id = extdata + (extend_up + kyres/2 + i)*xsize;
                    const gdouble *kdata = kernel[kno]->data;
                    gdouble v = 0.0;
                    for (ik = 0; ik < kyres; ik++, id -= xsize) {
                        const gdouble *jd = id + extend_left + kxres/2 + j;
                        const gdouble *krow = kdata + ik*kxres;
                        for (jk = 0; jk < kxres; jk++, jd--)
                            v += krow[jk] * *jd;
                    }
                    results[kno] = v;
                }
                trow[j] = combine_results(results);
            }
        }
    }

    g_free(extdata);
}

static void
convolve_fft(GwyDataField *field,
             guint col, guint row,
             guint width, guint height,
             GwyDataField *target,
             guint targetcol, guint targetrow,
             GwyDataField *kernel,
             RectExtendFunc extend_rect,
             gdouble fill_value)
{
    guint xres = field->xres, yres = field->yres,
          kxres = kernel->xres, kyres = kernel->yres;
    guint xsize = gwy_fft_find_nice_size(width + kxres - 1);
    guint ysize = gwy_fft_find_nice_size(height + kyres - 1);
    /* The innermost (contiguous) dimension of R2C the complex output is slightly larger than the real input.  If the
     * transform is in-place the input array needs to be padded.  Note @cstride is measured in fftw_complex, multiply
     * it by 2 for doubles. */
    guint cstride = xsize/2 + 1;
    /* Use in-place transforms.  Let FFTW figure out whether allocating temporary buffers is worth it or not. */
    fftw_complex *datac = gwy_fftw_new_complex(cstride*ysize);
    gdouble *extdata = (gdouble*)datac;
    fftw_complex *kernelc = gwy_fftw_new_complex(cstride*ysize);
    guint extend_left, extend_right, extend_up, extend_down;
    guint i, k;
    fftw_plan kplan, dplan, cplan;
    gdouble q;

    make_symmetrical_extension(width, xsize, &extend_left, &extend_right);
    make_symmetrical_extension(height, ysize, &extend_up, &extend_down);

    /* The R2C plan for transforming the extended kernel.  The input is in extdata to make it an out-of-place
     * transform (this also means the input data row stride is just xsize, not 2*cstride). */
    kplan = gwy_fftw_plan_dft_r2c_2d(ysize, xsize, extdata, kernelc, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    // The R2C plan for transforming the extended data.  This one is in-place.
    dplan = gwy_fftw_plan_dft_r2c_2d(ysize, xsize, extdata, datac, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    // The C2R plan the backward transform of the convolution.  The input is in fact in kernelc to make it an
    // out-of-place transform.  So, again, the output has cstride of only xsize.
    cplan = gwy_fftw_plan_dft_c2r_2d(ysize, xsize, kernelc, extdata, FFTW_ESTIMATE);

    // Transform the kernel.
    extend_kernel_rect(kernel->data, kxres, kyres, extdata, xsize, ysize, xsize);
    gwy_fftw_execute(kplan);

    // Convolve
    extend_rect(field->data, xres, extdata, 2*cstride,
                col, row, width, height, xres, yres,
                extend_left, extend_right, extend_up, extend_down, fill_value);
    gwy_fftw_execute(dplan);

    q = 1.0/(xsize*ysize);
    for (k = 0; k < cstride*ysize; k++) {
        complex_multiply_with(kernelc + k, datac + k);
        kernelc[k][0] *= q;
        kernelc[k][1] *= q;
    }
    gwy_fftw_execute(cplan);

    for (i = 0; i < height; i++) {
        gwy_assign(target->data + (targetrow + i)*target->xres + targetcol,
                   extdata + (extend_up + i)*xsize + extend_left,
                   width);
    }

    fftw_destroy_plan(kplan);
    fftw_destroy_plan(cplan);
    fftw_destroy_plan(dplan);
    fftw_free(kernelc);
    fftw_free(datac);
}

/**
 * gwy_data_field_area_ext_convolve:
 * @field: A two-dimensional data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @target: A two-dimensional data field where the result will be placed. It may be @field for an in-place
 *          modification.
 * @kernel: Kernel to convolve @field with.
 * @exterior: Exterior pixels handling.
 * @fill_value: The value to use with %GWY_EXTERIOR_FIXED_VALUE exterior.
 * @as_integral: %TRUE for normalisation and units as a convolution integral, %FALSE as a sum.
 *
 * Convolve a field with a two-dimensional kernel.
 *
 * Pixel dimensions of @target may match either @field or just the rectangular area.  In the former case the result is
 * written in the same rectangular area; in the latter case the result fills the entire @target.
 *
 * The convolution is performed with the kernel centred on the respective field pixels.  For directions in which the
 * kernel has an odd size this holds precisely.  For an even-sized kernel this means the kernel centre is placed 0.5
 * pixel left or up (towards lower indices) from the respective field pixel.
 *
 * See gwy_data_field_extend() for what constitutes the exterior and how it is handled.
 *
 * If @as_integral is %FALSE the function performs a simple discrete convolution sum and the value units of @target
 * are set to product of @field and @kernel units.
 *
 * If @as_integral is %TRUE the function approximates a convolution integral. In this case @kernel should be a sampled
 * continuous transfer function. The units of value @target are set to product of @field and @kernel value units and
 * @field lateral units squared.  Furthermore, the discrete sum is multiplied by the pixel size (i.e. d@x d@y in the
 * integral).
 *
 * In either case, the lateral units and pixel size of @kernel are assumed to be the same as for @field (albeit not
 * checked), because the convolution does not make sense otherwise.
 *
 * Since: 2.49
 **/
void
gwy_data_field_area_ext_convolve(GwyDataField *field,
                                 guint col, guint row,
                                 guint width, guint height,
                                 GwyDataField *target,
                                 GwyDataField *kernel,
                                 GwyExteriorType exterior,
                                 gdouble fill_value,
                                 gboolean as_integral)
{
    guint xres, yres, size;
    guint targetcol, targetrow;
    GwySIUnit *funit, *kunit, *tunit;
    RectExtendFunc extend_rect;
    gdouble dx, dy;

    if (!_gwy_data_field_check_area(field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(target));
    g_return_if_fail(GWY_IS_DATA_FIELD(kernel));
    xres = field->xres;
    yres = field->yres;
    g_return_if_fail((target->xres == xres && target->yres == yres)
                     || (target->xres == width && target->yres == height));
    targetcol = (target->xres == xres) ? col : 0;
    targetrow = (target->yres == yres) ? row : 0;

    ensure_defined_exterior(&exterior, &fill_value);
    if (!(extend_rect =_gwy_get_rect_extend_func(exterior)))
        return;

    size = height*width;
    if (size <= 25) {
        multiconvolve_direct(field, col, row, width, height,
                             target, targetcol, targetrow,
                             &kernel, 1, NULL, extend_rect, fill_value);
    }
    else
        convolve_fft(field, col, row, width, height, target, targetcol, targetrow, kernel, extend_rect, fill_value);

    dx = field->xreal/field->xres;
    dy = field->yreal/field->yres;
    if (target != field) {
        funit = gwy_data_field_get_si_unit_xy(field);
        tunit = gwy_data_field_get_si_unit_xy(target);
        gwy_si_unit_assign(tunit, funit);
        target->xreal = dx*target->xres;
        target->yreal = dy*target->yres;
    }

    funit = gwy_data_field_get_si_unit_z(field);
    kunit = gwy_data_field_get_si_unit_z(kernel);
    tunit = gwy_data_field_get_si_unit_z(target);
    gwy_si_unit_multiply(funit, kunit, tunit);
    if (as_integral) {
        funit = gwy_data_field_get_si_unit_xy(field);
        gwy_si_unit_power_multiply(tunit, 1, funit, 2, tunit);
        gwy_data_field_multiply(target, dx*dy);
    }

    gwy_data_field_invalidate(target);
}

/**
 * gwy_data_field_area_filter_laplacian:
 * @data_field: A data field to apply the filter to.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with Laplacian filter.
 **/
void
gwy_data_field_area_filter_laplacian(GwyDataField *data_field,
                                     gint col, gint row,
                                     gint width, gint height)
{
    const gdouble laplace[] = {
        0,  1, 0,
        1, -4, 1,
        0,  1, 0,
    };

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_convolve_3x3(data_field, laplace, col, row, width, height);
}

/**
 * gwy_data_field_filter_laplacian:
 * @data_field: A data field to apply the filter to.
 *
 * Filters a data field with Laplacian filter.
 **/
void
gwy_data_field_filter_laplacian(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_laplacian(data_field, 0, 0, data_field->xres, data_field->yres);
}

 /**
 * gwy_data_field_area_filter_laplacian_of_gaussians:
 * @data_field: A data field to apply the filter to.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with Laplacian of Gaussians filter.
 *
 * Since: 2.23
 **/
void
gwy_data_field_area_filter_laplacian_of_gaussians(GwyDataField *data_field,
                                                  gint col, gint row,
                                                  gint width,
                                                  gint height)
{
    /* optimized mexican hat from Scharr's works */
    enum { size = 5 };
    const gdouble laplacian_of_gaussians_data[size*size] = {
          1, -12,    3, -12,   1,
        -12,  78,  167,  78, -12,
          3, 167, -902, 167,   3,
        -12,  78,  167,  78, -12,
          1, -12,    3, -12,   1,
    };
    GwyDataField *laplacian_of_gaussians;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    laplacian_of_gaussians = gwy_data_field_new(size, size, 1.0, 1.0, TRUE);
    gwy_assign(laplacian_of_gaussians->data, laplacian_of_gaussians_data, size*size);
    gwy_data_field_area_convolve(data_field, laplacian_of_gaussians, col, row, width, height);
    g_object_unref(laplacian_of_gaussians);
}

/**
 * gwy_data_field_filter_laplacian_of_gaussians:
 * @data_field: A data field to apply the filter to.
 *
 * Filters a data field with Laplacian of Gaussians filter.
 *
 * Since: 2.23
 **/
void
gwy_data_field_filter_laplacian_of_gaussians(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_laplacian_of_gaussians(data_field, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_area_filter_sobel:
 * @data_field: A data field to apply the filter to.
 * @orientation: Filter orientation.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with a directional Sobel filter.
 **/
void
gwy_data_field_area_filter_sobel(GwyDataField *data_field,
                                 GwyOrientation orientation,
                                 gint col, gint row,
                                 gint width, gint height)
{
    static const gdouble hsobel[] = {
        0.25, 0, -0.25,
        0.5,  0, -0.5,
        0.25, 0, -0.25,
    };
    static const gdouble vsobel[] = {
         0.25,  0.5,  0.25,
         0,     0,    0,
        -0.25, -0.5, -0.25,
    };

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    if (orientation == GWY_ORIENTATION_HORIZONTAL)
        gwy_data_field_area_convolve_3x3(data_field, hsobel, col, row, width, height);
    else
        gwy_data_field_area_convolve_3x3(data_field, vsobel, col, row, width, height);
}

/**
 * gwy_data_field_filter_sobel:
 * @data_field: A data field to apply the filter to.
 * @orientation: Filter orientation.
 *
 * Filters a data field with a directional Sobel filter.
 **/
void
gwy_data_field_filter_sobel(GwyDataField *data_field,
                            GwyOrientation orientation)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_sobel(data_field, orientation, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_filter_sobel_total:
 * @data_field: A data field to apply the filter to.
 *
 * Filters a data field with total Sobel filter.
 *
 * Since: 2.31
 **/
void
gwy_data_field_filter_sobel_total(GwyDataField *data_field)
{
    GwyDataField *workspace;
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    workspace = gwy_data_field_duplicate(data_field);
    gwy_data_field_area_filter_sobel(data_field, GWY_ORIENTATION_HORIZONTAL, 0, 0, data_field->xres, data_field->yres);
    gwy_data_field_area_filter_sobel(workspace, GWY_ORIENTATION_VERTICAL, 0, 0, data_field->xres, data_field->yres);
    gwy_data_field_hypot_of_fields(data_field, data_field, workspace);
    g_object_unref(workspace);
}

/**
 * gwy_data_field_area_filter_prewitt:
 * @data_field: A data field to apply the filter to.
 * @orientation: Filter orientation.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with a directional Prewitt filter.
 **/
void
gwy_data_field_area_filter_prewitt(GwyDataField *data_field,
                                   GwyOrientation orientation,
                                   gint col, gint row,
                                   gint width, gint height)
{
    static const gdouble hprewitt[] = {
        1.0/3.0, 0, -1.0/3.0,
        1.0/3.0, 0, -1.0/3.0,
        1.0/3.0, 0, -1.0/3.0,
    };
    static const gdouble vprewitt[] = {
         1.0/3.0,  1.0/3.0,  1.0/3.0,
         0,        0,        0,
        -1.0/3.0, -1.0/3.0, -1.0/3.0,
    };

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    if (orientation == GWY_ORIENTATION_HORIZONTAL)
        gwy_data_field_area_convolve_3x3(data_field, hprewitt, col, row, width, height);
    else
        gwy_data_field_area_convolve_3x3(data_field, vprewitt, col, row, width, height);
}

/**
 * gwy_data_field_filter_prewitt:
 * @data_field: A data field to apply the filter to.
 * @orientation: Filter orientation.
 *
 * Filters a data field with Prewitt filter.
 **/
void
gwy_data_field_filter_prewitt(GwyDataField *data_field,
                              GwyOrientation orientation)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_prewitt(data_field, orientation, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_filter_prewitt_total:
 * @data_field: A data field to apply the filter to.
 *
 * Filters a data field with total Prewitt filter.
 *
 * Since: 2.31
 **/
void
gwy_data_field_filter_prewitt_total(GwyDataField *data_field)
{
    GwyDataField *workspace;
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    workspace = gwy_data_field_duplicate(data_field);
    gwy_data_field_area_filter_prewitt(data_field, GWY_ORIENTATION_HORIZONTAL,
                                       0, 0, data_field->xres, data_field->yres);
    gwy_data_field_area_filter_prewitt(workspace, GWY_ORIENTATION_VERTICAL,
                                       0, 0, data_field->xres, data_field->yres);
    gwy_data_field_hypot_of_fields(data_field, data_field, workspace);
    g_object_unref(workspace);
}

/**
 * gwy_data_field_area_filter_dechecker:
 * @data_field: A data field to apply the filter to.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with 5x5 checker pattern removal
 * filter.
 *
 * Since: 2.1
 **/
void
gwy_data_field_area_filter_dechecker(GwyDataField *data_field,
                                     gint col, gint row,
                                     gint width, gint height)
{
    enum { size = 5 };
    static const gdouble dechecker[size*size] = {
         0.0,        1.0/144.0, -1.0/72.0,  1.0/144.0,  0.0,
         1.0/144.0, -1.0/18.0,   1.0/9.0,  -1.0/18.0,   1.0/144.0,
        -1.0/72.0,   1.0/9.0,    7.0/9.0,   1.0/9.0,   -1.0/72.0,
         1.0/144.0, -1.0/18.0,   1.0/9.0,  -1.0/18.0,   1.0/144.0,
         0.0,        1.0/144.0, -1.0/72.0,  1.0/144.0,  0.0,
    };
    GwyDataField *kernel;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    kernel = gwy_data_field_new(size, size, 1.0, 1.0, FALSE);
    gwy_assign(kernel->data, dechecker, size*size);
    gwy_data_field_area_convolve(data_field, kernel, col, row, width, height);
    g_object_unref(kernel);
}

/**
 * gwy_data_field_filter_dechecker:
 * @data_field: A data field to apply the filter to.
 *
 * Filters a data field with 5x5 checker pattern removal filter.
 *
 * Since: 2.1
 **/
void
gwy_data_field_filter_dechecker(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_dechecker(data_field, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_area_filter_gaussian:
 * @data_field: A data field to apply the filter to.
 * @sigma: The sigma parameter of the Gaussian.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with a Gaussian filter.
 *
 * The Gausian is normalized, i.e. it is sum-preserving.
 *
 * Since: 2.4
 **/
void
gwy_data_field_area_filter_gaussian(GwyDataField *data_field,
                                    gdouble sigma,
                                    gint col, gint row,
                                    gint width, gint height)
{
    GwyDataLine *kernel;
    gdouble x;
    gint res, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(sigma >= 0.0);
    if (sigma == 0.0)
        return;

    res = (gint)ceil(5.0*sigma);
    res = 2*res + 1;
    /* FIXME */
    i = 3*MIN(data_field->xres, data_field->yres);
    if (res > i) {
        res = i;
        if (res % 2 == 0)
            res--;
    }

    kernel = gwy_data_line_new(res, 1.0, FALSE);
    for (i = 0; i < res; i++) {
        x = i - (res - 1)/2.0;
        x /= sigma;
        kernel->data[i] = exp(-x*x/2.0);
    }
    gwy_data_line_multiply(kernel, 1.0/gwy_data_line_get_sum(kernel));
    gwy_data_field_area_convolve_1d(data_field, kernel, GWY_ORIENTATION_HORIZONTAL, col, row, width, height);
    gwy_data_field_area_convolve_1d(data_field, kernel, GWY_ORIENTATION_VERTICAL, col, row, width, height);
    g_object_unref(kernel);
}

/**
 * gwy_data_field_filter_gaussian:
 * @data_field: A data field to apply the filter to.
 * @sigma: The sigma parameter of the Gaussian.
 *
 * Filters a data field with a Gaussian filter.
 *
 * Since: 2.4
 **/
void
gwy_data_field_filter_gaussian(GwyDataField *data_field,
                               gdouble sigma)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_gaussian(data_field, sigma, 0, 0, data_field->xres, data_field->yres);
}

/**
 * gwy_data_field_row_gaussian:
 * @data_field: A data field to apply the filter to.
 * @sigma: The sigma parameter of the Gaussian.
 *
 * Filters a data field with a Gaussian filter in horizontal direction.
 *
 * The Gausian is normalized, i.e. it is sum-preserving.
 *
 * Since: 2.54
 **/
void
gwy_data_field_row_gaussian(GwyDataField *data_field,
                            gdouble sigma)
{
    GwyDataLine *kernel;
    gdouble x;
    gint res, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(sigma >= 0.0);
    if (sigma == 0.0)
        return;

    res = (gint)ceil(5.0*sigma);
    res = 2*res + 1;
    i = 3*data_field->xres;
    if (res > i) {
        res = i;
        if (res % 2 == 0)
            res--;
    }

    kernel = gwy_data_line_new(res, 1.0, FALSE);
    for (i = 0; i < res; i++) {
        x = i - (res - 1)/2.0;
        x /= sigma;
        kernel->data[i] = exp(-x*x/2.0);
    }
    gwy_data_line_multiply(kernel, 1.0/gwy_data_line_get_sum(kernel));

    /* FIXME: use row_convolve_fft instead */
    gwy_data_field_area_convolve_1d(data_field, kernel, GWY_ORIENTATION_HORIZONTAL,
                                    0, 0, data_field->xres, data_field->yres);
    g_object_unref(kernel);
}

/**
 * gwy_data_field_column_gaussian:
 * @data_field: A data field to apply the filter to.
 * @sigma: The sigma parameter of the Gaussian.
 *
 * Filters a data field with a Gaussian filter in vertical direction.
 *
 * The Gausian is normalized, i.e. it is sum-preserving.
 *
 * Since: 2.54
 **/
void
gwy_data_field_column_gaussian(GwyDataField *data_field, gdouble sigma)
{
    GwyDataField *flipped = gwy_data_field_new_alike(data_field, FALSE);

    gwy_data_field_flip_xy(data_field, flipped, FALSE);
    gwy_data_field_row_gaussian(flipped, sigma);
    gwy_data_field_flip_xy(flipped, data_field, FALSE);
    g_object_unref(flipped);
}

static void
set_transfer_function_units(GwyDataField *ideal, GwyDataField *measured,
                            GwyDataField *transferfunc)
{
    GwySIUnit *sunit, *iunit, *tunit, *xyunit;

    xyunit = gwy_data_field_get_si_unit_xy(measured);
    sunit = gwy_data_field_get_si_unit_z(ideal);
    iunit = gwy_data_field_get_si_unit_z(measured);
    tunit = gwy_data_field_get_si_unit_z(transferfunc);
    gwy_si_unit_divide(iunit, sunit, tunit);
    gwy_si_unit_power_multiply(tunit, 1, xyunit, -2, tunit);
}

static inline void
deconvolve_do(const fftw_complex *fproduct,
              const fftw_complex *foperand,
              fftw_complex *fresult,
              guint n,
              gdouble lambda)
{
    guint k;

    for (k = 0; k < n; k++) {
        gdouble pre = gwycreal(fproduct[k]), pim = gwycimag(fproduct[k]);
        gdouble ore = gwycreal(foperand[k]), oim = gwycimag(foperand[k]);
        gdouble opnorm = ore*ore + oim*oim;
        gdouble denom = opnorm + lambda;

        gwycreal(fresult[k]) = (pre*ore + pim*oim)/denom;
        gwycimag(fresult[k]) = (-pre*oim + pim*ore)/denom;
    }
}

/**
 * gwy_data_field_deconvolve_regularized:
 * @dfield: A data field.
 * @operand: One of the factors entering the convolution resulting in @dfield. It must have the same dimensions as
 *           @dfield and it is assumed it has also the same physical size.
 * @out: Data field where to put the result into.  It will be resized to match @dfield.  It can also be @dfield
 *       itself.
 * @sigma: Regularization parameter.
 *
 * Performs deconvolution of a data field using a simple regularization.
 *
 * The operation can be used to deblur an image or conversely recover the point spread function from ideal response
 * image.
 *
 * Convolving the result with the operand using gwy_data_field_area_ext_convolve() with @as_integral=%TRUE will
 * recover (approximately) the image.  This means the deconvolution assumes continous convolution, not discrete sums.
 * Note that for the latter case this means the point spread function will be centered in @out.
 *
 * For recovery of transfer function, @dfield and @operand should be windowed beforehand if they are not periodic.
 *
 * Since: 2.51
 **/
void
gwy_data_field_deconvolve_regularized(GwyDataField *dfield,
                                      GwyDataField *operand,
                                      GwyDataField *out,
                                      gdouble sigma)
{
    gint xres, yres, cstride;
    gdouble lambda, msq;
    fftw_complex *ffield, *foper;
    fftw_plan fplan, bplan;

    g_return_if_fail(GWY_IS_DATA_FIELD(dfield));
    g_return_if_fail(GWY_IS_DATA_FIELD(operand));
    g_return_if_fail(GWY_IS_DATA_FIELD(out));
    xres = dfield->xres;
    yres = dfield->yres;
    g_return_if_fail(operand->xres == xres);
    g_return_if_fail(operand->yres == yres);

    cstride = xres/2 + 1;
    gwy_data_field_resample(out, xres, yres, GWY_INTERPOLATION_NONE);

    msq = gwy_data_field_get_mean_square(operand);
    if (!msq) {
        g_warning("Deconvolution by zero.");
        gwy_data_field_clear(out);
        return;
    }

    ffield = gwy_fftw_new_complex(cstride*yres);
    foper = gwy_fftw_new_complex(cstride*yres);
    fplan = gwy_fftw_plan_dft_r2c_2d(yres, xres, out->data, ffield, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    bplan = gwy_fftw_plan_dft_c2r_2d(yres, xres, ffield, out->data, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    gwy_data_field_copy(operand, out, FALSE);
    gwy_fftw_execute(fplan);
    gwy_assign(foper, ffield, cstride*yres);

    gwy_data_field_copy(dfield, out, FALSE);
    gwy_fftw_execute(fplan);
    fftw_destroy_plan(fplan);

    /* This seems wrong, but we just compensate the FFT.  Now the RMS of
     * fridata+i*fiidata is the same as of operand. */
    lambda = sigma * msq * xres*yres;
    deconvolve_do(ffield, foper, ffield, cstride*yres, lambda);
    fftw_free(foper);
    //gwycreal(ffield[0]) = gwycimag(ffield[0]) = 0.0; //FIXME
    gwy_fftw_execute(bplan);
    fftw_destroy_plan(bplan);
    fftw_free(ffield);

    /* NB: We normalize it as an integral.  So one recovers the convolution
     * with TRUE in ext-convolve! */
    gwy_data_field_multiply(out, 1.0/(dfield->xreal * dfield->yreal));
    gwy_data_field_2dfft_humanize(out);

    out->xreal = dfield->xreal;
    out->yreal = dfield->yreal;
    out->xoff = dfield->xoff;
    out->yoff = dfield->yoff;

    gwy_data_field_invalidate(out);
    set_transfer_function_units(operand, dfield, out);
}

static void
psf_sigmaopt_estimate_size(PSFSigmaOptData *sodata)
{
    const fftw_complex *fideal = sodata->fideal;
    const fftw_complex *fmeas = sodata->fmeas;
    fftw_complex *fpsf = sodata->fpsf;
    fftw_complex *cbuf = sodata->cbuf;
    GwyDataField *psf = sodata->psf;
    gint xres = psf->xres, yres = psf->yres;
    gint i, j, imin, jmin, imax, jmax, ext;
    const gdouble *d = sodata->psf->data;
    guint n;
    gdouble lambda, m;

    n = yres*sodata->cstride;
    sodata->col = xres/3;
    sodata->row = yres/3;
    sodata->width = xres - 2*sodata->col;
    sodata->height = yres - 2*sodata->row;

    /* Use a fairly large but not yet insane sigma value 4.0 to estimate the width.  We want to err on the side of
     * size overestimation here.
     * XXX: We might want to use a proportional to 1/sqrt(xres*yres) here. */
    lambda = 4.0 * sodata->sigma_scale;
    deconvolve_do(fmeas, fideal, fpsf, n, lambda);
    gwy_assign(cbuf, fpsf, n);
    gwy_fftw_execute(sodata->pplan);
    gwy_data_field_2dfft_humanize(psf);

    imax = yres/2;
    jmax = xres/2;
    m = 0.0;
    for (i = sodata->row; i < sodata->row + sodata->height; i++) {
        for (j = sodata->col; j < sodata->col + sodata->width; j++) {
            if (fabs(d[i*xres + j]) > m) {
                m = fabs(d[i*xres + j]);
                imax = i;
                jmax = j;
            }
        }
    }
    gwy_debug("maximum at (%d,%d)", imax, jmax);
    gwy_data_field_threshold(psf, 0.05*m, 0.0, 1.0);
    g_return_if_fail(d[imax*xres + jmax] > 0.0);
    gwy_data_field_grains_extract_grain(psf, jmax, imax);

    imin = imax;
    jmin = jmax;
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++) {
            if (d[i*xres + j] > 0.0) {
                if (i < imin)
                    imin = i;
                if (i > imax)
                    imax = i;
                if (j < jmin)
                    jmin = j;
                if (j > jmax)
                    jmax = j;
            }
        }
    }

    ext = GWY_ROUND(0.5*log(xres*yres)) + 1;
    sodata->col = jmin - ext;
    sodata->row = imin - ext;
    sodata->width = jmax+1 - jmin + 2*ext;
    sodata->height = imax+1 - imin + 2*ext;
    if (sodata->col < 0) {
        sodata->width += sodata->col;
        sodata->col = 0;
    }
    if (sodata->row < 0) {
        sodata->height += sodata->row;
        sodata->row = 0;
    }
    if (sodata->col + sodata->width > xres)
        sodata->width = xres - sodata->col;
    if (sodata->row + sodata->height > yres)
        sodata->height = yres - sodata->row;

    gwy_debug("estimated region: %dx%d centered at (%d,%d)",
              sodata->width, sodata->height, sodata->col + sodata->width/2, sodata->row + sodata->height/2);
}

static void
psf_sigmaopt_prepare(GwyDataField *measured, GwyDataField *ideal,
                     PSFSigmaOptData *sodata)
{
    GwyDataField *buf;
    guint n, xres = measured->xres, yres = measured->yres;
    fftw_plan fplan;

    sodata->cstride = xres/2 + 1;
    n = yres*sodata->cstride;

    sodata->psf = buf = gwy_data_field_new_alike(measured, FALSE);
    sodata->fideal = g_new(fftw_complex, n);
    sodata->fmeas = g_new(fftw_complex, n);
    sodata->fpsf = g_new(fftw_complex, n);
    sodata->cbuf = gwy_fftw_new_complex(n);

    fplan = gwy_fftw_plan_dft_r2c_2d(yres, xres, buf->data, sodata->cbuf, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    /* FFT(ideal) */
    gwy_data_field_copy(ideal, buf, FALSE);
    sodata->sigma_scale = gwy_data_field_get_mean_square(buf);
    sodata->sigma_scale *= xres*yres;
    gwy_fftw_execute(fplan);
    gwy_assign(sodata->fideal, sodata->cbuf, n);

    /* FFT(measured) */
    gwy_data_field_copy(measured, buf, FALSE);
    gwy_fftw_execute(fplan);
    gwy_assign(sodata->fmeas, sodata->cbuf, n);

    fftw_destroy_plan(fplan);
    sodata->pplan = gwy_fftw_plan_dft_c2r_2d(yres, xres, sodata->cbuf, buf->data, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    psf_sigmaopt_estimate_size(sodata);
}

static void
psf_sigmaopt_free(PSFSigmaOptData *sodata)
{
    fftw_destroy_plan(sodata->pplan);
    fftw_free(sodata->cbuf);
    g_object_unref(sodata->psf);
    g_free(sodata->fpsf);
    g_free(sodata->fideal);
    g_free(sodata->fmeas);
}

static gdouble
psf_sigmaopt_residuum(gdouble logsigma, gpointer user_data)
{
    PSFSigmaOptData *sodata = (PSFSigmaOptData*)user_data;
    GwyDataField *psf = sodata->psf;
    const fftw_complex *fideal = sodata->fideal;
    const fftw_complex *fmeas = sodata->fmeas;
    fftw_complex *fpsf = sodata->fpsf;
    fftw_complex *cbuf = sodata->cbuf;
    gdouble lambda, w, sigma = exp(logsigma);
    guint n;

    /* Calculate frequency-space PSF sodata->fpsf and real-space sodata->psf. This essentially reproduces
     * gwy_data_field_deconvolve_regularized(). */
    n = psf->yres * sodata->cstride;
    lambda = sigma * sodata->sigma_scale;
    deconvolve_do(fmeas, fideal, fpsf, n, lambda);
    gwy_assign(cbuf, fpsf, n);
    gwy_fftw_execute(sodata->pplan);
    gwy_data_field_2dfft_humanize(psf);
    gwy_data_field_area_abs(psf, sodata->row, sodata->col, sodata->width, sodata->height);
    w = gwy_data_field_area_get_dispersion(psf, NULL, GWY_MASK_IGNORE,
                                           sodata->row, sodata->col, sodata->width, sodata->height,
                                           NULL, NULL);
    w = sqrt(w);
    return w;
}

/**
 * gwy_data_field_find_regularization_sigma_for_psf:
 * @dfield: A data field with convolved noisy data.
 * @ideal: A data field with ideal sharp data.
 *
 * Finds regularization parameter for point spread function calculation using regularized deconvolution.
 *
 * The estimated value should be suitable for reconstruction of the point spread function using
 * gwy_data_field_deconvolve_regularized().  The estimate is only suitable for PSF, it does not work for image
 * sharpening using a known PSF.
 *
 * Returns: Estimated regularization parameter.
 *
 * Since: 2.51
 **/
gdouble
gwy_data_field_find_regularization_sigma_for_psf(GwyDataField *dfield,
                                                 GwyDataField *ideal)
{
    PSFSigmaOptData sodata;
    gdouble logsigma;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(dfield), 0.0);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(ideal), 0.0);
    g_return_val_if_fail(!gwy_data_field_check_compatibility(dfield, ideal,
                                                             GWY_DATA_COMPATIBILITY_RES
                                                             | GWY_DATA_COMPATIBILITY_REAL
                                                             | GWY_DATA_COMPATIBILITY_LATERAL),
                         0.0);

    psf_sigmaopt_prepare(dfield, ideal, &sodata);
    logsigma = gwy_math_find_minimum_1d(psf_sigmaopt_residuum, log(1e-8), log(1e3), &sodata);
    psf_sigmaopt_free(&sodata);
    /* Experimentally determined fudge factor from large-scale simulations. */
    return 0.276*exp(logsigma);
}

/* xlen and ylen are dimensions of segments from the beginning. We assume the end part is one pixel shorter, as
 * expected for usual dehumanized images of odd sizes. */
static void
copy_corners(GwyDataField *source, GwyDataField *dest,
             gint xlen, gint ylen)
{
    gint sxres = source->xres;
    gint syres = source->yres;
    gint dxres = dest->xres;
    gint dyres = dest->yres;

    g_return_if_fail(2*xlen - 1 <= sxres);
    g_return_if_fail(2*ylen - 1 <= syres);
    g_return_if_fail(2*xlen - 1 <= dxres);
    g_return_if_fail(2*ylen - 1 <= dyres);

    gwy_data_field_clear(dest);

    /* Top left. */
    gwy_data_field_area_copy(source, dest, 0, 0, xlen, ylen, 0, 0);
    /* Top right. */
    gwy_data_field_area_copy(source, dest, sxres - (xlen-1), 0, xlen-1, ylen, dxres - (xlen-1), 0);
    /* Bottom left. */
    gwy_data_field_area_copy(source, dest, 0, syres - (ylen-1), xlen, ylen-1, 0, dyres - (ylen-1));
    /* Bottom right. */
    gwy_data_field_area_copy(source, dest, sxres - (xlen-1), syres - (ylen-1), xlen-1, ylen-1,
                             dxres - (xlen-1), dyres - (ylen-1));
}

static gdouble
dotprod_fields(GwyDataField *a, GwyDataField *b)
{
    gint xres = a->xres, yres = a->yres;
    const gdouble *da = a->data;
    const gdouble *db = b->data;
    gdouble s = 0.0;
    gint k;

    for (k = 0; k < xres*yres; k++)
        s += da[k]*db[k];

    return s;
}

static void
conjgrad_update_field(GwyDataField *res, GwyDataField *a, gdouble q, GwyDataField *b)
{
    gint xres = a->xres, yres = a->yres;
    gdouble *dr = res->data;
    const gdouble *da = a->data, *db = b->data;
    gint k;

    for (k = 0; k < xres*yres; k++)
        dr[k] = da[k] - q*db[k];
}

static void
conjgrad_matrix_multiply(const fftw_complex *fmat,
                         fftw_complex *fvec,
                         GwyDataField *vec,  /* vector-sized */
                         GwyDataField *buf,  /* matrix-sized */
                         fftw_plan fplan,
                         fftw_plan bplan,
                         GwyDataField *result)
{
    gint mxres = buf->xres, myres = buf->yres;
    gint vxres = vec->xres, vyres = vec->yres;
    gint cstride = mxres/2 + 1;
    gint k;

    /* Zero-extend the vector in the middle. */
    copy_corners(vec, buf, vxres/2+1, vyres/2+1);
    /* FFT(buf) -> fvec */
    gwy_fftw_execute(fplan);

    /* Multiply in frequency domain.  The matrix is symmetrical so ignore
     * the imaginary part. */
    for (k = 0; k < cstride*myres; k++) {
        gdouble mre = gwycreal(fmat[k]);
        gdouble vre = gwycreal(fvec[k]), vim = gwycimag(fvec[k]);
        gwycreal(fvec[k]) = mre*vre;
        gwycimag(fvec[k]) = mre*vim;
    }
    /* IFFT(fvec) -> buf */
    gwy_fftw_execute(bplan);
    /* Cut. */
    copy_corners(buf, result, vxres/2+1, vyres/2+1);
}

static gdouble
conjgrad_make_equations(GwyDataField *ideal, GwyDataField *measured,
                        gint txres, gint tyres,
                        GwyDataField *matrix, GwyDataField *rhs,
                        GwyDataField *buf)
{
    GwyDataField *autocor, *product;
    gint xres = ideal->xres;
    gint yres = ideal->yres;
    gint cstride = xres/2 + 1;
    fftw_complex *f = g_new(fftw_complex, cstride*yres);
    fftw_complex *g = g_new(fftw_complex, cstride*yres);
    gdouble *b = gwy_data_field_get_data(buf);
    fftw_plan fplan, bplan;
    gdouble fnorm;
    gint k;

    fplan = gwy_fftw_plan_dft_r2c_2d(yres, xres, b, f, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    bplan = gwy_fftw_plan_dft_c2r_2d(yres, xres, f, b, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

    product = gwy_data_field_new_alike(ideal, FALSE);
    autocor = gwy_data_field_new_alike(ideal, FALSE);

    gwy_data_field_copy(ideal, buf, FALSE);
    gwy_fftw_execute(fplan);
    gwy_assign(g, f, cstride*yres);
    for (k = 0; k < cstride*yres; k++) {
        gdouble re = gwycreal(f[k]), im = gwycimag(f[k]);
        gwycreal(f[k]) = re*re + im*im;
        gwycimag(f[k]) = 0.0;
    }
    gwy_fftw_execute(bplan);
    gwy_data_field_copy(buf, autocor, FALSE);
    gwy_data_field_multiply(autocor, 1.0/(xres*yres));

    gwy_data_field_copy(measured, buf, FALSE);
    gwy_fftw_execute(fplan);
    fftw_destroy_plan(fplan);
    for (k = 0; k < cstride*yres; k++) {
        gdouble fre = gwycreal(f[k]), fim = gwycimag(f[k]);
        gdouble gre = gwycreal(g[k]), gim = gwycimag(g[k]);
        gwycreal(f[k]) = gre*fre + gim*fim;
        gwycimag(f[k]) = gre*fim - gim*fre;
    }
    g_free(g);
    gwy_fftw_execute(bplan);
    fftw_destroy_plan(bplan);
    g_free(f);
    gwy_data_field_copy(buf, product, FALSE);
    gwy_data_field_multiply(product, 1.0/(xres*yres));

    /* Calculate the left hand side matrix.  This is the autocorrelation ideal⊛ideal.  We subsequently cut it to the
     * central part with dimensions (2txres-1)×(2tyres-1) and then possibly zero-extend slightly to implement matrix
     * multiplication by convolution.  This is done in single step by copying the four corners. */
    copy_corners(autocor, matrix, txres, tyres);
    g_object_unref(autocor);

    /* Calculate the right hand side vector.  This is the cross-correlation ideal⊛measured. Here we cut the result
     * just to txres×tyres size and zero-extend it. */
    copy_corners(product, rhs, txres/2+1, tyres/2+1);
    g_object_unref(product);

    /* The matrix and rhs data can be of wild orders of magnitude.  Scale them in sync to some reasonable range,
     * preferably a bit larger to have many orders of magnitude to go down. */
    fnorm = sqrt(gwy_data_field_get_mean_square(rhs));
    if (fnorm > 0.0) {
        fnorm = 1e24/fnorm;
        gwy_data_field_multiply(rhs, fnorm);
        gwy_data_field_multiply(matrix, fnorm);
    }
    else
        fnorm = 1.0;

    return fnorm;
}

static void
psf_least_squares_conjgrad_iterate(GwyDataField *ffield,
                                   GwyDataField *vfield,
                                   GwyDataField *wfield,
                                   GwyDataField *tf,
                                   GwyDataField *rhs,
                                   GwyDataField *buf,
                                   fftw_complex *fmat,
                                   fftw_complex *fvec,
                                   fftw_plan fplan,
                                   fftw_plan bplan,
                                   guint maxiter,
                                   gdouble eps)
{
    gdouble fnorm0, fnorm;
    guint iter;

    /* Initialise conjugated gradients.
     * f0 = v0 = matrix tf0 - rhs
     * w0 = matrix v0
     */
    conjgrad_matrix_multiply(fmat, fvec, tf, buf, fplan, bplan, vfield);
    gwy_data_field_subtract_fields(vfield, vfield, rhs);
    gwy_data_field_copy(vfield, ffield, FALSE);

    conjgrad_matrix_multiply(fmat, fvec, vfield, buf, fplan, bplan, wfield);
    fnorm0 = gwy_data_field_get_mean_square(ffield);

    /* Iterate. */
    for (iter = 0; iter < maxiter; iter++) {
        gdouble numer = dotprod_fields(ffield, vfield);
        gdouble denom = dotprod_fields(vfield, wfield);

        if (numer == 0.0 || denom == 0.0)
            break;

        /* Update solution p(j+1) = p(j) - n/d v(j) */
        conjgrad_update_field(tf, tf, numer/denom, vfield);
        /* Update residual f(j+1) = f(j) - n/d w(j) */
        conjgrad_update_field(ffield, ffield, numer/denom, wfield);
        /* Check norm of ffield and maybe stop iterations... */
        fnorm = gwy_data_field_get_mean_square(ffield);
        if (fnorm <= eps*fnorm0)
            break;
        /* Update gradient v(j+1) = f(j+1) - n/d v(j) */
        conjgrad_update_field(vfield, ffield, dotprod_fields(ffield, wfield)/denom, vfield);
        conjgrad_matrix_multiply(fmat, fvec, vfield, buf, fplan, bplan, wfield);
    }
}

/* Set up plans for FFT between buf (real) and fvec (complex). */
static gint
conjgrad_make_fft_plans(GwyDataField *buf,
                        fftw_complex **fmat, fftw_complex **fvec,
                        fftw_plan *fplan, fftw_plan *bplan,
                        gint xsize, gint ysize)
{
    gint cstride;

    /* Set up the frequency domain operations. */
    gwy_data_field_resample(buf, xsize, ysize, GWY_INTERPOLATION_NONE);
    cstride = xsize/2 + 1;
    *fmat = gwy_fftw_new_complex(cstride*ysize);
    *fvec = gwy_fftw_new_complex(cstride*ysize);
    *fplan = gwy_fftw_plan_dft_r2c_2d(ysize, xsize, buf->data, *fvec, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    *bplan = gwy_fftw_plan_dft_c2r_2d(ysize, xsize, *fvec, buf->data, FFTW_DESTROY_INPUT | FFTW_ESTIMATE);
    return cstride;
}

static gboolean
psf_least_squares_warn(gint *want_txres, gint *want_tyres, gint border,
                       gint xres, gint yres)
{
    gboolean ok = TRUE;

    if (!(*want_txres & 1)) {
        g_warning("Transfer function x-res is even.  Fixing to odd.");
        *want_txres |= 1;
        ok = FALSE;
    }
    if (!(*want_tyres & 1)) {
        g_warning("Transfer function y-res is even.  Fixing to odd.");
        *want_tyres |= 1;
        ok = FALSE;
    }

    if (*want_txres + 2*border > xres/2 || *want_tyres + 2*border > yres/2) {
        g_warning("Transfer function size is more than half image size. "
                  "I will proceed but the result will be rubbish.");
        ok = FALSE;
    }
    return ok;
}

/**
 * gwy_data_field_deconvolve_psf_leastsq:
 * @dfield: A data field.
 * @operand: Ideal sharp measurement (before convolution). It must have the same dimensions as @dfield and it is
 *           assumed it has also the same physical size.
 * @out: Output field for the transfer function.  Its dimensions are preserved and determine the transfer function
 *       support.  It must be smaller than half of @dfield.
 * @sigma: Regularization parameter.
 * @border: Number of pixel to extend and cut off the transfer function.
 *
 * Performs reconstruction of transfer function from convolved and ideal sharp images.
 *
 * The transfer function is reconstructed by solving the corresponding least squares problem.  This method is suitable
 * when the dimensions of @out are much smaller than the images.
 *
 * Since the method accumulates errors close to edges, they can be removed within the procedure by reconstructing
 * a slightly larger transfer function and then cutting the result.  The extension is given by @border, typical
 * suitable values are 2 or 3.
 *
 * Convolving the result with the operand using gwy_data_field_area_ext_convolve() with @as_integral=%TRUE will
 * recover (approximately) the image.  This means the deconvolution assumes continous convolution, not discrete sums.
 * Note that for the latter case this means the point spread function will be centered in @out.
 *
 * Fields @dfield and @operand should be windowed beforehand if they are not periodic.
 *
 * Since: 2.52
 **/
void
gwy_data_field_deconvolve_psf_leastsq(GwyDataField *measured,
                                      GwyDataField *ideal,
                                      GwyDataField *tf,
                                      gdouble sigma,
                                      gint border)
{
    gint xres, yres, want_txres, want_tyres;
    gint txres, tyres, xsize, ysize, cstride;
    gdouble dx, dy, q, mnorm, fnorm, lambda;
    GwyDataField *buf, *matrix, *rhs, *ffield, *vfield, *wfield;
    fftw_complex *fmat, *fvec;
    fftw_plan fplan, bplan;

    g_return_if_fail(GWY_IS_DATA_FIELD(measured));
    g_return_if_fail(GWY_IS_DATA_FIELD(ideal));
    g_return_if_fail(GWY_IS_DATA_FIELD(tf));
    xres = measured->xres;
    yres = measured->yres;
    g_return_if_fail(ideal->xres == xres);
    g_return_if_fail(ideal->yres == yres);
    g_return_if_fail(border >= 0);

    want_txres = tf->xres;
    want_tyres = tf->yres;
    dx = gwy_data_field_get_dx(measured);
    dy = gwy_data_field_get_dy(measured);
    q = dx*dy;

    want_txres = tf->xres;
    want_tyres = tf->yres;
    psf_least_squares_warn(&want_txres, &want_tyres, border, xres, yres);
    txres = want_txres + 2*border;
    tyres = want_tyres + 2*border;
    xsize = gwy_fft_find_nice_size(2*txres-1);
    ysize = gwy_fft_find_nice_size(2*tyres-1);

    gwy_data_field_resample(tf, txres, tyres, GWY_INTERPOLATION_NONE);
    gwy_data_field_clear(tf);

    buf = gwy_data_field_new_alike(ideal, FALSE);
    matrix = gwy_data_field_new(xsize, ysize, xsize, ysize, FALSE);
    rhs = gwy_data_field_new_alike(tf, FALSE);
    fnorm = conjgrad_make_equations(ideal, measured, txres, tyres, matrix, rhs, buf);

    /* Normalise sigma relatively to data. */
    /* Element matrix_00 is (unnormalised) ACF at zero distance, i.e. sum of squared data.  This is a good scaling
     * norm for the matrix we can use instead of actually calculating some matrix norm -- its elements are kind of
     * proportional to matrix_00 around the origin anyway, so this would add just some odd multiplicative constant. */
    mnorm = gwy_data_field_get_mean_square(ideal);
    /* We should multiply with the number of columns (not rows!) of matrix implementing the convolution operator
     * ideal*, see Golub's book. However, we then want to divide it by the square root of that to make sigma more or
     * less resolution-independent.  */
    mnorm *= xres*yres;
    /* This compensates different scaling of LSM with transfer function support size than the other two methods.  */
    mnorm *= cbrt(want_txres*want_tyres);
    lambda = sigma * mnorm * fnorm;
    gwy_data_field_set_val(matrix, 0, 0, gwy_data_field_get_val(matrix, 0, 0) + lambda);

    /* Set up the frequency domain operations. */
    gwy_data_field_resample(buf, xsize, ysize, GWY_INTERPOLATION_NONE);
    cstride = conjgrad_make_fft_plans(buf, &fmat, &fvec, &fplan, &bplan, xsize, ysize);

    gwy_data_field_copy(matrix, buf, FALSE);
    GWY_OBJECT_UNREF(matrix);
    /* Compensate unnormalised FFTs here and in conjgrad_matrix_multiply() */
    gwy_data_field_multiply(buf, 1.0/(xsize*ysize));
    gwy_fftw_execute(fplan);
    gwy_assign(fmat, fvec, cstride*ysize);

    vfield = gwy_data_field_new(txres, tyres, txres, tyres, FALSE);
    ffield = gwy_data_field_new_alike(vfield, FALSE);
    wfield = gwy_data_field_new_alike(vfield, FALSE);

    psf_least_squares_conjgrad_iterate(ffield, vfield, wfield, tf, rhs, buf, fmat, fvec, fplan, bplan, 150, 1e-40);
    fftw_free(fvec);
    fftw_free(fmat);

    gwy_data_field_2dfft_humanize(tf);
    if (border > 0)
        gwy_data_field_resize(tf, border, border, want_txres + border, want_tyres + border);
    gwy_data_field_multiply(tf, 1.0/q);
    gwy_data_field_set_xreal(tf, want_txres*dx);
    gwy_data_field_set_yreal(tf, want_tyres*dy);
    gwy_data_field_set_xoffset(tf, -0.5*want_txres*dx);
    gwy_data_field_set_yoffset(tf, -0.5*want_tyres*dy);
    set_transfer_function_units(ideal, measured, tf);

    g_object_unref(buf);
    g_object_unref(rhs);
    g_object_unref(ffield);
    g_object_unref(vfield);
    g_object_unref(wfield);
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);
}

static gdouble
psf_conjgrad_residuum(gdouble logsigma, gpointer user_data)
{
    PSFConjGradData *cgdata = (PSFConjGradData*)user_data;
    gdouble lambda, sigma, width;
    gint want_txres = cgdata->want_txres, want_tyres = cgdata->want_tyres;
    gint xsize = cgdata->xsize, ysize = cgdata->ysize;
    gint border = cgdata->border;

    sigma = exp(logsigma);

    /* cgdata->mnorm already includes all the scaling factors, so just multiply
     * it with sigma. */
    lambda = sigma * cgdata->mnorm;
    gwy_data_field_set_val(cgdata->matrix, 0, 0, cgdata->matrix00 + lambda);

    gwy_data_field_copy(cgdata->matrix, cgdata->buf, FALSE);
    /* Compensate unnormalised FFTs here and in conjgrad_matrix_multiply() */
    gwy_data_field_multiply(cgdata->buf, 1.0/(xsize*ysize));
    gwy_fftw_execute(cgdata->fplan);
    gwy_assign(cgdata->fmat, cgdata->fvec, cgdata->cstride*ysize);

    gwy_data_field_clear(cgdata->tf);
    psf_least_squares_conjgrad_iterate(cgdata->ffield, cgdata->vfield, cgdata->wfield, cgdata->tf,
                                       cgdata->rhs, cgdata->buf, cgdata->fmat, cgdata->fvec,
                                       cgdata->fplan, cgdata->bplan,
                                       60, 1e-24);
    gwy_data_field_2dfft_humanize(cgdata->tf);

    gwy_data_field_abs(cgdata->tf);
    width = gwy_data_field_area_get_dispersion(cgdata->tf, NULL, GWY_MASK_IGNORE,
                                               border, border, want_txres, want_tyres, NULL, NULL);
    return sqrt(width);
}

/**
 * gwy_data_field_find_regularization_sigma_leastsq:
 * @dfield: A data field with convolved noisy data.
 * @ideal: A data field with ideal sharp data.
 * @width: Horizontal size of transfer function support.
 * @height: Vertical size of transfer function support.
 * @border: Number of pixel to extend and cut off the transfer function.
 *
 * Finds regularization parameter for point spread function calculation using least squares method.
 *
 * The estimated value should be suitable for reconstruction of the point spread function using
 * gwy_data_field_deconvolve_psf_leastsq().
 *
 * Returns: Estimated regularization parameter.
 *
 * Since: 2.52
 **/
gdouble
gwy_data_field_find_regularization_sigma_leastsq(GwyDataField *measured,
                                                 GwyDataField *ideal,
                                                 gint want_txres,
                                                 gint want_tyres,
                                                 gint border)
{
    gdouble logsigma, dx, dy, fnorm;
    gint xres, yres, txres, tyres, xsize, ysize;
    PSFConjGradData cgdata;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(measured), 0.001);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(ideal), 0.001);
    xres = measured->xres;
    yres = measured->yres;
    g_return_val_if_fail(ideal->xres == xres, 0.001);
    g_return_val_if_fail(ideal->yres == yres, 0.001);
    g_return_val_if_fail(want_txres > 0, 0.001);
    g_return_val_if_fail(want_tyres > 0, 0.001);
    g_return_val_if_fail(border >= 0, 0.001);
    dx = measured->xreal/xres;
    dy = measured->yreal/yres;

    psf_least_squares_warn(&want_txres, &want_tyres, border, xres, yres);
    cgdata.want_txres = want_txres;
    cgdata.want_tyres = want_tyres;
    cgdata.border = border;
    cgdata.txres = txres = want_txres + 2*border;
    cgdata.tyres = tyres = want_tyres + 2*border;
    cgdata.xsize = xsize = gwy_fft_find_nice_size(2*cgdata.txres-1);
    cgdata.ysize = ysize = gwy_fft_find_nice_size(2*cgdata.tyres-1);
    cgdata.tf = gwy_data_field_new(txres, tyres, dx*txres, dy*tyres, TRUE);
    cgdata.buf = gwy_data_field_new_alike(ideal, FALSE);
    cgdata.matrix = gwy_data_field_new(xsize, ysize, xsize, ysize, FALSE);
    cgdata.rhs = gwy_data_field_new_alike(cgdata.tf, FALSE);
    fnorm = conjgrad_make_equations(ideal, measured, txres, tyres, cgdata.matrix, cgdata.rhs, cgdata.buf);
    /* See gwy_data_field_deconvolve_psf_leastsq() for explanation. */
    cgdata.mnorm = gwy_data_field_get_mean_square(ideal);
    cgdata.mnorm *= xres*yres;
    cgdata.mnorm *= cbrt(want_txres*want_tyres);
    cgdata.mnorm *= fnorm;
    cgdata.matrix00 = gwy_data_field_get_val(cgdata.matrix, 0, 0);

    /* Set up the frequency domain operations. */
    gwy_data_field_resample(cgdata.buf, cgdata.xsize, cgdata.ysize,
                            GWY_INTERPOLATION_NONE);
    cgdata.cstride = conjgrad_make_fft_plans(cgdata.buf, &cgdata.fmat, &cgdata.fvec, &cgdata.fplan, &cgdata.bplan,
                                             xsize, ysize);

    cgdata.vfield = gwy_data_field_new(txres, tyres, txres, tyres, FALSE);
    cgdata.ffield = gwy_data_field_new_alike(cgdata.vfield, FALSE);
    cgdata.wfield = gwy_data_field_new_alike(cgdata.vfield, FALSE);

    logsigma = gwy_math_find_minimum_1d(psf_conjgrad_residuum, log(1e-8), log(1e3), &cgdata);

    fftw_free(cgdata.fvec);
    fftw_free(cgdata.fmat);
    g_object_unref(cgdata.tf);
    g_object_unref(cgdata.buf);
    g_object_unref(cgdata.rhs);
    g_object_unref(cgdata.matrix);
    g_object_unref(cgdata.ffield);
    g_object_unref(cgdata.vfield);
    g_object_unref(cgdata.wfield);
    fftw_destroy_plan(cgdata.fplan);
    fftw_destroy_plan(cgdata.bplan);

    /* Experimentally determined fudge factor from large-scale simulations. */
    return 0.298*exp(logsigma);
}

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
