/*
 *  $Id: mfm.c 22422 2019-08-20 12:17:42Z klapetek $
 *  Copyright (C) 2016-2018 David Necas (Yeti), Petr Klapetek.
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
#include <libgwyddion/gwymath.h>
#include <libprocess/arithmetic.h>
#include <libprocess/stats.h>
#include <libprocess/filters.h>
#include <libprocess/inttrans.h>
#include <libprocess/mfm.h>

#define MU_0 1.256637061435917295e-6
#define EPSILON_0 8.854187817620389850e-12

typedef struct {
    GwyDataField *ztf;
    GwyDataField *freq;
    GwyDataField *rea;
    GwyDataField *ima;
    GwyDataField *reb;
    GwyDataField *imb;
} GwyMFMZShiftData;

/* Calculate a data field with frequency magnitudes in FFT (not humanized)
 * arrangement. */
static void
precaulcate_frequency_field(GwyDataField *model,
                            GwyDataField *freq)
{
    gint i, j, xres = model->xres, yres = model->yres;
    gdouble kx, ky, sx, sy, r;
    gdouble *f = freq->data;

    sx = 1.0/model->xreal;
    sy = 1.0/model->yreal;
    f[0] = 0.0;

    for (j = 1; j <= xres/2; j++) {
        kx = j*sx;
        r = kx;
        f[j] = r;
        f[xres-j] = r;
    }

    for (i = 1; i <= yres/2; i++) {
        ky = i*sy;
        r = ky;
        f[(yres-i)*xres] = r;
        f[i*xres] = r;
    }

    for (i = 1; i <= yres/2; i++) {
        ky = i*sy;
        for (j = 1; j <= xres/2; j++) {
            kx = j*sx;
            r = sqrt(kx*kx + ky*ky);
            f[(yres-i)*xres + xres-j] = r;
            f[i*xres + xres-j] = r;
            f[(yres-i)*xres + j] = r;
            f[i*xres + j] = r;
        }
    }

    gwy_data_field_invalidate(freq);
}

/* Calculate two data fields with frequency components in FFT (not humanized)
 * arrangement.
 * XXX: Duplicate with psf-fit.c.  Make public?  */
static void
precaulcate_xy_frequency_fields(GwyDataField *model,
                                GwyDataField *freq_x,
                                GwyDataField *freq_y)
{
    guint xres = model->xres, yres = model->yres;
    gdouble sx = 1.0/model->xreal, sy = 1.0/model->yreal;
    gdouble *fx = freq_x->data, *fy = freq_y->data;
    gdouble vx, vy;
    guint i, j;

    fx[0] = fy[0] = 0.0;

    for (j = 1; j <= xres/2; j++) {
        vx = j*sx;
        fx[xres-j] = -vx;
        fx[j] = vx;
        fy[j] = fy[xres-j] = 0.0;
    }

    for (i = 1; i <= yres/2; i++) {
        vy = i*sy;
        fx[i*xres] = fx[(yres-i)*xres] = 0.0;
        fy[(yres-i)*xres] = -vy;
        fy[i*xres] = vy;
    }

    for (i = 1; i <= yres/2; i++) {
        vy = i*sy;
        for (j = 1; j <= xres/2; j++) {
            vx = j*sx;
            fx[(yres-i)*xres + xres-j] = -vx;
            fx[i*xres + xres-j] = -vx;
            fx[(yres-i)*xres + j] = vx;
            fx[i*xres + j] = vx;
            fy[(yres-i)*xres + xres-j] = -vy;
            fy[(yres-i)*xres + j] = -vy;
            fy[i*xres + xres-j] = vy;
            fy[i*xres + j] = vy;
        }
    }

    gwy_data_field_invalidate(freq_x);
    gwy_data_field_invalidate(freq_y);
}

/* Petr TODO: Use mirror extension for domain walls */
#if 0
static void
mirror_extend_for_filtering(GwyDataField *source,
                            gint col, gint row,
                            gint width, gint height,
                            GwyDataField *target)
{
    gint xres, yres, i, j;
    const gdouble *srow;
    gdouble *trow;

    g_return_if_fail(GWY_IS_DATA_FIELD(source));
    g_return_if_fail(GWY_IS_DATA_FIELD(target));
    xres = source->xres;
    yres = source->yres;
    g_assert(col >= 0 && width > 0 && col + width <= xres);
    g_assert(row >= 0 && height > 0 && row + height <= yres);
    /* Maybe we should resize automatically?  But the caller probably does not
     * expect @target size to change, so better fail noisily. */
    g_assert(target->xres == 2*width);
    g_assert(target->yres == 2*height);

    for (i = 0; i < height; i++) {
        /* Direct copy of the area. */
        srow = source->data + (row + i)*xres + col;
        trow = target->data + i*2*width;
        gwy_assign(trow, srow, width);

        /* L-R mirrored image to the right. */
        srow = trow;
        trow += 2*width-1;
        for (j = width; j; j--, srow++, trow--)
            *trow = *srow;

        /* And the entire thing mirrored upside down in the bottom part. */
        srow = target->data + i*2*width;
        trow = target->data + (2*height-1 - i)*2*width;
        gwy_assign(trow, srow, 2*width);
    }
}
#endif

static void
mfm_perpendicular_create_wall_mask(GwyDataField *wm, gdouble delta)
{
    gint xres, yres, i, j;
    gdouble dx, dy, x, y, s, sum, *data;

    g_return_if_fail(GWY_IS_DATA_FIELD(wm));
    xres = wm->xres;
    yres = wm->yres;
    dx = wm->xreal/xres;
    dy = wm->yreal/yres;
    data = wm->data;

    sum = 0.0;
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++, data++) {
            x = (i - xres/2)*dx;
            y = (j - yres/2)*dy;
            s = sqrt(x*x + y*y);
            sum += *data = 1.0/(cosh(G_PI*s/delta)*cosh(G_PI*s/delta));
        }
    }
    gwy_data_field_invalidate(wm);
    gwy_data_field_multiply(wm, 1.0/sum);
}

static void
mfm_perpendicular_create_ftf(GwyDataField *ftf, GwyDataField *buf,
                             gdouble mtip,
                             gdouble bx, gdouble by, gdouble length,
                             GwyMFMProbeType type)
{
    guint i, n;
    gdouble c, k;
    gdouble *data, *bdata;

    g_return_if_fail(GWY_IS_DATA_FIELD(ftf));

    c = -MU_0*mtip*bx*by;
    if (type == GWY_MFM_PROBE_CHARGE) {
        gwy_data_field_fill(ftf, c);
        return;
    }
    precaulcate_xy_frequency_fields(ftf, ftf, buf);
    n = ftf->xres * ftf->yres;
    data = ftf->data;
    bdata = buf->data;
    for (i = 0; i < n; i++) {
        k = sqrt(data[i]*data[i] + bdata[i]*bdata[i]);
        data[i] = c * gwy_sinc(data[i]*bx/2) * gwy_sinc(bdata[i]*by/2)
                  * (1.0 - exp(-k*length));
    }
    gwy_data_field_invalidate(ftf);
}

static void
mfm_create_ztf_from_frequencies(GwyDataField *ztf, GwyDataField *freq,
                                gdouble zdiff)
{
    guint n, i;
    gdouble *data, *fdata;

    g_return_if_fail(GWY_IS_DATA_FIELD(ztf));
    g_return_if_fail(GWY_IS_DATA_FIELD(freq));
    data = ztf->data;
    fdata = freq->data;
    n = ztf->xres * ztf->yres;
    for (i = 0; i < n; i++)
        data[i] = exp(2*M_PI*fdata[i]*zdiff);
    gwy_data_field_invalidate(ztf);
}

static void
mfm_create_ztf(GwyDataField *ztf, gdouble zdiff)
{
    precaulcate_frequency_field(ztf, ztf);
    mfm_create_ztf_from_frequencies(ztf, ztf, zdiff);
}

static void
mfm_perpendicular_create_field_mask(GwyDataField *fieldmask,
                                    gdouble height,
                                    gdouble thickness)
{
    guint n, i;
    gdouble *data;

    g_return_if_fail(GWY_IS_DATA_FIELD(fieldmask));
    precaulcate_frequency_field(fieldmask, fieldmask);
    data = fieldmask->data;
    n = fieldmask->xres*fieldmask->yres;
    for (i = 0; i < n; i++)
        data[i] = 0.5*exp(-2*M_PI*data[i]*height)*(1.0 - exp(-2*M_PI*data[i]*thickness)); //doubts about 0.5 before meff/hz split
    gwy_data_field_invalidate(fieldmask);
}

/**
 * gwy_data_field_mfm_perpendicular_stray_field:
 * @mfield: Mask representing the magnetisation orientation.
 * @out: Target data field to put the result to.  It will be resized to match
 *       @mfield.
 * @height: Height above the surface.
 * @thickness: Film thickness.
 * @sigma: Magnetic charge.
 * @walls: Include domain walls.
 * @wall_delta: Domain wall thickness
 *
 * Calculates stray field for perpendicular media, based on a mask showing the
 * magnetisation orientation.
 *
 * Since: 2.51
 **/
void
gwy_data_field_mfm_perpendicular_stray_field(GwyDataField *mfield,
                                             GwyDataField *out,
                                             gdouble height,
                                             gdouble thickness,
                                             gdouble sigma,
                                             gboolean walls,
                                             gdouble wall_delta)
{
    GwyDataField *rea, *ima, *reb, *imb, *fieldmask, *wallmask = NULL;

    g_return_if_fail(GWY_IS_DATA_FIELD(mfield));
    g_return_if_fail(GWY_IS_DATA_FIELD(out));

    rea = gwy_data_field_new_alike(mfield, TRUE);
    reb = gwy_data_field_new_alike(mfield, TRUE);
    ima = gwy_data_field_new_alike(mfield, TRUE);
    imb = gwy_data_field_new_alike(mfield, TRUE);
    fieldmask = gwy_data_field_new_alike(mfield, TRUE);

    gwy_data_field_copy(mfield, rea, FALSE);
    gwy_data_field_multiply(rea, 2*sigma);
    gwy_data_field_add(rea, -sigma);

    if (walls) {
        wallmask = gwy_data_field_new_alike(mfield, TRUE);
        mfm_perpendicular_create_wall_mask(wallmask, wall_delta);
        gwy_data_field_area_ext_convolve(rea,
                                         0, 0,
                                         gwy_data_field_get_xres(rea),
                                         gwy_data_field_get_yres(rea),
                                         rea, wallmask,
                                         GWY_EXTERIOR_MIRROR_EXTEND, 0.0,
                                         FALSE);
    }


    gwy_data_field_2dfft_raw(rea, NULL, reb, imb,
                             GWY_TRANSFORM_DIRECTION_FORWARD);
    mfm_perpendicular_create_field_mask(fieldmask, height, thickness);
    gwy_data_field_multiply_fields(reb, reb, fieldmask);
    gwy_data_field_multiply_fields(imb, imb, fieldmask);
    gwy_data_field_2dfft_raw(reb, imb, rea, ima,
                             GWY_TRANSFORM_DIRECTION_BACKWARD);
    gwy_data_field_resample(out, mfield->xres, mfield->yres,
                            GWY_INTERPOLATION_NONE);

    gwy_data_field_copy(rea, out, TRUE);

    gwy_si_unit_set_from_string(gwy_data_field_get_si_unit_z(out), "A/m");

    g_object_unref(rea);
    g_object_unref(reb);
    g_object_unref(ima);
    g_object_unref(imb);
    g_object_unref(fieldmask);
    GWY_OBJECT_UNREF(wallmask);
}


/**
 * gwy_data_field_mfm_perpendicular_stray_field_angle_correction:
 * @field: Field to be processed. It will be changed by the correction.
 * @out: Cantilever angle in degrees.
 * @orientation: Cantilever orientation with respect of the data.
 *
 * Performs correction of magnetic data for cantilever tilt.
 *
 * Since: 2.54
 **/
void
gwy_data_field_mfm_perpendicular_stray_field_angle_correction(GwyDataField *field,
                                                      gdouble angle,
                                                      GwyOrientation orientation)
{
    GwyDataField *rea, *ima, *reb, *imb, *kx, *ky;
    gdouble *rd, *id, *kxdata, *kydata, kval, kdirval, theta, ctheta, stheta;
    gint n, i;
    gdouble a, b, c, d;

    g_return_if_fail(GWY_IS_DATA_FIELD(field));

    rea = gwy_data_field_new_alike(field, TRUE);
    reb = gwy_data_field_new_alike(field, TRUE);
    ima = gwy_data_field_new_alike(field, TRUE);
    imb = gwy_data_field_new_alike(field, TRUE);
    kx = gwy_data_field_new_alike(field, FALSE);
    ky = gwy_data_field_new_alike(field, FALSE);

    gwy_data_field_copy(field, rea, FALSE);
    gwy_data_field_2dfft_raw(rea, NULL, reb, imb,
                             GWY_TRANSFORM_DIRECTION_FORWARD);

    precaulcate_xy_frequency_fields(rea, kx, ky);
    kxdata = kx->data;
    kydata = ky->data;
    rd = reb->data;
    id = imb->data;

    theta = angle*M_PI/180;
    ctheta = cos(theta);
    stheta = sin(theta);
    n = reb->xres*reb->yres;
    for (i = 0; i < n; i++) {

        a = rd[i];
        b = id[i];
        kval = sqrt(kxdata[i]*kxdata[i] + kydata[i]*kydata[i]);

        if (orientation == GWY_ORIENTATION_HORIZONTAL)
            kdirval = kxdata[i];
        else
            kdirval = kydata[i];

        if (kval == 0)
            c = ctheta*ctheta;
        else
            c = ctheta*ctheta - kdirval*kdirval*stheta*stheta/kval/kval;

        if (kval == 0)
            d = 0;
        else
            d = 2*kdirval*stheta*ctheta/kval;

        rd[i] = a*c - b*d;
        id[i] = a*d + b*c;
    }

    gwy_data_field_2dfft_raw(reb, imb, rea, ima,
                             GWY_TRANSFORM_DIRECTION_BACKWARD);
    gwy_data_field_copy(rea, field, TRUE);

    g_object_unref(rea);
    g_object_unref(reb);
    g_object_unref(ima);
    g_object_unref(imb);
    g_object_unref(kx);
    g_object_unref(ky);
}


/**
 * gwy_data_field_mfm_perpendicular_medium_force:
 * @hz: Data field contaning the Z-component of the magnetic H field.
 * @fz: Target data field to put the result to.  It will be resized to match
 *       @hz.
 * @type: Probe type.
 * @mtip: Probe magnetic moment.
 * @bx: x size for parallelpiped probe.
 * @by: y size for parallelpiped probe.
 * @length: length (z size) for parallelpiped probe.
 *
 * Calculates force as evaluated from z-component of the magnetic field for a given probe type.
 *
 * Since: 2.51
 **/
void
gwy_data_field_mfm_perpendicular_medium_force(GwyDataField *hz,
                                              GwyDataField *fz,
                                              GwyMFMProbeType type,
                                              gdouble mtip,
                                              gdouble bx,
                                              gdouble by,
                                              gdouble length)
{
    GwyDataField *rea, *ima, *reb, *imb, *ftf;
    GwySIUnit *unit, *unit2;

    g_return_if_fail(GWY_IS_DATA_FIELD(hz));
    g_return_if_fail(GWY_IS_DATA_FIELD(fz));

    rea = gwy_data_field_new_alike(hz, TRUE);
    reb = gwy_data_field_new_alike(hz, TRUE);
    ima = gwy_data_field_new_alike(hz, TRUE);
    imb = gwy_data_field_new_alike(hz, TRUE);
    ftf = gwy_data_field_new_alike(hz, TRUE);

    gwy_data_field_copy(hz, rea, FALSE);
    gwy_data_field_2dfft_raw(rea, NULL, reb, imb,
                             GWY_TRANSFORM_DIRECTION_FORWARD);
    gwy_data_field_resample(fz, hz->xres, hz->yres, GWY_INTERPOLATION_NONE);
    mfm_perpendicular_create_ftf(ftf, fz, mtip, bx, by, length, type);
    gwy_data_field_multiply_fields(reb, reb, ftf);
    gwy_data_field_multiply_fields(imb, imb, ftf);
    gwy_data_field_2dfft_raw(reb, imb, rea, ima,
                             GWY_TRANSFORM_DIRECTION_BACKWARD);
    gwy_data_field_copy(rea, fz, TRUE);

    g_object_unref(rea);
    g_object_unref(reb);
    g_object_unref(ima);
    g_object_unref(imb);
    g_object_unref(ftf);

    /* Set the units by trying the exact expected units first.  When we fail,
     * simply multiply the units by J/A, which is the dimension factor between
     * N and A/m. */
    unit = gwy_data_field_get_si_unit_z(fz);
    unit2 = gwy_si_unit_new("A/m");
    if (gwy_si_unit_equal(unit, unit2))
        gwy_si_unit_set_from_string(unit, "N");
    else {
        gwy_si_unit_set_from_string(unit2, "A/m^2");
        if (gwy_si_unit_equal(unit, unit2))
            gwy_si_unit_set_from_string(unit, "N/m");
        else {
            gwy_si_unit_set_from_string(unit2, "A/m^3");
            if (gwy_si_unit_equal(unit, unit2))
                gwy_si_unit_set_from_string(unit, "N/m^2");
            else {
                gwy_si_unit_set_from_string(unit2, "J/A");
                gwy_si_unit_multiply(unit, unit2, unit);
            }
        }
    }
    g_object_unref(unit2);
}

/**
 * gwy_data_field_mfm_shift_z:
 * @dfield: Data field containing magnetic field component.
 * @out: Target data field to put the result to.
 * @zdiff: The shift distance in physical units.
 *
 * Shifts magnetic field to a different lift height above the surface.
 *
 * Positive @zdiff means away from the measured surface and blurring the data.
 * Negative @zdiff means shifting towards (or within) the measured surface
 * and sharpening the data.  For negative @zdiff the result grows exponentially
 * and is generally not very useful.
 *
 * Since: 2.51
 **/
void
gwy_data_field_mfm_shift_z(GwyDataField *dfield,
                           GwyDataField *out,
                           gdouble zdiff)
{
    GwyDataField *rea, *ima, *reb, *imb, *ztf;

    g_return_if_fail(GWY_IS_DATA_FIELD(dfield));
    g_return_if_fail(GWY_IS_DATA_FIELD(out));

    rea = gwy_data_field_new_alike(dfield, TRUE);
    reb = gwy_data_field_new_alike(dfield, TRUE);
    ima = gwy_data_field_new_alike(dfield, TRUE);
    imb = gwy_data_field_new_alike(dfield, TRUE);
    ztf = gwy_data_field_new_alike(dfield, TRUE);

    gwy_data_field_copy(dfield, rea, FALSE);

    gwy_data_field_2dfft_raw(rea, NULL, reb, imb,
                             GWY_TRANSFORM_DIRECTION_FORWARD);
    mfm_create_ztf(ztf, zdiff);
    gwy_data_field_multiply_fields(reb, reb, ztf);
    gwy_data_field_multiply_fields(imb, imb, ztf);
    gwy_data_field_2dfft_raw(reb, imb, rea, ima,
                             GWY_TRANSFORM_DIRECTION_BACKWARD);
    gwy_data_field_resample(out, dfield->xres, dfield->yres,
                            GWY_INTERPOLATION_NONE);
    gwy_data_field_copy(rea, out, TRUE);

    g_object_unref(rea);
    g_object_unref(reb);
    g_object_unref(ima);
    g_object_unref(imb);
    g_object_unref(ztf);
}

static gdouble
square_sum_diff_from_product(GwyDataField *pfield,
                             GwyDataField *arefield,
                             GwyDataField *aimfield,
                             GwyDataField *brefield,
                             GwyDataField *bimfield)
{
    const gdouble *p = pfield->data;
    const gdouble *are = arefield->data, *aim = aimfield->data;
    const gdouble *bre = brefield->data, *bim = bimfield->data;
    guint i, n = pfield->xres * pfield->yres;
    gdouble dre, dim, s = 0.0;

    for (i = 0; i < n; i++) {
        dre = p[i]*are[i] - bre[i];
        dim = p[i]*aim[i] - bim[i];
        s += dre*dre + dim*dim;
    }

    return s;
}

static gdouble
zshift_residuum(gdouble zshift, gpointer user_data)
{
    const GwyMFMZShiftData *zsdata = (const GwyMFMZShiftData*)user_data;

    mfm_create_ztf_from_frequencies(zsdata->ztf, zsdata->freq, zshift);
    return square_sum_diff_from_product(zsdata->ztf,
                                        zsdata->rea, zsdata->ima,
                                        zsdata->reb, zsdata->imb);
}

/**
 * gwy_data_field_mfm_find_shift_z:
 * @dfield: Data field containing magnetic field component.
 * @shifted: Data field containing magnetic field component measured at
 *           a different lift height.
 * @zdiffmin: Start of shift scan range.
 * @zdiffmax: Start of shift scan range.
 *
 * Estimates the height difference between two magnetic field images.
 *
 * See gwy_data_field_mfm_shift_z() for the sign convention.  It is generally
 * only meaningful to estimate the shift whe @shifted was measured at larger
 * lift height than @dfield.
 *
 * Returns: The estimated shift between @shifted and @dfield.
 *
 * Since: 2.51
 **/
gdouble
gwy_data_field_mfm_find_shift_z(GwyDataField *dfield,
                                GwyDataField *shifted,
                                gdouble zdiffmin,
                                gdouble zdiffmax)
{
    GwyMFMZShiftData zsdata;
    gdouble zshift;

    zshift = 0.5*(zdiffmin + zdiffmax);

    g_return_val_if_fail(GWY_IS_DATA_FIELD(dfield), zshift);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(shifted), zshift);
    g_return_val_if_fail(!gwy_data_field_check_compatibility
                                                 (dfield, shifted,
                                                  GWY_DATA_COMPATIBILITY_RES),
                         zshift);

    zsdata.rea = gwy_data_field_new_alike(dfield, TRUE);
    zsdata.reb = gwy_data_field_new_alike(dfield, TRUE);
    zsdata.ima = gwy_data_field_new_alike(dfield, TRUE);
    zsdata.imb = gwy_data_field_new_alike(dfield, TRUE);
    zsdata.ztf = gwy_data_field_new_alike(dfield, TRUE);
    zsdata.freq = gwy_data_field_new_alike(dfield, TRUE);

    precaulcate_frequency_field(dfield, zsdata.freq);
    gwy_data_field_2dfft_raw(dfield, NULL, zsdata.rea, zsdata.ima,
                             GWY_TRANSFORM_DIRECTION_FORWARD);
    gwy_data_field_2dfft_raw(shifted, NULL, zsdata.reb, zsdata.imb,
                             GWY_TRANSFORM_DIRECTION_FORWARD);

    zshift = gwy_math_find_minimum_1d(zshift_residuum, zdiffmin, zdiffmax,
                                      &zsdata);

    g_object_unref(zsdata.rea);
    g_object_unref(zsdata.reb);
    g_object_unref(zsdata.ima);
    g_object_unref(zsdata.imb);
    g_object_unref(zsdata.ztf);
    g_object_unref(zsdata.freq);

    return zshift;
}

/**
 * gwy_data_field_mfm_parallel_medium:
 * @hfield: Resulting array.
 * @height: Height above surface.
 * @size_a: Left direction oriented area width.
 * @size_b: Right direction orientated area width.
 * @size_c: Gap size.
 * @magnetisation: Remanent magnetisation.
 * @thickness: Film thickness.
 * @component: Component to output.
 *
 * Calculates magnetic field or its derivatives above a simple medium
 * consisting of stripes of left and right direction magnetisation.
 * Results are added to the @hfield array, so it should be cleared if
 * function is run only once.
 *
 * Since: 2.51
 **/
void
gwy_data_field_mfm_parallel_medium(GwyDataField *hfield,
                                   gdouble height,
                                   gdouble size_a,
                                   gdouble size_b,
                                   gdouble size_c,
                                   gdouble magnetisation,
                                   gdouble thickness,
                                   GwyMFMComponentType component)
{
    gint xres, yres, i, j, k, n, listsize;
    gdouble d, xreal, x, pos, val, u, v, m;
    gdouble *data, *row;
    gint *xlist, *dirlist;
    const gchar *unitstr;

    g_return_if_fail(GWY_IS_DATA_FIELD(hfield));

    xres = hfield->xres;
    yres = hfield->xres;
    xreal = hfield->xreal;

    listsize = 10*xres;
    xlist = g_new(gint, listsize);
    dirlist = g_new(gint, listsize);

    m = MU_0*magnetisation/G_PI;
    d = 20.0*(size_a + size_b + thickness + height);
    pos = -d;

    n = 0;
    xlist[n] = pos*xres/xreal;
    dirlist[n] = 1;
    n++;

    do {
       pos += size_a + size_c/2;
       xlist[n] = pos*xres/xreal;
       dirlist[n] = -1;
       n++;

       pos += size_b + size_c/2;
       xlist[n] = pos*xres/xreal;
       dirlist[n] = 1;
       n++;
    } while (pos < xreal + d && n < listsize);

    row = g_new0(gdouble, xres);
    for (k = 0; k < n; k++) {
        for (j = 0; j < xres; j++) {
            x = ((gdouble)j - xlist[k])*xreal/xres;
            if (component == GWY_MFM_COMPONENT_HX) {
                u = x*x + size_c*size_c + size_c*(thickness+height);
                v = x*x + size_c*size_c + size_c*height;
                val = -m*dirlist[k]*(atan(x*(thickness+height)/u)
                                     - atan(x*height/v));
            }
            else if (component == GWY_MFM_COMPONENT_HY)
                val = 0;
            else if (component == GWY_MFM_COMPONENT_HZ) {
                u = size_c + height + thickness;
                v = size_c + height;
                val = m*dirlist[k]/2*log((x*x + u*u)/(x*x + v*v));
            }
            else if (component == GWY_MFM_COMPONENT_DHZ_DZ) {
                u = size_c + height + thickness;
                v = size_c + height;
                val = m*dirlist[k]*(u/(x*x + u*u) - v/(x*x + v*v));
            }
            else if (component == GWY_MFM_COMPONENT_D2HZ_DZ2) {
                u = size_c + height + thickness;
                v = size_c + height;
                val = m*dirlist[k]*((x*x - u*u)/(x*x + u*u)
                                    - (x*x - v*v)/(x*x + v*v));
            }
            else {
                g_return_if_reached();
            }
            row[j] += val;
        }
    }
    g_free(xlist);
    g_free(dirlist);

    data = hfield->data;
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++, data++)
            *data += row[j];
    }
    g_free(row);
    gwy_data_field_invalidate(hfield);

    if (component == GWY_MFM_COMPONENT_DHZ_DZ)
        unitstr = "A/m^2";
    else if (component == GWY_MFM_COMPONENT_D2HZ_DZ2)
        unitstr = "A/m^3";
    else
        unitstr = "A/m";

    gwy_si_unit_set_from_string(gwy_data_field_get_si_unit_z(hfield), unitstr);
}

/**
 * gwy_data_field_mfm_current_line:
 * @hfield: Resulting array.
 * @height: Height above surface.
 * @width: Current line width.
 * @position: Current line x position in the resulting array.
 * @current: Curent passing through the line.
 * @component: Component to output.
 *
 * Calculates magnetic field or its derivatives above a flat current line (stripe).
 * Results are added to the @hfield array, so it should be cleared if
 * function is run only once.
 *
 * Since: 2.51
 **/
void
gwy_data_field_mfm_current_line(GwyDataField *hfield,
                                gdouble height,
                                gdouble width,
                                gdouble position,
                                gdouble current,
                                GwyMFMComponentType component)
{
    gint xres, yres, i, j;
    gdouble x, xreal, val, m, w2, hh, xpw2h2, xmw2h2, t;
    gdouble *data, *row;
    const gchar *unitstr;

    g_return_if_fail(GWY_IS_DATA_FIELD(hfield));
    xres = hfield->xres;
    yres = hfield->xres;
    xreal = hfield->xreal;

    m = current/(2*G_PI*width);
    w2 = 0.5*width;
    hh = height*height;

    row = g_new0(gdouble, xres);
    for (j = 0; j < xres; j++) {
        x = j*xreal/xres - position;
        if (component == GWY_MFM_COMPONENT_HX)
            val = m*atan(height*width/(hh + x*x - w2*w2));
        else if (component == GWY_MFM_COMPONENT_HY)
            val = 0;
        else if (component == GWY_MFM_COMPONENT_HZ) {
            xmw2h2 = (x - w2)*(x - w2) + hh;
            xpw2h2 = (x + w2)*(x + w2) + hh;
            val = 0.5*m*log(xmw2h2/xpw2h2);
        }
        else if (component == GWY_MFM_COMPONENT_DHZ_DZ) {
            xmw2h2 = (x - w2)*(x - w2) + hh;
            xpw2h2 = (x + w2)*(x + w2) + hh;
            t = 1.0/(xmw2h2*xpw2h2);
            val = m*x*height*width*t;
        }
        else if (component == GWY_MFM_COMPONENT_D2HZ_DZ2) {
            xmw2h2 = (x - w2)*(x - w2) + hh;
            xpw2h2 = (x + w2)*(x + w2) + hh;
            t = 1.0/(xmw2h2*xpw2h2);
            val = m*x*width*t*(1.0 - 2.0*hh*(xmw2h2 + xpw2h2)*t);
        }
        else {
            g_return_if_reached();
        }
        row[j] += val;
    }

    data = hfield->data;
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++, data++)
            *data += row[j];
    }
    g_free(row);
    gwy_data_field_invalidate(hfield);

    if (component == GWY_MFM_COMPONENT_DHZ_DZ)
        unitstr = "A/m^2";
    else if (component == GWY_MFM_COMPONENT_D2HZ_DZ2)
        unitstr = "A/m^3";
    else
        unitstr = "A/m";

    gwy_si_unit_set_from_string(gwy_data_field_get_si_unit_z(hfield), unitstr);
}

/************************** Documentation ****************************/

/**
 * GwyMFMProbeType:
 * @GWY_MFM_PROBE_CHARGE: Magnetic point charge probe.
 * @GWY_MFM_PROBE_BAR: Finite rectangular bar.
 *
 * Type of probe for calculation of force in magnetic field microscopy.
 *
 * Since: 2.51
 **/

/**
 * GwyMFMComponentType:
 * @GWY_MFM_COMPONENT_HX: X-component of magnetic field H.
 * @GWY_MFM_COMPONENT_HY: Y-component of magnetic field H.
 * @GWY_MFM_COMPONENT_HZ: Z-component of magnetic field H.
 * @GWY_MFM_COMPONENT_DHZ_DZ: Z-derivative of Z-component of magnetic field H.
 * @GWY_MFM_COMPONENT_D2HZ_DZ2: Second Z-derivative of Z-component of magnetic
 *                              field H.
 *
 * Type of field component calculated in magnetic field microscopy.
 *
 * Since: 2.51
 **/

/**
 * SECTION:mfm
 * @title: MFM
 * @short_description: Magnetic force microscopy
 **/

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
