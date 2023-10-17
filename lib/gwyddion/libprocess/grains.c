/*
 *  $Id: grains.c 24720 2022-03-22 13:04:22Z yeti-dn $
 *  Copyright (C) 2003-2022 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
 *
 *  Copyright (C) 2013 Brazilian Nanotechnology National Laboratory
 *  E-mail: Vinicius Barboza <vinicius.barboza@lnnano.cnpem.br>
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
#include <libprocess/linestats.h>
#include <libprocess/filters.h>
#include <libprocess/arithmetic.h>
#include <libprocess/stats.h>
#include <libprocess/grains.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"

enum {
    FOREGROUND_FLAG = 1,
    BACKGROUND_FLAG = 0,
};

static gdouble  class_weight                 (GwyDataLine *hist,
                                              gint t,
                                              gint flag);
static gdouble  class_mean                   (GwyDataLine *hist,
                                              gdouble min,
                                              gdouble max,
                                              gint t,
                                              gint flag);
static gint*    gwy_data_field_fill_grain    (GwyDataField *data_field,
                                              gint col,
                                              gint row,
                                              gint *nindices);
static gint     gwy_data_field_fill_one_grain(gint xres,
                                              gint yres,
                                              const gint *data,
                                              gint col,
                                              gint row,
                                              gint *visited,
                                              gint grain_no,
                                              IntList *listv,
                                              IntList *listh);

/**
 * gwy_data_field_grains_mark_height:
 * @data_field: Data to be used for marking.
 * @grain_field: Data field to store the resulting mask to.
 * @threshval: Relative height threshold, in percents.
 * @below: If %TRUE, data below threshold are marked, otherwise data above threshold are marked.
 *
 * Marks data that are above/below height threshold.
 **/
void
gwy_data_field_grains_mark_height(GwyDataField *data_field,
                                  GwyDataField *grain_field,
                                  gdouble threshval,
                                  gboolean below)
{
    gdouble min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    gwy_data_field_copy(data_field, grain_field, FALSE);
    gwy_data_field_get_min_max(grain_field, &min, &max);
    if (below)
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 1, 0);
    else
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 0, 1);

    gwy_data_field_invalidate(grain_field);
}

/**
 * gwy_data_field_grains_mark_curvature:
 * @data_field: Data to be used for marking.
 * @grain_field: Data field to store the resulting mask to.
 * @threshval: Relative curvature threshold, in percents.
 * @below: If %TRUE, data below threshold are marked, otherwise data above threshold are marked.
 *
 * Marks data that are above/below curvature threshold.
 **/
void
gwy_data_field_grains_mark_curvature(GwyDataField *data_field,
                                     GwyDataField *grain_field,
                                     gdouble threshval,
                                     gboolean below)
{
    gdouble min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    gwy_data_field_copy(data_field, grain_field, FALSE);
    gwy_data_field_filter_laplacian(grain_field);

    gwy_data_field_get_min_max(grain_field, &min, &max);
    if (below)
        gwy_data_field_threshold(grain_field, min + threshval*(max - min)/100.0, 1, 0);
    else
        gwy_data_field_threshold(grain_field, min + threshval*(max - min)/100.0, 0, 1);

    gwy_data_field_invalidate(grain_field);
}

/**
 * gwy_data_field_grains_mark_slope:
 * @data_field: Data to be used for marking.
 * @grain_field: Data field to store the resulting mask to.
 * @threshval: Relative slope threshold, in percents.
 * @below: If %TRUE, data below threshold are marked, otherwise data above threshold are marked.
 *
 * Marks data that are above/below slope threshold.
 **/
void
gwy_data_field_grains_mark_slope(GwyDataField *data_field,
                                 GwyDataField *grain_field,
                                 gdouble threshval,
                                 gboolean below)
{
    GwyDataField *masky;
    const gdouble *mdata;
    gdouble *gdata;
    gint i, n;
    gdouble min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    n = data_field->xres * data_field->yres;

    masky = gwy_data_field_duplicate(data_field);
    gwy_data_field_copy(data_field, grain_field, FALSE);
    gwy_data_field_filter_sobel(grain_field, GWY_ORIENTATION_HORIZONTAL);
    gwy_data_field_filter_sobel(masky, GWY_ORIENTATION_VERTICAL);

    gdata = grain_field->data;
    mdata = masky->data;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i) \
            shared(gdata,mdata,n)
#endif
    for (i = 0; i < n; i++)
        gdata[i] = sqrt(gdata[i]*gdata[i] + mdata[i]*mdata[i]);

    gwy_data_field_get_min_max(grain_field, &min, &max);
    if (below)
        gwy_data_field_threshold(grain_field, min + threshval*(max - min)/100.0, 1, 0);
    else
        gwy_data_field_threshold(grain_field, min + threshval*(max - min)/100.0, 0, 1);

    g_object_unref(masky);
    gwy_data_field_invalidate(grain_field);
}

/**
 * gwy_data_field_otsu_threshold:
 * @data_field: A data field.
 *
 * Finds Otsu's height threshold for a data field.
 *
 * The Otsu's threshold is optimal in the sense that it minimises the inter-class variances of two classes of pixels:
 * above and below theshold.
 *
 * Since: 2.37
 **/
gdouble
gwy_data_field_otsu_threshold(GwyDataField *data_field)
{
    guint nstats, nn, t;
    GwyDataLine *hist;
    gdouble min, max, thresh, max_var;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);

    /* Getting histogram length and max/min values for the data field */
    nn = data_field->xres * data_field->yres;
    gwy_data_field_get_min_max(data_field, &min, &max);
    if (min == max)
        return min;

    nstats = floor(13.49*cbrt(nn) + 0.5);
    hist = gwy_data_line_new(nstats, 1, FALSE);
    gwy_data_field_dh(data_field, hist, nstats);

    /* Calculating the threshold */
    thresh = max_var = 0.0;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(t) \
            shared(hist,nstats,min,max,thresh,max_var)
#endif
    for (t = 0; t < nstats; t++) {
        /* Background statistics*/
        gdouble weight_0, mean_0;
        /* Foreground statistics*/
        gdouble  weight_1, mean_1;
        gdouble bin, var;

        /* Getting histogram statistics for fg and bg classes */
        weight_0 = class_weight(hist, t, BACKGROUND_FLAG);
        weight_1 = class_weight(hist, t, FOREGROUND_FLAG);
        mean_0 = class_mean(hist, min, max, t, BACKGROUND_FLAG);
        mean_1 = class_mean(hist, min, max, t, FOREGROUND_FLAG);

        /* Interclass variance */
        var = weight_0 * weight_1 * (mean_0 - mean_1) * (mean_0 - mean_1);
        bin = min + (t + 0.5)*(max - min)/nstats;

#ifdef _OPENMP
#pragma omp critical
#endif
        /* Check for greater interclass variance */
        if (var > max_var) {
            max_var = var;
            thresh = bin;
        }
    }

    g_object_unref(hist);

    return thresh;
}

/*
 * class_weight:
 * @hist: A #GwyDataLine histogram from where to calculate the class probability for @t.
 * @t: A threshold value.
 * @flag: Alternates between %BACKGROUND or %FOREGROUND class probabilities.
 *
 * Returns: The probability for a class given a threshold value.
 **/
static gdouble
class_weight(GwyDataLine *hist,
             gint t,
             gint flag)
{
    gint i;
    gint len;
    gdouble roi, total, weight;
    gdouble *data;

    len = hist->res;
    data = hist->data;
    roi = 0;
    total = 0;

    for (i = 0; i < len; i++)
        total += data[i];

    for (i = flag*t; i < (1 - flag)*t + flag*len; i++)
        roi += data[i];

    weight = roi/total;

    return weight;
}

/*
 * class_mean:
 * @hist: A #GwyDataLine histogram from where to calculate the class mean for @t.
 * @t: A threshold value.
 * @flag: Alternates between %BACKGROUND or %FOREGROUND class mean.
 *
 * Returns: The mean value for a class given a threshold value.
 **/
static gdouble
class_mean(GwyDataLine *hist,
           gdouble min, gdouble max,
           gint t,
           gint flag)
{
    gint len, i;
    gdouble bin, val;
    gdouble roi, total, mean;
    gdouble *data;

    len = hist->res;
    data = hist->data;
    roi = 0;
    total = 0;

    if (t == 0 && flag == 0) {
        return 0.0;
    }

    for (i = flag*t; i < (1 - flag)*t + flag*len; i++) {
        val = data[i];
        bin = min + ( (i + 0.5) * (max-min) / len );
        roi += bin * val;
        total += val;
    }

    mean = roi/total;

    return mean;
}

/**
 * gwy_data_field_grains_remove_grain:
 * @grain_field: Field of marked grains (mask).
 * @col: Column inside a grain.
 * @row: Row inside a grain.
 *
 * Removes one grain at given position.
 *
 * Returns: %TRUE if a grain was actually removed, i.e. (@col,@row) was inside a grain.
 **/
gboolean
gwy_data_field_grains_remove_grain(GwyDataField *grain_field,
                                   gint col,
                                   gint row)
{
    gint *points;
    gint npoints = 0;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(grain_field), FALSE);
    g_return_val_if_fail(col >= 0 && col < grain_field->xres, FALSE);
    g_return_val_if_fail(row >= 0 && row < grain_field->yres, FALSE);

    if (!grain_field->data[grain_field->xres*row + col])
        return FALSE;

    points = gwy_data_field_fill_grain(grain_field, col, row, &npoints);
    while (npoints) {
        npoints--;
        grain_field->data[points[npoints]] = 0.0;
    }
    g_free(points);
    gwy_data_field_invalidate(grain_field);

    return TRUE;
}

/**
 * gwy_data_field_grains_extract_grain:
 * @grain_field: Field of marked grains (mask).
 * @col: Column inside a grain.
 * @row: Row inside a grain.
 *
 * Removes all grains except that one at given position.
 *
 * If there is no grain at (@col, @row), all grains are removed.
 *
 * Returns: %TRUE if a grain remained (i.e., (@col,@row) was inside a grain).
 **/
gboolean
gwy_data_field_grains_extract_grain(GwyDataField *grain_field,
                                    gint col,
                                    gint row)
{
    gint *points;
    gint npoints = 0;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(grain_field), FALSE);
    g_return_val_if_fail(col >= 0 && col < grain_field->xres, FALSE);
    g_return_val_if_fail(row >= 0 && row < grain_field->yres, FALSE);

    if (!grain_field->data[grain_field->xres*row + col]) {
        gwy_data_field_clear(grain_field);
        return FALSE;
    }

    points = gwy_data_field_fill_grain(grain_field, col, row, &npoints);
    gwy_data_field_clear(grain_field);
    while (npoints) {
        npoints--;
        grain_field->data[points[npoints]] = 1.0;
    }
    g_free(points);
    gwy_data_field_invalidate(grain_field);

    return TRUE;
}

/**
 * gwy_data_field_grains_remove_by_number:
 * @grain_field: Field of marked grains (mask).
 * @number: Grain number was filled by gwy_data_field_number_grains().
 *
 * Removes grain identified by @number.
 *
 * Since: 2.35
 **/
void
gwy_data_field_grains_remove_by_number(GwyDataField *grain_field,
                                       gint number)
{
    gint i, xres, yres;
    gdouble *data;
    gint *grains;

    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    xres = grain_field->xres;
    yres = grain_field->yres;
    data = grain_field->data;

    grains = g_new0(gint, xres*yres);
    gwy_data_field_number_grains(grain_field, grains);

    for (i = 0; i < xres*yres; i++) {
        if (grains[i] == number) {
            data[i] = 0;
        }
    }

    g_free(grains);
}

/**
 * gwy_data_field_grains_remove_by_size:
 * @grain_field: Field of marked grains (mask).
 * @size: Grain area threshold, in square pixels.
 *
 * Removes all grains below specified area.
 **/
void
gwy_data_field_grains_remove_by_size(GwyDataField *grain_field,
                                     gint size)
{
    gint i, n, ngrains;
    gdouble *data;
    gint *grain_size;
    gint *grains;

    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    n = grain_field->xres * grain_field->yres;
    data = grain_field->data;

    grains = g_new0(gint, n);
    ngrains = gwy_data_field_number_grains(grain_field, grains);
    grain_size = gwy_data_field_get_grain_sizes(grain_field, ngrains, grains, NULL);
    /* Avoid some no-op work below. */
    grain_size[0] = size;

    /* Too trivial to parallelise. */
    /* remove grains */
    for (i = 0; i < n; i++) {
        if (grain_size[grains[i]] < size)
            data[i] = 0;
    }
    g_free(grains);

    for (i = 1; i <= ngrains; i++) {
        if (grain_size[i] < size) {
            gwy_data_field_invalidate(grain_field);
            break;
        }
    }
    g_free(grain_size);
}

/**
 * gwy_data_field_grains_remove_by_height:
 * @data_field: Data to be used for marking
 * @grain_field: Field of marked grains (mask)
 * @threshval: Relative height threshold, in percents.
 * @below: If %TRUE, grains below threshold are removed, otherwise grains above threshold are removed.
 *
 * Removes grains that are higher/lower than given threshold value.
 **/
void
gwy_data_field_grains_remove_by_height(GwyDataField *data_field,
                                       GwyDataField *grain_field,
                                       gdouble threshval,
                                       gboolean below)
{
    gint i, n, ngrains;
    gdouble *data;
    gboolean *grain_kill;
    gdouble min, max;
    gint *grains;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    n = grain_field->xres * grain_field->yres;
    data = grain_field->data;

    gwy_data_field_get_min_max(data_field, &min, &max);
    threshval = min + threshval*(max - min)/100.0;

    grains = g_new0(gint, n);
    ngrains = gwy_data_field_number_grains(grain_field, grains);

    /* find grains to remove */
    grain_kill = g_new0(gboolean, ngrains + 1);
    if (below) {
        /* Too trivial to parallelise. */
        for (i = 0; i < n; i++) {
            if (grains[i] && data[i] < threshval)
                grain_kill[grains[i]] = TRUE;
        }
    }
    else {
        /* Too trivial to parallelise. */
        for (i = 0; i < n; i++) {
            if (grains[i] && data[i] > threshval)
                grain_kill[grains[i]] = TRUE;
        }
    }
    /* Avoid some no-op work below. */
    grain_kill[0] = FALSE;

    /* remove them */
    /* Too trivial to parallelise. */
    for (i = 0; i < n; i++) {
        if (grain_kill[grains[i]])
            data[i] = 0;
    }
    for (i = 1; i <= ngrains; i++) {
        if (grain_kill[i]) {
            gwy_data_field_invalidate(grain_field);
            break;
        }
    }

    g_free(grains);
    g_free(grain_kill);
}

/**
 * gwy_data_field_grains_remove_touching_border:
 * @grain_field: Field of marked grains (mask).
 *
 * Removes all grains that touch field borders.
 *
 * Since: 2.30
 **/
void
gwy_data_field_grains_remove_touching_border(GwyDataField *grain_field)
{
    gint i, xres, yres, ngrains;
    gdouble *data;
    gint *grains;
    gboolean *touching;

    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    xres = grain_field->xres;
    yres = grain_field->yres;
    data = grain_field->data;

    grains = g_new0(gint, xres*yres);
    ngrains = gwy_data_field_number_grains(grain_field, grains);

    /* Remember grains that touch any border. */
    touching = g_new0(gboolean, ngrains + 1);
    for (i = 0; i < xres; i++)
        touching[grains[i]] = TRUE;
    for (i = 1; i < yres-1; i++) {
        touching[grains[i*xres]] = TRUE;
        touching[grains[i*xres + xres-1]] = TRUE;
    }
    for (i = 0; i < xres; i++)
        touching[grains[(yres-1)*xres + i]] = TRUE;
    /* Avoid some no-op work below. */
    touching[0] = FALSE;

    /* Remove grains. */
    for (i = 0; i < xres*yres; i++) {
        if (touching[grains[i]])
            data[i] = 0;
    }
    for (i = 1; i <= ngrains; i++) {
        if (touching[i]) {
            gwy_data_field_invalidate(grain_field);
            break;
        }
    }

    g_free(grains);
    g_free(touching);
}

/**
 * gwy_data_field_grains_add:
 * @grain_field: Field of marked grains (mask).
 * @add_field: Field of marked grains (mask) to be added.
 *
 * Adds @add_field grains to @grain_field.
 *
 * Note: This function is equivalent to
 * |[
 * gwy_data_field_max_of_fields(grain_field, grain_field, add_field);
 * ]|
 **/
void
gwy_data_field_grains_add(GwyDataField *grain_field, GwyDataField *add_field)
{
    gwy_data_field_max_of_fields(grain_field, grain_field, add_field);
}

/**
 * gwy_data_field_grains_intersect:
 * @grain_field: Field of marked grains (mask).
 * @intersect_field: Field of marked grains (mask).
 *
 * Performs intersection betweet two grain fields, result is stored in @grain_field.
 *
 * Note: This function is equivalent to
 * |[
 * gwy_data_field_min_of_fields(grain_field, grain_field, intersect_field);
 * ]|
 **/
void
gwy_data_field_grains_intersect(GwyDataField *grain_field,
                                GwyDataField *intersect_field)
{
    gwy_data_field_min_of_fields(grain_field, grain_field, intersect_field);
}

/**
 * gwy_data_field_grains_invert:
 * @grain_field: Data field (mask) of marked grains.
 *
 * Inverts a data field representing a mask.
 *
 * All non-positive values are transformed to 1.0.  All positive values are transformed to 0.0.
 *
 * Since: 2.43
 **/
void
gwy_data_field_grains_invert(GwyDataField *grain_field)
{
    guint xres, yres, k;

    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));
    xres = grain_field->xres;
    yres = grain_field->yres;
    for (k = 0; k < xres*yres; k++) {
        if (grain_field->data[k] > 0.0)
            grain_field->data[k] = 0.0;
        else
            grain_field->data[k] = 1.0;
    }
}

/**
 * gwy_data_field_grains_autocrop:
 * @mask_field: Data field representing a mask.
 * @symmetrically: %TRUE to remove borders symmetrically, i.e the same number of pixels from left and right, and also
 *                 top and bottom. %FALSE to remove as many empty rows and columns as possible.
 * @left: Location to store how many column were removed from the left, or %NULL.
 * @right: Location to store how many column were removed from the right, or %NULL.
 * @up: Location to store how many row were removed from the top, or %NULL.
 * @down: Location to store how many row were removed from the bottom, or %NULL.
 *
 * Removes empty border rows and columns from a data field representing a mask.
 *
 * If there are border rows and columns filled completely with non-positive values the size of the data field is
 * reduced, removing these rows.  The parameter @symmetrically controls whether the size reduction is maximum possible
 * or symmetrical.
 *
 * When there is no positive value in the field the field size is reduced to the smallest possible.  This means 1x1
 * for @symmetrical being %FALSE and even original dimensions to 2 for @symmetrical being %TRUE.
 *
 * Returns: %TRUE if the field size was reduced at all.  Detailed information about the reduction can be obtained from
 *          @left, @right, @up and @down.
 *
 * Since: 2.43
 **/
gboolean
gwy_data_field_grains_autocrop(GwyDataField *mask_field,
                               gboolean symmetrically,
                               guint *left,
                               guint *right,
                               guint *up,
                               guint *down)
{
    gint xres, yres, i, j, firstcol, firstrow, lastcol, lastrow;
    const gdouble *d;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(mask_field), FALSE);
    xres = mask_field->xres;
    yres = mask_field->yres;
    firstcol = xres;
    firstrow = yres;
    lastcol = lastrow = -1;
    d = mask_field->data;
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++, d++) {
            if (*d > 0.0) {
                if (G_UNLIKELY(i < firstrow))
                    firstrow = i;
                if (G_UNLIKELY(j < firstcol))
                    firstcol = j;
                if (G_UNLIKELY(i > lastrow))
                    lastrow = i;
                if (G_UNLIKELY(j > lastcol))
                    lastcol = j;
            }
        }
    }
    gwy_debug("first (%d,%d) last (%d,%d)",
              firstcol, firstrow, lastcol, lastrow);
    if (firstcol > lastcol) {
        g_assert(firstrow > lastrow);
        /* Anticipate the reduction to 2 for even-sized dimensions. */
        lastcol = (xres - 1)/2;
        firstcol = xres - lastcol;
        lastrow = (yres - 1)/2;
        firstrow = yres - lastrow;
    }
    if (symmetrically) {
        firstcol = MIN(firstcol, xres-1 - lastcol);
        lastcol = xres-1 - firstcol;
        firstrow = MIN(firstrow, yres-1 - lastrow);
        lastrow = yres-1 - firstrow;
    }
    lastcol++;
    lastrow++;

    if (left)
        *left = firstcol;
    if (right)
        *right = xres - lastcol;
    if (up)
        *up = firstrow;
    if (down)
        *down = yres - lastrow;

    gwy_debug("%dx%d at (%d,%d) of %dx%d",
              lastcol-firstcol, lastrow-firstrow, firstcol, firstrow, xres, yres);
    if (firstcol == 0 && firstrow == 0 && lastcol == xres && lastrow == yres)
        return FALSE;

    gwy_data_field_resize(mask_field, firstcol, firstrow, lastcol, lastrow);
    return TRUE;
}

/* Merge grains i and j in map with full resolution */
static inline void
resolve_grain_map(gint *m, gint i, gint j)
{
    gint ii, jj, k;

    /* Find what i and j fully resolve to */
    for (ii = i; m[ii] != ii; ii = m[ii])
        ;
    for (jj = j; m[jj] != jj; jj = m[jj])
        ;
    k = MIN(ii, jj);

    /* Fix partial resultions to full */
    for (ii = m[i]; m[ii] != ii; ii = m[ii]) {
        m[i] = k;
        i = ii;
    }
    m[ii] = k;
    for (jj = m[j]; m[jj] != jj; jj = m[jj]) {
        m[j] = k;
        j = jj;
    }
    m[jj] = k;
}

static inline gint
finalise_grain_numbering(gint *m, gint max_id,
                         gint *grains, gint n)
{
    gint *mm;
    gint i, id;

    /* Resolve remianing grain number links in map */
    for (i = 1; i <= max_id; i++)
        m[i] = m[m[i]];

    /* Compactify grain numbers */
    mm = g_new0(gint, max_id + 1);
    id = 0;
    for (i = 1; i <= max_id; i++) {
        if (!mm[m[i]]) {
            id++;
            mm[m[i]] = id;
        }
        m[i] = mm[m[i]];
    }
    g_free(mm);

    /* Renumber grains (we make use of the fact m[0] = 0) */
    for (i = 0; i < n; i++)
        grains[i] = m[grains[i]];

    return id;
}

/* Given an image of integers, with zeros corresponding to non-grains and any non-zero values to grains, number
 * grains.  In other words, the non-zero values become grain numbers while zeros are untouched. */
static gint
renumber_grains(gint *grains, gint xres, gint yres, IntList *m)
{
    gint i, j, k, grain_id, max_id, id;
    gboolean must_free = FALSE;

    g_return_val_if_fail(grains, 0);

    if (!m) {
        m = int_list_new(xres + yres);
        must_free = TRUE;
    }
    else
        m->len = 0;

    int_list_add(m, 0);

    /* Number grains with simple unidirectional grain number propagation, updating map m for later full grain join */
    /* OpenMP: This seems too simple and fast to parallelise reasonably. */
    max_id = 0;
    grain_id = 0;
    k = 0;
    for (i = 0; i < yres; i++) {
        grain_id = 0;
        for (j = 0; j < xres; j++, k++) {
            if (grains[k]) {
                /* Grain number is kept from left neighbour unless it does not exist (a new number is assigned) or
                 * a join with top neighbour occurs (m is updated) */
                if (i > 0 && (id = grains[k-xres])) {
                    if (!grain_id)
                        grain_id = id;
                    else if (id != grain_id) {
                        resolve_grain_map(m->data, id, grain_id);
                        grain_id = m->data[id];
                    }
                }
                if (!grain_id) {
                    grain_id = ++max_id;
                    int_list_add(m, grain_id);
                }
                grains[k] = grain_id;
            }
            else
                grain_id = 0;
        }
    }

    max_id = finalise_grain_numbering(m->data, max_id, grains, xres*yres);
    if (must_free)
        int_list_free(m);

    return max_id;
}

/**
 * gwy_data_field_number_grains:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.
 * @grains: Zero-filled array of integers of equal size to @mask_field to put grain numbers to.  Empty space will be
 *          left 0, pixels inside a grain will be set to grain number.  Grains are numbered sequentially 1, 2, 3, ...
 *
 * Numbers grains in a mask data field.
 *
 * Returns: The number of last grain (note they are numbered from 1).
 **/
/* Since 2.53 the caller does not need to pre-fill grains[] with zeros, but do not advertise it much to preserve
 * compatibility... */
gint
gwy_data_field_number_grains(GwyDataField *mask_field,
                             gint *grains)
{
    const gdouble *data;
    gint xres, yres, i;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(mask_field), 0);
    g_return_val_if_fail(grains, 0);

    xres = mask_field->xres;
    yres = mask_field->yres;
    data = mask_field->data;
    /* Too trivial to parallelise. */
    for (i = 0; i < xres*yres; i++)
        grains[i] = (data[i] > 0.0);

    return renumber_grains(grains, xres, yres, NULL);
}

/**
 * gwy_data_field_number_grains_periodic:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.
 * @grains: Zero-filled array of integers of equal size to @mask_field to put grain numbers to.  Empty space will be
 *          left 0, pixels inside a grain will be set to grain number.  Grains are numbered sequentially 1, 2, 3, ...
 *
 * Numbers grains in a periodic mask data field.
 *
 * This function differs from gwy_data_field_number_grains() by the assumption of periodicity, i.e. grains can touch
 * across the opposite field edges.
 *
 * Note that some grain operations assume that grains are contiguous within the image and do not work with periodic
 * grains.
 *
 * You can use gwy_data_field_get_grain_sizes() and there is gwy_data_field_get_grain_bounding_boxes_periodic() for
 * bouding boxes of periodic grains.  As for grain quantities, simple quantities that do not depend on the shape
 * (mean, median, etc.) are evaluated correctly even for periodic grains.  You cannot evaluate quantities depending on
 * grain boundaries though.
 *
 * Returns: The number of last grain (note they are numbered from 1).
 *
 * Since: 2.38
 **/
gint
gwy_data_field_number_grains_periodic(GwyDataField *mask_field,
                                      gint *grains)
{
    gint xres, yres, i, j, ngrains;
    gboolean merged_anything = FALSE;
    gint *m;

    ngrains = gwy_data_field_number_grains(mask_field, grains);
    if (ngrains < 2)
        return ngrains;

    xres = mask_field->xres;
    yres = mask_field->yres;

    m = g_new0(gint, ngrains+1);
    for (j = 0; j <= ngrains; j++)
        m[j] = j;

    for (i = 0; i < xres; i++) {
        gint gno1 = grains[i], gno2 = grains[i + (yres-1)*xres];
        if (gno1 && gno2 && gno1 != gno2) {
            resolve_grain_map(m, gno1, gno2);
            merged_anything = TRUE;
        }
    }

    for (i = 0; i < yres; i++) {
        gint gno1 = grains[i*xres], gno2 = grains[i*xres + xres-1];
        if (gno1 && gno2 && gno1 != gno2) {
            resolve_grain_map(m, gno1, gno2);
            merged_anything = TRUE;
        }
    }

    if (merged_anything)
        ngrains = finalise_grain_numbering(m, ngrains, grains, xres*yres);

    g_free(m);

    return ngrains;
}

static inline void
init_bbox_ranges(gint *bboxes, gint ngrains)
{
    gint i;

    for (i = 1; i <= ngrains; i++) {
        bboxes[4*i] = bboxes[4*i + 1] = G_MAXINT;
        bboxes[4*i + 2] = bboxes[4*i + 3] = -1;
    }
}

static inline void
extend_bbox_range(gint *bboxes, gint id, gint j, gint i)
{
    id *= 4;
    if (j < bboxes[id])
        bboxes[id] = j;
    if (i < bboxes[id + 1])
        bboxes[id + 1] = i;
    if (j > bboxes[id + 2])
        bboxes[id + 2] = j;
    if (i > bboxes[id + 3])
        bboxes[id + 3] = i;
}

static inline void
transform_range_to_bbox(gint *bbox)
{
    bbox[2] = bbox[2] + 1 - bbox[0];
    bbox[3] = bbox[3] + 1 - bbox[1];
}

/**
 * gwy_data_field_get_grain_bounding_boxes:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.  However its contents is
 *              ignored as all grain information is taken from @grains (its dimensions determine the dimensions of
 *              @grains).
 * @ngrains: The number of grains as returned by gwy_data_field_number_grains().
 * @grains: Grain numbers filled with gwy_data_field_number_grains().
 * @bboxes: Array of size at least 4*(@ngrains+1) to fill with grain bounding boxes.  It can be %NULL to allocate
 *          a new array.
 *
 * Find bounding boxes of all grains.
 *
 * As usual the zeroth element of @bboxes does not correspond to any grain; grain numbers start from 1. The bounding
 * boxes are stored as quadruples of indices: (column, row, width, height).
 *
 * Returns: Either @bboxes (if it was not %NULL), or a newly allocated array of size 4(@ngrains + 1).
 *
 * Since: 2.3
 **/
gint*
gwy_data_field_get_grain_bounding_boxes(GwyDataField *mask_field,
                                        gint ngrains,
                                        const gint *grains,
                                        gint *bboxes)
{
    gint xres, yres, i, j, id;
    const gint *grow;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(mask_field), NULL);
    g_return_val_if_fail(grains, NULL);

    xres = mask_field->xres;
    yres = mask_field->yres;
    if (!bboxes)
        bboxes = g_new(gint, 4*(ngrains + 1));

    init_bbox_ranges(bboxes, ngrains);
    for (i = 0; i < yres; i++) {
        grow = grains + i*xres;
        for (j = 0; j < xres; j++) {
            if ((id = grow[j]))
                extend_bbox_range(bboxes, id, j, i);
        }
    }
    for (i = 1; i <= ngrains; i++)
        transform_range_to_bbox(bboxes + 4*i);

    return bboxes;
}

/* Find bounding box from projection. */
static void
make_bbox_from_projection(const gboolean *projection, gint len,
                          gint *pos, gint *size)
{
    gint bl, br;

    if (projection[0] && projection[len-1]) {
        /* At most one gap in the middle, possibly entire width filled. */
        for (bl = 0; bl < len && projection[bl]; bl++)
            ;
        /* If all projection pixels are covered no smaller bbox is possible. */
        if (bl == len) {
            *pos = 0;
            *size = len;
        }
        else {
            /* There must be a gap. */
            for (br = len-1; projection[br]; br--)
                ;
            *pos = br+1;
            *size = bl + len-1-br;
        }
    }
    else {
        /* One contiguous segment, not covering the entire line.  So neither cycle needs an explicit length
         * stop-condition. */
        for (bl = 0; !projection[bl]; bl++)
            ;
        for (br = len-1; !projection[br]; br--)
            ;
        *pos = bl;
        *size = br+1 - bl;
    }
}

/**
 * gwy_data_field_get_grain_bounding_boxes_periodic:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.  However its contents is
 *              ignored as all grain information is taken from @grains (its dimensions determine the dimensions of
 *              @grains).
 * @ngrains: The number of grains as returned by gwy_data_field_number_grains_periodic().
 * @grains: Grain numbers filled with gwy_data_field_number_grains_periodic().
 * @bboxes: Array of size at least 4*(@ngrains+1) to fill with grain bounding boxes.  It can be %NULL to allocate
 *          a new array.
 *
 * Find bounding boxes of all grains.
 *
 * As usual zeroth element of @bboxes does not correspond to any grain; grain numbers start from 1. The bounding boxes
 * are stored as quadruples of indices: (column, row, width, height).
 *
 * The row and column always lie inside the the image.  However, width and height may specify an area which sticks
 * outside.  In this case periodicity needs to be taken into account.
 *
 * Returns: Either @bboxes (if it was not %NULL), or a newly allocated array of size 4(@ngrains + 1).
 *
 * Since: 2.51
 **/
gint*
gwy_data_field_get_grain_bounding_boxes_periodic(GwyDataField *mask_field,
                                                 gint ngrains,
                                                 const gint *grains,
                                                 gint *bboxes)
{
    gint xres, yres, i, j, ii, k, id, ntodo;
    gint *sbboxes;
    gboolean *todograins, *proj_storage;
    gboolean **xprojection, **yprojection;
    const gint *grow;

    bboxes = gwy_data_field_get_grain_bounding_boxes(mask_field, ngrains, grains, bboxes);
    g_return_val_if_fail(bboxes, NULL);

    xres = mask_field->xres;
    yres = mask_field->yres;

    /* Grains passing through the boundary will have full image width or height according to the normal bounding box
     * algorithm.  If there are none we can avoid the slow path. */
    todograins = g_new0(gboolean, ngrains+1);
    ntodo = 0;
    for (id = 1; id <= ngrains; id++) {
        if (bboxes[4*id + 2] == xres || bboxes[4*id + 3] == yres) {
            todograins[id] = TRUE;
            ntodo++;
        }
        else {
            gwy_debug("normal bbox[%d] %dx%d at (%d,%d)",
                      id,
                      bboxes[4*id + 2], bboxes[4*id + 3],
                      bboxes[4*id + 0], bboxes[4*id + 1]);
        }
    }
    gwy_debug("after normal bboxing %d possibly cross-border grains remaining",
              ntodo);
    if (!ntodo) {
        g_free(todograins);
        return bboxes;
    }

    /* Now try virtually shifting the image periodically by half (horizontally (0), vertically (1) and both (2)) and
     * find the bounding boxes again. This catches all grains that simply cross the boundary but are not actually too
     * large. */
    sbboxes = g_new(gint, 4*(ngrains + 1));
    for (k = 0; k < 3 && ntodo; k++) {
        init_bbox_ranges(sbboxes, ngrains);

        for (i = 0; i < yres; i++) {
            grow = grains + i*xres;
            ii = (k == 0) ? i : (i >= yres/2 ? i : i + yres);
            if (k == 1) {
                /* Horizontal lines are preserved. */
                for (j = 0; j < xres; j++) {
                    if (todograins[id = grow[j]])
                        extend_bbox_range(sbboxes, id, j, ii);
                }
            }
            else {
                /* Horizontal lines are split and shifted. */
                for (j = 0; j < xres/2; j++) {
                    if (todograins[id = grow[j]])
                        extend_bbox_range(sbboxes, id, j + xres, ii);
                }
                for (j = xres/2; j < xres; j++) {
                    if (todograins[id = grow[j]])
                        extend_bbox_range(sbboxes, id, j, ii);
                }
            }
        }

        for (id = 1; id <= ngrains; id++) {
            if (!todograins[id])
                continue;

            transform_range_to_bbox(sbboxes + 4*id);
            if (sbboxes[4*id + 2] == xres || sbboxes[4*id + 3] == yres)
                continue;

            /* We were able to find shifted finite bbox for grain @i. Furthermore if the bbox is finite now it must
             * begin inside the unshifted half -- otherwise it would lie *entirely* within the shifted half and be
             * found by the normal algorithm. */
            g_assert(sbboxes[4*id + 0] < xres);
            g_assert(sbboxes[4*id + 1] < yres);
            gwy_assign(bboxes + 4*id, sbboxes + 4*id, 4);
            gwy_debug("shifted-%d bbox[%d] %dx%d at (%d,%d)",
                      k, id, bboxes[4*id + 2], bboxes[4*id + 3], bboxes[4*id + 0], bboxes[4*id + 1]);
            todograins[id] = FALSE;
            ntodo--;
        }
    }
    g_free(sbboxes);

    gwy_debug("after shifted bboxing %d possibly cross-border grains remaining", ntodo);
    if (!ntodo) {
        g_free(todograins);
        return bboxes;
    }

    /* For any remaining grains, hopefully relatively few, find horizontal and vertical projections. */
    proj_storage = g_new0(gboolean, ntodo*(xres + yres));
    xprojection = g_new0(gboolean*, 2*(ngrains + 1));
    yprojection = xprojection + ngrains+1;
    j = 0;
    for (id = 1; id <= ngrains; id++) {
        if (!todograins[id])
            continue;

        xprojection[id] = proj_storage + j*(xres + yres);
        yprojection[id] = proj_storage + j*(xres + yres) + xres;
        j++;
    }
    g_free(todograins);

    for (i = 0; i < yres; i++) {
        grow = grains + i*xres;
        for (j = 0; j < xres; j++) {
            if (xprojection[id = grow[j]]) {
                xprojection[id][j] = TRUE;
                yprojection[id][i] = TRUE;
            }
        }
    }

    for (id = 1; id <= ngrains; id++) {
        if (!xprojection[id])
            continue;

        make_bbox_from_projection(xprojection[id], xres, bboxes + 4*id + 0, bboxes + 4*id + 2);
        make_bbox_from_projection(yprojection[id], yres, bboxes + 4*id + 1, bboxes + 4*id + 3);
    }
    g_free(xprojection);
    g_free(proj_storage);

    return bboxes;
}

static void
close_rectangles(gint *rectwidths, gint from, gint to, gint *wmax, gint *hmax)
{
    gint hm = 0, wm = 0, areamax = 0, h;

    /* Ignore the zeroth element. */
    from = MAX(from, 1);
    for (h = from; h <= to; h++) {
        if (rectwidths[h] * h > areamax) {
            areamax = rectwidths[h] * h;
            hm = h;
            wm = rectwidths[h];
        }
        rectwidths[h] = 0;
    }

    *hmax = hm;
    *wmax = wm;
}

static inline void
remember_largest_rectangle(gint *iboxes, gint gno,
                           gint j, gint i, gint wmax, gint hmax)
{
    if (!gno)
        return;

    iboxes += 4*gno;
    if (hmax*wmax > iboxes[2]*iboxes[3]) {
        iboxes[0] = j - wmax;
        iboxes[1] = i;
        iboxes[2] = wmax;
        iboxes[3] = hmax;
    }
}

/**
 * gwy_data_field_get_grain_inscribed_boxes:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.  However its contents is
 *              ignored as all grain information is taken from @grains (its dimensions determine the dimensions of
 *              @grains).
 * @ngrains: The number of grains as returned by gwy_data_field_number_grains().
 * @grains: Grain numbers filled with gwy_data_field_number_grains().
 * @iboxes: Array of size at least 4*(@ngrains+1) to fill with grain bounding boxes.  It can be %NULL to allocate
 *          a new array.
 *
 * Find maximum-area inscribed boxes of all grains.
 *
 * As usual zeroth element of @iboxes does not correspond to any grain; grain numbers start from 1. The bounding boxes
 * are stored as quadruples of indices: (column, row, width, height).
 *
 * Returns: Either @iboxes (if it was not %NULL), or a newly allocated array of size 4(@ngrains + 1).
 *
 * Since: 2.53
 **/
gint*
gwy_data_field_get_grain_inscribed_boxes(GwyDataField *mask_field,
                                         gint ngrains,
                                         const gint *grains,
                                         gint *iboxes)
{
    gint xres, yres, i, j, k, h, hmax, wmax, gno, maxopen;
    gint *colsizes, *csrow, *rectwidths;
    const gint *grow;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(mask_field), NULL);
    g_return_val_if_fail(grains, NULL);

    xres = mask_field->xres;
    yres = mask_field->yres;

    /* Pass 1.  Each value in colsizes[] will be filled with how may pixels we we can extend this one while keeping in
     * the same grain.  One means just the single pixels.
     *
     * Unlike in the standard presentation of the algorithm, we go by row for good data locality. */
    colsizes = g_new0(gint, xres*yres);
    csrow = colsizes + (yres - 1)*xres;
    grow = grains + (yres - 1)*xres;
    for (j = 0; j < xres; j++)
        csrow[j] = !!grow[j];

    for (i = yres-1; i; i--) {
        /* These rows we set, reading from the row below. */
        csrow = colsizes + (i-1)*xres;
        grow = grains + (i-1)*xres;
        for (j = 0; j < xres; j++) {
            if (grow[j]) {
                if (grow[j] == grow[xres + j])
                    csrow[j] = csrow[xres + j] + 1;
                else
                    csrow[j] = 1;
            }
        }
    }

    if (!iboxes)
        iboxes = g_new(gint, 4*(ngrains + 1));
    gwy_clear(iboxes, 4*(ngrains + 1));

    /* Pass 2.  Go from the top and close rectangles, remembering the maximum rectangle for each grain number. */
    rectwidths = g_new0(gint, yres);

    maxopen = 0;
    for (i = 0; i < yres; i++) {
        csrow = colsizes + i*xres;
        grow = grains + i*xres;
        gwy_clear(rectwidths, maxopen);
        maxopen = 0;
        gno = 0;
        for (j = 0; j < xres; j++) {
            h = csrow[j];
            /* We hit a new grain -- the algorithm works even with touching grains -- or empty space. */
            if (grow[j] != gno) {
                close_rectangles(rectwidths, 0, maxopen, &wmax, &hmax);
                remember_largest_rectangle(iboxes, gno, j, i, wmax, hmax);
            }
            else {
                /* Grain continuation.  Close rectangles higher than h (up to maxopen). */
                close_rectangles(rectwidths, h+1, maxopen, &wmax, &hmax);
                remember_largest_rectangle(iboxes, gno, j, i, wmax, hmax);
            }
            gno = grow[j];
            if (gno) {
                /* Open or extend all rectangles at most h high.  Both is realised by adding 1 because closed ones are
                 * zeros. */
                for (k = 1; k <= h; k++)
                    rectwidths[k]++;
                maxopen = h;
            }
            else
                maxopen = 0;
        }

        /* Hitting the end of row is the same as hitting empty space. */
        close_rectangles(rectwidths, 0, maxopen, &wmax, &hmax);
        remember_largest_rectangle(iboxes, gno, j, i, wmax, hmax);
    }

    g_free(rectwidths);
    g_free(colsizes);

    return iboxes;
}

/**
 * gwy_data_field_get_grain_sizes:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.  However its contents is
 *              ignored as all grain information is taken from @grains (its dimensions determine the dimensions of
 *              @grains).
 * @ngrains: The number of grains as returned by gwy_data_field_number_grains().
 * @grains: Grain numbers filled with gwy_data_field_number_grains().
 * @sizes: Array of size at least @ngrains+1 to fill with grain sizes (as usual zero does not correspond to any grain,
 *         grains start from 1). It can be %NULL to allocate a new array.
 *
 * Find sizes of all grains in a mask data field.
 *
 * Size is the number of pixels in the grain.
 *
 * The zeroth element of @sizes is filled with the number of pixels not covered by the mask.
 *
 * Returns: Either @sizes (if it was not %NULL), or a newly allocated array of size @ngrains+1.
 *
 * Since: 2.47
 **/
gint*
gwy_data_field_get_grain_sizes(GwyDataField *mask_field,
                               gint ngrains,
                               const gint *grains,
                               gint *sizes)
{
    gint xres, yres, k;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(mask_field), NULL);
    g_return_val_if_fail(grains, NULL);

    xres = mask_field->xres;
    yres = mask_field->yres;
    if (!sizes)
        sizes = g_new(gint, ngrains + 1);

    gwy_clear(sizes, ngrains + 1);
    for (k = 0; k < xres*yres; k++)
        sizes[grains[k]]++;

    return sizes;
}

/**
 * gwy_data_field_area_grains_tgnd:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to the requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @below: If %TRUE, valleys are marked, otherwise mountains are marked.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, a suitable resolution is
 *          determined automatically.
 *
 * Calculates threshold grain number distribution.
 *
 * This function is a simple gwy_data_field_area_grains_tgnd_range() that calculates the distribution in the full
 * range.
 **/
void
gwy_data_field_area_grains_tgnd(GwyDataField *data_field,
                                GwyDataLine *target_line,
                                gint col, gint row,
                                gint width, gint height,
                                gboolean below,
                                gint nstats)
{
    gdouble min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_get_min_max(data_field, NULL, col, row, width, height, &min, &max);
    gwy_data_field_area_grains_tgnd_range(data_field, target_line, col, row, width, height, min, max, below, nstats);
}

static void
tgnd_discrete_range(const gint *heights, gint width, gint height,
                    const gint *hcdist, const gint *hindex,
                    gdouble *line, gint hfrom, gint hto)
{
    gint *grains, *m, *mm;
    gint n, h, i, j, k;
    gint grain_no, last_grain_no;
    IntList *listv, *listh;
    guint msize;

    if (hto == hfrom)
        return;

    n = width*height;

    /* Mark all grains that correspond to discrete level 0, i.e. before we start counting. */
    grains = g_new(gint, n);
    for (i = 0; i < n; i++)
        grains[i] = (heights[i] <= hfrom);

    /* Reuse listv for grain numbering. */
    listv = int_list_new(64);
    last_grain_no = renumber_grains(grains, width, height, listv);
    line[hfrom] = last_grain_no;

    /* When the number of steps is small, just number grains directly at each level. */
    if (hto - hfrom < 20) {
        for (h = hfrom+1; h < hto; h++) {
            if (hcdist[h+1] == hcdist[h]) {
                line[h] = last_grain_no;
                continue;
            }

            for (i = hcdist[h]; i < hcdist[h+1]; i++)
                grains[hindex[i]] = 1;
            last_grain_no = renumber_grains(grains, width, height, listv);
            line[h] = last_grain_no;
        }

        g_free(grains);
        int_list_free(listv);

        return;
    }

    listv->len = 0;
    listh = int_list_new(64);

    m = g_new(gint, 2);
    mm = m;
    msize = 1;

    /* Main iteration */
    for (h = hfrom+1; h < hto; h++) {
        /* Mark new subgrains corresponding just to height h */
        grain_no = last_grain_no;
        gwy_debug("Height %d, number of old grains: %d", h, grain_no);
        for (i = hcdist[h]; i < hcdist[h+1]; i++) {
            j = hindex[i];
            if (!grains[j]) {
                grain_no++;
                gwy_data_field_fill_one_grain(width, height, heights, j % width, j/width,
                                              grains, grain_no, listv, listh);
            }
        }
        gwy_debug("new subgrains: %d", grain_no-last_grain_no);

        if (grain_no == last_grain_no) {
            gwy_debug("skipping empty height level");
            line[h] = line[h-1];
            continue;
        }

        /* Initialize grains number maps for merge scan */
        if (grain_no+1 > msize) {
            g_free(m);
            msize = 2*(grain_no + 1);
            m = g_new(gint, 2*msize);
            mm = m + msize;
        }
        for (i = 0; i <= grain_no; i++) {
            m[i] = i;
            mm[i] = 0;
        }

        /* Find grains that touch each other for merge.
         *
         * Grains that did not touch before don't touch now.  So we are only interested in neighbours of pixels of new
         * subgrains. */
        for (i = hcdist[h]; i < hcdist[h+1]; i++) {
            j = hindex[i];
            /* Left */
            if (j % width && grains[j-1] && m[grains[j]] != m[grains[j-1]])
                resolve_grain_map(m, grains[j], grains[j-1]);
            /* Right */
            if ((j+1) % width && grains[j+1] && m[grains[j]] != m[grains[j+1]])
                resolve_grain_map(m, grains[j], grains[j+1]);
            /* Up */
            if (j/width && grains[j-width] && m[grains[j]] != m[grains[j-width]])
                resolve_grain_map(m, grains[j], grains[j-width]);
            /* Down */
            if (j/width < height-1 && grains[j+width] && m[grains[j]] != m[grains[j+width]])
                resolve_grain_map(m, grains[j], grains[j+width]);
        }

        /* Resolve remianing grain number links in m */
        for (i = 1; i <= grain_no; i++)
            m[i] = m[m[i]];

        /* Compactify grain numbers */
        k = 0;
        for (i = 1; i <= grain_no; i++) {
            if (!mm[m[i]]) {
                k++;
                mm[m[i]] = k;
            }
            m[i] = mm[m[i]];
        }

        /* Renumber grains (we make use of the fact m[0] = 0).
         *
         * This is the only place where we have to scan complete data field. Since grain numbers usually vary wildly
         * and globally, we probably can't avoid it. */
        for (i = 0; i < n; i++)
            grains[i] = m[grains[i]];

        /* The number of grains for this h */
        line[h] = k;
        last_grain_no = k;
    }

    g_free(m);
    g_free(grains);
    int_list_free(listv);
    int_list_free(listh);
}

/**
 * gwy_data_field_area_grains_tgnd_range:
 * @data_field: A data field.
 * @target_line: A data line to store the distribution to.  It will be resampled to the requested width.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @min: Minimum threshold value.
 * @max: Maximum threshold value.
 * @below: If %TRUE, valleys are marked, otherwise mountains are marked.
 * @nstats: The number of samples to take on the distribution function.  If nonpositive, a suitable resolution is
 *          determined automatically.
 *
 * Calculates threshold grain number distribution in given height range.
 *
 * This is the number of grains for each of @nstats equidistant height threshold levels.  For large @nstats this
 * function is much faster than the equivalent number of gwy_data_field_grains_mark_height() calls.
 **/
void
gwy_data_field_area_grains_tgnd_range(GwyDataField *data_field,
                                      GwyDataLine *target_line,
                                      gint col, gint row,
                                      gint width, gint height,
                                      gdouble min, gdouble max,
                                      gboolean below,
                                      gint nstats)
{
    gint *heights, *hindex, *hcdist, *htmp;
    gdouble *data;
    gdouble q;
    gint i, j, h, n, xres;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    if (nstats < 1) {
        nstats = floor(3.49*cbrt(width*height) + 0.5);
        nstats = MAX(nstats, 2);
    }

    gwy_data_line_resample(target_line, nstats, GWY_INTERPOLATION_NONE);
    target_line->off = min;

    n = width*height;
    if (max == min || n == 0) {
        gwy_data_line_clear(target_line);
        target_line->real = 1.0;
        return;
    }
    target_line->real = max - min;

    xres = data_field->xres;
    data = data_field->data + row*xres + col;

    /* Calculate discrete heights.
     *
     * The buckets are a bit unusual.  We want the distribution symmetric from above- and below- viewpoint.  That
     * means if we invert the input data field, negate @below, and reverse the contents of @hcdist we want to get the
     * same distribution (except points that are excatly equal to some bucket edge value, they get rounded to the
     * other bucker).
     *
     * Since min and max may not be the full range of data, we define one more height than nstats and reserve
     * heights[i] == 0 for all values below min and heights[i] == nstats for all values above max (or the other way
     * round if below). */
    heights = g_new(gint, n);
    q = (nstats - 1.0)/(max - min);
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(i,j) \
            shared(data,heights,xres,width,height,min,max,nstats,q,below)
#endif
    for (i = 0; i < height; i++) {
        const gdouble *drow = data + i*xres;
        gint *hrow = heights + i*width;

        if (below) {
            for (j = 0; j < width; j++) {
                gdouble z = drow[j];

                if (z < min)
                    hrow[j] = 0;
                else if (z > max)
                    hrow[j] = nstats;
                else
                    hrow[j] = (gint)floor((z - min)*q + 1);
            }
        }
        else {
            for (j = 0; j < width; j++) {
                gdouble z = drow[j];

                if (z > max)
                    hrow[j] = 0;
                else if (z < min)
                    hrow[j] = nstats;
                else
                    hrow[j] = (gint)floor((max - z)*q + 1);
            }
        }
    }

    /* Calculate cumulative discrete distribution function hcdist.  In other words, group pixels of the same discrete
     * height.  hcdist then holds indices in hindex where each height starts.
     *
     * The purpose of all this is to scale well with nstats.  For small nstats we use more memory and thus trash CPU
     * more badly, but the big advantage is that we do not scan complete @heights each iteration and only touch the
     * pixels that actually constitute the grains (and their neighbours).
     */
    hcdist = g_new0(gint, nstats+1);
    for (i = 0; i < n; i++)
        hcdist[heights[i]]++;

    for (i = 1; i <= nstats; i++)
        hcdist[i] += hcdist[i-1];
    g_assert(hcdist[nstats] == n);

    for (i = nstats; i; i--)
        hcdist[i] = hcdist[i-1];
    hcdist[0] = 0;

    hindex = g_new(gint, n);
    htmp = g_new0(gint, nstats+1);
    for (i = 0; i < n; i++) {
        h = heights[i];
        hindex[hcdist[h] + htmp[h]] = i;
        htmp[h]++;
    }
    g_free(htmp);

#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            shared(heights,width,height,hcdist,hindex,target_line,nstats)
#endif
    {
        gint hfrom = gwy_omp_chunk_start(nstats);
        gint hto = gwy_omp_chunk_end(nstats);

        tgnd_discrete_range(heights, width, height, hcdist, hindex, target_line->data, hfrom, hto);
    }

    g_free(hindex);
    g_free(hcdist);
    g_free(heights);
}

/**
 * gwy_data_field_fill_grain:
 * @data_field: A data field with zeroes in empty space and nonzeroes in grains.
 * @col: Column inside a grain.
 * @row: Row inside a grain.
 * @nindices: Where the number of points in the grain at (@col, @row) should be stored.
 *
 * Finds all the points belonging to the grain at (@col, @row).
 *
 * Returns: A newly allocated array of indices of grain points in @dfield's data, the size of the list is returned in
 *          @nindices.
 **/
static gint*
gwy_data_field_fill_grain(GwyDataField *data_field,
                          gint col, gint row, gint *nindices)
{
    IntList *listv, *listh;
    gint *data, *visited;
    gint *indices;
    gint xres, yres, n, count;
    gint i, j;
    gint initial;

    xres = data_field->xres;
    yres = data_field->yres;
    initial = row*xres + col;
    g_return_val_if_fail(data_field->data[initial], NULL);

    /* check for a single point */
    if ((!col || data_field->data[initial - 1] <= 0)
        && (!row || data_field->data[initial - xres] <= 0)
        && (col + 1 == xres || data_field->data[initial + 1] <= 0)
        && (row + 1 == yres || data_field->data[initial + xres] <= 0)) {
        indices = g_new(gint, 1);

        indices[0] = initial;
        *nindices = 1;

        return indices;
    }

    n = xres*yres;
    visited = g_new0(gint, n);
    data = g_new(gint, n);
    listv = int_list_new(64);
    listh = int_list_new(64);

    for (i = 0; i < n; i++)
        data[i] = data_field->data[i] > 0;

    count = gwy_data_field_fill_one_grain(xres, yres, data, col, row, visited, 1, listv, listh);

    int_list_free(listh);
    int_list_free(listv);
    g_free(data);

    indices = g_new(gint, count);

    j = 0;
    for (i = 0; i < n; i++) {
        if (visited[i])
            indices[j++] = i;
    }
    g_free(visited);

    *nindices = count;
    return indices;
}

/**
 * gwy_data_field_fill_one_grain:
 * @xres: The number of columns in @data.
 * @yres: The number of rows in @data.
 * @data: Arbitrary integer data.  Grain is formed by values equal to the value at (@col, @row).
 * @col: Column inside a grain.
 * @row: Row inside a grain.
 * @visited: An array @col x @row that contain zeroes in empty space and yet unvisited grains.  Current grain will be
 *           filled with @grain_no.
 * @grain_no: Value to fill current grain with.
 * @listv: A working buffer of size at least @col x @row/2 + 2, its content is owerwritten.
 * @listh: A working buffer of size at least @col x @row/2 + 2, its content is owerwritten.
 *
 * Internal function to fill/number a one grain.
 *
 * The @visited, @listv, and @listh buffers are recyclable between calls so they don't have to be allocated and freed
 * for each grain, speeding up sequential grain processing.  Generally, this function itself does not allocate or free
 * any memory.
 *
 * Returns: The number of pixels in the grain.
 **/
static gint
gwy_data_field_fill_one_grain(gint xres,
                              gint yres,
                              const gint *data,
                              gint col, gint row,
                              gint *visited,
                              gint grain_no,
                              IntList *listv,
                              IntList *listh)
{
    gint n, count;
    gint i, p, j;
    gint initial;
    gint look_for;

    g_return_val_if_fail(grain_no, 0);
    initial = row*xres + col;
    look_for = data[initial];

    /* check for a single point */
    visited[initial] = grain_no;
    count = 1;
    if ((!col || data[initial - 1] != look_for)
        && (!row || data[initial - xres] != look_for)
        && (col + 1 == xres || data[initial + 1] != look_for)
        && (row + 1 == yres || data[initial + xres] != look_for)) {

        return count;
    }

    n = xres*yres;
    listv->len = listh->len = 0;
    int_list_add(listv, initial);
    int_list_add(listv, initial);
    int_list_add(listh, initial);
    int_list_add(listh, initial);

    while (listv->len) {
        /* go through vertical lines and expand them horizontally */
        for (i = 0; i < listv->len; i += 2) {
            for (p = listv->data[i]; p <= listv->data[i + 1]; p += xres) {
                gint start, stop;

                /* scan left */
                start = p - 1;
                stop = (p/xres)*xres;
                for (j = start; j >= stop; j--) {
                    if (visited[j] || data[j] != look_for)
                        break;
                    visited[j] = grain_no;
                    count++;
                }
                if (j < start) {
                    int_list_add(listh, j + 1);
                    int_list_add(listh, start);
                }

                /* scan right */
                start = p + 1;
                stop = (p/xres + 1)*xres;
                for (j = start; j < stop; j++) {
                    if (visited[j] || data[j] != look_for)
                        break;
                    visited[j] = grain_no;
                    count++;
                }
                if (j > start) {
                    int_list_add(listh, start);
                    int_list_add(listh, j - 1);
                }
            }
        }
        listv->len = 0;

        /* go through horizontal lines and expand them vertically */
        for (i = 0; i < listh->len; i += 2) {
            for (p = listh->data[i]; p <= listh->data[i + 1]; p++) {
                gint start, stop;

                /* scan up */
                start = p - xres;
                stop = p % xres;
                for (j = start; j >= stop; j -= xres) {
                    if (visited[j] || data[j] != look_for)
                        break;
                    visited[j] = grain_no;
                    count++;
                }
                if (j < start) {
                    int_list_add(listv, j + xres);
                    int_list_add(listv, start);
                }

                /* scan down */
                start = p + xres;
                stop = p % xres + n;
                for (j = start; j < stop; j += xres) {
                    if (visited[j] || data[j] != look_for)
                        break;
                    visited[j] = grain_no;
                    count++;
                }
                if (j > start) {
                    int_list_add(listv, start);
                    int_list_add(listv, j - xres);
                }
            }
        }
        listh->len = 0;
    }

    return count;
}

/**
 * gwy_data_field_fill_voids:
 * @data_field: A data field with zeroes in empty space and nonzeroes in grains.
 * @nonsimple: Pass %TRUE to fill also voids that are not simple-connected (e.g. ring-like).  This can result in grain
 *             merging if a small grain is contained within a void.  Pass %FALSE to fill only simple-connected grains.
 *
 * Fills voids in grains in a data field representing a mask.
 *
 * Voids in grains are zero pixels in @data_field from which no path exists through other zero pixels to the field
 * boundary.  The paths are considered in 8-connectivity because grains themselves are considered in 4-connectivity.
 *
 * Returns: %TRUE if any voids were filled at all, %FALSE if no change was made.
 *
 * Since: 2.37
 **/
gboolean
gwy_data_field_fill_voids(GwyDataField *data_field,
                          gboolean nonsimple)
{
    GwyDataField *voids;
    gdouble *data;
    gint xres, yres;
    gint i, j, k, gno, gno2;
    gint *grains = NULL, *vgrains = NULL, *grain_neighbours = NULL;
    gboolean *unbound_vgrains;
    IntList **vgrain_neighbours;
    guint nvgrains;
    gboolean changed, changed_ever = FALSE, onlysimple = !nonsimple;

    g_return_val_if_fail(data_field, FALSE);

    xres = gwy_data_field_get_xres(data_field);
    yres = gwy_data_field_get_yres(data_field);
    data = gwy_data_field_get_data(data_field);

    voids = gwy_data_field_duplicate(data_field);
    gwy_data_field_grains_invert(voids);
    vgrains = g_new0(gint, xres*yres);
    nvgrains = gwy_data_field_number_grains(voids, vgrains);
    g_object_unref(voids);

    unbound_vgrains = g_new0(gboolean, nvgrains+1);
    for (i = 0; i < xres; i++) {
        unbound_vgrains[vgrains[i]] = TRUE;
        unbound_vgrains[vgrains[xres*(yres - 1) + i]] = TRUE;
    }
    for (j = 0; j < yres; j++) {
        unbound_vgrains[vgrains[j*xres]] = TRUE;
        unbound_vgrains[vgrains[j*xres + xres-1]] = TRUE;
    }

    if (onlysimple) {
        grains = g_new0(gint, xres*yres);
        grain_neighbours = g_new0(gint, nvgrains+1);
        gwy_data_field_number_grains(data_field, grains);
    }

    vgrain_neighbours = g_new0(IntList*, nvgrains+1);
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++) {
            k = i*xres + j;
            if (!(gno = vgrains[k]))
                continue;

            /* We must take into account grain separators (vgrains) have 8-connectivity while all grain functions work
             * with 4-connectivity.  So construct a map of diagonally touching grains and spread the unboundness
             * through it. */
            if (i && j && (gno2 = vgrains[k-xres-1]) && gno2 != gno)
                int_list_add_unique(vgrain_neighbours + gno, gno2);
            if (i && j < xres-1 && (gno2 = vgrains[k-xres+1]) && gno2 != gno)
                int_list_add_unique(vgrain_neighbours + gno, gno2);
            if (i < yres-1 && j && (gno2 = vgrains[k+xres-1]) && gno2 != gno)
                int_list_add_unique(vgrain_neighbours + gno, gno2);
            if (i < yres-1 && j < xres-1 && (gno2 = vgrains[k+xres+1]) && gno2 != gno)
                int_list_add_unique(vgrain_neighbours + gno, gno2);

            /* To detect non-simple-connected grains, we remember the number of the first normal grain the vgrain
             * touches.  If we find later that the void grain also touches a normal grain with a different number, it
             * means it is not simple-connected and we indicate by setting the neighbour number to G_MAXINT. */
            if (!onlysimple || grain_neighbours[gno] == G_MAXINT)
                continue;

            if (i && (gno2 = grains[k-xres])) {
                if (!grain_neighbours[gno])
                    grain_neighbours[gno] = gno2;
                else if (grain_neighbours[gno] != gno2)
                    grain_neighbours[gno] = G_MAXINT;
            }
            if (j && (gno2 = grains[k-1])) {
                if (!grain_neighbours[gno])
                    grain_neighbours[gno] = gno2;
                else if (grain_neighbours[gno] != gno2)
                    grain_neighbours[gno] = G_MAXINT;
            }
            if (j < xres-1 && (gno2 = grains[k+1])) {
                if (!grain_neighbours[gno])
                    grain_neighbours[gno] = gno2;
                else if (grain_neighbours[gno] != gno2)
                    grain_neighbours[gno] = G_MAXINT;
            }
            if (i < yres-1 && (gno2 = grains[k+xres])) {
                if (!grain_neighbours[gno])
                    grain_neighbours[gno] = gno2;
                else if (grain_neighbours[gno] != gno2)
                    grain_neighbours[gno] = G_MAXINT;
            }
        }
    }

    do {
        changed = FALSE;
        for (gno = 1; gno <= nvgrains; gno++) {
            IntList *list = vgrain_neighbours[gno];

            if (!list || !unbound_vgrains[gno])
                continue;

            changed = TRUE;
            for (k = 0; k < list->len; k++)
                unbound_vgrains[list->data[k]] = TRUE;

            /* We have propagated everything we can from @list, so avoid doing that again and again. */
            int_list_free(list);
            vgrain_neighbours[gno] = NULL;
        }

        changed_ever |= changed;
    } while (changed);

    k = 0;
    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++) {
            if (!data[k] && !unbound_vgrains[vgrains[k]] && (!onlysimple || grain_neighbours[vgrains[k]] != G_MAXINT))
                data[k] = 1.0;
            k++;
        }
    }

    for (gno = 1; gno <= nvgrains; gno++) {
        IntList *list = vgrain_neighbours[gno];
        if (list)
            int_list_free(list);
    }
    g_free(vgrain_neighbours);
    g_free(grain_neighbours);
    g_free(unbound_vgrains);
    g_free(vgrains);
    g_free(grains);

    if (changed_ever)
        gwy_data_field_invalidate(data_field);

    return changed_ever;
}

static inline void
add_boundary_grain(IntList *bpixels, int *bindex, gint gno, gint k)
{
    int_list_add(bpixels, gno);
    int_list_add(bpixels, k);
    bindex[gno]++;
}

/**
 * gwy_data_field_grains_find_boundaries:
 * @mask_field: Data field containing positive values in grains, nonpositive in free space.  However its contents is
 *              ignored as all grain information is taken from @grains (its dimensions determine the dimensions of
 *              @grains).
 * @grains: Grain numbers filled with gwy_data_field_number_grains().
 * @ngrains: The number of grains as returned by gwy_data_field_number_grains().
 * @bindex: Array of size at least @ngrains+2.  It will be filled block starts in the returned array.
 * @from_border: %TRUE to consider image edges to be grain boundaries.
 *
 * Find boundary pixels of all grains in a data field.
 *
 * The returned array contains pixel indices of grain boundaries, concatenated all one after another.  The block for
 * grain with number @i starts at position @bindex[@i] and ends one before position @bindex[@i+1] where the next block
 * starts.  The indices in each block are stored in ascending order.
 *
 * Boundary pixels are considered in the usual 4-connected metric.
 *
 * The boundary of the no-grain area is not computed. Therefore, the first two elements of @bindex will be always
 * zeros.
 *
 * Returns: A newly allocated array of size given by @bindex[@ngrains+1].
 *
 * Since: 2.61
 **/
gint*
gwy_data_field_grains_find_boundaries(GwyDataField *data_field,
                                      const gint *grains,
                                      gint ngrains,
                                      gint *bindex,
                                      gboolean from_border)
{
    gint xres, yres, i, j, k, g;
    const gint *grow, *gprev, *gnext;
    IntList *bpixels;
    gint *retval;

    g_return_val_if_fail(data_field, NULL);
    g_return_val_if_fail(grains, NULL);
    g_return_val_if_fail(bindex, NULL);

    xres = data_field->xres;
    yres = data_field->yres;
    gwy_clear(bindex, ngrains+2);
    bpixels = int_list_new(ngrains+1);

    /* First row. */
    i = 0;
    grow = grains + i*xres;
    gnext = grow + xres;
    j = 0;
    if ((g = grow[j]) && (from_border | (g ^ grow[j+1]) | (g ^ gnext[j])))
        add_boundary_grain(bpixels, bindex, g, i*xres + j);
    for (j = 1; j < xres-1; j++) {
        if ((g = grow[j]) && (from_border | (g ^ grow[j-1]) | (g ^ grow[j+1]) | (g ^ gnext[j])))
            add_boundary_grain(bpixels, bindex, g, i*xres + j);
    }
    j = xres-1;
    if ((g = grow[j]) && (from_border | (g ^ grow[j-1]) | (g ^ gnext[j])))
        add_boundary_grain(bpixels, bindex, g, i*xres + j);

    /* Middle rows - we do not have to care abound borders, except for the first and last pixels. */
    for (i = 1; i < yres-1; i++) {
        grow = grains + i*xres;
        gprev = grow - xres;
        gnext = grow + xres;
        j = 0;
        if ((g = grow[j]) && (from_border | (g ^ grow[j+1]) | (g ^ gprev[j]) | (g ^ gnext[j])))
            add_boundary_grain(bpixels, bindex, g, i*xres + j);
        for (j = 1; j < xres-1; j++) {
            if ((g = grow[j]) && ((g ^ grow[j-1]) | (g ^ grow[j+1]) | (g ^ gprev[j]) | (g ^ gnext[j])))
                add_boundary_grain(bpixels, bindex, g, i*xres + j);
        }
        j = xres-1;
        if ((g = grow[j]) && (from_border | (g ^ grow[j-1]) | (g ^ gprev[j]) | (g ^ gnext[j])))
            add_boundary_grain(bpixels, bindex, g, i*xres + j);
    }

    /* Last row. */
    i = yres-1;
    grow = grains + i*xres;
    gprev = grow - xres;
    j = 0;
    if ((g = grow[j]) && (from_border | (g ^ grow[j+1]) | (g ^ gprev[j])))
        add_boundary_grain(bpixels, bindex, g, i*xres + j);
    for (j = 1; j < xres-1; j++) {
        if ((g = grow[j]) && (from_border | (g ^ grow[j-1]) | (g ^ grow[j+1]) | (g ^ gprev[j])))
            add_boundary_grain(bpixels, bindex, g, i*xres + j);
    }
    j = xres-1;
    if ((g = grow[j]) && (from_border | (g ^ grow[j-1]) | (g ^ gprev[j])))
        add_boundary_grain(bpixels, bindex, g, i*xres + j);

    /* Sum and rewind the index. */
    for (i = 1; i <= ngrains+1; i++)
        bindex[i] += bindex[i-1];
    retval = g_new(gint, bindex[ngrains+1]);
    for (i = ngrains+1; i > 0; i--)
        bindex[i] = bindex[i-1];
    bindex[0] = 0;

    /* Fill the result using remembered bpixels. */
    for (i = 0; i < bpixels->len; i += 2) {
        g = bpixels->data[i];
        k = bpixels->data[i+1];
        retval[bindex[g]++] = k;
    }
    int_list_free(bpixels);
    for (i = ngrains+1; i > 0; i--)
        bindex[i] = bindex[i-1];
    bindex[0] = 0;

    return retval;
}

/************************** Documentation ****************************/

/**
 * SECTION:grains
 * @title: grains
 * @short_description: Grain detection and processing
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
