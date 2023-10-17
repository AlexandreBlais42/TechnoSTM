/*
 *  $Id: grains-wshed.c 22335 2019-07-25 08:23:38Z yeti-dn $
 *  Copyright (C) 2003-2017 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
 *
 *  The quicksort algorithm was copied from GNU C library,
 *  Copyright (C) 1991, 1992, 1996, 1997, 1999 Free Software Foundation, Inc.
 *  See below.
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

#include "config.h"
#include <string.h>
#include <libgwyddion/gwymacros.h>
#include <libprocess/filters.h>
#include <libprocess/stats.h>
#include <libprocess/grains.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"

#define GRAIN_BARRIER G_MAXINT

/* Watershed iterator */
typedef struct {
    GwyComputationState cs;
    GwyDataField *data_field;
    GwyDataField *grain_field;
    gint locate_steps;
    gint locate_thresh;
    gdouble locate_dropsize;
    gint wshed_steps;
    gdouble wshed_dropsize;
    gboolean prefilter;
    gboolean below;
    gint internal_i;
    GwyDataField *min;
    GwyDataField *water;
    GwyDataField *mark_dfield;
} GwyWatershedState;

static gboolean step_by_one                  (GwyDataField *data_field,
                                              gint *rcol,
                                              gint *rrow);
static void     drop_step                    (GwyDataField *data_field,
                                              GwyDataField *water_field,
                                              gdouble dropsize);
static void     drop_minima                  (GwyDataField *water_field,
                                              GwyDataField *min_field,
                                              gint threshval);
static void     process_mask                 (GwyDataField *grain_field,
                                              gint col,
                                              gint row);
static void     wdrop_step                   (GwyDataField *data_field,
                                              GwyDataField *min_field,
                                              GwyDataField *water_field,
                                              GwyDataField *grain_field,
                                              gdouble dropsize);
static void     mark_grain_boundaries        (GwyDataField *grain_field);
static void     waterpour_sort               (const gdouble *d,
                                              gint *idx,
                                              gint n);


/**
 * gwy_data_field_grains_mark_watershed:
 * @data_field: Data to be used for marking.
 * @grain_field: Result of marking (mask).
 * @locate_steps: Locating algorithm steps.
 * @locate_thresh: Locating algorithm threshold.
 * @locate_dropsize: Locating drop size.
 * @wshed_steps: Watershed steps.
 * @wshed_dropsize: Watershed drop size.
 * @prefilter: Use prefiltering.
 * @below: If %TRUE, valleys are marked, otherwise mountains are marked.
 *
 * Performs watershed algorithm.
 **/
void
gwy_data_field_grains_mark_watershed(GwyDataField *data_field,
                                     GwyDataField *grain_field,
                                     gint locate_steps,
                                     gint locate_thresh,
                                     gdouble locate_dropsize,
                                     gint wshed_steps,
                                     gdouble wshed_dropsize,
                                     gboolean prefilter,
                                     gboolean below)
{
    GwyDataField *min, *water, *mark_dfield;
    gint xres, yres, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    xres = data_field->xres;
    yres = data_field->yres;

    min = gwy_data_field_new_alike(data_field, TRUE);
    water = gwy_data_field_new_alike(data_field, TRUE);
    mark_dfield = gwy_data_field_duplicate(data_field);
    if (below)
        gwy_data_field_multiply(mark_dfield, -1.0);
    if (prefilter)
        gwy_data_field_filter_median(mark_dfield, 6);

    gwy_data_field_resample(grain_field, xres, yres, GWY_INTERPOLATION_NONE);
    gwy_data_field_clear(grain_field);

    /* odrop */
    for (i = 0; i < locate_steps; i++)
        drop_step(mark_dfield, water, locate_dropsize);
    drop_minima(water, min, locate_thresh);

    /* owatershed */
    gwy_data_field_copy(data_field, mark_dfield, FALSE);
    if (below)
        gwy_data_field_multiply(mark_dfield, -1.0);
    for (i = 0; i < wshed_steps; i++)
        wdrop_step(mark_dfield, min, water, grain_field, wshed_dropsize);

    mark_grain_boundaries(grain_field);

    g_object_unref(min);
    g_object_unref(water);
    g_object_unref(mark_dfield);
    gwy_data_field_invalidate(grain_field);
}

/**
 * gwy_data_field_grains_watershed_init:
 * @data_field: Data to be used for marking.
 * @grain_field: Result of marking (mask).
 * @locate_steps: Locating algorithm steps.
 * @locate_thresh: Locating algorithm threshold.
 * @locate_dropsize: Locating drop size.
 * @wshed_steps: Watershed steps.
 * @wshed_dropsize: Watershed drop size.
 * @prefilter: Use prefiltering.
 * @below: If %TRUE, valleys are marked, otherwise mountains are marked.
 *
 * Initializes the watershed algorithm.
 *
 * This iterator reports its state as #GwyWatershedStateType.
 *
 * Returns: A new watershed iterator.
 **/
GwyComputationState*
gwy_data_field_grains_watershed_init(GwyDataField *data_field,
                                     GwyDataField *grain_field,
                                     gint locate_steps,
                                     gint locate_thresh,
                                     gdouble locate_dropsize,
                                     gint wshed_steps,
                                     gdouble wshed_dropsize,
                                     gboolean prefilter,
                                     gboolean below)
{
    GwyWatershedState *state;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(grain_field), NULL);

    state = g_new0(GwyWatershedState, 1);

    state->cs.state = GWY_WATERSHED_STATE_INIT;
    state->cs.fraction = 0.0;
    state->data_field = g_object_ref(data_field);
    state->grain_field = g_object_ref(grain_field);
    state->locate_steps = locate_steps;
    state->locate_thresh = locate_thresh;
    state->locate_dropsize = locate_dropsize;
    state->wshed_steps = wshed_steps;
    state->wshed_dropsize = wshed_dropsize;
    state->prefilter = prefilter;
    state->below = below;
    state->internal_i = 0;

    return (GwyComputationState*)state;
}

/**
 * gwy_data_field_grains_watershed_iteration:
 * @state: Watershed iterator.
 *
 * Performs one iteration of the watershed algorithm.
 *
 * Fields @state and progress @fraction of watershed state are updated
 * (fraction is calculated for each phase individually).  Once @state
 * becomes %GWY_WATERSHED_STATE_FINISHED, the calculation is finised.
 *
 * A watershed iterator can be created with
 * gwy_data_field_grains_watershed_init().  When iteration ends, either
 * by finishing or being aborted, gwy_data_field_grains_watershed_finalize()
 * must be called to release allocated resources.
 **/
void
gwy_data_field_grains_watershed_iteration(GwyComputationState *cstate)
{
    GwyWatershedState *state = (GwyWatershedState*)cstate;

    if (state->cs.state == GWY_WATERSHED_STATE_INIT) {
        state->min = gwy_data_field_new_alike(state->data_field, TRUE);
        state->water = gwy_data_field_new_alike(state->data_field, TRUE);
        state->mark_dfield = gwy_data_field_duplicate(state->data_field);
        if (state->below)
            gwy_data_field_multiply(state->mark_dfield, -1.0);
        if (state->prefilter)
            gwy_data_field_filter_median(state->mark_dfield, 6);

        gwy_data_field_resample(state->grain_field,
                                state->data_field->xres,
                                state->data_field->yres,
                                GWY_INTERPOLATION_NONE);
        gwy_data_field_clear(state->grain_field);

        state->cs.state = GWY_WATERSHED_STATE_LOCATE;
        state->internal_i = 0;
        state->cs.fraction = 0.0;
    }
    else if (state->cs.state == GWY_WATERSHED_STATE_LOCATE) {
        if (state->internal_i < state->locate_steps) {
            drop_step(state->mark_dfield, state->water, state->locate_dropsize);
            state->internal_i += 1;
            state->cs.fraction = (gdouble)state->internal_i/state->locate_steps;
        }
        else {
            state->cs.state = GWY_WATERSHED_STATE_MIN;
            state->internal_i = 0;
            state->cs.fraction = 0.0;
        }
    }
    else if (state->cs.state == GWY_WATERSHED_STATE_MIN) {
        drop_minima(state->water, state->min, state->locate_thresh);
        state->cs.state = GWY_WATERSHED_STATE_WATERSHED;
        state->internal_i = 0;
        state->cs.fraction = 0.0;
    }
    else if (state->cs.state == GWY_WATERSHED_STATE_WATERSHED) {
        if (state->internal_i == 0) {
            gwy_data_field_copy(state->data_field, state->mark_dfield, FALSE);
            if (state->below)
                gwy_data_field_multiply(state->mark_dfield, -1.0);
        }
        if (state->internal_i < state->wshed_steps) {
            wdrop_step(state->mark_dfield, state->min, state->water,
                       state->grain_field, state->wshed_dropsize);
            state->internal_i += 1;
            state->cs.fraction = (gdouble)state->internal_i/state->wshed_steps;
        }
        else {
            state->cs.state = GWY_WATERSHED_STATE_MARK;
            state->internal_i = 0;
            state->cs.fraction = 0.0;
        }
    }
    else if (state->cs.state == GWY_WATERSHED_STATE_MARK) {
        mark_grain_boundaries(state->grain_field);
        state->cs.state = GWY_WATERSHED_STATE_FINISHED;
        state->cs.fraction = 1.0;
    }
    else if (state->cs.state == GWY_WATERSHED_STATE_FINISHED)
        return;

    gwy_data_field_invalidate(state->grain_field);
}

/**
 * gwy_data_field_grains_watershed_finalize:
 * @state: Watershed iterator.
 *
 * Destroys a watershed iterator, freeing all resources.
 **/
void
gwy_data_field_grains_watershed_finalize(GwyComputationState *cstate)
{
    GwyWatershedState *state = (GwyWatershedState*)cstate;

    GWY_OBJECT_UNREF(state->min);
    GWY_OBJECT_UNREF(state->water);
    GWY_OBJECT_UNREF(state->mark_dfield);
    GWY_OBJECT_UNREF(state->data_field);
    GWY_OBJECT_UNREF(state->grain_field);
    g_free(state);
}

void
gwy_data_field_grains_splash_water(GwyDataField *data_field,
                                   GwyDataField *water,
                                   gint locate_steps,
                                   gdouble locate_dropsize)
{
    GwyDataField *mark_dfield;
    gint i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    mark_dfield = gwy_data_field_duplicate(data_field);

    /* odrop */
    gwy_data_field_clear(water);
    for (i = 0; i < locate_steps; i++)
        drop_step(mark_dfield, water, locate_dropsize);

    gwy_data_field_invalidate(water);
    g_object_unref(mark_dfield);
}

static gboolean
step_by_one(GwyDataField *data_field, gint *rcol, gint *rrow)
{
    gint xres, yres;
    gdouble a, b, c, d, v;
    const gdouble *data;

    xres = data_field->xres;
    yres = data_field->yres;
    data = data_field->data;

    if (*rcol < (xres - 1))
        a = data[*rcol + 1 + xres*(*rrow)];
    else
        a = -G_MAXDOUBLE;

    if (*rcol > 0)
        b = data[*rcol - 1 + xres*(*rrow)];
    else
        b = -G_MAXDOUBLE;

    if (*rrow < (yres - 1))
        c = data[*rcol + xres*(*rrow + 1)];
    else
        c = -G_MAXDOUBLE;

    if (*rrow > 0)
        d = data[*rcol + xres*(*rrow - 1)];
    else
        d = -G_MAXDOUBLE;

    v = data[*rcol + xres*(*rrow)];

    if (v >= a && v >= b && v >= c && v >= d) {
        return TRUE;
    }
    else if (a >= v && a >= b && a >= c && a >= d) {
        *rcol += 1;
        return FALSE;
    }
    else if (b >= v && b >= a && b >= c && b >= d) {
        *rcol -= 1;
        return FALSE;
    }
    else if (c >= v && c >= b && c >= a && c >= d) {
        *rrow += 1;
        return FALSE;
    }
    else {
        *rrow -= 1;
        return FALSE;
    }

    return FALSE;
}

static void
drop_step(GwyDataField *data_field, GwyDataField *water_field, gdouble dropsize)
{
    gint xres, yres, i;
    gint col, row;
    gboolean retval;

    xres = data_field->xres;
    yres = data_field->yres;

    for (i = 0; i < xres*yres; i++) {
        row = (gint)floor((gdouble)i/(gdouble)xres);
        col = i - xres*row;
        if (col == 0 || row == 0 || col == (xres - 1) || row == (yres - 1))
            continue;

        do {
            retval = step_by_one(data_field, &col, &row);
        } while (!retval);

        water_field->data[col + xres*row] += 1;
        data_field->data[col + xres*row] -= dropsize;

    }
    gwy_data_field_invalidate(water_field);
    gwy_data_field_invalidate(data_field);
}

static void
drop_minima(GwyDataField *water_field, GwyDataField *min_field, gint threshval)
{
    gint xres, yres, i, j, ngrains;
    gint *grain_maxima, *grain_size;
    gdouble *data;
    gint *grains;

    xres = water_field->xres;
    yres = water_field->yres;
    data = water_field->data;

    grains = g_new0(gint, xres*yres);
    ngrains = gwy_data_field_number_grains(water_field, grains);
    grain_size = g_new0(gint, ngrains + 1);
    grain_maxima = g_new(gint, ngrains + 1);
    for (i = 1; i <= ngrains; i++)
        grain_maxima[i] = -1;

    /* sum grain sizes and find maxima */
    for (i = 0; i < xres*yres; i++) {
        j = grains[i];
        if (!j)
            continue;

        grain_size[j]++;
        if (grain_maxima[j] < 0
            || data[grain_maxima[j]] < data[i])
            grain_maxima[j] = i;
    }
    g_free(grains);

    /* mark maxima */
    for (i = 1; i <= ngrains; i++) {
        if (grain_size[i] <= threshval)
            continue;

        min_field->data[grain_maxima[i]] = i;
    }

    g_free(grain_maxima);
    g_free(grain_size);
}

static void
process_mask(GwyDataField *grain_field, gint col, gint row)
{
    gint xres, yres, ival[4], val, i;
    gboolean stat;
    gdouble *data;

    xres = grain_field->xres;
    yres = grain_field->yres;
    data = grain_field->data;

    if (col == 0 || row == 0 || col == (xres - 1) || row == (yres - 1)) {
        data[col + xres*row] = -1;
        return;
    }

    /* if this is grain or boundary, keep it */
    if (data[col + xres*row] != 0)
        return;

    /* if there is nothing around, do nothing */
    if ((fabs(data[col + 1 + xres*row]) + fabs(data[col - 1 + xres*row])
         + fabs(data[col + xres*(row + 1)]) + fabs(data[col + xres*(row - 1)]))
        == 0)
        return;

    /* now count the grain values around */
    ival[0] = data[col - 1 + xres*row];
    ival[1] = data[col + xres*(row - 1)];
    ival[2] = data[col + 1 + xres*row];
    ival[3] = data[col + xres*(row + 1)];

    val = 0;
    stat = FALSE;
    for (i = 0; i < 4; i++) {
        if (val > 0 && ival[i] > 0 && ival[i] != val) {
            /* if some value already was there and the now one is different */
            stat = TRUE;
            break;
        }
        else {
            /* if there is some value */
            if (ival[i] > 0) {
                val = ival[i];
            }
        }
    }

    /* it will be boundary or grain */
    data[col + xres*row] = stat ? -1 : val;
}

static void
wdrop_step(GwyDataField *data_field, GwyDataField *min_field,
           GwyDataField *water_field, GwyDataField *grain_field,
           gdouble dropsize)
{
    gint xres, yres, vcol, vrow, col, row, grain;
    gboolean retval;

    xres = data_field->xres;
    yres = data_field->yres;

    grain = 0;
    for (col = 0; col < xres; col++) {
        for (row = 0; row < yres; row++) {
            if (min_field->data[col + xres*row] > 0)
                grain_field->data[col + xres*row] = grain++;
        }
    }
    for (col = 1; col < xres - 1; col++) {
        for (row = 1; row < yres - 1; row++) {

            vcol = col;
            vrow = row;
            do {
                retval = step_by_one(data_field, &vcol, &vrow);
            } while (!retval);

            /*now, distinguish what to change at point vi, vj */
            process_mask(grain_field, vcol, vrow);
            water_field->data[vcol + xres*(vrow)] += 1;
            data_field->data[vcol + xres*(vrow)] -= dropsize;

        }
    }
}

static void
mark_grain_boundaries(GwyDataField *grain_field)
{
    gint xres, yres, col, row;
    GwyDataField *buffer;
    gdouble *data;

    xres = grain_field->xres;
    yres = grain_field->yres;
    /* FIXME: it is not necessary to duplicate complete data field to check
     * a few boundary pixels. */
    buffer = gwy_data_field_duplicate(grain_field);
    data = buffer->data;

    for (col = 1; col < xres - 1; col++) {
        for (row = 1; row < yres - 1; row++) {
            if (data[col + xres*row] != data[col + 1 + xres*row]
                || data[col + xres*row] != data[col + xres*(row + 1)])
                grain_field->data[col + xres*row] = 0;
        }
    }
    g_object_unref(buffer);
}

static gint
waterpour_decide(const gint *assigned, gint xres, gint yres, gint km)
{
    gint i = km/xres, j = km % xres;
    gint idu = (i ? assigned[km-xres] : GRAIN_BARRIER),
         idl = (j ? assigned[km-1] : GRAIN_BARRIER),
         idr = (j < xres-1 ? assigned[km+1] : GRAIN_BARRIER),
         idd = (i < yres-1 ? assigned[km+xres] : GRAIN_BARRIER);

    if (idu == GRAIN_BARRIER || idu == 0) {
        if (idl == GRAIN_BARRIER || idl == 0) {
            if (idr == GRAIN_BARRIER || idr == 0) {
                if (idd == GRAIN_BARRIER || idd == 0)
                    return 0;
                return idd;
            }
            if (idd == GRAIN_BARRIER || idd == 0|| idd == idr)
                return idr;
            return GRAIN_BARRIER;
        }
        if (idr == GRAIN_BARRIER || idr == 0 || idr == idl) {
            if (idd == GRAIN_BARRIER || idd == 0|| idd == idl)
                return idl;
        }
        return GRAIN_BARRIER;
    }
    if (idl == GRAIN_BARRIER || idl == 0 || idl == idu) {
        if (idr == GRAIN_BARRIER || idr == 0 || idr == idu) {
            if (idd == GRAIN_BARRIER || idd == 0 || idd == idu)
                return idu;
        }
    }
    return GRAIN_BARRIER;
}

static gint
mark_one_grain(const gdouble *d, gint *assigned,
               gint xres, gint yres,
               gint km, gint gno,
               IntList *inqueue, IntList *outqueue)
{
    gint m, i, j, count = 1;
    gdouble z = d[km];

    inqueue->len = 0;
    int_list_add(inqueue, km);
    assigned[km] = gno;

    while (inqueue->len) {
        outqueue->len = 0;
        for (m = 0; m < inqueue->len; m++) {
            km = inqueue->data[m];
            i = km/xres;
            j = km % xres;

            if (i > 0 && !assigned[km-xres] && d[km-xres] == z)
                int_list_add(outqueue, km-xres);
            if (j > 0 && !assigned[km-1] && d[km-1] == z)
                int_list_add(outqueue, km-1);
            if (j < xres-1 && !assigned[km+1] && d[km+1] == z)
                int_list_add(outqueue, km+1);
            if (i < yres-1 && !assigned[km+xres] && d[km+xres] == z)
                int_list_add(outqueue, km+xres);
        }

        inqueue->len = 0;
        for (m = 0; m < outqueue->len; m++) {
            km = outqueue->data[m];
            if (!assigned[km]) {
                assigned[km] = gno;
                int_list_add(inqueue, km);
                count++;
            }
        }
    }

    return count;
}

static void
fix_grain_numbers(gint *grains, gint *buf, gint n)
{
    gint gno = 0, k;

    for (k = 0; k < n; k++) {
        gint gnok = grains[k];
        if (gnok && !buf[gnok]) {
            buf[gnok] = ++gno;
        }
    }
    for (k = 0; k < n; k++)
        grains[k] = buf[grains[k]];
}

/**
 * gwy_data_field_waterpour:
 * @data_field: A data field to segmentate.
 * @result: Data field that will be filled with the resulting mask.  It will be
 *          resized to the dimensions of @data_field and its properties set
 *          accordingly.
 * @grains: Optionally, an array with the same number of items as @data_field.
 *          If non-%NULL, it will be filled with grain numbers in the same
 *          manner as gwy_data_field_number_grains().  Pass %NULL to ignore.
 *
 * Performs the classical Vincent watershed segmentation of a data field.
 *
 * The segmentation always results in the entire field being masked with the
 * exception of thin (8-connectivity) lines separating the segments (grains).
 *
 * Compared to gwy_data_field_grains_mark_watershed(), this algorithm is very
 * fast.  However, when used alone, it typically results in a serious
 * oversegmentation as each local minimum gives raise to a grain.  Furthermore,
 * the full segmentation means that also pixels which would be considered
 * outside any grain in the topographical sense will be assigned to some
 * catchment basin.  Therefore, pre- or postprocessing is usually necessary,
 * using the gradient image or a more sophisticated method.
 *
 * The function does not assign pixels with value %HUGE_VAL or larger to any
 * segment.  This can be used to pre-mark certain areas explicitly as
 * boundaries.
 *
 * Since the algorithm numbers the grains as a side effect, you can pass a
 * @grains array and get the grain numbers immediatelly, avoiding the
 * relatively (although not drastically) expensive
 * gwy_data_field_number_grains() call.
 *
 * Returns: The number of segments (grains) in the result, excluding the
 *          separators, i.e. the same convention as in
 *          gwy_data_field_number_grains() is used.
 *
 * Since: 2.37
 **/
gint
gwy_data_field_waterpour(GwyDataField *data_field,
                         GwyDataField *result,
                         gint *grains)
{
    IntList *flatqueue, *fillqueue;
    gint xres, yres, n, k, kq, gno;
    gint *queue, *assigned;
    const gdouble *d;
    gdouble *rd;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(result), 0);

    xres = data_field->xres;
    yres = data_field->yres;
    n = xres*yres;
    gwy_data_field_resample(result, xres, yres, GWY_INTERPOLATION_NONE);

    queue = g_new(gint, n);
    for (k = 0; k < n; k++)
        queue[k] = k;

    d = data_field->data;
    waterpour_sort(d, queue, n);

    assigned = grains ? grains : g_new0(gint, n);
    flatqueue = int_list_new(0);
    fillqueue = int_list_new(0);
    kq = gno = 0;
    while (kq < n) {
        gint len = 1, todo, um;
        gdouble z;

        k = queue[kq];
        z = d[k];
        if (z >= HUGE_VAL)
            break;

        while (kq + len < n && d[queue[kq + len]] == z)
            len++;

        todo = len;
        while (todo) {
            gint m;

            um = GRAIN_BARRIER;
            flatqueue->len = 0;
            for (m = 0; m < len; m++) {
                gint km = queue[kq + m];
                gint id = assigned[km];
                if (id)
                    continue;

                if ((id = waterpour_decide(assigned, xres, yres, km))) {
                    /* Must postpone the grain number assignment. There may
                     * be conflicts later so only queue the position; we have
                     * to waterpour_decide() again. */
                    int_list_add(flatqueue, km);
                }
                else if (um == GRAIN_BARRIER)
                    um = m;
            }

            if (flatqueue->len) {
                /* We have some modifications to commit. */
                for (m = 0; m < flatqueue->len; m++) {
                    gint km = flatqueue->data[m];
                    gint id = waterpour_decide(assigned, xres, yres, km);
                    g_assert(id);
                    assigned[km] = id;
                }
                todo -= flatqueue->len;
            }
            else {
                /* We do not have any modifications.  All unassigned pixels
                 * of this height belong to new grains. */
                break;
            }
        }

        /* Create new grains from remaining pixels. */
        while (todo) {
            gint km = GRAIN_BARRIER;
            while (um < len) {
                k = queue[kq + um];
                um++;
                if (!assigned[k]) {
                    km = k;
                    break;
                }
            }
            g_assert(km != GRAIN_BARRIER);
            todo -= mark_one_grain(d, assigned, xres, yres,
                                   km, ++gno, flatqueue, fillqueue);
        }

        kq += len;
    }

    while (kq < n) {
        k = queue[kq++];
        assigned[k] = GRAIN_BARRIER;
    }

    rd = result->data;
    for (k = 0; k < n; k++) {
        gint gnok = assigned[k];
        gnok = (gnok == GRAIN_BARRIER) ? 0 : gnok;
        assigned[k] = gnok;
        rd[k] = !!gnok;
    }

    /* The grain numbering differs from gwy_data_field_number_grains() which
     * performs the numbering from top left to bottom right.  Since we
     * guarantee stable grain numbers, renumber the grains to match that.
     * Recycle @queue as a scratch buffer.  */
    if (grains) {
        gwy_clear(queue, gno+1);
        fix_grain_numbers(grains, queue, n);
    }

    int_list_free(fillqueue);
    int_list_free(flatqueue);
    g_free(queue);
    if (!grains)
        g_free(assigned);

    result->xreal = data_field->xreal;
    result->yreal = data_field->yreal;
    result->xoff = data_field->xoff;
    result->yoff = data_field->yoff;

    _gwy_copy_si_unit(data_field->si_unit_xy, &result->si_unit_xy);
    if (result->si_unit_z)
        gwy_si_unit_set_from_string(result->si_unit_z, NULL);

    gwy_data_field_invalidate(result);

    return gno;
}

/* Mark sharp maxima with 2, known non-maxima with 1. */
static guint
mark_maxima(GwyDataField *field,
            guint *types)
{
    guint xres = field->xres, yres = field->yres;
    const gdouble *d = field->data;
    guint i, j, marked = 0;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:marked) \
            private(i,j) \
            shared(d,types,xres,yres)
#endif
    for (i = 0; i < yres; i++) {
        gint k = i*xres;
        for (j = 0; j < xres; j++, k++) {
            /* Mark non-maxima. */
            if ((i && d[k] < d[k-xres])
                || (j && d[k] < d[k-1])
                || (j < xres-1 && d[k] < d[k+1])
                || (i < yres-1 && d[k] < d[k+xres])) {
                types[k] = 1;
                marked++;
            }
            /* Mark maxima. */
            else if ((!i || d[k] > d[k-xres])
                     && (!j || d[k] > d[k-1])
                     && (j == xres-1 || d[k] > d[k+1])
                     && (i == yres-1 || d[k] > d[k+xres])) {
                types[k] = 2;
                marked++;
            }
        }
    }

    return xres*yres - marked;
}

/* Mark sharp minima with 2, known non-minima with 1. */
static guint
mark_minima(GwyDataField *field,
            guint *types)
{
    guint xres = field->xres, yres = field->yres;
    const gdouble *d = field->data;
    guint i, j, marked = 0;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:marked) \
            private(i,j) \
            shared(d,types,xres,yres)
#endif
    for (i = 0; i < yres; i++) {
        gint k = i*xres;
        for (j = 0; j < xres; j++, k++) {
            /* Mark non-minima. */
            if ((i && d[k] > d[k-xres])
                || (j && d[k] > d[k-1])
                || (j < xres-1 && d[k] > d[k+1])
                || (i < yres-1 && d[k] > d[k+xres])) {
                types[k] = 1;
                marked++;
            }
            /* Mark minima. */
            else if ((!i || d[k] < d[k-xres])
                     && (!j || d[k] < d[k-1])
                     && (j == xres-1 || d[k] < d[k+1])
                     && (i == yres-1 || d[k] < d[k+xres])) {
                types[k] = 2;
                marked++;
            }
        }
    }

    return xres*yres - marked;
}

/* Propagate non-maxima type to all pixels of the same value.  Or minima. This
 * alogorithm no longer depends on how the states was marked, it just
 * propagates the 1 state though identical values. */
static void
propagate_non_extrema_marking(guint *types, const gdouble *d,
                              guint xres, guint yres)
{
    IntList *inqueue = int_list_new(16);
    IntList *outqueue = int_list_new(16);
    guint i, j, m, k = 0;

    for (i = 0; i < yres; i++) {
        for (j = 0; j < xres; j++, k++) {
            if (types[k])
                continue;
            /* If the value is equal to some neighbour which is a known
             * non-maximum then this pixel is also non-maximum.  (If the
             * neighbour is a known maximum this pixel cannot be unknown.) */
            if ((i && types[k-xres] && d[k] == d[k-xres])
                || (j && types[k-1] && d[k] == d[k-1])
                || (j < xres-1 && types[k+1] && d[k] == d[k+1])
                || (i < yres-1 && types[k+xres] && d[k] == d[k+xres])) {
                types[k] = 1;
                int_list_add(outqueue, k);
            }
        }
    }
    GWY_SWAP(IntList*, inqueue, outqueue);

    while (inqueue->len) {
        for (m = 0; m < inqueue->len; m++) {
            k = inqueue->data[m];
            i = k/xres;
            j = k % xres;

            /* Propagate the non-maximum type to all still unknown
             * neighbours.  Since we set them to known immediately, double
             * queuing is avoided. */
            if (i && !types[k-xres]) {
                types[k-xres] = 1;
                int_list_add(outqueue, k-xres);
            }
            if (j && !types[k-1]) {
                types[k-1] = 1;
                int_list_add(outqueue, k-1);
            }
            if (j < xres-1 && !types[k+1]) {
                types[k+1] = 1;
                int_list_add(outqueue, k+1);
            }
            if (i < yres-1 && !types[k+xres]) {
                types[k+xres] = 1;
                int_list_add(outqueue, k+xres);
            }
        }

        inqueue->len = 0;
        GWY_SWAP(IntList*, inqueue, outqueue);
    }

    int_list_free(inqueue);
    int_list_free(outqueue);
}

/**
 * gwy_data_field_mark_extrema:
 * @dfield: A two-dimensional data field.
 * @extrema: Target field for the extrema mask.
 * @maxima: %TRUE to mark maxima, %FALSE to mark minima.
 *
 * Marks local maxima or minima in a two-dimensional field.
 *
 * Local (or regional) maximum is a contiguous set of pixels that have the same
 * value and this value is sharply greater than the value of any pixel touching
 * the set.  A minimum is defined analogously.  A field filled with a single
 * value is considered to have neither minimum nor maximum.
 *
 * Since: 2.37
 **/
void
gwy_data_field_mark_extrema(GwyDataField *dfield,
                            GwyDataField *extrema,
                            gboolean maxima)
{
    gdouble min, max;
    guint xres, yres, unmarked, k;
    guint *types;

    g_return_if_fail(GWY_IS_DATA_FIELD(dfield));
    g_return_if_fail(GWY_IS_DATA_FIELD(extrema));
    xres = dfield->xres;
    yres = dfield->yres;
    g_return_if_fail(extrema->xres == xres);
    g_return_if_fail(extrema->yres == yres);

    gwy_data_field_clear(extrema);

    gwy_data_field_get_min_max(dfield, &min, &max);
    /* This takes care of 1×1 fields too. */
    if (min == max)
        return;

    types = g_new0(guint, xres*yres);
    unmarked = (maxima ? mark_maxima : mark_minima)(dfield, types);

    if (unmarked)
        propagate_non_extrema_marking(types, dfield->data, xres, yres);

    /* Mark 1 as 0 (non-extremum); mark 0 and 2 as 1 (extremum).  The remaining
       0s are exactly those flat areas which cannot be made non-maximum, i.e.
       they must be maxima.  Assume extrema are relatively sparse so prefer
       fast iteration compared to fast mask bit setting. */
    for (k = 0; k < xres*yres; k++) {
        if (!(types[k] & 1))
            extrema->data[k] = 1.0;
    }

    g_free(types);
    gwy_data_field_invalidate(extrema);
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

typedef struct {
    gdouble z;
    guint k;
} Pair;

#define is_smaller_pair(a, b) \
    ((a)->z < (b)->z || ((a)->z == (b)->z && (a)->k < (b)->k))

#define PSWAP(x, y) GWY_SWAP(Pair, x, y)

#define STACK_SIZE      (CHAR_BIT * sizeof(gsize))
#define PUSH(low, high) ((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define POP(low, high)  ((void) (--top, (low = top->lo), (high = top->hi)))
#define STACK_NOT_EMPTY (stack < top)

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
   next array partition to sort.  To save time, this maximum amount
   of space required to store an array of SIZE_MAX is allocated on the
   stack.  Assuming a 32-bit (64 bit) integer for size_t, this needs
   only 32 * sizeof(stack_node) == 256 bytes (for 64 bit: 1024 bytes).
   Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
   This reduces the probability of selecting a bad pivot value and
   eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
   insertion sort to order the MAX_THRESH items within each partition.
   This is a big win, since insertion sort is faster for small, mostly
   sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the
   stack first, with the algorithm then concentrating on the
   smaller partition.  This *guarantees* no more than log(n)
   stack size is needed (actually O(1) in this case)!  */

static void
sort_pairs(Pair *array,
           gsize n)
{
    /* Note: Specialization makes the insertion sort part relatively more
     * efficient, after some benchmarking this seems be about the best value
     * on Athlon 64. */
    enum { MAX_THRESH = 20 };

    // Stack node declarations used to store unfulfilled partition obligations.
    typedef struct {
        Pair *lo;
        Pair *hi;
    } stack_node;

    if (n < 2)
        /* Avoid lossage with unsigned arithmetic below.  */
        return;

    if (n > MAX_THRESH) {
        Pair *lo = array;
        Pair *hi = lo + (n - 1);
        stack_node stack[STACK_SIZE];
        stack_node *top = stack + 1;

        while (STACK_NOT_EMPTY) {
            Pair *left_ptr;
            Pair *right_ptr;

            /* Select median value from among LO, MID, and HI. Rearrange
               LO and HI so the three values are sorted. This lowers the
               probability of picking a pathological pivot value and
               skips a comparison for both the LEFT_PTR and RIGHT_PTR in
               the while loops. */

            Pair *mid = lo + ((hi - lo) >> 1);

            if (is_smaller_pair(mid, lo))
                PSWAP(*mid, *lo);
            if (is_smaller_pair(hi, mid))
                PSWAP(*mid, *hi);
            else
                goto jump_over;
            if (is_smaller_pair(mid, lo))
                PSWAP(*mid, *lo);

jump_over:
          left_ptr  = lo + 1;
          right_ptr = hi - 1;

          /* Here's the famous ``collapse the walls'' section of quicksort.
             Gotta like those tight inner loops!  They are the main reason
             that this algorithm runs much faster than others. */
          do {
              while (is_smaller_pair(left_ptr, mid))
                  left_ptr++;

              while (is_smaller_pair(mid, right_ptr))
                  right_ptr--;

              if (left_ptr < right_ptr) {
                  PSWAP(*left_ptr, *right_ptr);
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
        Pair *const end_ptr = array + (n - 1);
        Pair *tmp_ptr = array;
        Pair *thresh = MIN(end_ptr, array + MAX_THRESH);
        Pair *run_ptr;

        /* Find smallest element in first threshold and place it at the
           array's beginning.  This is the smallest array element,
           and the operation speeds up insertion sort's inner loop. */

        for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++) {
            if (is_smaller_pair(run_ptr, tmp_ptr))
                tmp_ptr = run_ptr;
        }

        if (tmp_ptr != array)
            PSWAP(*tmp_ptr, *array);

        /* Insertion sort, running from left-hand-side up to right-hand-side.
         */

        run_ptr = array + 1;
        while (++run_ptr <= end_ptr) {
            tmp_ptr = run_ptr - 1;
            while (is_smaller_pair(run_ptr, tmp_ptr))
                tmp_ptr--;

            tmp_ptr++;
            if (tmp_ptr != run_ptr) {
                Pair *hi, *lo;
                Pair d;

                d = *run_ptr;
                for (hi = lo = run_ptr; --lo >= tmp_ptr; hi = lo)
                    *hi = *lo;
                *hi = d;
            }
        }
    }
}

static void
waterpour_sort(const gdouble *d, gint *idx, gint n)
{
    Pair *pairs = g_new(Pair, n);
    gint k;

    for (k = 0; k < n; k++)
        pairs[k] = (Pair){ d[k], k };

    sort_pairs(pairs, n);

    for (k = 0; k < n; k++)
        idx[k] = pairs[k].k;

    g_free(pairs);
}

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
