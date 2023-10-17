/*
 *  $Id: linestats.c 25308 2023-04-21 13:16:28Z yeti-dn $
 *  Copyright (C) 2003-2019 David Necas (Yeti), Petr Klapetek.
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
#include <libgwyddion/gwymath.h>
#include <libprocess/inttrans.h>
#include <libprocess/simplefft.h>
#include <libprocess/linestats.h>
#include "libgwyddion/gwyomp.h"
#include "gwyfftw.h"

/**
 * gwy_data_line_get_max:
 * @data_line: A data line.
 *
 * Finds the maximum value of a data line.
 *
 * Returns: The maximum value.
 **/
gdouble
gwy_data_line_get_max(GwyDataLine *data_line)
{
    gint i;
    gdouble max;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), -G_MAXDOUBLE);

    max = data_line->data[0];
    for (i = 1; i < data_line->res; i++) {
        if (G_UNLIKELY(data_line->data[i] > max))
            max = data_line->data[i];
    }
    return max;
}

/**
 * gwy_data_line_get_min:
 * @data_line: A data line.
 *
 * Finds the minimum value of a data line.
 *
 * Returns: The minimum value.
 **/
gdouble
gwy_data_line_get_min(GwyDataLine *data_line)
{
    gint i;
    gdouble min;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), G_MAXDOUBLE);

    min = data_line->data[0];
    for (i = 1; i < data_line->res; i++) {
        if (G_UNLIKELY(data_line->data[i] < min))
            min = data_line->data[i];
    }
    return min;
}

/**
 * gwy_data_line_min_pos_i:
 * @data_line: A data line.
 *
 * Finds the minimum pixel position of a data line.
 *
 * For historical reasons the value is returned as double, but it is always an integer.
 *
 * Returns: The minimum pixel position.
 *
 * Since 2.48
 **/
gdouble
gwy_data_line_min_pos_i(GwyDataLine *data_line)
{
    gint i, res, minpos = 0;
    gdouble min = G_MAXDOUBLE, *val;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    res = data_line->res;
    val = data_line->data;
    for (i = 0; i < res; i++) {
        if (G_UNLIKELY(*(val++) < min)) {
            min = *val;
            minpos = i;
        }
    }
    return minpos;
}

/**
 * gwy_data_line_max_pos_i:
 * @data_line: A data line.
 *
 * Finds the maximum pixel position of a data line.
 *
 * For historical reasons the value is returned as double, but it is always an integer.
 *
 * Returns: The maximum pixel position.
 *
 * Since 2.48
 **/
gdouble
gwy_data_line_max_pos_i(GwyDataLine *data_line)
{
    gint i, res, maxpos = 0;
    gdouble max = G_MINDOUBLE, *val;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    res = data_line->res;;
    val = data_line->data;
    for (i = 0; i < res; i++) {
        if (G_UNLIKELY(*(val++) > max)) {
            max = *val;
            maxpos = i;
        }
    }
    return maxpos;
}

/**
 * gwy_data_line_min_pos_r:
 * @data_line: A data line.
 *
 * Finds the real minimum position in a data line.
 *
 * Returns: Real position for the minimum.
 *
 * Since 2.48
 **/
gdouble
gwy_data_line_min_pos_r(GwyDataLine *data_line)
{
    gint res;
    gdouble real, offset, pos;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    res = data_line->res;
    if (G_UNLIKELY(res == 0.0)) {
        return 0.0;
    }
    real = data_line->real;
    offset = data_line->off;

    pos = gwy_data_line_min_pos_i(data_line);

    return pos * real / res + offset;
}

/**
 * gwy_data_line_max_pos_r:
 * @data_line: A data line.
 *
 * Finds the real maximum position in a data line.
 *
 * Returns: Real position for the maximum.
 *
 * Since 2.48
 **/
gdouble
gwy_data_line_max_pos_r(GwyDataLine *data_line)
{
    gint res;
    gdouble real, offset, pos;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    res = data_line->res;
    if (G_UNLIKELY(res == 0.0)) {
        return 0.0;
    }
    real = data_line->real;
    offset = data_line->off;

    pos = gwy_data_line_max_pos_i(data_line);

    return pos * real / res + offset;
}

/**
 * gwy_data_line_get_avg:
 * @data_line: A data line.
 *
 * Computes average value of a data line.
 *
 * Returns: Average value
 **/
gdouble
gwy_data_line_get_avg(GwyDataLine *data_line)
{
    gint i;
    gdouble avg = 0;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    for (i = 0; i < data_line->res; i++)
        avg += data_line->data[i];

    return avg/(gdouble)data_line->res;
}

/**
 * gwy_data_line_get_rms:
 * @data_line: A data line.
 *
 * Computes root mean square value of a data line.
 *
 * Returns: Root mean square deviation of values.
 **/
gdouble
gwy_data_line_get_rms(GwyDataLine *data_line)
{
    gint i;
    gdouble avg, sum2 = 0;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    avg = gwy_data_line_get_avg(data_line);
    for (i = 0; i < data_line->res; i++)
        sum2 += (data_line->data[i] - avg)*(data_line->data[i] - avg);

    return sqrt(sum2/data_line->res);
}

/**
 * gwy_data_line_get_sum:
 * @data_line: A data line.
 *
 * Computes sum of all values in a data line.
 *
 * Returns: sum of all the values.
 **/
gdouble
gwy_data_line_get_sum(GwyDataLine *data_line)
{
    gint i;
    gdouble sum = 0;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    for (i = 0; i < data_line->res; i++)
        sum += data_line->data[i];

    return sum;
}

/**
 * gwy_data_line_part_get_max:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Finds the maximum value of a part of a data line.
 *
 * Returns: Maximum within given interval.
 **/
gdouble
gwy_data_line_part_get_max(GwyDataLine *a,
                           gint from, gint to)
{
    gint i;
    gdouble max = -G_MAXDOUBLE;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), max);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, max);

    for (i = from; i < to; i++) {
        if (max < a->data[i])
            max = a->data[i];
    }
    return max;
}

/**
 * gwy_data_line_part_get_min:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Finds the minimum value of a part of a data line.
 *
 * Returns: Minimum within given interval.
 **/
gdouble
gwy_data_line_part_get_min(GwyDataLine *a,
                           gint from, gint to)
{
    gint i;
    gdouble min = G_MAXDOUBLE;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), min);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, min);

    for (i = from; i < to; i++) {
        if (min > a->data[i])
            min = a->data[i];
    }

    return min;
}

/**
 * gwy_data_line_part_get_min_max:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 * @min: Location to store minimum to.
 * @max: Location to store maximum to.
 *
 * Finds the minimum and maximum values of a part of a data line.
 *
 * Since 2.42
 **/
void
gwy_data_line_part_get_min_max(GwyDataLine *a,
                               gint from, gint to,
                               gdouble *min, gdouble *max)
{
    gint i;
    gdouble min1 = G_MAXDOUBLE, max1 = -G_MAXDOUBLE;

    g_return_if_fail(GWY_IS_DATA_LINE(a));
    g_return_if_fail(from >= 0 && to <= a->res);
    /* If both are not requested do not bother finding them. */
    if (!max && min)
        *min = gwy_data_line_part_get_min(a, from, to);
    if (max && !min)
        *max = gwy_data_line_part_get_max(a, from, to);
    if (!max || !min)
        return;

    GWY_ORDER(gint, from, to);

    for (i = from; i < to; i++) {
        if (min1 > a->data[i])
            min1 = a->data[i];
        if (max1 < a->data[i])
            max1 = a->data[i];
    }

    *min = min1;
    *max = max1;
}

/**
 * gwy_data_line_get_min_max:
 * @data_line: A data line.
 * @min: Location to store minimum to.
 * @max: Location to store maximum to.
 *
 * Finds the minimum and maximum values of a data line.
 *
 * Since 2.42
 **/
void
gwy_data_line_get_min_max(GwyDataLine *a,
                          gdouble *min, gdouble *max)
{
    g_return_if_fail(GWY_IS_DATA_LINE(a));
    gwy_data_line_part_get_min_max(a, 0, a->res, min, max);
}

/**
 * gwy_data_line_part_get_avg:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes mean value of all values in a part of a data line.
 *
 * Returns: Average value within given interval.
 **/
gdouble
gwy_data_line_part_get_avg(GwyDataLine *a, gint from, gint to)
{
    return gwy_data_line_part_get_sum(a, from, to)/(gdouble)(to-from);
}

/**
 * gwy_data_line_part_get_rms:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes root mean square value of a part of a data line.
 *
 * Returns: Root mean square deviation of heights within a given interval
 **/
gdouble
gwy_data_line_part_get_rms(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble rms = 0.0;
    gdouble avg;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), rms);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, rms);

    avg = gwy_data_line_part_get_avg(a, from, to);
    for (i = from; i < to; i++)
        rms += (avg - a->data[i])*(avg - a->data[i]);

    return sqrt(rms)/(to-from);
}

/**
 * gwy_data_line_get_ra:
 * @data_line: A data line.
 *
 * Computes the mean absolute deviation of a data line.
 *
 * Returns: The mean absolute deviation of height values.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_ra(GwyDataLine *a)
{
    return gwy_data_line_part_get_ra(a, 0, a->res);
}

/**
 * gwy_data_line_get_skew:
 * @data_line: A data line.
 *
 * Computes the skew of a data line.
 *
 * Returns: The skew of height values.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_skew(GwyDataLine *a)
{
    return gwy_data_line_part_get_skew(a, 0, a->res);
}

/**
 * gwy_data_line_get_kurtosis:
 * @data_line: A data line.
 *
 * Computes the kurtosis of a data line.
 *
 * Returns: The kurtosis of height values.
 *
 * Note the kurtosis returned by this function returns is the excess kurtosis which is zero for the Gaussian
 * distribution (not 3).
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_kurtosis(GwyDataLine *a)
{
    return gwy_data_line_part_get_kurtosis(a, 0, a->res);
}

/**
 * gwy_data_line_part_get_ra:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes mean absolute deviation value of a part of a data line.
 *
 * Returns: Mean absolute deviation of heights within a given interval.
 **/
gdouble
gwy_data_line_part_get_ra(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble ra = 0.0;
    gdouble avg;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), ra);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, ra);

    avg = gwy_data_line_part_get_avg(a, from, to);
    for (i = from; i < to; i++)
        ra += fabs(a->data[i] - avg);

    return ra/(to - from);
}

/**
 * gwy_data_line_part_get_skew:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes skew value of a part of a data line.
 *
 * Returns: Skew of heights within a given interval.
 **/
gdouble
gwy_data_line_part_get_skew(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble rms = 0.0, skew = 0.0;
    gdouble avg;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), skew);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, skew);

    avg = gwy_data_line_part_get_avg(a, from, to);
    for (i = from; i < to; i++) {
        gdouble d = a->data[i] - avg;
        rms += d*d;
        skew += d*d*d;
    }

    if (!rms)
        return 0.0;

    return skew*sqrt(to+1 - from)/pow(rms, 1.5);
}

/**
 * gwy_data_line_part_get_kurtosis:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes kurtosis value of a part of a data line.
 *
 * Note the kurtosis returned by this function returns is the excess kurtosis which is zero for the Gaussian
 * distribution (not 3).
 *
 * Returns: Kurtosis of heights within a given interval.
 **/
gdouble
gwy_data_line_part_get_kurtosis(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble rms = 0.0, kurtosis = 0.0;
    gdouble avg;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), kurtosis);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, kurtosis);

    avg = gwy_data_line_part_get_avg(a, from, to);
    for (i = from; i < to; i++) {
        gdouble d = a->data[i] - avg;
        d *= d;
        rms += d;
        kurtosis += d*d;
    }

    if (!rms)
        return 0.0;

    return kurtosis*(to+1 - from)/(rms*rms) - 3.0;
}

/**
 * gwy_data_line_get_tan_beta0:
 * @data_line: A data line.
 *
 * Computes root mean square slope in a data line.
 *
 * Returns: Root mean square slope within a given interval.
 *
 * Since: 2.2
 **/
gdouble
gwy_data_line_get_tan_beta0(GwyDataLine *a)
{
    return gwy_data_line_part_get_tan_beta0(a, 0, a->res);
}

/**
 * gwy_data_line_part_get_tan_beta0:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes root mean square slope in a part of a data line.
 *
 * This is the root mean square of value derivatives, it is also proportional to the second derivative of both HHCF
 * and ACF at zero.
 *
 * This roughness quantity is also known as Dq.
 *
 * Returns: Root mean square slope within a given interval.
 *
 * Since: 2.2
 **/
gdouble
gwy_data_line_part_get_tan_beta0(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble rms = 0.0;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), rms);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, rms);

    if (to - from < 2)
        return rms;

    for (i = from + 1; i < to; i++)
        rms += (a->data[i] - a->data[i-1])*(a->data[i] - a->data[i-1]);

    return sqrt(rms/(to-from - 1)) * a->res/a->real;
}

/**
 * gwy_data_line_part_get_variation:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes the total variation of a part of a data line.
 *
 * The total variation is estimated as the integral of the absolute value of local gradient.  For one dimensional
 * data, the variation reduces to the integral of absolute value of the derivative.  Its units are thus the same as
 * the value units of the line.  See also gwy_data_field_area_get_variation() for some more discussion.
 *
 * Returns: The total variation within a given interval.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_part_get_variation(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble var = 0.0;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), var);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, var);

    if (to - from < 2)
        return var;

    for (i = from + 1; i < to; i++)
        var += fabs(a->data[i] - a->data[i-1]);

    return var;
}

/**
 * gwy_data_line_get_variation:
 * @data_line: A data line.
 *
 * Computes the total variation of a data line.
 *
 * See gwy_data_line_part_get_variation() for definition and discussion.
 *
 * Returns: The total variation.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_variation(GwyDataLine *a)
{
    return gwy_data_line_part_get_variation(a, 0, a->res);
}

/**
 * gwy_data_line_part_get_sum:
 * @data_line: A data line.
 * @from: Index the line part starts at.
 * @to: Index the line part ends at + 1.
 *
 * Computes sum of all values in a part of a data line.
 *
 * Returns: Sum of all values within the interval.
 **/
gdouble
gwy_data_line_part_get_sum(GwyDataLine *a, gint from, gint to)
{
    gint i;
    gdouble sum = 0;

    g_return_val_if_fail(GWY_IS_DATA_LINE(a), sum);
    GWY_ORDER(gint, from, to);

    g_return_val_if_fail(from >= 0 && to <= a->res, sum);

    for (i = from; i < to; i++)
        sum += a->data[i];

    return sum;
}

static void
gwy_data_line_func_fft(GwyDataLine *data_line, GwyDataLine *target_line,
                       GwyFFTAreaFunc func)
{
    gint i, n, res;
    gdouble avg;
    fftw_plan plan;
    GwyDataLine *din, *dout;

    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    avg = gwy_data_line_get_avg(data_line);
    n = data_line->res;
    res = gwy_fft_find_nice_size(2*n);
    gwy_data_line_resample(target_line, n, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(target_line);

    din = gwy_data_line_new(res, 1.0, FALSE);
    dout = gwy_data_line_new(res, 1.0, FALSE);
    plan = gwy_fftw_plan_r2r_1d(res, din->data, dout->data, FFTW_R2HC, FFTW_ESTIMATE);

    /* In principle, we only need this for ACF, but removing constant offset
     * never hurts. */
    for (i = 0; i < n; i++)
        din->data[i] = data_line->data[i] - avg;

    func(plan, din, dout, target_line);
    g_object_unref(din);
    g_object_unref(dout);
    fftw_destroy_plan(plan);

    /* Set correct properties and units. */
    target_line->real = data_line->real;
    target_line->off = 0.0;
    gwy_data_line_multiply(target_line, 1.0/res);
    gwy_si_unit_assign(gwy_data_line_get_si_unit_x(target_line), gwy_data_line_get_si_unit_x(data_line));
    gwy_si_unit_power(gwy_data_line_get_si_unit_y(data_line), 2, gwy_data_line_get_si_unit_y(target_line));
}

/**
 * gwy_data_line_acf:
 * @data_line: A data line.
 * @target_line: Data line to store autocorrelation function to.  It will be
 *               resized to @data_line size.
 *
 * Coputes autocorrelation function and stores the values in @target_line
 *
 * Up to version 2.53 it did not set the output units properly.
 **/
void
gwy_data_line_acf(GwyDataLine *data_line, GwyDataLine *target_line)
{
    gwy_data_line_func_fft(data_line, target_line, do_fft_acf);
}

/**
 * gwy_data_line_hhcf:
 * @data_line: A data line.
 * @target_line: Data line to store height-height function to.  It will be resized to @data_line size.
 *
 * Computes height-height correlation function and stores results in @target_line.
 *
 * Up to version 2.53 it did not set the output units properly.
 **/
void
gwy_data_line_hhcf(GwyDataLine *data_line, GwyDataLine *target_line)
{
    gwy_data_line_func_fft(data_line, target_line, do_fft_hhcf);
}

/**
 * gwy_data_line_psdf:
 * @data_line: A data line.
 * @target_line: Data line to store power spectral density function to. It will be resized to @data_line size.
 * @windowing: Windowing method to use.
 * @interpolation: Interpolation type. Ignored since 2.8 as no resampling is performed.
 *
 * Calculates the power spectral density function of a data line.
 *
 * Up to version 2.45 it destroyed the input data and did not set the output units properly.
 **/
void
gwy_data_line_psdf(GwyDataLine *data_line,
                   GwyDataLine *target_line,
                   gint windowing,
                   gint interpolation)
{
    GwyDataLine *iin, *rout, *iout;
    GwySIUnit *xunit, *yunit, *lineunit;
    gdouble *data, *rdata, *idata;
    gdouble q;
    gint i, res;

    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    res = data_line->res;
    iin = gwy_data_line_new_alike(data_line, TRUE);
    rout = gwy_data_line_new_alike(data_line, FALSE);
    iout = gwy_data_line_new_alike(data_line, FALSE);
    gwy_data_line_resample(target_line, res/2, GWY_INTERPOLATION_NONE);
    gwy_data_line_fft(data_line, iin, rout, iout, windowing, GWY_TRANSFORM_DIRECTION_FORWARD, interpolation, TRUE, 2);

    data = target_line->data;
    rdata = rout->data;
    idata = iout->data;
    q = data_line->real/(res*2.0*G_PI);

    /* Calculate modulus */
    for (i = 0; i < res/2; i++)
        data[i] = q*(rdata[i]*rdata[i] + idata[i]*idata[i]);

    target_line->real = 2*G_PI*target_line->res/data_line->real;
    target_line->off = 0.0;

    g_object_unref(rout);
    g_object_unref(iin);
    g_object_unref(iout);

    /* Set proper units */
    xunit = gwy_data_line_get_si_unit_x(data_line);
    yunit = gwy_data_line_get_si_unit_y(data_line);
    lineunit = gwy_data_line_get_si_unit_x(target_line);
    gwy_si_unit_power(xunit, -1, lineunit);
    lineunit = gwy_data_line_get_si_unit_y(target_line);
    gwy_si_unit_power(yunit, 2, lineunit);
    gwy_si_unit_multiply(lineunit, xunit, lineunit);
}

/**
 * gwy_data_line_distribution:
 * @data_line: A data line.
 * @distribution: Data line to put the distribution of @data_line values to. It will be resampled to @nstats samples
 *                (or the automatically chosen number of bins).
 * @ymin: Start of value range, pass @ymin = @ymax = 0.0 for the full range.
 * @ymax: End of value range.
 * @normalize_to_unity: %TRUE to normalize the integral to unity (including setting y-units of output to the inverse
 *                      of x-units), %FALSE to keep plain counts in the output (and set y-units to none).
 * @nstats: The requested number of histogram bins, pass a non-positive number to automatically choose a suitable
 *          number of bins.
 *
 * Calculates the distribution of data line values.
 *
 * This function is quite similar to gwy_data_line_dh(), the differences are: output normalization (chosen to make the
 * integral unity), output units (again set to make the integral unity), automated binning.
 *
 * Note the @i-th bin is [@i*@dx+@off,(@i+1)*@dx+@off] so the central value you probably want to use for plotting is
 * (@i+0.5)*@dx+@off (where @dx is the @distribution data line pixel size, @off is its offset).
 *
 * If all values are equal and @ymin, @ymax are not explictly specified, the range is chosen as [@v-|@v|/2,@v+|@v/2]
 * where @v is the unique value, except when @v=0, in which case the range is set to [-1,1].
 *
 * Since: 2.8
 **/
void
gwy_data_line_distribution(GwyDataLine *data_line,
                           GwyDataLine *distribution,
                           gdouble ymin,
                           gdouble ymax,
                           gboolean normalize_to_unity,
                           gint nstats)
{
    GwySIUnit *yunit, *lineunit;
    const gdouble *data;
    gint i, res, ndata;
    guint *counts;
    gdouble s;

    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(distribution));

    /* Find reasonable binning */
    GWY_ORDER(gdouble, ymin, ymax);

    /* if ymin == ymax == 0 use the full range */
    if (!ymin && !ymax) {
        ymin = gwy_data_line_get_min(data_line);
        ymax = gwy_data_line_get_max(data_line);
        if (ymin > 0.0 && ymin <= 0.1*ymax)
            ymin = 0.0;
        else if (ymax < 0.0 && ymax >= 0.1*ymin)
            ymax = 0.0;
    }
    if (ymin == ymax) {
        if (ymax) {
            ymin -= 0.5*fabs(ymin);
            ymax += 0.5*fabs(ymax);
        }
        else {
            ymin = -1.0;
            ymax = 1.0;
        }
    }

    res = data_line->res;
    data = data_line->data;
    if (nstats < 1) {
        ndata = 0;
        for (i = 0; i < res; i++) {
            if (data[i] >= ymin && data[i] <= ymax)
                ndata++;
        }
        nstats = floor(3.49*cbrt(ndata) + 0.5);
        nstats = MAX(nstats, 2);
    }

    gwy_debug("min: %g, max: %g, nstats: %d", ymin, ymax, nstats);
    s = (ymax - ymin)/(nstats - 1e-9);

    /* Fill histogram */
    gwy_data_line_resample(distribution, nstats, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(distribution);

    counts = g_new0(guint, nstats);
    ndata = gwy_math_histogram(data, res, ymin, ymax, nstats, counts);
    for (i = 0; i < nstats; i++)
        distribution->data[i] = counts[i];
    g_free(counts);

    /* Set proper units and scales */
    distribution->real = ymax - ymin;
    distribution->off = ymin;

    yunit = gwy_data_line_get_si_unit_y(data_line);
    lineunit = gwy_data_line_get_si_unit_x(distribution);
    gwy_si_unit_power(yunit, 1, lineunit);
    lineunit = gwy_data_line_get_si_unit_y(distribution);
    if (normalize_to_unity) {
        gwy_data_line_multiply(distribution, 1.0/(ndata*s));
        gwy_si_unit_power(yunit, -1, lineunit);
    }
    else
        gwy_si_unit_set_from_string(lineunit, NULL);
}

/**
 * gwy_data_line_dh:
 * @data_line: A data line.
 * @target_line: Data line to store height distribution function to. It will be resized to @nsteps.
 * @ymin: Height distribution minimum value.
 * @ymax: Height distribution maximum value.
 * @nsteps: Number of histogram steps.
 *
 * Computes distribution of heights in interval [@ymin, @ymax).
 *
 * If the interval is (0, 0) it computes the distribution from real data minimum and maximum value.
 **/
void
gwy_data_line_dh(GwyDataLine *data_line,
                 GwyDataLine *target_line,
                 gdouble ymin, gdouble ymax,
                 gint nsteps)
{
    gint i, n;
    gdouble step;
    guint *counts;

    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    n = data_line->res;
    gwy_data_line_resample(target_line, nsteps, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(target_line);

    /* if ymin == ymax == 0 we want to set up histogram area */
    if (!ymin && !ymax) {
        ymin = gwy_data_line_get_min(data_line);
        ymax = gwy_data_line_get_max(data_line);
    }
    step = (ymax - ymin)/(nsteps - 1.0);

    counts = g_new0(guint, nsteps);
    n = gwy_math_histogram(data_line->data, n, ymin, ymax, nsteps, counts);
    for (i = 0; i < nsteps; i++)
        target_line->data[i] = counts[i];
    g_free(counts);

    gwy_data_line_multiply(target_line, 1.0/(n*step));
    target_line->off = ymin;
    target_line->real = ymax - ymin;
}

/**
 * gwy_data_line_cdh:
 * @data_line: A data line.
 * @target_line: Data line to store height distribution function to. It will be resized to @nsteps.
 * @ymin: Height distribution minimum value.
 * @ymax: Height distribution maximum value.
 * @nsteps: Number of histogram steps.
 *
 * Computes cumulative distribution of heighs in interval [@ymin, @ymax).
 *
 * If the interval is (0, 0) it computes the distribution from real data minimum and maximum value.
 **/
void
gwy_data_line_cdh(GwyDataLine *data_line,
                  GwyDataLine *target_line,
                  gdouble ymin, gdouble ymax,
                  gint nsteps)
{
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    gwy_data_line_dh(data_line, target_line, ymin, ymax, nsteps);
    gwy_data_line_cumulate(target_line);
    target_line->data[target_line->res-1] = 1.0;   /* Fix rounding errors. */
}

/**
 * gwy_data_line_da:
 * @data_line: A data line.
 * @target_line: Data line to store angle distribution function to.
 * @ymin: Angle distribution minimum value.
 * @ymax: Angle distribution maximum value.
 * @nsteps: Mumber of angular histogram steps.
 *
 * Computes distribution of angles in interval [@ymin, @ymax).
 *
 * If the interval is (0, 0) it computes the distribution from real data minimum and maximum angle value.
 **/
void
gwy_data_line_da(GwyDataLine *data_line,
                 GwyDataLine *target_line,
                 gdouble ymin, gdouble ymax,
                 gint nsteps)
{
    gint i, n, val;
    gdouble step, angle, imin;

    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    n = data_line->res;
    gwy_data_line_resample(target_line, nsteps, GWY_INTERPOLATION_NONE);
    gwy_data_line_clear(target_line);

    /* if ymin == ymax == 0 we want to set up histogram area */
    if (!ymin && !ymax) {
        ymin = G_MAXDOUBLE;
        ymax = -G_MAXDOUBLE;
        for (i = 0; i < n; i++) {
            angle = gwy_data_line_get_der(data_line, i);
            if (ymin > angle)
                ymin = angle;
            if (ymax < angle)
                ymax = angle;
        }
    }
    step = (ymax - ymin)/(nsteps - 1.0);
    imin = ymin/step;

    for (i = 0; i < n; i++) {
        val = (gint)(gwy_data_line_get_der(data_line, i)/step - imin);
        if (G_UNLIKELY(val < 0))
            val = 0; /* this should never happened */
        if (G_UNLIKELY(val >= nsteps))
            val = nsteps-1; /* this should never happened */
        target_line->data[val] += 1.0;
    }
    target_line->real = ymax - ymin;
    target_line->off = ymin;
}

/**
 * gwy_data_line_cda:
 * @data_line: A data line.
 * @target_line: Data line to store angle distribution function to. It will be resized to @nsteps.
 * @ymin: Angle distribution minimum value.
 * @ymax: Angle distribution maximum value.
 * @nsteps: Number of angular histogram steps.
 *
 * Computes cumulative distribution of angles in interval [@ymin, @ymax).
 *
 * If the interval is (0, 0) it computes the distribution from real data minimum and maximum angle value.
 **/
void
gwy_data_line_cda(GwyDataLine *data_line,
                  GwyDataLine *target_line,
                  gdouble ymin, gdouble ymax,
                  gint nsteps)
{
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_LINE(target_line));

    gwy_data_line_da(data_line, target_line, ymin, ymax, nsteps);
    gwy_data_line_cumulate(target_line);
    target_line->data[target_line->res-1] = 1.0;   /* Fix rounding errors. */
}

/**
 * gwy_data_line_get_length:
 * @data_line: A data line to compute length of.
 *
 * Calculates physical length of a data line.
 *
 * The length is calculated from approximation by straight segments between values.
 *
 * Returns: The line length.
 **/
gdouble
gwy_data_line_get_length(GwyDataLine *data_line)
{
    gdouble sum, q;
    gint i, n;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0.0);

    n = data_line->res;
    q = data_line->real/n;
    if (G_UNLIKELY(data_line->res == 1))
        return q;

    sum = 0.0;
    for (i = 1; i < n; i++)
        sum += hypot(q, data_line->data[i] - data_line->data[i-1]);

    /* We calculate length of inner part of a segment.  If we assume the average properties of border are the same as
     * of the inner part, we can simply multiply the sum with the total/inner length ratio */
    sum *= n/(n - 1.0);

    return sum;
}

/**
 * gwy_data_line_part_get_modus:
 * @data_line: A data line.
 * @from: The index in @data_line to start from (inclusive).
 * @to: The index in @data_line to stop (noninclusive).
 * @histogram_steps: Number of histogram steps used for modus searching, pass a nonpositive number to autosize.
 *
 * Finds approximate modus of a data line part.
 *
 * As each number in the data line is usually unique, this function does not return modus of the data itself, but
 * modus of a histogram.
 *
 * Returns: The modus.
 **/
gdouble
gwy_data_line_part_get_modus(GwyDataLine *data_line,
                             gint from, gint to,
                             gint histogram_steps)
{
    gint *histogram;
    gint n, i, j, m;
    gdouble min, max, sum;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0);
    g_return_val_if_fail(from >= 0 && to <= data_line->res, 0);
    g_return_val_if_fail(from != to, 0);

    GWY_ORDER(gint, from, to);
    n = to - from;

    if (n == 1)
        return data_line->data[from];

    if (histogram_steps < 1) {
        /*
        gdouble sigma = gwy_data_line_part_get_rms(data_line, from, to);
        histogram_steps = floor(0.49*sigma*cbrt(n) + 0.5);
        */
        histogram_steps = floor(3.49*cbrt(n) + 0.5);
        gwy_debug("histogram_steps = %d", histogram_steps);
    }

    min = gwy_data_line_part_get_min(data_line, from, to);
    max = gwy_data_line_part_get_max(data_line, from, to);
    if (min == max)
        return min;

    histogram = g_new0(gint, histogram_steps);
    for (i = from; i < to; i++) {
        j = (data_line->data[i] - min)/(max - min)*histogram_steps;
        j = CLAMP(j, 0, histogram_steps-1);
        histogram[j]++;
    }

    m = 0;
    for (i = 1; i < histogram_steps; i++) {
        if (histogram[i] > histogram[m])
            m = i;
    }

    n = 0;
    sum = 0.0;
    for (i = from; i < to; i++) {
        j = (data_line->data[i] - min)/(max - min)*histogram_steps;
        j = CLAMP(j, 0, histogram_steps-1);
        if (j == m) {
            sum += data_line->data[i];
            n++;
        }
    }

    g_free(histogram);
    gwy_debug("modus = %g", sum/n);

    return sum/n;
}

/**
 * gwy_data_line_get_modus:
 * @data_line: A data line.
 * @histogram_steps: Number of histogram steps used for modus searching, pass a nonpositive number to autosize.
 *
 * Finds approximate modus of a data line.
 *
 * See gwy_data_line_part_get_modus() for details and caveats.
 *
 * Returns: The modus.
 **/
gdouble
gwy_data_line_get_modus(GwyDataLine *data_line,
                        gint histogram_steps)
{
    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0);

    return gwy_data_line_part_get_modus(data_line, 0, data_line->res, histogram_steps);
}

/**
 * gwy_data_line_part_get_median:
 * @data_line: A data line.
 * @from: The index in @data_line to start from (inclusive).
 * @to: The index in @data_line to stop (noninclusive).
 *
 * Finds median of a data line part.
 *
 * Returns: The median.
 *
 * Since: 2.1
 **/
gdouble
gwy_data_line_part_get_median(GwyDataLine *data_line,
                              gint from, gint to)
{
    gdouble *d;
    gdouble med;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0);
    g_return_val_if_fail(from >= 0 && to <= data_line->res, 0);
    g_return_val_if_fail(from != to, 0);

    d = g_memdup(data_line->data + from, (to - from)*sizeof(gdouble));
    med = gwy_math_median(to - from, d);
    g_free(d);

    return med;
}

/**
 * gwy_data_line_get_median:
 * @data_line: A data line.
 *
 * Finds median of a data line.
 *
 * Returns: The median.
 *
 * Since: 2.1
 **/
gdouble
gwy_data_line_get_median(GwyDataLine *data_line)
{
    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0);

    return gwy_data_line_part_get_median(data_line, 0, data_line->res);
}

static inline void
add_sorted(gdouble *vals, gint len, gint *n, gdouble value)
{
    gint i, j, nn = *n;

    for (i = 0; i < nn; i++) {
        if (value > vals[i]) {
            for (j = MIN(len-1, nn); j > i; j--)
                vals[j] = vals[j-1];
            vals[i] = value;
            *n = MIN(nn+1, len);
            return;
        }
    }
    if (!nn && len) {
        vals[0] = value;
        *n = 1;
    }
}

/**
 * gwy_data_line_get_kth_peaks:
 * @data_line: A data line.
 * @m: Number of sampling lengths the line is split into.
 * @rank: Rank of the peak to find.  One means the highest peak, three the third highers, etc.
 * @peaks: %TRUE for peaks, %FALSE for valleys.  If you pass %FALSE, swap the meanings of peaks and valleys in the
 *         description.  Valley depths are positive.
 * @average: Calculate the average of the first @rank peaks instead of the height of @rank-th peak.
 * @pthreshold: Peak height threshold.  Peaks must stick above this threshold.
 * @vthreshold: Valley depth threshold.  Valleys must fall below this threshold.  The depth is a positive value.
 * @peakvalues: Array of length at least @m where the peak heights in each sampling length should be stored.
 *
 * Calculate k-th largers peaks or valleys in a data line split into given number of sampling lengths.
 *
 * This is a general function that can be used as the base for various standard roughness quantities such as Rp, Rpm,
 * Rv, Rvm or R3z.  It is assumed the line is already levelled, the form removed, etc.
 *
 * See gwy_data_line_count_peaks() for the description what is considered a peak.
 *
 * For larger thresholds and/or short lines some sampling lengths may not contain the requested number of peaks.  If
 * there are any peaks at all, the smallest peak height (even though it is not @rank-th) is used.  If there are no
 * peaks, a large negative value is stored in the corresponding @peakvalues item.
 *
 * Returns: The actual number of peaks found (i.e. number of positive values in @peakvalues).
 *
 * Since: 2.50
 **/
gint
gwy_data_line_get_kth_peaks(GwyDataLine *data_line,
                            gint m, gint rank,
                            gboolean peaks, gboolean average,
                            gdouble pthreshold, gdouble vthreshold,
                            gdouble *peakvalues)
{
    const gdouble *data;
    gdouble *pvals = NULL;
    gint res, len, i, mm, npeaks, n;
    gdouble currpeak, d;
    gboolean seen_valley;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0);
    g_return_val_if_fail(rank > 0, 0);
    g_return_val_if_fail(m > 0, 0);
    g_return_val_if_fail(pthreshold >= 0.0, 0);
    g_return_val_if_fail(vthreshold >= 0.0, 0);
    g_return_val_if_fail(peakvalues, 0);

    if (!peaks)
        GWY_SWAP(gdouble, pthreshold, vthreshold);

    pvals = g_new(gdouble, rank);
    res = data_line->res;
    npeaks = 0;
    for (mm = 0; mm < m; mm++) {
        i = mm*res/m;
        len = (mm + 1)*res/m - i;
        data = data_line->data + i;
        /* Peak is a segment of the line that goes above the positive threshold, separated segments that go below the
         * negative threshold. Valley is the opposite.  Between them there can be things that are neither. */
        seen_valley = FALSE;
        currpeak = 0.0;
        n = 0;
        for (i = 0; i < len; i++) {
            d = peaks ? data[i] : -data[i];
            if (d > pthreshold) {
                /* Finish the previous peak, if any occured. */
                if (seen_valley) {
                    if (currpeak > 0.0)
                        add_sorted(pvals, rank, &n, currpeak);
                    seen_valley = FALSE;
                    currpeak = 0.0;
                }
                if (d > currpeak)
                    currpeak = d;
            }
            else if (d < -vthreshold)
                seen_valley = TRUE;
        }
        if (currpeak > 0.0)
            add_sorted(pvals, rank, &n, currpeak);

        /* XXX: There are several reasonable things we can do when we do not
         * find enough peaks... */
        /*if (n == rank)*/
        if (n) {
            if (average) {
                d = 0.0;
                for (i = 0; i < n; i++)
                    d += pvals[i];
                peakvalues[mm] = d/n;
            }
            else
                peakvalues[mm] = pvals[n-1];
            npeaks++;
        }
        else
            peakvalues[mm] = -G_MAXDOUBLE;
    }

    g_free(pvals);

    return npeaks;
}

/**
 * gwy_data_line_count_peaks:
 * @data_line: A data line.
 * @peaks: %TRUE for peaks, %FALSE for valleys.  If you pass %FALSE, swap the meanings of peaks and valleys in the
 *         description.  Valley depths are positive.
 * @pthreshold: Peak height threshold.  Peaks must stick above this threshold.
 * @vthreshold: Valley depth threshold.  Valleys must fall below this threshold.
 *
 * Counts peaks or valleys defined by thresholds in a data line.
 *
 * Peak is defined as a part of the profile that extends above the peak threshold and is separarted by valleys that
 * extend below the valley threshold.  For non-zero thresholds there may be parts between that are neither peaks not
 * valleys because the local maxima in them are insignificant.
 *
 * In either case, values of @pthreshold and @vthreshold must be non-negative. Usually one passes the same value for
 * both.
 *
 * Returns: The number of peaks found.
 *
 * Since: 2.50
 **/
gint
gwy_data_line_count_peaks(GwyDataLine *data_line, gboolean peaks,
                          gdouble pthreshold, gdouble vthreshold)
{
    gint res, i, peakcount = 0;
    gboolean seen_interruption = TRUE;
    const gdouble *data;
    gdouble d;

    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), 0);
    g_return_val_if_fail(pthreshold >= 0.0, 0);
    g_return_val_if_fail(vthreshold >= 0.0, 0);

    res = data_line->res;
    data = data_line->data;
    if (peaks) {
        for (i = 0; i < res; i++) {
            d = data[i];
            if (seen_interruption && d > pthreshold) {
                peakcount++;
                seen_interruption = FALSE;
            }
            else if (d < -vthreshold)
                seen_interruption = TRUE;
        }
    }
    else {
        for (i = 0; i < res; i++) {
            d = data[i];
            if (seen_interruption && d < -vthreshold) {
                peakcount++;
                seen_interruption = FALSE;
            }
            else if (d > pthreshold)
                seen_interruption = TRUE;
        }
    }

    return peakcount;
}

/**
 * gwy_data_line_get_xpm:
 * @data_line: A data line.
 * @m: Number of sampling lengths.
 * @k: Number of peaks to consider.
 *
 * Calculates a peak roughness quantity for a data line.
 *
 * Depending on @m and @k, the function can calculate Average Maximum Profile Peak Height @Rpm or Maximum Profile Peak
 * Height @Rp, @Pp, @Wp.
 *
 * Returns: The peak roughness quantity defined by @m and @k.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_xpm(GwyDataLine *data_line, gint m, gint k)
{
    gdouble *peaks;
    gdouble Xpm = 0.0;
    gint i, n;

    g_return_val_if_fail(m >= 1, Xpm);
    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), Xpm);
    g_return_val_if_fail(m >= 1, Xpm);

    peaks = g_new(gdouble, m);
    gwy_data_line_get_kth_peaks(data_line, m, k, TRUE, FALSE, 0.0, 0.0, peaks);
    n = 0;
    for (i = 0; i < m; i++) {
        if (peaks[i] > 0.0) {
            Xpm += peaks[i];
            n++;
        }
    }

    return n ? Xpm/n : 0.0;
}

/**
 * gwy_data_line_get_xvm:
 * @data_line: A data line.
 * @m: Number of sampling lengths.
 * @k: Number of valleys to consider.
 *
 * Calculates a valley roughness quantity for a data line.
 *
 * Depending on @m and @k, the function can calculate Average Maximum Profile Valley Depth @Rvm or Maximum Profile
 * Peak Depth @Rv, @Pv, @Wv.
 *
 * Returns: The valley roughness quantity defined by @m and @k.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_xvm(GwyDataLine *data_line, gint m, gint k)
{
    gdouble *peaks;
    gdouble Xpm = 0.0;
    gint i, n;

    g_return_val_if_fail(m >= 1, Xpm);
    g_return_val_if_fail(GWY_IS_DATA_LINE(data_line), Xpm);
    g_return_val_if_fail(m >= 1, Xpm);

    peaks = g_new(gdouble, m);
    gwy_data_line_get_kth_peaks(data_line, m, k, FALSE, FALSE, 0.0, 0.0, peaks);
    n = 0;
    for (i = 0; i < m; i++) {
        if (peaks[i] > 0.0) {
            Xpm += peaks[i];
            n++;
        }
    }

    return n ? Xpm/n : 0.0;
}

/**
 * gwy_data_line_get_xtm:
 * @data_line: A data line.
 * @m: Number of sampling lengths.
 * @k: Number of peaks and valleys to consider.
 *
 * Calculates a total roughness quantity for a data line.
 *
 * The total quantity is just the sum of the corresponding quantities obtained by gwy_data_line_get_xpm() and
 * gwy_data_line_get_xvm().
 *
 * Returns: The total roughness quantity defined by @m and @k.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_line_get_xtm(GwyDataLine *data_line, gint m, gint k)
{
    return gwy_data_line_get_xpm(data_line, m, k) + gwy_data_line_get_xvm(data_line, m, k);
}

/************************** Documentation ****************************/

/**
 * SECTION:linestats
 * @title: linestats
 * @short_description: One-dimensional statistical functions
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
