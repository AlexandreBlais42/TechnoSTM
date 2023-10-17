/*
 *  $Id: peaks.c 24837 2022-05-27 12:26:05Z yeti-dn $
 *  Copyright (C) 2016-2021 David Necas (Yeti).
 *  E-mail: yeti@gwyddion.net.
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

#include <libgwyddion/gwymacros.h>
#include <libgwyddion/gwymath.h>
#include "peaks.h"

typedef struct {
    gdouble prominence;
    gdouble x;
    gdouble height;
    gdouble area;
    gdouble width;
    gint i;
} Peak;

struct _GwyPeaks {
    GArray *peaks;
    GwyPeakBackgroundType background;
    GwyPeakOrderType order;
};

static gint compare_prominence_descending(gconstpointer a,
                                          gconstpointer b);
static gint compare_abscissa_ascending   (gconstpointer a,
                                          gconstpointer b);
GType
gwy_peaks_get_type(void)
{
    static GType peaks_type = 0;

    if (G_UNLIKELY(!peaks_type)) {
        peaks_type = g_boxed_type_register_static("GwyPeaks",
                                                  (GBoxedCopyFunc)gwy_peaks_copy,
                                                  (GBoxedFreeFunc)gwy_peaks_free);
    }

    return peaks_type;
}

/**
 * gwy_peaks_new:
 *
 * Creates a new empty peak analyser.
 *
 * Returns: A new peak analyser.
 *
 * Since: 2.46
 **/
GwyPeaks*
gwy_peaks_new(void)
{
    GwyPeaks *peaks;

    peaks = g_slice_new(GwyPeaks);
    peaks->peaks = g_array_new(FALSE, FALSE, sizeof(Peak));
    peaks->background = GWY_PEAK_BACKGROUND_MMSTEP;
    peaks->order = GWY_PEAK_ORDER_ABSCISSA;

    return peaks;
}

/**
 * gwy_peaks_copy:
 * @peaks: A peak analyser.
 *
 * Creates a copy of a peak analyser.
 *
 * This is mostly useful for language bindings.
 *
 * Returns: A newly created peak analyser.
 *
 * Since: 2.47
 **/
GwyPeaks*
gwy_peaks_copy(GwyPeaks *peaks)
{
    GwyPeaks *copy = g_slice_dup(GwyPeaks, peaks);

    copy->peaks = g_array_new(FALSE, FALSE, sizeof(Peak));
    g_array_append_vals(copy->peaks, peaks->peaks->data, peaks->peaks->len);

    return copy;
}

/**
 * gwy_peaks_free:
 * @peaks: A peak analyser.
 *
 * Frees a peak analyser and all associated data.
 *
 * Since: 2.46
 **/
void
gwy_peaks_free(GwyPeaks *peaks)
{
    g_return_if_fail(peaks);
    g_array_free(peaks->peaks, TRUE);
    g_slice_free(GwyPeaks, peaks);
}

/**
 * gwy_peaks_set_background:
 * @peaks: A peak analyser.
 * @background: Background type to use in future analyses.
 *
 * Sets the background type a peak analyser will use.
 *
 * The default background is %GWY_PEAK_BACKGROUND_MMSTEP.  Note that the new background type will only be used in
 * future analyses; it does not change the results of the already performed analysis.
 *
 * Since: 2.46
 **/
void
gwy_peaks_set_background(GwyPeaks *peaks,
                         GwyPeakBackgroundType background)
{
    g_return_if_fail(peaks);
    peaks->background = background;
}

/**
 * gwy_peaks_set_order:
 * @peaks: A peak analyser.
 * @order: Order type to use in future analyses.
 *
 * Sets the order type a peak analyser will use.
 *
 * The default order is %GWY_PEAK_ORDER_ABSCISSA.  Note that the new order will only be effective in future analyses;
 * it does not change the results of the already performed analysis.
 *
 * Since: 2.46
 **/
void
gwy_peaks_set_order(GwyPeaks *peaks,
                    GwyPeakOrderType order)
{
    g_return_if_fail(peaks);
    peaks->order = order;
}

/**
 * gwy_peaks_analyze_xy:
 * @peaks: A peak analyser.
 * @xydata: Curve points (array with @n items) that must be ordered by @x values in ascending order.
 * @n: Number of data points in the curve.
 * @maxpeaks: Maximum number of the most prominent peaks to locate.
 *
 * Finds peaks a graph curve given as GwyXY data.
 *
 * The peaks are remembered by the analyser and their properties can be subsequently requested using
 * gwy_peaks_get_quantity().
 *
 * Returns: The number of peaks found.
 *
 * Since: 2.46
 **/
guint
gwy_peaks_analyze_xy(GwyPeaks *peaks,
                     const GwyXY *xydata,
                     guint n,
                     guint maxpeaks)
{
    gdouble *data;
    guint i, retval;

    g_return_val_if_fail(xydata, 0);

    data = g_new(gdouble, 2*n);
    for (i = 0; i < n; i++) {
        data[i] = xydata[i].x;
        data[n+i] = xydata[i].y;
    }
    retval = gwy_peaks_analyze(peaks, data, data+n, n, maxpeaks);
    g_free(data);

    return retval;
}

/**
 * gwy_peaks_analyze_dataline:
 * @peaks: A peak analyser.
 * @dline: Curve data as a data line.
 * @maxpeaks: Maximum number of the most prominent peaks to locate.
 *
 * Finds peaks a graph curve given as GwyDataLine.
 *
 * The peaks are remembered by the analyser and their properties can be subsequently requested using
 * gwy_peaks_get_quantity().
 *
 * Returns: The number of peaks found.
 *
 * Since: 2.46
 **/
guint
gwy_peaks_analyze_dataline(GwyPeaks *peaks,
                           GwyDataLine *dline,
                           guint maxpeaks)
{
    gdouble *xdata;
    gdouble xoff, dx;
    guint i, n, retval;

    g_return_val_if_fail(GWY_IS_DATA_LINE(dline), 0);
    g_return_val_if_fail(peaks, 0);

    n = gwy_data_line_get_res(dline);
    dx = gwy_data_line_get_real(dline)/n;
    xoff = gwy_data_line_get_offset(dline);
    xdata = g_new(gdouble, n);
    for (i = 0; i < n; i++)
        xdata[i] = (i + 0.5)*dx + xoff;
    retval = gwy_peaks_analyze(peaks, xdata, gwy_data_line_get_data(dline), n, maxpeaks);
    g_free(xdata);

    return retval;
}


/**
 * gwy_peaks_analyze:
 * @peaks: A peak analyser.
 * @xdata: Abscissa values (array with @n items), must be ordered in ascending order.
 * @ydata: Ordinate values corresponding to @xdata.
 * @n: Number of data points in the curve.
 * @maxpeaks: Maximum number of the most prominent peaks to locate.
 *
 * Finds peaks a graph curve given as separated @x and @y data.
 *
 * The peaks are remembered by the analyser and their properties can be subsequently requested using
 * gwy_peaks_get_quantity().
 *
 * Returns: The number of peaks found.
 *
 * Since: 2.46
 **/
guint
gwy_peaks_analyze(GwyPeaks *peaks,
                  const gdouble *xdata,
                  const gdouble *ydata,
                  guint n,
                  guint maxpeaks)
{
    GArray *p;
    gint i, k, flatsize;
    gdouble *ydata_filtered, *ydata_filtered2;

    g_return_val_if_fail(peaks, 0);
    g_return_val_if_fail(xdata, 0);
    g_return_val_if_fail(ydata, 0);

    p = peaks->peaks;
    g_array_set_size(p, 0);
    if (!n || !maxpeaks)
        return 0;

    /* Perform simple closing. */
    ydata_filtered = g_new(gdouble, n);
    ydata_filtered2 = g_new(gdouble, n);
    gwy_assign(ydata_filtered, ydata, n);
    for (k = 1; k < log(n) - 1.4; k++) {
        ydata_filtered2[0] = ydata_filtered[0];
        for (i = 1; i+1 < n; i++) {
            gdouble y = ydata_filtered[i];
            gdouble yl = 0.5*(ydata_filtered[i+1] + ydata_filtered[i-1]);
            ydata_filtered2[i] = MAX(y, yl);
        }
        ydata_filtered2[n-1] = ydata_filtered[n-1];
        GWY_SWAP(gdouble*, ydata_filtered2, ydata_filtered);
    }
    g_free(ydata_filtered2);

    /* Find local maxima. */
    flatsize = 0;
    for (i = 1; i+1 < n; i++) {
        gdouble y = ydata_filtered[i];
        gdouble yp = ydata_filtered[i-1];
        gdouble yn = ydata_filtered[i+1];

        /* The normal cases. */
        if (y < yp || y < yn)
            continue;
        if (y > yp && y > yn) {
            Peak peak;
            peak.i = i;
            g_array_append_val(p, peak);
            continue;
        }

        /* Flat tops. */
        if (y == yn && y > yp)
            flatsize = 0;
        else if (y == yn && y == yp)
            flatsize++;
        else if (y == yp && y > yn) {
            Peak peak;
            peak.i = i - flatsize/2;
            g_array_append_val(p, peak);
        }
    }

    /* Analyse prominence. */
    for (k = 0; k < p->len; k++) {
        Peak *peak = &g_array_index(p, Peak, k);
        gint ileft, iright;
        gdouble yleft, yright, arealeft, arearight, disp2left, disp2right;

        /* Find the peak extents. */
        for (ileft = peak->i - 1; ileft && ydata_filtered[ileft] == ydata_filtered[ileft+1]; ileft--)
            ;

        if (peaks->background == GWY_PEAK_BACKGROUND_ZERO) {
            while (ileft && ydata_filtered[ileft] > ydata_filtered[ileft-1] && ydata[ileft] > 0.0)
                ileft--;
            if (ydata[ileft] < 0.0)
                ileft++;
            yleft = 0.0;
        }
        else {
            while (ileft && ydata_filtered[ileft] > ydata_filtered[ileft-1])
                ileft--;
            yleft = ydata[ileft];
        }

        for (iright = peak->i + 1; iright+1 < n && ydata_filtered[iright] == ydata_filtered[iright-1]; iright++)
            ;

        if (peaks->background == GWY_PEAK_BACKGROUND_ZERO) {
            while (iright+1 < n && ydata_filtered[iright] > ydata_filtered[iright+1] && ydata[iright] > 0.0)
                iright++;
            if (ydata[iright] < 0.0)
                iright--;
            yright = 0.0;
        }
        else {
            while (iright+1 < n && ydata_filtered[iright] > ydata_filtered[iright+1])
                iright++;
            yright = ydata[iright];
        }

        /* Calculate height, area, etc. */
        arealeft = arearight = 0.0;
        disp2left = disp2right = 0.0;
        peak->x = xdata[peak->i];
        for (i = ileft; i < peak->i; i++) {
            gdouble xl = xdata[i] - peak->x, xr = xdata[i+1] - peak->x,
                    yl = fmax(ydata[i] - yleft, 0.0), yr = fmax(ydata[i+1] - yleft, 0.0);
            arealeft += (xr - xl)*(yl + yr)/2.0;
            disp2left += (xr - xl)*((3.0*yr + yl)*xr*xr + 2.0*(yl + yl)*xr*xl + (yr + 3.0*yl)*xl*xl)/12.0;
        }
        for (i = iright; i > peak->i; i--) {
            gdouble xl = xdata[i-1] - peak->x, xr = xdata[i] - peak->x,
                    yl = fmax(ydata[i-1] - yright, 0.0), yr = fmax(ydata[i] - yright, 0.0);
            arearight += (xr - xl)*(yl + yr)/2.0;
            disp2right += (xr - xl)*((3.0*yr + yl)*xr*xr + 2.0*(yl + yl)*xr*xl + (yr + 3.0*yl)*xl*xl)/12.0;
        }

        peak->area = arealeft + arearight;
        if (arealeft > 0.0 && arearight > 0.0)
            peak->width = sqrt(0.5*(disp2left/arealeft + disp2right/arearight));
        else if (arealeft > 0.0)
            peak->width = sqrt(disp2left/arealeft);
        else if (arearight > 0.0)
            peak->width = sqrt(disp2right/arearight);
        else
            peak->width = 0.0;

        i = peak->i;
        peak->height = ydata[i] - 0.5*(yleft + yright);
        if (ydata[i] > ydata[i-1] || ydata[i] > ydata[i+1]) {
            gdouble epsp = ydata[i] - ydata[i+1];
            gdouble epsm = ydata[i] - ydata[i-1];
            gdouble dp = xdata[i+1] - xdata[i];
            gdouble dm = xdata[i] - xdata[i-1];
            gdouble xdiff = 0.5*(epsm*dp*dp - epsp*dm*dm)/(epsm*dp + epsp*dm);
            if (peak->x + xdiff < xdata[i+1] && peak->x + xdiff > xdata[i-1])
                peak->x += xdiff;
        }
    }
    g_free(ydata_filtered);

    for (k = 0; k < p->len; ) {
        Peak *peak = &g_array_index(p, Peak, k);
        gdouble xleft = (k > 0 ? g_array_index(p, Peak, k-1).x : xdata[0]);
        gdouble xright = (k+1 < p->len ? g_array_index(p, Peak, k+1).x : xdata[n-1]);

        if (peak->height <= 0.0 || peak->area <= 0.0 || peak->x >= xright || peak->x <= xleft)
            g_array_remove_index(p, k);
        else {
            peak->prominence = log(peak->height * peak->area * (xright - peak->x) * (peak->x - xleft));
            k++;
        }
    }

    g_array_sort(p, compare_prominence_descending);
    if (p->len > maxpeaks)
        g_array_set_size(p, maxpeaks);

    if (peaks->order == GWY_PEAK_ORDER_ABSCISSA)
        g_array_sort(p, compare_abscissa_ascending);

    return p->len;
}

/**
 * gwy_peaks_n_peaks:
 * @peaks: A peak analyser.
 *
 * Gets the current number of peaks of a peak analyser.
 *
 * Returns: The currently remembered number of peaks.
 *
 * Since: 2.46
 **/
guint
gwy_peaks_n_peaks(GwyPeaks *peaks)
{
    g_return_val_if_fail(peaks, 0);
    return peaks->peaks->len;
}

/**
 * gwy_peaks_get_quantity:
 * @peaks: A peak analyser.
 * @quantity: Peak property to return.
 * @data: Array of sufficient length to hold values for all peaks (their number is returned by gwy_peaks_n_peaks()).
 *
 * Obtaines values of a given quantity for all found peaks.
 *
 * Since: 2.46
 **/
void
gwy_peaks_get_quantity(GwyPeaks *peaks,
                       GwyPeakQuantity quantity,
                       gdouble *data)
{
    GArray *p;
    guint i;

    g_return_if_fail(peaks);
    g_return_if_fail(data);
    g_return_if_fail(quantity <= GWY_PEAK_WIDTH);

    p = peaks->peaks;
    for (i = 0; i < p->len; i++) {
        const Peak *peak = &g_array_index(p, Peak, i);
        if (quantity == GWY_PEAK_PROMINENCE)
            data[i] = peak->prominence;
        else if (quantity == GWY_PEAK_ABSCISSA)
            data[i] = peak->x;
        else if (quantity == GWY_PEAK_HEIGHT)
            data[i] = peak->height;
        else if (quantity == GWY_PEAK_AREA)
            data[i] = peak->area;
        else if (quantity == GWY_PEAK_WIDTH)
            data[i] = peak->width;
    }
}

static gint
compare_prominence_descending(gconstpointer a, gconstpointer b)
{
    const gdouble pa = ((const Peak*)a)->prominence;
    const gdouble pb = ((const Peak*)b)->prominence;

    if (pa > pb)
        return -1;
    if (pa < pb)
        return 1;
    return 0;
}

static gint
compare_abscissa_ascending(gconstpointer a, gconstpointer b)
{
    const gdouble xa = ((const Peak*)a)->x;
    const gdouble xb = ((const Peak*)b)->x;

    if (xb > xa)
        return -1;
    if (xb < xa)
        return 1;
    return 0;
}

/************************** Documentation ****************************/

/**
 * SECTION:peaks
 * @title: GwyPeaks
 * @short_description: Graph peak analyser
 *
 * One #GwyPeaks analyser can be used repeatedly to find and characterise peaks in several data sets.
 *
 * Feed the curve data to the analyser using functions such as gwy_peaks_analyze(), gwy_peaks_analyze_dataline() or
 * gwy_peaks_analyze_xy(). It will process the data, locate peaks and remember their various characteristics that you
 * can subsequently obtain using gwy_peaks_get_quantity().  The results are static; if you change the analyser
 * settings (or the curve data) you need to re-run the analysis function.
 **/

/**
 * GwyPeaks:
 *
 * #GwyPeaks is an opaque data structure and should be only manipulated
 * with the functions below.
 *
 * Since: 2.46
 **/

/**
 * GwyPeakBackgroundType:
 * @GWY_PEAK_BACKGROUND_ZERO: The background is fixed at zero value.
 * @GWY_PEAK_BACKGROUND_MMSTEP: The background is a step function connecting the nearest minima on the left and right
 *                              side.
 *
 * Type of background available in graph peak analysers.
 *
 * Since: 2.46
 **/

/**
 * GwyPeakOrderType:
 * @GWY_PEAK_ORDER_ABSCISSA: Peaks are ordered by abscissa values, from left to right.
 * @GWY_PEAK_ORDER_PROMINENCE: Peaks are ordered by prominence, from most to least prominent.
 *
 * Type of peak ordering by in the graph peak analyser results.
 *
 * Since: 2.46
 **/

/**
 * GwyPeakQuantity:
 * @GWY_PEAK_PROMINENCE: Compound quantity characteristing the overall peak prominence (taking into account height,
 *                       area, distance from other peaks, ...).
 * @GWY_PEAK_ABSCISSA: Position of peak maximum.
 * @GWY_PEAK_HEIGHT: Peak height (with respect to the chosen background function).
 * @GWY_PEAK_AREA: Peak area (with respect to the chosen background function).
 * @GWY_PEAK_WIDTH: Peak width, more or less corresponding to standard deviation.
 *
 * Type of characteristics graph peak analysers can provide.
 *
 * Since: 2.46
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
