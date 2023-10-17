/*
 *  $Id: stats-entropy.c 25308 2023-04-21 13:16:28Z yeti-dn $
 *  Copyright (C) 2003-2017 David Necas (Yeti), Petr Klapetek.
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
#include <libprocess/datafield.h>
#include <libprocess/stats.h>
#include "gwyprocessinternal.h"

typedef struct _BinTreeNode BinTreeNode;
typedef struct _QuadTreeNode QuadTreeNode;

struct _BinTreeNode {
    /* This optimally uses memory on 64bit architectures where pt and children have the same size (16 bytes). */
    union {
        /* The at most two points inside for non-max-depth leaves. */
        struct {
            gdouble a;
            gdouble b;
        } pt;
        /* Children for non-max-depth non-leaves. */
        BinTreeNode *children[2];
    } u;
    /* Always set; for max-depth leaves it is the only meaningful field. */
    guint count;
};

typedef struct {
    gdouble min;
    gdouble max;
    BinTreeNode *root;
    guint maxdepth;
    gboolean degenerate;
    gdouble degenerateS;
} BinTree;

struct _QuadTreeNode {
    /* This optimally uses memory on 64bit architectures where pt and children have the same size (32 bytes). */
    union {
        /* The at most two points inside for non-max-depth leaves. */
        struct {
            GwyXY a;
            GwyXY b;
        } pt;
        /* Children for non-max-depth non-leaves. */
        QuadTreeNode *children[4];
    } u;
    /* Always set; for max-depth leaves it is the only meaningful field. */
    guint count;
};

typedef struct {
    GwyXY min;
    GwyXY max;
    QuadTreeNode *root;
    guint maxdepth;
    gboolean degenerate;
    gdouble degenerateS;
} QuadTree;

/* Find the flattest part of the curve representing scaling histogram-based entropy on scale and use the value there
 * as the entropy estimate.  Handle the too-few-pixels cases gracefully.
 *
 * NB: We assume
 * (1) ecurve beings from large scales.  This is important only when it has lots of points because we may skip a few
 * at the beginning then to avoid mistaking the flat part of the curve there for the inflexion point.
 * (2) ecurve goes by powers of 2 scales, this is for the mindiff filtering.
 */
static gdouble
calculate_entropy_from_scaling(const gdouble *ecurve, guint maxdiv)
{
    /* Initialise S to the δ-function entropy and mindiff to the half of the asymptotic value for distribution that is
     * sum of δ-functions. This means only if the differences drops substantially from this asymptotic value we will
     * consider is as potential inflexion point. If we get ecurve[] essentially corresponding to a set of δ-functions
     * then we return -G_MAXDOUBLE. */
    gdouble S = -G_MAXDOUBLE, mindiff = 0.6*G_LN2;
    guint i, from = (maxdiv >= 12) + (maxdiv >= 36);

    if (maxdiv < 1)
        return ecurve[0];

    if (maxdiv < 5) {
        for (i = from; i <= maxdiv-2; i++) {
            gdouble diff = 0.5*(fabs(ecurve[i+1] - ecurve[i]) + fabs(ecurve[i+2] - ecurve[i+1]))/G_LN2;
            gdouble diff2 = 0.5*(fabs(ecurve[i] + ecurve[i+2] - 2.0*ecurve[i+1]))/(G_LN2*G_LN2);
            if (diff + diff2 < mindiff) {
                S = ecurve[i+1];
                mindiff = diff + diff2;
            }
        }
    }
    else {
        for (i = from; i <= maxdiv-4; i++) {
            gdouble diff = 0.25*(fabs(ecurve[i+1] - ecurve[i]) + fabs(ecurve[i+2] - ecurve[i+1])
                                 + fabs(ecurve[i+3] - ecurve[i+2]) + fabs(ecurve[i+4] - ecurve[i+3]));
            gdouble diff2 = 0.5*(fabs(ecurve[i+1] + ecurve[i+4] - 2.0*ecurve[i+2]))/(G_LN2*G_LN2);
            if (diff + diff2 < mindiff) {
                S = (ecurve[i+1] + ecurve[i+2] + ecurve[i+3])/3.0;
                mindiff = diff + diff2;
            }
        }
    }

    return S;
}

/* This is what we get on average from all possible two-point configurations if they are randomly distributed.
 * A fairly good estimate that in practice seems to result in some deviation on the 5th significant digit, which is
 * hardly significant at all.  The contribution is the same in 1D and 2D. */
static void
add_estimated_unsplit_node_entropy(gdouble *S, guint maxdepth, gdouble w)
{
    gdouble q = 2.0*G_LN2*w;
    guint i;

    for (i = 0; i <= maxdepth; i++, S++) {
        S[0] += q;
        q *= 0.5;
    }
}

static BinTreeNode*
bin_tree_node_new(const gdouble pt)
{
    BinTreeNode *btnode = g_slice_new(BinTreeNode);
    btnode->u.pt.a = pt;
    btnode->count = 1;
    return btnode;
}

static void
bin_tree_add_node(BinTreeNode *btnode, const gdouble pt,
                  gdouble min, gdouble max, guint maxdepth)
{
    BinTreeNode *child;
    gdouble centre;
    guint i;

    /* We reached maximum allowed subdivision.  Just increase the count. */
    if (!maxdepth) {
        if (btnode->count <= 2)
            gwy_clear(&btnode->u, 1);
        btnode->count++;
        return;
    }

    /* We will descend into subtrees. */
    centre = 0.5*(min + max);

    /* If this node has just one point add the other there and we are done. */
    if (btnode->count == 1) {
        btnode->u.pt.b = pt;
        btnode->count++;
        return;
    }

    /* We will be recursing.  So if this node is a leaf start by making it non-leaf. */
    if (btnode->count == 2) {
        gdouble pta = btnode->u.pt.a;
        gdouble ptb = btnode->u.pt.b;
        guint ia = (pta > centre);
        guint ib = (ptb > centre);

        gwy_clear(&btnode->u, 1);
        child = btnode->u.children[ia] = bin_tree_node_new(pta);
        /* Must distinguish between creating two child nodes and creating one two-point child node. */
        if (ia == ib) {
            child->u.pt.b = ptb;
            child->count = 2;
        }
        else
            btnode->u.children[ib] = bin_tree_node_new(ptb);
    }

    /* Add the new point to the appropriate child. */
    i = (pt > centre);
    maxdepth--;
    btnode->count++;

    if ((child = btnode->u.children[i])) {
        /* Recurse.  This will end either by reaching maxdepth=0 or by successful separation in the other branch of
         * this conditon. */
        if (i == 0)
            bin_tree_add_node(child, pt, min, centre, maxdepth);
        else
            bin_tree_add_node(child, pt, centre, max, maxdepth);
    }
    else {
        /* There is nothing here yet.  Add the point as a new leaf. */
        btnode->u.children[i] = bin_tree_node_new(pt);
    }
}

static void
bin_tree_add(BinTree *btree, const gdouble pt)
{
    if (G_LIKELY(btree->root))
        bin_tree_add_node(btree->root, pt, btree->min, btree->max, btree->maxdepth);
    else
        btree->root = bin_tree_node_new(pt);
}

static void
bin_tree_find_range(BinTree *btree, const gdouble *xdata, guint n)
{
    gdouble min = G_MAXDOUBLE;
    gdouble max = -G_MAXDOUBLE;
    guint i;

    for (i = 0; i < n; i++) {
        gdouble x = xdata[i];

        if (x < min)
            min = x;
        if (x > max)
            max = x;
    }

    btree->min = min;
    btree->max = max;
}

static void
bin_tree_node_free(BinTreeNode *btnode)
{
    guint i;

    if (btnode->count > 2) {
        for (i = 0; i < G_N_ELEMENTS(btnode->u.children); i++) {
            if (btnode->u.children[i])
                bin_tree_node_free(btnode->u.children[i]);
        }
    }
    g_slice_free(BinTreeNode, btnode);
}

static void
bin_tree_free(BinTree *btree)
{
    if (!btree->degenerate)
        bin_tree_node_free(btree->root);
    g_free(btree);
}

static BinTree*
bin_tree_new(const gdouble *xdata, guint n, guint maxdepth)
{
    BinTree *btree;
    guint i;

    btree = g_new0(BinTree, 1);

    if (!maxdepth)
        maxdepth = 24;
    btree->maxdepth = maxdepth;

    bin_tree_find_range(btree, xdata, n);
    if (!(btree->min < btree->max)) {
        btree->degenerate = TRUE;
        btree->degenerateS = G_MAXDOUBLE;
        return btree;
    }

    /* Return explicit estimates for n < 4, making maxdiv at least 1 (with half-scales included, ecurve will have at
     * least 3 points then). */
    if (n == 2) {
        btree->degenerate = TRUE;
        btree->degenerateS = log(btree->max - btree->min);
        return btree;
    }
    if (n == 3) {
        btree->degenerate = TRUE;
        btree->degenerateS = (log(btree->max - btree->min) + 0.5*log(1.5) - G_LN2/3.0);
        return btree;
    }

    for (i = 0; i < n; i++) {
        gdouble pt = xdata[i];
        bin_tree_add(btree, pt);
    }

    return btree;
}

static void
bin_tree_node_entropies_at_scales(BinTreeNode *btnode, guint maxdepth,
                                  gdouble *S, guint *unsplit)
{
    BinTreeNode *child;
    guint i;

    /* Singular points contribute to p*ln(p) always with zero.  So we can stop recursion to finer subdivisions when
     * count == 1. */
    if (btnode->count <= 1)
        return;

    if (!maxdepth) {
        S[0] += gwy_xlnx_int(btnode->count);
        return;
    }

    if (btnode->count == 2) {
        unsplit[0]++;
        return;
    }

    S[0] += gwy_xlnx_int(btnode->count);
    S++;

    maxdepth--;
    unsplit++;
    for (i = 0; i < G_N_ELEMENTS(btnode->u.children); i++) {
        if ((child = btnode->u.children[i]))
            bin_tree_node_entropies_at_scales(child, maxdepth, S, unsplit);
    }
}

static gdouble*
bin_tree_entropies_at_scales(BinTree *btree, guint maxdepth)
{
    gdouble *S;
    guint *unsplit;
    guint i, n, npts;
    gdouble Sscale;

    if (!maxdepth)
        maxdepth = btree->maxdepth;

    n = maxdepth + 1;
    S = g_new0(gdouble, n);

    if (btree->degenerate) {
        S[0] = btree->degenerateS;
        for (i = 1; i < n; i++)
            S[i] = S[i-1] - G_LN2;
        return S;
    }

    unsplit = g_new0(guint, maxdepth);
    bin_tree_node_entropies_at_scales(btree->root, MIN(maxdepth, btree->maxdepth), S, unsplit);

    for (i = 0; i < maxdepth; i++) {
        if (unsplit[i])
            add_estimated_unsplit_node_entropy(S + i, maxdepth - i, unsplit[i]);
    }
    g_free(unsplit);

    npts = btree->root->count;
    Sscale = log(npts*(btree->max - btree->min));
    for (i = 0; i < n; i++)
        S[i] = Sscale - i*G_LN2 - S[i]/npts;

    return S;
}

static gdouble*
calculate_entropy_at_scales(GwyDataField *dfield,
                            GwyDataField *mask,
                            GwyMaskingType mode,
                            gint col, gint row,
                            gint width, gint height,
                            guint *maxdiv,
                            gdouble *S)
{
    gint xres;
    guint i, j, n;
    gdouble *xdata;
    const gdouble *base;
    gboolean must_free_xdata = TRUE;
    gdouble *ecurve;
    BinTree *btree;

    if (mask) {
        gwy_data_field_area_count_in_range(mask, NULL, col, row, width, height, G_MAXDOUBLE, 1.0, NULL, &n);
        if (mode == GWY_MASK_EXCLUDE)
            n = width*height - n;
    }
    else
        n = width*height;

    if (!*maxdiv) {
        if (n >= 2)
            *maxdiv = (guint)floor(3.0*log(n)/G_LN2 + 1e-12);
        else
            *maxdiv = 2;

        /* We will run out of significant digits in coordinates after that. */
        *maxdiv = MIN(*maxdiv, 50);
    }

    if (n < 2) {
        ecurve = g_new(gdouble, *maxdiv+1);
        for (i = 0; i <= *maxdiv; i++)
            ecurve[i] = -G_MAXDOUBLE;
        if (S)
            *S = -G_MAXDOUBLE;
        return ecurve;
    }

    xres = dfield->xres;
    base = dfield->data + row*xres + col;
    if (n == xres*dfield->yres) {
        /* Handle the full-field case without allocating anything. */
        xdata = dfield->data;
        must_free_xdata = FALSE;
    }
    else {
        xdata = g_new(gdouble, n);
        if (mask) {
            const gdouble *mbase = mask->data + row*xres + col;
            const gboolean invert = (mode == GWY_MASK_EXCLUDE);
            guint k = 0;

            for (i = 0; i < height; i++) {
                const gdouble *d = base + i*xres;
                const gdouble *m = mbase + i*xres;
                for (j = width; j; j--, d++, m++) {
                    if ((*m < 1.0) == invert)
                        xdata[k++] = *d;
                }
            }
            g_assert(k == n);
        }
        else {
            for (i = 0; i < height; i++)
                gwy_assign(xdata + i*width, base + i*xres, width);
        }
    }

    /* FIXME: How can we parallelise this mess?
     * Trees covering overlapping (or the same) intervals are lots of work to merge.  Processing them separately is
     * possible, but then we hit different number of add_estimated_unsplit_node_entropy() depending on the number of
     * threads.
     *
     * So one option is add a gwy_math_kth_ranks() preprocessing.  It finds the values that split the data to equal
     * chunks.  And not only that, it also physically splits them, so each chunk becomes a continuous block in the
     * array.  Neat!
     *
     * So each thread can construct its own, independent tree.  However, the result still depends on the number of
     * threads because the bin boundaries are not the same as for serial processing.
     *
     * For deterministic results each thread must essentially construct exactly one branch of the tree.  We must
     * assign an interval of the input of length (max-min)/2^n for some n and also at position corresponding to
     * a branch.  Then we just merge the branches to form a common root.
     *
     * Splitting the work reasonably requires a fair amount of initial analysis.   The result could be for instance
     * that if we have 6 threads they should process 1/4, 1/4, 1/16, 1/16, 1/8 and 1/4 of the entire interval of
     * values.
     *
     * For a key commonly used function it may be woth implementing this. Here, it is hardly a priority. */
    btree = bin_tree_new(xdata, n, *maxdiv);
    if (must_free_xdata)
        g_free(xdata);

    ecurve = bin_tree_entropies_at_scales(btree, *maxdiv);
    if (S) {
        if (btree->degenerate)
            *S = btree->degenerateS;
        else
            *S = calculate_entropy_from_scaling(ecurve, *maxdiv);
    }
    bin_tree_free(btree);

    return ecurve;
}

/**
 * gwy_data_field_area_get_entropy_at_scales:
 * @data_field: A data field.
 * @target_line: A data line to store the result to.  It will be resampled to @maxdiv+1 items.
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @mode: Masking mode to use.  See the introduction for description of masking modes.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @maxdiv: Maximum number of divisions of the value range.  Pass zero to choose it automatically.
 *
 * Calculates estimates of value distribution entropy at various scales.
 *
 * Returns: The best estimate, as gwy_data_field_area_get_entropy().
 *
 * Since: 2.44
 **/
gdouble
gwy_data_field_area_get_entropy_at_scales(GwyDataField *data_field,
                                          GwyDataLine *target_line,
                                          GwyDataField *mask,
                                          GwyMaskingType mode,
                                          gint col, gint row,
                                          gint width, gint height,
                                          gint maxdiv)
{
    GwySIUnit *lineunit;
    guint umaxdiv = (maxdiv > 0 ? maxdiv : 0);
    gdouble *ecurve;
    gdouble min, max, S = -G_MAXDOUBLE;
    gint i;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height)
        || !_gwy_data_field_check_mask(data_field, &mask, &mode))
        return S;
    g_return_val_if_fail(GWY_IS_DATA_LINE(target_line), S);

    ecurve = calculate_entropy_at_scales(data_field, mask, mode, col, row, width, height, &umaxdiv, &S);
    maxdiv = maxdiv ? maxdiv : umaxdiv + 1;
    gwy_data_line_resample(target_line, maxdiv, GWY_INTERPOLATION_NONE);
    target_line->real = maxdiv*G_LN2;
    for (i = 0; i < maxdiv; i++)
        target_line->data[maxdiv-1 - i] = ecurve[i];
    g_free(ecurve);

    gwy_data_field_area_get_min_max_mask(data_field, mask, mode, col, row, width, height, &min, &max);
    if (max > min)
        target_line->off = log(max - min) - (maxdiv - 0.5)*G_LN2;

    lineunit = gwy_data_line_get_si_unit_x(target_line);
    gwy_si_unit_set_from_string(lineunit, NULL);
    lineunit = gwy_data_line_get_si_unit_y(target_line);
    gwy_si_unit_set_from_string(lineunit, NULL);

    return S;
}

/**
 * gwy_data_field_get_entropy:
 * @data_field: A data field.
 *
 * Computes the entropy of a data field.
 *
 * See gwy_data_field_area_get_entropy() for the definition.
 *
 * This quantity is cached.
 *
 * Returns: The value distribution entropy.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_field_get_entropy(GwyDataField *data_field)
{
    gdouble S = -G_MAXDOUBLE;
    gdouble *ecurve;
    guint maxdiv = 0;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), S);

    gwy_debug("%s", CTEST(data_field, ENT) ? "cache" : "lame");
    if (CTEST(data_field, ENT))
        return CVAL(data_field, ENT);

    ecurve = calculate_entropy_at_scales(data_field, NULL, GWY_MASK_IGNORE, 0, 0, data_field->xres, data_field->yres,
                                         &maxdiv, &S);
    g_free(ecurve);

    CVAL(data_field, ENT) = S;
    data_field->cached |= CBIT(ENT);

    return S;
}

/**
 * gwy_data_field_area_get_entropy:
 * @data_field: A data field.
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @mode: Masking mode to use.  See the introduction for description of masking modes.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Estimates the entropy of field data distribution.
 *
 * The estimate is calculated as @S = ln(@n Δ) − 1/@n ∑ @n_i ln(@n_i), where @n is the number of pixels considered,
 * Δ the bin size and @n_i the count in the @i-th bin.  If @S is plotted as a function of the bin size Δ, it is,
 * generally, a growing function with a plateau for ‘reasonable’ bin sizes. The estimate is taken at the plateau.  If
 * no plateau is found, which means the distribution is effectively a sum of δ-functions, -%G_MAXDOUBLE is returned.
 *
 * It should be noted that this estimate may be biased.
 *
 * Returns: The estimated entropy of the data values.  The entropy of no data or a single single is returned as
 *          -%G_MAXDOUBLE.
 *
 * Since: 2.42
 **/
gdouble
gwy_data_field_area_get_entropy(GwyDataField *data_field,
                                GwyDataField *mask,
                                GwyMaskingType mode,
                                gint col, gint row,
                                gint width, gint height)
{
    gdouble S = -G_MAXDOUBLE;
    gdouble *ecurve;
    guint maxdiv = 0;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height)
        || !_gwy_data_field_check_mask(data_field, &mask, &mode))
        return S;

    /* The result is the same, but it can be cached. */
    if (!mask && row == 0 && col == 0 && width == data_field->xres && height == data_field->yres)
        return gwy_data_field_get_entropy(data_field);

    ecurve = calculate_entropy_at_scales(data_field, mask, mode, col, row, width, height, &maxdiv, &S);
    g_free(ecurve);
    return S;
}

static QuadTreeNode*
quad_tree_node_new(const GwyXY *pt)
{
    QuadTreeNode *qtnode = g_slice_new(QuadTreeNode);
    qtnode->u.pt.a = *pt;
    qtnode->count = 1;
    return qtnode;
}

static void
quad_tree_add_node(QuadTreeNode *qtnode, const GwyXY *pt,
                   GwyXY min, GwyXY max, guint maxdepth)
{
    QuadTreeNode *child;
    GwyXY centre;
    guint i;

    /* We reached maximum allowed subdivision.  Just increase the count. */
    if (!maxdepth) {
        if (qtnode->count <= 2)
            gwy_clear(&qtnode->u, 1);
        qtnode->count++;
        return;
    }

    /* We will descend into subtrees. */
    centre.x = 0.5*(min.x + max.x);
    centre.y = 0.5*(min.y + max.y);

    /* If this node has just one point add the other there and we are done. */
    if (qtnode->count == 1) {
        qtnode->u.pt.b = *pt;
        qtnode->count++;
        return;
    }

    /* We will be recursing.  So if this node is a leaf start by making it
     * non-leaf. */
    if (qtnode->count == 2) {
        GwyXY pta = qtnode->u.pt.a;
        GwyXY ptb = qtnode->u.pt.b;
        guint ia = (pta.x > centre.x) + 2*(pta.y > centre.y);
        guint ib = (ptb.x > centre.x) + 2*(ptb.y > centre.y);

        gwy_clear(&qtnode->u, 1);
        child = qtnode->u.children[ia] = quad_tree_node_new(&pta);
        /* Must distinguish between creating two child nodes and creating one
         * two-point child node. */
        if (ia == ib) {
            child->u.pt.b = ptb;
            child->count = 2;
        }
        else
            qtnode->u.children[ib] = quad_tree_node_new(&ptb);
    }

    /* Add the new point to the appropriate child. */
    i = (pt->x > centre.x) + 2*(pt->y > centre.y);
    maxdepth--;
    qtnode->count++;

    if ((child = qtnode->u.children[i])) {
        /* Recurse.  This will end either by reaching maxdepth=0 or by successful separation in the other branch of
         * this conditon. */
        if (i == 0)
            quad_tree_add_node(child, pt, min, centre, maxdepth);
        else if (i == 1) {
            min.x = centre.x;
            max.y = centre.y;
            quad_tree_add_node(child, pt, min, max, maxdepth);
        }
        else if (i == 2) {
            max.x = centre.x;
            min.y = centre.y;
            quad_tree_add_node(child, pt, min, max, maxdepth);
        }
        else
            quad_tree_add_node(child, pt, centre, max, maxdepth);
    }
    else {
        /* There is nothing here yet.  Add the point as a new leaf. */
        qtnode->u.children[i] = quad_tree_node_new(pt);
    }
}

static void
quad_tree_add(QuadTree *qtree, const GwyXY *pt)
{
    if (G_LIKELY(qtree->root))
        quad_tree_add_node(qtree->root, pt, qtree->min, qtree->max, qtree->maxdepth);
    else
        qtree->root = quad_tree_node_new(pt);
}

static void
quad_tree_find_range(QuadTree *qtree,
                     const gdouble *xdata, const gdouble *ydata, guint n)
{
    GwyXY min = { G_MAXDOUBLE, G_MAXDOUBLE };
    GwyXY max = { -G_MAXDOUBLE, -G_MAXDOUBLE };
    guint i;

    for (i = 0; i < n; i++) {
        gdouble x = xdata[i];
        gdouble y = ydata[i];

        if (x < min.x)
            min.x = x;
        if (x > max.x)
            max.x = x;
        if (y < min.y)
            min.y = y;
        if (y > max.y)
            max.y = y;
    }

    qtree->min = min;
    qtree->max = max;
}

static void
quad_tree_node_free(QuadTreeNode *qtnode)
{
    guint i;

    if (qtnode->count > 2) {
        for (i = 0; i < G_N_ELEMENTS(qtnode->u.children); i++) {
            if (qtnode->u.children[i])
                quad_tree_node_free(qtnode->u.children[i]);
        }
    }
    g_slice_free(QuadTreeNode, qtnode);
}

static void
quad_tree_free(QuadTree *qtree)
{
    quad_tree_node_free(qtree->root);
    g_free(qtree);
}

static QuadTree*
quad_tree_new(const gdouble *xdata, const gdouble *ydata, guint n,
              guint maxdepth)
{
    QuadTree *qtree;
    guint i;

    qtree = g_new0(QuadTree, 1);

    if (!maxdepth)
        maxdepth = 16;
    qtree->maxdepth = maxdepth;

    quad_tree_find_range(qtree, xdata, ydata, n);
    if (!(qtree->min.x < qtree->max.x) || !(qtree->min.y < qtree->max.y)) {
        qtree->degenerate = TRUE;
        qtree->degenerateS = G_MAXDOUBLE;
        return qtree;
    }

    /* Return explicit estimates for n < 4, making maxdiv at least 1 (with half-scales included, ecurve will have at
     * least 3 points then). */
    if (n == 2) {
        qtree->degenerate = TRUE;
        qtree->degenerateS = (log(qtree->max.x - qtree->min.x) + log(qtree->max.y - qtree->min.y));
        return qtree;
    }
    if (n == 3) {
        qtree->degenerate = TRUE;
        qtree->degenerateS = (log(qtree->max.x - qtree->min.x) + log(qtree->max.y - qtree->min.y)
                              + 0.5*log(1.5) - 2.0*G_LN2/3.0);
        return qtree;
    }

    for (i = 0; i < n; i++) {
        GwyXY pt = { xdata[i], ydata[i] };
        quad_tree_add(qtree, &pt);
    }

    return qtree;
}

static gdouble
quad_tree_node_half_scale_entropy(QuadTreeNode *qtnode)
{
    QuadTreeNode *child;
    guint cnt[G_N_ELEMENTS(qtnode->u.children)] = { 0, 0, 0, 0 };
    guint i;

    for (i = 0; i < G_N_ELEMENTS(qtnode->u.children); i++) {
        if ((child = qtnode->u.children[i]))
            cnt[i] = child->count;
    }
    return 0.5*(gwy_xlnx_int(cnt[0] + cnt[1]) + gwy_xlnx_int(cnt[2] + cnt[3])
                + gwy_xlnx_int(cnt[0] + cnt[2]) + gwy_xlnx_int(cnt[1] + cnt[3]));
}

static void
quad_tree_node_entropies_at_scales(QuadTreeNode *qtnode, guint maxdepth,
                                   gdouble *S, guint *unsplit)
{
    QuadTreeNode *child;
    guint i;

    /* Singular points contribute to p*ln(p) always with zero.  So we can stop recursion to finer subdivisions when
     * count == 1. */
    if (qtnode->count <= 1)
        return;

    if (!maxdepth) {
        S[0] += gwy_xlnx_int(qtnode->count);
        return;
    }

    if (qtnode->count == 2) {
        unsplit[0]++;
        return;
    }

    S[0] += gwy_xlnx_int(qtnode->count);
    S++;

    /* Half-scale entropies we estimate as averages of horizontal and vertical binning. */
    S[0] += quad_tree_node_half_scale_entropy(qtnode);
    S++;

    maxdepth--;
    unsplit++;
    for (i = 0; i < G_N_ELEMENTS(qtnode->u.children); i++) {
        if ((child = qtnode->u.children[i]))
            quad_tree_node_entropies_at_scales(child, maxdepth, S, unsplit);
    }
}

static gdouble*
quad_tree_entropies_at_scales(QuadTree *qtree, guint maxdepth)
{
    gdouble *S;
    guint *unsplit;
    guint i, n, npts;
    gdouble Sscale;

    if (!maxdepth)
        maxdepth = qtree->maxdepth;

    n = 2*maxdepth + 1;
    S = g_new0(gdouble, n);
    unsplit = g_new0(guint, maxdepth);
    quad_tree_node_entropies_at_scales(qtree->root, MIN(maxdepth, qtree->maxdepth), S, unsplit);

    for (i = 0; i < maxdepth; i++) {
        if (unsplit[i])
            add_estimated_unsplit_node_entropy(S + 2*i, 2*(maxdepth - i), unsplit[i]);
    }
    g_free(unsplit);

    npts = qtree->root->count;
    Sscale = log(npts*(qtree->max.x - qtree->min.x)*(qtree->max.y - qtree->min.y));
    for (i = 0; i < n; i++)
        S[i] = Sscale - i*G_LN2 - S[i]/npts;

    return S;
}

static gdouble*
calculate_entropy_2d_at_scales(GwyDataField *xfield,
                               GwyDataField *yfield,
                               guint *maxdiv,
                               gdouble *S)
{
    guint xres, yres, n, i;
    gdouble *ecurve;
    QuadTree *qtree;

    xres = xfield->xres;
    yres = xfield->yres;
    n = xres*yres;

    if (!*maxdiv) {
        if (n >= 2)
            *maxdiv = (guint)floor(1.5*log(n)/G_LN2 + 1e-12);
        else
            *maxdiv = 1;

        /* We will run out of significant digits in coordinates after that. */
        *maxdiv = MIN(*maxdiv, 50);
    }

    if (n < 2) {
        ecurve = g_new(gdouble, *maxdiv+1);
        for (i = 0; i <= *maxdiv; i++)
            ecurve[i] = -G_MAXDOUBLE;
        if (S)
            *S = -G_MAXDOUBLE;
        return ecurve;
    }

    qtree = quad_tree_new(xfield->data, yfield->data, n, *maxdiv);
    ecurve = quad_tree_entropies_at_scales(qtree, *maxdiv);
    if (S) {
        if (qtree->degenerate)
            *S = qtree->degenerateS;
        else
            *S = calculate_entropy_from_scaling(ecurve, 2*(*maxdiv));
    }
    quad_tree_free(qtree);

    return ecurve;
}

/**
 * gwy_data_field_get_entropy_2d_at_scales:
 * @xfield: A data field containing the @x-coordinates.
 * @yfield: A data field containing the @y-coordinates.
 * @target_line: A data line to store the result to.  It will be resampled to @maxdiv+1 items.
 * @maxdiv: Maximum number of divisions of the value range.  Pass zero to choose it automatically.
 *
 * Calculates estimates of entropy of two-dimensional point cloud at various scales.
 *
 * Returns: The best estimate, as gwy_data_field_get_entropy_2d().
 *
 * Since: 2.44
 **/
gdouble
gwy_data_field_get_entropy_2d_at_scales(GwyDataField *xfield,
                                        GwyDataField *yfield,
                                        GwyDataLine *target_line,
                                        gint maxdiv)
{
    GwySIUnit *lineunit;
    guint umaxdiv = (maxdiv > 0 ? maxdiv/2 : 0);
    gdouble *ecurve;
    gdouble xmin, xmax, ymin, ymax, S = -G_MAXDOUBLE;
    gint i;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(xfield), S);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(yfield), S);
    g_return_val_if_fail(GWY_IS_DATA_LINE(target_line), S);
    g_return_val_if_fail(xfield->xres == yfield->xres, S);
    g_return_val_if_fail(xfield->yres == yfield->yres, S);

    ecurve = calculate_entropy_2d_at_scales(xfield, yfield, &umaxdiv, &S);
    maxdiv = maxdiv ? maxdiv : 2*umaxdiv + 1;
    gwy_data_line_resample(target_line, maxdiv, GWY_INTERPOLATION_NONE);
    target_line->real = maxdiv*G_LN2;
    for (i = 0; i < maxdiv; i++)
        target_line->data[maxdiv-1 - i] = ecurve[i];
    g_free(ecurve);

    gwy_data_field_get_min_max(xfield, &xmin, &xmax);
    gwy_data_field_get_min_max(xfield, &ymin, &ymax);
    if ((xmax > xmin) && (ymax > ymin))
        target_line->off = (log((xmax - xmin)*(ymax - ymin)) - (maxdiv - 0.5)*G_LN2);

    lineunit = gwy_data_line_get_si_unit_x(target_line);
    gwy_si_unit_set_from_string(lineunit, NULL);
    lineunit = gwy_data_line_get_si_unit_y(target_line);
    gwy_si_unit_set_from_string(lineunit, NULL);

    return S;
}

/**
 * gwy_data_field_get_entropy_2d:
 * @xfield: A data field containing the @x-coordinates.
 * @yfield: A data field containing the @y-coordinates.
 *
 * Computes the entropy of a two-dimensional point cloud.
 *
 * Each pair of corresponding @xfield and @yfield pixels is assumed to represent the coordinates (@x,@y) of a point in
 * plane.  Hence they must have the same dimensions.
 *
 * Returns: The two-dimensional distribution entropy.
 *
 * Since: 2.44
 **/
gdouble
gwy_data_field_get_entropy_2d(GwyDataField *xfield,
                              GwyDataField *yfield)
{
    gdouble *ecurve;
    guint maxdiv = 0;
    gdouble S = -G_MAXDOUBLE;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(xfield), S);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(yfield), S);
    g_return_val_if_fail(xfield->xres == yfield->xres, S);
    g_return_val_if_fail(xfield->yres == yfield->yres, S);

    ecurve = calculate_entropy_2d_at_scales(xfield, yfield, &maxdiv, &S);
    g_free(ecurve);

    return S;
}

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
