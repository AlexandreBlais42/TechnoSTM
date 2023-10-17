/*
 *  $Id: datafield.c 25308 2023-04-21 13:16:28Z yeti-dn $
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

#include "config.h"
#include <string.h>
#include <stdlib.h>
#include <libgwyddion/gwymacros.h>
#include <libprocess/datafield.h>
#include <libprocess/interpolation.h>
#include <libprocess/stats.h>
#include <libprocess/grains.h>
#include "libgwyddion/gwyomp.h"
#include "gwyprocessinternal.h"

#define GWY_DATA_FIELD_TYPE_NAME "GwyDataField"

enum { BLOCK_SIZE = 64 };

enum {
    DATA_CHANGED,
    LAST_SIGNAL
};

typedef struct {
    gdouble dist;
    gint i;
    gint j;
} MaskedPoint;

static void        gwy_data_field_finalize         (GObject *object);
static void        gwy_data_field_serializable_init(GwySerializableIface *iface);
static GByteArray* gwy_data_field_serialize        (GObject *obj,
                                                    GByteArray *buffer);
static gsize       gwy_data_field_get_size         (GObject *obj);
static GObject*    gwy_data_field_deserialize      (const guchar *buffer,
                                                    gsize size,
                                                    gsize *position);
static GObject*    gwy_data_field_duplicate_real   (GObject *object);
static void        gwy_data_field_clone_real       (GObject *source,
                                                    GObject *copy);
static void        set_cache_for_constant_field    (GwyDataField *data_field,
                                                    gdouble value);
static gboolean    data_field_is_constant          (GwyDataField *dfield,
                                                    gdouble *z);

static guint data_field_signals[LAST_SIGNAL] = { 0 };

G_DEFINE_TYPE_EXTENDED(GwyDataField, gwy_data_field, G_TYPE_OBJECT, 0,
                       GWY_IMPLEMENT_SERIALIZABLE(gwy_data_field_serializable_init))

static void
gwy_data_field_serializable_init(GwySerializableIface *iface)
{
    iface->serialize = gwy_data_field_serialize;
    iface->deserialize = gwy_data_field_deserialize;
    iface->get_size = gwy_data_field_get_size;
    iface->duplicate = gwy_data_field_duplicate_real;
    iface->clone = gwy_data_field_clone_real;
}

static void
gwy_data_field_class_init(GwyDataFieldClass *klass)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(klass);

    gobject_class->finalize = gwy_data_field_finalize;

    /**
    * GwyDataField::data-changed:
    * @gwydatafield: The #GwyDataField which received the signal.
    *
    * The ::data-changed signal is never emitted by data field itself.  It is intended as a means to notify others
    * data field users they should update themselves.
    */
    data_field_signals[DATA_CHANGED] = g_signal_new("data-changed",
                                                    G_OBJECT_CLASS_TYPE(gobject_class),
                                                    G_SIGNAL_RUN_FIRST,
                                                    G_STRUCT_OFFSET(GwyDataFieldClass, data_changed),
                                                    NULL, NULL,
                                                    g_cclosure_marshal_VOID__VOID,
                                                    G_TYPE_NONE, 0);
}

static void
gwy_data_field_init(G_GNUC_UNUSED GwyDataField *data_field)
{
}

static void
gwy_data_field_finalize(GObject *object)
{
    GwyDataField *data_field = (GwyDataField*)object;

    GWY_OBJECT_UNREF(data_field->si_unit_xy);
    GWY_OBJECT_UNREF(data_field->si_unit_z);
    g_free(data_field->data);

    G_OBJECT_CLASS(gwy_data_field_parent_class)->finalize(object);
}

/**
 * gwy_data_field_new:
 * @xres: X-resolution, i.e., the number of columns.
 * @yres: Y-resolution, i.e., the number of rows.
 * @xreal: Real horizontal physical dimension.
 * @yreal: Real vertical physical dimension.
 * @nullme: Whether the data field should be initialized to zeroes. If %FALSE, the data will not be initialized.
 *
 * Creates a new data field.
 *
 * Returns: A newly created data field.
 **/
GwyDataField*
gwy_data_field_new(gint xres, gint yres,
                   gdouble xreal, gdouble yreal,
                   gboolean nullme)
{
    GwyDataField *data_field;

    data_field = g_object_new(GWY_TYPE_DATA_FIELD, NULL);

    data_field->xreal = xreal;
    data_field->yreal = yreal;
    data_field->xres = xres;
    data_field->yres = yres;
    if (nullme) {
        data_field->data = g_new0(gdouble, data_field->xres*data_field->yres);
        set_cache_for_constant_field(data_field, 0.0);
    }
    else
        data_field->data = g_new(gdouble, data_field->xres*data_field->yres);

    return data_field;
}

/**
 * gwy_data_field_new_alike:
 * @model: A data field to take resolutions and units from.
 * @nullme: Whether the data field should be initialized to zeroes. If %FALSE, the data will not be initialized.
 *
 * Creates a new data field similar to an existing one.
 *
 * Use gwy_data_field_duplicate() if you want to copy a data field including data.
 *
 * Returns: A newly created data field.
 **/
GwyDataField*
gwy_data_field_new_alike(GwyDataField *model,
                         gboolean nullme)
{
    GwyDataField *data_field;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(model), NULL);
    data_field = g_object_new(GWY_TYPE_DATA_FIELD, NULL);

    data_field->xreal = model->xreal;
    data_field->yreal = model->yreal;
    data_field->xres = model->xres;
    data_field->yres = model->yres;
    data_field->xoff = model->xoff;
    data_field->yoff = model->yoff;
    if (nullme) {
        data_field->data = g_new0(gdouble, data_field->xres*data_field->yres);
        set_cache_for_constant_field(data_field, 0.0);
    }
    else
        data_field->data = g_new(gdouble, data_field->xres*data_field->yres);

    gwy_data_field_copy_units(model, data_field);

    return data_field;
}

/**
 * gwy_data_field_new_resampled:
 * @data_field: A data field.
 * @xres: Desired X resolution.
 * @yres: Desired Y resolution.
 * @interpolation: Interpolation method to use.
 *
 * Creates a new data field by resampling an existing one.
 *
 * This method is equivalent to gwy_data_field_duplicate() followed by gwy_data_field_resample(), but it is more
 * efficient.
 *
 * Returns: A newly created data field.
 **/
GwyDataField*
gwy_data_field_new_resampled(GwyDataField *data_field,
                             gint xres, gint yres,
                             GwyInterpolationType interpolation)
{
    GwyDataField *result;
    gdouble z;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    if (data_field->xres == xres && data_field->yres == yres)
        return gwy_data_field_duplicate(data_field);

    g_return_val_if_fail(xres > 0 && yres > 0, NULL);

    result = gwy_data_field_new(xres, yres, data_field->xreal, data_field->yreal, FALSE);
    result->xoff = data_field->xoff;
    result->yoff = data_field->yoff;
    gwy_data_field_copy_units(data_field, result);

    /* Prevent rounding errors from introducing different values in constants
     * field during resampling. */
    if (data_field_is_constant(data_field, &z)) {
        gwy_data_field_fill(result, z);
        return result;
    }

    gwy_interpolation_resample_block_2d(data_field->xres, data_field->yres,
                                        data_field->xres, data_field->data,
                                        result->xres, result->yres,
                                        result->xres, result->data,
                                        interpolation, TRUE);

    return result;
}

static GByteArray*
gwy_data_field_serialize(GObject *obj,
                         GByteArray *buffer)
{
    GwyDataField *data_field;
    guint32 datasize;
    gpointer pxoff, pyoff, pxyunit, pzunit;

    data_field = GWY_DATA_FIELD(obj);
    datasize = data_field->xres*data_field->yres;
    pxoff = data_field->xoff ? &data_field->xoff : NULL;
    pyoff = data_field->yoff ? &data_field->yoff : NULL;
    pxyunit = unit_pointer_if_nonempty(data_field->si_unit_xy);
    pzunit = unit_pointer_if_nonempty(data_field->si_unit_z);
    {
        GwySerializeSpec spec[] = {
            { 'i', "xres", &data_field->xres, NULL, },
            { 'i', "yres", &data_field->yres, NULL, },
            { 'd', "xreal", &data_field->xreal, NULL, },
            { 'd', "yreal", &data_field->yreal, NULL, },
            { 'd', "xoff", pxoff, NULL, },
            { 'd', "yoff", pyoff, NULL, },
            { 'o', "si_unit_xy", pxyunit, NULL, },
            { 'o', "si_unit_z", pzunit, NULL, },
            { 'D', "data", &data_field->data, &datasize, },
        };
        return gwy_serialize_pack_object_struct(buffer, GWY_DATA_FIELD_TYPE_NAME, G_N_ELEMENTS(spec), spec);
    }
}

static gsize
gwy_data_field_get_size(GObject *obj)
{
    GwyDataField *data_field;
    guint32 datasize;
    gpointer pxoff, pyoff, pxyunit, pzunit;

    data_field = GWY_DATA_FIELD(obj);
    datasize = data_field->xres*data_field->yres;
    pxoff = data_field->xoff ? &data_field->xoff : NULL;
    pyoff = data_field->yoff ? &data_field->yoff : NULL;
    pxyunit = unit_pointer_if_nonempty(data_field->si_unit_xy);
    pzunit = unit_pointer_if_nonempty(data_field->si_unit_z);
    {
        GwySerializeSpec spec[] = {
            { 'i', "xres", &data_field->xres, NULL, },
            { 'i', "yres", &data_field->yres, NULL, },
            { 'd', "xreal", &data_field->xreal, NULL, },
            { 'd', "yreal", &data_field->yreal, NULL, },
            { 'd', "xoff", pxoff, NULL, },
            { 'd', "yoff", pyoff, NULL, },
            { 'o', "si_unit_xy", pxyunit, NULL, },
            { 'o', "si_unit_z", pzunit, NULL, },
            { 'D', "data", &data_field->data, &datasize, },
        };
        return gwy_serialize_get_struct_size(GWY_DATA_FIELD_TYPE_NAME, G_N_ELEMENTS(spec), spec);
    }
}

static GObject*
gwy_data_field_deserialize(const guchar *buffer,
                           gsize size,
                           gsize *position)
{
    guint32 datasize = 0;
    gint xres, yres;
    gdouble xreal, yreal, xoff = 0.0, yoff = 0.0, *data = NULL;
    GwySIUnit *si_unit_xy = NULL, *si_unit_z = NULL;
    GwyDataField *data_field;
    GwySerializeSpec spec[] = {
        { 'i', "xres", &xres, NULL, },
        { 'i', "yres", &yres, NULL, },
        { 'd', "xreal", &xreal, NULL, },
        { 'd', "yreal", &yreal, NULL, },
        { 'd', "xoff", &xoff, NULL, },
        { 'd', "yoff", &yoff, NULL, },
        { 'o', "si_unit_xy", &si_unit_xy, NULL, },
        { 'o', "si_unit_z", &si_unit_z, NULL, },
        { 'D', "data", &data, &datasize, },
        /* Ignored */
        { 'i', "cache_bits", NULL, NULL, },
        { 'D', "cache_data", NULL, NULL, },
    };

    gwy_debug("");

    g_return_val_if_fail(buffer, NULL);

    if (!gwy_serialize_unpack_object_struct(buffer, size, position, GWY_DATA_FIELD_TYPE_NAME,
                                            G_N_ELEMENTS(spec), spec)) {
        g_free(data);
        GWY_OBJECT_UNREF(si_unit_xy);
        GWY_OBJECT_UNREF(si_unit_z);
        return NULL;
    }
    if (datasize != (gsize)(xres*yres)) {
        g_critical("Serialized %s size mismatch %u != %u", GWY_DATA_FIELD_TYPE_NAME, datasize, xres*yres);
        g_free(data);
        GWY_OBJECT_UNREF(si_unit_xy);
        GWY_OBJECT_UNREF(si_unit_z);
        return NULL;
    }

    if (xreal <= 0.0) {
        g_warning("Non-positive xreal (%g)", xreal);
        if (xreal)
            xreal = fabs(xreal);
        else
            xreal = 1.0;  /* Does not matter, just make it sane. */
    }
    if (yreal <= 0.0) {
        g_warning("Non-positive xreal (%g)", xreal);
        if (yreal)
            yreal = fabs(yreal);
        else
            yreal = 1.0;  /* Does not matter, just make it sane. */
    }

    /* don't allocate large amount of memory just to immediately free it */
    data_field = gwy_data_field_new(1, 1, xreal, yreal, FALSE);
    g_free(data_field->data);
    data_field->data = data;
    data_field->xres = xres;
    data_field->yres = yres;
    data_field->xoff = xoff;
    data_field->yoff = yoff;
    if (si_unit_z) {
        GWY_OBJECT_UNREF(data_field->si_unit_z);
        data_field->si_unit_z = si_unit_z;
    }
    if (si_unit_xy) {
        GWY_OBJECT_UNREF(data_field->si_unit_xy);
        data_field->si_unit_xy = si_unit_xy;
    }

    return (GObject*)data_field;
}

static GObject*
gwy_data_field_duplicate_real(GObject *object)
{
    GwyDataField *data_field, *duplicate;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(object), NULL);
    data_field = GWY_DATA_FIELD(object);
    duplicate = gwy_data_field_new_alike(data_field, FALSE);
    gwy_assign(duplicate->data, data_field->data, data_field->xres*data_field->yres);
    duplicate->cached = data_field->cached;
    gwy_assign(duplicate->cache, data_field->cache, GWY_DATA_FIELD_CACHE_SIZE);

    return (GObject*)duplicate;
}

static void
gwy_data_field_clone_real(GObject *source, GObject *copy)
{
    GwyDataField *data_field, *clone;
    guint n;

    g_return_if_fail(GWY_IS_DATA_FIELD(source));
    g_return_if_fail(GWY_IS_DATA_FIELD(copy));

    data_field = GWY_DATA_FIELD(source);
    clone = GWY_DATA_FIELD(copy);

    n = data_field->xres*data_field->yres;
    if (clone->xres*clone->yres != n)
        clone->data = g_renew(gdouble, clone->data, n);
    clone->xres = data_field->xres;
    clone->yres = data_field->yres;

    gwy_data_field_copy(data_field, clone, TRUE);
    clone->xoff = data_field->xoff;
    clone->yoff = data_field->yoff;
}

/**
 * gwy_data_field_data_changed:
 * @data_field: A data field.
 *
 * Emits signal "data-changed" on a data field.
 **/
void
gwy_data_field_data_changed(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_signal_emit(data_field, data_field_signals[DATA_CHANGED], 0);
}

/**
 * gwy_data_field_copy:
 * @src: Source data field.
 * @dest: Destination data field.
 * @nondata_too: Whether non-data (units) should be copied too.
 *
 * Copies the contents of an already allocated data field to a data field of the same size.
 **/
void
gwy_data_field_copy(GwyDataField *src,
                    GwyDataField *dest,
                    gboolean nondata_too)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(src));
    g_return_if_fail(GWY_IS_DATA_FIELD(dest));
    g_return_if_fail(src->xres == dest->xres && src->yres == dest->yres);

    if (src == dest)
        return;

    gwy_assign(dest->data, src->data, src->xres*src->yres);

    dest->xreal = src->xreal;
    dest->yreal = src->yreal;

    dest->cached = src->cached;
    gwy_assign(dest->cache, src->cache, GWY_DATA_FIELD_CACHE_SIZE);

    if (!nondata_too)
        return;

    gwy_data_field_copy_units(src, dest);
}

/**
 * gwy_data_field_area_copy:
 * @src: Source data field.
 * @dest: Destination data field.
 * @col: Area upper-left column coordinate in @src.
 * @row: Area upper-left row coordinate @src.
 * @width: Area width (number of columns), pass -1 for full @src widdth.
 * @height: Area height (number of rows), pass -1 for full @src height.
 * @destcol: Destination column in @dest.
 * @destrow: Destination row in @dest.
 *
 * Copies a rectangular area from one data field to another.
 *
 * The area starts at (@col, @row) in @src and its dimension is @width*@height. It is copied to @dest starting from
 * (@destcol, @destrow).
 *
 * The source area has to be completely contained in @src.  No assumptions are made about destination position,
 * however, parts of the source area sticking out the destination data field @dest are cut off.
 *
 * If @src is equal to @dest, the areas may not overlap.
 **/
void
gwy_data_field_area_copy(GwyDataField *src,
                         GwyDataField *dest,
                         gint col, gint row,
                         gint width, gint height,
                         gint destcol, gint destrow)
{
    gint i;

    g_return_if_fail(GWY_IS_DATA_FIELD(src));
    g_return_if_fail(GWY_IS_DATA_FIELD(dest));
    if (width == -1)
        width = src->xres;
    if (height == -1)
        height = src->yres;
    g_return_if_fail(col >= 0 && row >= 0 && width >= 0 && height >= 0
                     && col + width <= src->xres && row + height <= src->yres);

    if (destcol + width > dest->xres)
        width = dest->xres - destcol;
    if (destrow + height > dest->yres)
        height = dest->yres - destrow;
    if (destcol < 0) {
        col -= destcol;
        width += destcol;
        destcol = 0;
    }
    if (destrow < 0) {
        row -= destrow;
        height += destrow;
        destrow = 0;
    }
    if (width <= 0 || height <= 0)
        return;

    gwy_data_field_invalidate(dest);
    if (width == src->xres && width == dest->xres) {
        /* make it as fast as gwy_data_field_copy() whenever possible (and maybe faster, as we don't play with units */
        g_assert(col == 0 && destcol == 0);
        gwy_assign(dest->data + width*destrow, src->data + width*row, width*height);
    }
    else {
        for (i = 0; i < height; i++)
            gwy_assign(dest->data + dest->xres*(destrow + i) + destcol, src->data + src->xres*(row + i) + col, width);
    }
}

/**
 * gwy_data_field_resample:
 * @data_field: A data field to be resampled.
 * @xres: Desired X resolution.
 * @yres: Desired Y resolution.
 * @interpolation: Interpolation method to use.
 *
 * Resamples a data field using given interpolation method
 *
 * This method may invalidate raw data buffer returned by gwy_data_field_get_data().
 **/
void
gwy_data_field_resample(GwyDataField *data_field,
                        gint xres, gint yres,
                        GwyInterpolationType interpolation)
{
    gdouble *bdata;
    gdouble z;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    if (data_field->xres == xres && data_field->yres == yres)
        return;
    g_return_if_fail(xres > 0 && yres > 0);

    if (interpolation == GWY_INTERPOLATION_NONE) {
        gwy_data_field_invalidate(data_field);
        data_field->xres = xres;
        data_field->yres = yres;
        data_field->data = g_renew(gdouble, data_field->data, data_field->xres*data_field->yres);
        return;
    }

    /* Prevent rounding errors from introducing different values in constants
     * field during resampling. */
    if (data_field_is_constant(data_field, &z)) {
        data_field->xres = xres;
        data_field->yres = yres;
        data_field->data = g_renew(gdouble, data_field->data, data_field->xres*data_field->yres);
        gwy_data_field_fill(data_field, z);
        return;
    }

    gwy_data_field_invalidate(data_field);
    bdata = g_new(gdouble, xres*yres);
    gwy_interpolation_resample_block_2d(data_field->xres, data_field->yres,
                                        data_field->xres, data_field->data,
                                        xres, yres, xres, bdata,
                                        interpolation, FALSE);
    g_free(data_field->data);
    data_field->data = bdata;
    data_field->xres = xres;
    data_field->yres = yres;
}

static gboolean
data_field_is_constant(GwyDataField *data_field, gdouble *z)
{
    const gdouble *d = data_field->data;
    gint k, n = data_field->xres * data_field->yres;

    if (CTEST(data_field, MIN) && CTEST(data_field, MAX)) {
        gdouble min = CVAL(data_field, MIN);
        gdouble max = CVAL(data_field, MAX);
        if (min == max) {
            *z = min;
            set_cache_for_constant_field(data_field, min);
            return TRUE;
        }
    }

    /* This check normally either finishes fast (there are different values) or goes to the end of the field, but then
     * the caller should save a lot of time by knowing the field is constant-valued. */
    for (k = 0; k < n-1; k++) {
        if (d[k] != d[k+1])
            return FALSE;
    }

    *z = d[0];
    set_cache_for_constant_field(data_field, d[0]);
    return TRUE;
}

/**
 * gwy_data_field_resize:
 * @data_field: A data field to be resized
 * @ulcol: Upper-left column coordinate.
 * @ulrow: Upper-left row coordinate.
 * @brcol: Bottom-right column coordinate + 1.
 * @brrow: Bottom-right row coordinate + 1.
 *
 * Resizes (crops) a data field.
 *
 * Crops a data field to a rectangle between upper-left and bottom-right points, recomputing real size.
 *
 * This method may invalidate raw data buffer returned by gwy_data_field_get_data().
 **/
void
gwy_data_field_resize(GwyDataField *data_field,
                      gint ulcol, gint ulrow,
                      gint brcol, gint brrow)
{
    GwyDataField *b;
    gint i, xres, yres;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    GWY_ORDER(gint, ulcol, brcol);
    GWY_ORDER(gint, ulrow, brrow);

    g_return_if_fail(ulcol >= 0 && ulrow >= 0 && brcol <= data_field->xres && brrow <= data_field->yres);

    yres = brrow - ulrow;
    xres = brcol - ulcol;
    if (xres == data_field->xres && yres == data_field->yres)
        return;

    /* FIXME: don't allocate second field, use memmove */
    b = gwy_data_field_new(xres, yres, 1.0, 1.0, FALSE);

    for (i = ulrow; i < brrow; i++)
        gwy_assign(b->data + (i-ulrow)*xres, data_field->data + i*data_field->xres + ulcol, xres);
    data_field->xreal *= (gdouble)xres/data_field->xres;
    data_field->yreal *= (gdouble)yres/data_field->yres;
    data_field->xres = xres;
    data_field->yres = yres;
    GWY_SWAP(gdouble*, data_field->data, b->data);
    g_object_unref(b);

    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_area_extract:
 * @data_field: A data field to be resized
 * @row: Upper-left row coordinate.
 * @col: Upper-left column coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Extracts a rectangular part of a data field to a new data field.
 *
 * Returns: The extracted area as a newly created data field.
 **/
GwyDataField*
gwy_data_field_area_extract(GwyDataField *data_field,
                            gint col, gint row,
                            gint width, gint height)
{
    GwyDataField *result;
    gint i;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return NULL;

    if (col == 0
        && row == 0
        && width == data_field->xres
        && height == data_field->yres)
        return gwy_data_field_duplicate(data_field);

    result = gwy_data_field_new(width, height,
                                data_field->xreal*width/data_field->xres,
                                data_field->yreal*height/data_field->yres,
                                FALSE);
    for (i = 0; i < height; i++)
        gwy_assign(result->data + i*width, data_field->data + (i + row)*data_field->xres + col, width);
    gwy_data_field_copy_units(data_field, result);

    return result;
}

/**
 * gwy_data_field_get_dval:
 * @data_field: A data field
 * @x: Horizontal position in pixel units, in range [0, x-resolution].
 * @y: Vertical postition in pixel units, in range [0, y-resolution].
 * @interpolation: Interpolation method to be used.
 *
 * Gets interpolated value at arbitrary data field point indexed by pixel coordinates.
 *
 * Note pixel values are centered in pixels, so to get the same value as gwy_data_field_get_val(@data_field, @j, @i)
 * returns, it's necessary to add 0.5: gwy_data_field_get_dval(@data_field, @j+0.5, @i+0.5, @interpolation).
 *
 * See also gwy_data_field_get_dval_real() that does the same, but takes real coordinates.
 *
 * Returns: Interpolated value at position (@x,@y).
 **/
gdouble
gwy_data_field_get_dval(GwyDataField *a,
                        gdouble x, gdouble y,
                        GwyInterpolationType interp)
{
    gint ix, iy, ixp, iyp;
    gint floorx, floory;
    gdouble restx, resty, valxy, valpy, valxp, valpp, va, vb, vc, vd;
    gdouble *data;
    gdouble intline[4];

    g_return_val_if_fail(GWY_IS_DATA_FIELD(a), 0.0);

    if (G_UNLIKELY(interp == GWY_INTERPOLATION_NONE))
        return 0.0;

    if (interp == GWY_INTERPOLATION_ROUND) {
        /* floor() centers pixel value */
        floorx = floor(x);
        floory = floor(y);
        ix = CLAMP(floorx, 0, a->xres - 1);
        iy = CLAMP(floory, 0, a->yres - 1);
        return a->data[ix + a->xres*iy];
    }
    if (interp == GWY_INTERPOLATION_LINEAR) {
        /* To centered pixel value */
        x -= 0.5;
        y -= 0.5;
        floorx = floor(x);
        floory = floor(y);
        restx = x - floorx;
        resty = y - floory;
        ix = CLAMP(floorx, 0, a->xres - 1);
        iy = CLAMP(floory, 0, a->yres - 1);
        ixp = CLAMP(floorx + 1, 0, a->xres - 1);
        iyp = CLAMP(floory + 1, 0, a->yres - 1);

        valxy = (1.0 - restx)*(1.0 - resty)*a->data[ix + a->xres*iy];
        valxp = (1.0 - restx)*resty*a->data[ix + a->xres*iyp];
        valpy = restx*(1.0 - resty)*a->data[ixp + a->xres*iy];
        valpp = restx*resty*a->data[ixp + a->xres*iyp];
        return valxy + valpy + valxp + valpp;
    }

    /* To centered pixel value */
    x -= 0.5;
    y -= 0.5;
    floorx = floor(x);
    floory = floor(y);
    restx = x - floorx;
    resty = y - floory;

    /* fall back to bilinear for border pixels. */
    if (floorx < 1 || floory < 1 || floorx >= a->xres-2 || floory >= a->yres-2) {
        ix = CLAMP(floorx, 0, a->xres - 1);
        iy = CLAMP(floory, 0, a->yres - 1);
        ixp = CLAMP(floorx + 1, 0, a->xres - 1);
        iyp = CLAMP(floory + 1, 0, a->yres - 1);

        valxy = (1.0 - restx)*(1.0 - resty)*a->data[ix + a->xres*iy];
        valxp = (1.0 - restx)*resty*a->data[ix + a->xres*iyp];
        valpy = restx*(1.0 - resty)*a->data[ixp + a->xres*iy];
        valpp = restx*resty*a->data[ixp + a->xres*iyp];
        return valxy + valpy + valxp + valpp;
    }

    /* interpolation in x direction */
    data = a->data + floorx-1 + a->xres*(floory-1);
    gwy_assign(intline, data, 4);
    va = gwy_interpolation_get_dval_of_equidists(restx, intline, interp);
    gwy_assign(intline, data + a->xres, 4);
    vb = gwy_interpolation_get_dval_of_equidists(restx, intline, interp);
    gwy_assign(intline, data + 2*a->xres, 4);
    vc = gwy_interpolation_get_dval_of_equidists(restx, intline, interp);
    gwy_assign(intline, data + 3*a->xres, 4);
    vd = gwy_interpolation_get_dval_of_equidists(restx, intline, interp);

    /*interpolation in y direction*/
    intline[0] = va;
    intline[1] = vb;
    intline[2] = vc;
    intline[3] = vd;
    return gwy_interpolation_get_dval_of_equidists(resty, intline, interp);
}

/**
 * gwy_data_field_get_data:
 * @data_field: A data field
 *
 * Gets the raw data buffer of a data field.
 *
 * The returned buffer is not guaranteed to be valid through whole data field life time.  Some function may change it,
 * most notably gwy_data_field_resize() and gwy_data_field_resample().
 *
 * This function invalidates any cached information, use gwy_data_field_get_data_const() if you are not going to
 * change the data.
 *
 * See gwy_data_field_invalidate() for some discussion.
 *
 * Returns: The data field as a pointer to an array of gwy_data_field_get_xres()*gwy_data_field_get_yres() #gdouble's,
 *          ordered by lines.  I.e., they are to be accessed as data[row*xres + column].
 **/
gdouble*
gwy_data_field_get_data(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    gwy_data_field_invalidate(data_field);
    return data_field->data;
}

/**
 * gwy_data_field_get_data_const:
 * @data_field: A data field.
 *
 * Gets the raw data buffer of a data field, read-only.
 *
 * The returned buffer is not guaranteed to be valid through whole data field life time.  Some function may change it,
 * most notably gwy_data_field_resize() and gwy_data_field_resample().
 *
 * Use gwy_data_field_get_data() if you want to change the data.
 *
 * See gwy_data_field_invalidate() for some discussion.
 *
 * Returns: The data field as a pointer to an array of gwy_data_field_get_xres()*gwy_data_field_get_yres() #gdouble's,
 *          ordered by lines.  I.e., they are to be accessed as data[row*xres + column].
 **/
const gdouble*
gwy_data_field_get_data_const(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    return (const gdouble*)data_field->data;
}

/**
 * gwy_data_field_get_xres:
 * @data_field: A data field.
 *
 * Gets X resolution (number of columns) of a data field.
 *
 * Returns: X resolution.
 **/
gint
gwy_data_field_get_xres(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0);
    return data_field->xres;
}

/**
 * gwy_data_field_get_yres:
 * @data_field: A data field.
 *
 * Gets Y resolution (number of rows) of the field.
 *
 * Returns: Y resolution.
 **/
gint
gwy_data_field_get_yres(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0);
    return data_field->yres;
}

/**
 * gwy_data_field_get_xreal:
 * @data_field: A data field.
 *
 * Gets the X real (physical) size of a data field.
 *
 * Returns: X real size value.
 **/
gdouble
gwy_data_field_get_xreal(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return data_field->xreal;
}

/**
 * gwy_data_field_get_yreal:
 * @data_field: A data field
 *
 * Gets the Y real (physical) size of a data field.
 *
 * Returns: Y real size value.
 **/
gdouble
gwy_data_field_get_yreal(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return data_field->yreal;
}

/**
 * gwy_data_field_set_xreal:
 * @data_field: A data field.
 * @xreal: New X real size value.
 *
 * Sets X real (physical) size value of a data field.
 **/
void
gwy_data_field_set_xreal(GwyDataField *data_field, gdouble xreal)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(xreal > 0.0);
    if (xreal != data_field->xreal) {
        data_field->cached &= ~(CBIT(ARE) | CBIT(VAR));
        data_field->xreal = xreal;
    }
}

/**
 * gwy_data_field_set_yreal:
 * @data_field: A data field.
 * @yreal: New Y real size value.
 *
 * Sets Y real (physical) size value of a data field.
 **/
void
gwy_data_field_set_yreal(GwyDataField *data_field, gdouble yreal)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(yreal > 0.0);
    if (yreal != data_field->yreal) {
        data_field->cached &= ~(CBIT(ARE) | CBIT(VAR));
        data_field->yreal = yreal;
    }
}

/**
 * gwy_data_field_get_xoffset:
 * @data_field: A data field.
 *
 * Gets the X offset of data field origin.
 *
 * Returns: X offset value.
 **/
gdouble
gwy_data_field_get_xoffset(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return data_field->xoff;
}

/**
 * gwy_data_field_get_yoffset:
 * @data_field: A data field
 *
 * Gets the Y offset of data field origin.
 *
 * Returns: Y offset value.
 **/
gdouble
gwy_data_field_get_yoffset(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return data_field->yoff;
}

/**
 * gwy_data_field_set_xoffset:
 * @data_field: A data field.
 * @xoff: New X offset value.
 *
 * Sets the X offset of a data field origin.
 *
 * Note offsets don't affect any calculation, nor functions like gwy_data_field_rtoj().
 **/
void
gwy_data_field_set_xoffset(GwyDataField *data_field, gdouble xoff)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    data_field->xoff = xoff;
}

/**
 * gwy_data_field_set_yoffset:
 * @data_field: A data field.
 * @yoff: New Y offset value.
 *
 * Sets the Y offset of a data field origin.
 *
 * Note offsets don't affect any calculation, nor functions like gwy_data_field_rtoi().
 **/
void
gwy_data_field_set_yoffset(GwyDataField *data_field, gdouble yoff)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    data_field->yoff = yoff;
}

/**
 * gwy_data_field_get_dx:
 * @data_field: A data field.
 *
 * Gets the horizontal pixel size of a data field in real units.
 *
 * The result is the same as gwy_data_field_get_xreal(data_field)/gwy_data_field_get_xres(data_field).
 *
 * Returns: Horizontal pixel size.
 *
 * Since: 2.52
 **/
gdouble
gwy_data_field_get_dx(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return data_field->xreal/data_field->xres;
}

/**
 * gwy_data_field_get_dy:
 * @data_field: A data field.
 *
 * Gets the vertical pixel size of a data field in real units.
 *
 * The result is the same as gwy_data_field_get_yreal(data_field)/gwy_data_field_get_yres(data_field).
 *
 * Returns: Vertical pixel size.
 *
 * Since: 2.52
 **/
gdouble
gwy_data_field_get_dy(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return data_field->yreal/data_field->yres;
}

/**
 * gwy_data_field_get_si_unit_xy:
 * @data_field: A data field.
 *
 * Returns lateral SI unit of a data field.
 *
 * Returns: SI unit corresponding to the lateral (XY) dimensions of the data field.  Its reference count is not
 *          incremented.
 **/
GwySIUnit*
gwy_data_field_get_si_unit_xy(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);

    if (!data_field->si_unit_xy)
        data_field->si_unit_xy = gwy_si_unit_new(NULL);

    return data_field->si_unit_xy;
}

/**
 * gwy_data_field_get_si_unit_z:
 * @data_field: A data field.
 *
 * Returns value SI unit of a data field.
 *
 * Returns: SI unit corresponding to the "height" (Z) dimension of the data field.  Its reference count is not
 *          incremented.
 **/
GwySIUnit*
gwy_data_field_get_si_unit_z(GwyDataField *data_field)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);

    if (!data_field->si_unit_z)
        data_field->si_unit_z = gwy_si_unit_new(NULL);

    return data_field->si_unit_z;
}

/**
 * gwy_data_field_set_si_unit_xy:
 * @data_field: A data field.
 * @si_unit: SI unit to be set.
 *
 * Sets the SI unit corresponding to the lateral (XY) dimensions of a data field.
 *
 * It does not assume a reference on @si_unit, instead it adds its own reference.
 **/
void
gwy_data_field_set_si_unit_xy(GwyDataField *data_field,
                              GwySIUnit *si_unit)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    _gwy_set_object_si_unit(si_unit, &data_field->si_unit_xy);
}

/**
 * gwy_data_field_set_si_unit_z:
 * @data_field: A data field.
 * @si_unit: SI unit to be set.
 *
 * Sets the SI unit corresponding to the "height" (Z) dimension of a data field.
 *
 * It does not assume a reference on @si_unit, instead it adds its own reference.
 **/
void
gwy_data_field_set_si_unit_z(GwyDataField *data_field,
                             GwySIUnit *si_unit)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    _gwy_set_object_si_unit(si_unit, &data_field->si_unit_z);
}

/**
 * gwy_data_field_get_value_format_xy:
 * @data_field: A data field.
 * @style: Unit format style.
 * @format: A SI value format to modify, or %NULL to allocate a new one.
 *
 * Finds value format good for displaying coordinates of a data field.
 *
 * Returns: The value format.  If @format is %NULL, a newly allocated format is returned, otherwise (modified) @format
 *          itself is returned.
 **/
GwySIValueFormat*
gwy_data_field_get_value_format_xy(GwyDataField *data_field,
                                   GwySIUnitFormatStyle style,
                                   GwySIValueFormat *format)
{
    gdouble max, unit;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);

    max = MAX(data_field->xreal, data_field->yreal);
    unit = MIN(data_field->xreal/data_field->xres,
               data_field->yreal/data_field->yres);
    return gwy_si_unit_get_format_with_resolution(gwy_data_field_get_si_unit_xy(data_field), style, max, unit, format);
}

/**
 * gwy_data_field_get_value_format_z:
 * @data_field: A data field.
 * @style: Unit format style.
 * @format: A SI value format to modify, or %NULL to allocate a new one.
 *
 * Finds value format good for displaying values of a data field.
 *
 * Returns: The value format.  If @format is %NULL, a newly allocated format is returned, otherwise (modified) @format
 *          itself is returned.
 **/
GwySIValueFormat*
gwy_data_field_get_value_format_z(GwyDataField *data_field,
                                  GwySIUnitFormatStyle style,
                                  GwySIValueFormat *format)
{
    GwySIUnit *siunit;
    gdouble max, min;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);

    gwy_data_field_get_min_max(data_field, &min, &max);
    gwy_data_field_get_autorange(data_field, &min, &max);
    if (max == min) {
        max = ABS(max);
        min = 0.0;
    }
    siunit = gwy_data_field_get_si_unit_z(data_field);

    return gwy_si_unit_get_format_with_digits(siunit, style, max - min, 3, format);
}

/**
 * gwy_data_field_itor:
 * @data_field: A data field.
 * @row: Vertical pixel coordinate.
 *
 * Transforms vertical pixel coordinate to real (physical) Y coordinate.
 *
 * That is it maps range [0..y-resolution] to range [0..real-y-size]. It is not suitable for conversion of matrix
 * indices to physical coordinates, you have to use gwy_data_field_itor(@data_field, @row + 0.5) for that.
 *
 * Returns: Real Y coordinate.
 **/
gdouble
gwy_data_field_itor(GwyDataField *data_field, gdouble row)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return row * data_field->yreal/data_field->yres;
}

/**
 * gwy_data_field_jtor:
 * @data_field: A data field.
 * @col: Horizontal pixel coordinate.
 *
 * Transforms horizontal pixel coordinate to real (physical) X coordinate.
 *
 * That is it maps range [0..x-resolution] to range [0..real-x-size]. It is not suitable for conversion of matrix
 * indices to physical coordinates, you have to use gwy_data_field_jtor(@data_field, @col + 0.5) for that.
 *
 * Returns: Real X coordinate.
 **/
gdouble
gwy_data_field_jtor(GwyDataField *data_field, gdouble col)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return col * data_field->xreal/data_field->xres;
}


/**
 * gwy_data_field_rtoi:
 * @data_field: A data field.
 * @realy: Real (physical) Y coordinate.
 *
 * Transforms real (physical) Y coordinate to row.
 *
 * That is it maps range [0..real-y-size] to range [0..y-resolution].
 *
 * Returns: Vertical pixel coodinate.
 **/
gdouble
gwy_data_field_rtoi(GwyDataField *data_field, gdouble realy)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return realy * data_field->yres/data_field->yreal;
}


/**
 * gwy_data_field_rtoj:
 * @data_field: A data field.
 * @realx: Real (physical) X coodinate.
 *
 * Transforms real (physical) X coordinate to column.
 *
 * That is it maps range [0..real-x-size] to range [0..x-resolution].
 *
 * Returns: Horizontal pixel coordinate.
 **/
gdouble
gwy_data_field_rtoj(GwyDataField *data_field, gdouble realx)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return realx * data_field->xres/data_field->xreal;
}

static inline gboolean
gwy_data_field_inside(GwyDataField *data_field, gint i, gint j)
{
    if (i >= 0 && j >= 0 && i < data_field->xres && j < data_field->yres)
        return TRUE;
    else
        return FALSE;
}

/**
 * gwy_data_field_get_val:
 * @data_field: A data field.
 * @col: Column index.
 * @row: Row index.
 *
 * Gets value at given position in a data field.
 *
 * Do not access data with this function inside inner loops, it's slow. Get the raw data buffer with
 * gwy_data_field_get_data_const() and access it directly instead.
 *
 * Returns: Value at (@col, @row).
 **/
gdouble
gwy_data_field_get_val(GwyDataField *data_field, gint col, gint row)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    g_return_val_if_fail(gwy_data_field_inside(data_field, col, row), 0.0);
    return data_field->data[col + data_field->xres*row];
}

/**
 * gwy_data_field_set_val:
 * @data_field: A data field.
 * @col: Column index.
 * @row: Row index.
 * @value: Value to set.
 *
 * Sets value at given position in a data field.
 *
 * Do not set data with this function inside inner loops, it's slow.  Get the raw data buffer with
 * gwy_data_field_get_data() and write to it directly instead.
 **/
void
gwy_data_field_set_val(GwyDataField *data_field,
                       gint col, gint row,
                       gdouble value)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(gwy_data_field_inside(data_field, col, row));
    gwy_data_field_invalidate(data_field);
    data_field->data[col + data_field->xres*row] = value;
}

/**
 * gwy_data_field_get_dval_real:
 * @data_field: A data field.
 * @x: X postion in real coordinates.
 * @y: Y postition in real coordinates.
 * @interpolation: Interpolation method to use.
 *
 * Gets interpolated value at arbitrary data field point indexed by real coordinates.
 *
 * See also gwy_data_field_get_dval() that does the same, but takes pixel coordinates.
 *
 * Returns: Value at position (@x,@y).
 **/
gdouble
gwy_data_field_get_dval_real(GwyDataField *data_field, gdouble x, gdouble y,
                             GwyInterpolationType interpolation)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    return gwy_data_field_get_dval(data_field,
                                   gwy_data_field_rtoj(data_field, x), gwy_data_field_rtoi(data_field, y),
                                   interpolation);
}

/**
 * gwy_data_field_rotate:
 * @data_field: A data field.
 * @angle: Rotation angle (in radians).
 * @interpolation: Interpolation method to use.
 *
 * Rotates a data field by a given angle.
 *
 * This function is mostly obsolete.  See gwy_data_field_new_rotated() and gwy_data_field_new_rotated_90().
 *
 * Values that get outside of data field by the rotation are lost. Undefined values from outside of data field that
 * get inside are set to data field minimum value.
 *
 * The rotation is performed in pixel space, i.e. it can be in fact a more general affine transform in the real
 * coordinates when pixels are not square.
 **/
void
gwy_data_field_rotate(GwyDataField *a,
                      gdouble angle,
                      GwyInterpolationType interpolation)
{
    GwyDataField *b;
    gdouble icor, jcor, sn, cs, val, x, y, v;
    gdouble *coeff;
    gint xres, yres, newi, newj, oldi, oldj, i, j, ii, jj, suplen, sf, st;

    g_return_if_fail(GWY_IS_DATA_FIELD(a));

    suplen = gwy_interpolation_get_support_size(interpolation);
    if (suplen <= 0)
        return;

    angle = gwy_canonicalize_angle(angle, TRUE, TRUE);
    if (fabs(angle) < 1e-15)
        return;
    if (fabs(angle - G_PI) < 2e-15) {
        gwy_data_field_invert(a, TRUE, TRUE, FALSE);
        return;
    }

    if (fabs(angle - G_PI/2) < 1e-15) {
        sn = 1.0;
        cs = 0.0;
    }
    else if (fabs(angle - 3*G_PI/4) < 3e-15) {
        sn = -1.0;
        cs = 0.0;
    }
    else {
        sn = sin(angle);
        cs = cos(angle);
    }

    xres = a->xres;
    yres = a->yres;
    icor = ((yres - 1.0)*(1.0 - cs) - (xres - 1.0)*sn)/2.0;
    jcor = ((xres - 1.0)*(1.0 - cs) + (yres - 1.0)*sn)/2.0;

    coeff = g_new(gdouble, suplen*suplen);
    sf = -((suplen - 1)/2);
    st = suplen/2;

    /* FIXME: Shouldn't we implement this in terms of gwy_data_field_distort()? */
    val = gwy_data_field_get_min(a);
    b = gwy_data_field_duplicate(a);
    gwy_interpolation_resolve_coeffs_2d(xres, yres, xres, b->data, interpolation);

    for (newi = 0; newi < yres; newi++) {
        for (newj = 0; newj < xres; newj++) {
            y = newi*cs + newj*sn + icor;
            x = -newi*sn + newj*cs + jcor;
            if (y > yres || x > xres || y < 0.0 || x < 0.0)
                v = val;
            else {
                oldi = (gint)floor(y);
                y -= oldi;
                oldj = (gint)floor(x);
                x -= oldj;
                for (i = sf; i <= st; i++) {
                    ii = (oldi + i + 2*st*yres) % (2*yres);
                    if (G_UNLIKELY(ii >= yres))
                        ii = 2*yres-1 - ii;
                    for (j = sf; j <= st; j++) {
                        jj = (oldj + j + 2*st*xres) % (2*xres);
                        if (G_UNLIKELY(jj >= xres))
                            jj = 2*xres-1 - jj;
                        coeff[(i - sf)*suplen + j - sf] = b->data[ii*xres + jj];
                    }
                }
                v = gwy_interpolation_interpolate_2d(x, y, suplen, coeff, interpolation);
            }
            a->data[newj + xres*newi] = v;
        }
    }

    g_free(coeff);
    g_object_unref(b);
    gwy_data_field_invalidate(a);
}

/**
 * gwy_data_field_new_rotated_90:
 * @data_field: A data field.
 * @clockwise: %TRUE to rotate clocwise, %FALSE to rotate anti-clockwise.
 *
 * Creates a new data field by rotating a data field by 90 degrees.
 *
 * Returns: A newly created data field.
 *
 * Since: 2.46
 **/
GwyDataField*
gwy_data_field_new_rotated_90(GwyDataField *data_field,
                              gboolean clockwise)
{
    GwyDataField *result;

    result = gwy_data_field_new_alike(data_field, FALSE);
    gwy_data_field_flip_xy(data_field, result, FALSE);
    /* Clockwise = flip + rowinv; Anti-clockwise = flip + colinv. */
    gwy_data_field_invert(result, !clockwise, clockwise, FALSE);
    result->xoff = data_field->yoff;
    result->yoff = data_field->xoff;

    return result;
}

static gboolean
rotate_find_out_dimensions(GwyDataField *dfield,
                           gdouble phi, GwyRotateResizeType resize,
                           gdouble *newxreal, gdouble *newyreal)
{
    gdouble xreal, yreal, sphi, cphi;

    xreal = dfield->xreal;
    yreal = dfield->yreal;
    if (resize == GWY_ROTATE_RESIZE_SAME_SIZE) {
        gdouble q;

        /* FIXME: This should be same area or something like that.  We really do not want to cut the ~π/2-rotated
         * field to the original rectangle for non-square fields! */
        sphi = fabs(sin(phi));
        cphi = fabs(cos(phi));
        q = sqrt(1.0 + (xreal/yreal + yreal/xreal)*cphi*sphi);
        *newxreal = (xreal*cphi + yreal*sphi)/q;
        *newyreal = (yreal*cphi + xreal*sphi)/q;
    }
    else if (resize == GWY_ROTATE_RESIZE_CUT) {
        gdouble c2phi, s2phi;

        /* Make 0 ≤ φ ≤ π. */
        phi = gwy_canonicalize_angle(phi, TRUE, FALSE);
        /* Make 0 ≤ φ ≤ π/2. */
        if (phi > 0.5*G_PI)
            phi = G_PI - phi;

        sphi = sin(phi);
        cphi = cos(phi);
        s2phi = sin(2.0*phi);
        c2phi = cos(2.0*phi);

        if (yreal <= xreal*s2phi) {
            *newxreal = 0.5*yreal/sphi;
            *newyreal = 0.5*yreal/cphi;
        }
        else if (xreal <= yreal*s2phi) {
            *newxreal = 0.5*xreal/cphi;
            *newyreal = 0.5*xreal/sphi;
        }
        else {
            *newxreal = (xreal*cphi - yreal*sphi)/c2phi;
            *newyreal = (yreal*cphi - xreal*sphi)/c2phi;
        }
    }
    else if (resize == GWY_ROTATE_RESIZE_EXPAND) {
        sphi = fabs(sin(phi));
        cphi = fabs(cos(phi));
        *newxreal = xreal*cphi + yreal*sphi;
        *newyreal = yreal*cphi + xreal*sphi;
    }
    else
        return FALSE;

    return TRUE;
}

/**
 * gwy_data_field_new_rotated:
 * @dfield: A data field.
 * @exterior_mask: Optional data field where pixels corresponding to exterior will be set to 1.  It will be resized to
 *                 match the returned field.
 * @angle: Rotation angle (in radians).
 * @interp: Interpolation type to use.
 * @resize: Controls how the result size is determined.
 *
 * Creates a new data field by rotating a data field by an atribtrary angle.
 *
 * The returned data field can have pixel corresponding to exterior in @dfield (unless @resize is
 * %GWY_ROTATE_RESIZE_CUT).  They are filled with a neutral value; pass @exterior_mask and replace them as you wish if
 * you need more control.
 *
 * The rotation is performed in real space, i.e. it is a more general affine transform in the pixel space for data
 * field with non-square pixels. See gwy_data_field_rotate() which rotates in the pixel space.
 *
 * The returned data field has always square pixels.  If you want to rotate by a multiple of %G_PI/2 while preserving
 * non-square pixels, you must use explicitly a function such as gwy_data_field_new_rotated_90().
 *
 * Returns: A newly created data field.
 *
 * Since: 2.46
 **/
GwyDataField*
gwy_data_field_new_rotated(GwyDataField *dfield,
                           GwyDataField *exterior_mask,
                           gdouble angle,
                           GwyInterpolationType interp,
                           GwyRotateResizeType resize)
{
    GwyDataField *result, *coeffield;
    gint xres, yres, newxres, newyres, sf, st, suplen, n;
    gdouble xreal, yreal, newxreal, newyreal, sphi, cphi;
    gdouble dx, dy, h, q;
    gdouble axx, axy, ayx, ayy, bx, by, avg;
    gboolean nonsquare;
    gdouble *dest, *m = NULL;
    const gdouble *src;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(dfield), NULL);
    g_return_val_if_fail(!exterior_mask || GWY_IS_DATA_FIELD(exterior_mask), NULL);

    angle = gwy_canonicalize_angle(angle, TRUE, TRUE);
    if (!rotate_find_out_dimensions(dfield, angle, resize,
                                    &newxreal, &newyreal)) {
        g_return_val_if_reached(NULL);
    }

    suplen = gwy_interpolation_get_support_size(interp);
    g_return_val_if_fail(suplen > 0, NULL);

    xres = dfield->xres;
    yres = dfield->yres;
    if (gwy_interpolation_has_interpolating_basis(interp))
        coeffield = g_object_ref(dfield);
    else {
        coeffield = gwy_data_field_duplicate(dfield);
        gwy_interpolation_resolve_coeffs_2d(xres, yres, xres, gwy_data_field_get_data(coeffield), interp);
    }
    src = coeffield->data;

    xreal = dfield->xreal;
    yreal = dfield->yreal;
    dx = xreal/xres;
    dy = yreal/yres;
    nonsquare = !(fabs(log(dx/dy)) < 1e-9);

    if (nonsquare)
        h = fmin(dx, dy);
    else
        h = sqrt(dx*dy);

    newxres = GWY_ROUND(newxreal/h);
    newyres = GWY_ROUND(newyreal/h);
    newxres = CLAMP(newxres, 1, 32768);
    newyres = CLAMP(newyres, 1, 32768);
    q = (newxreal/newxres)/(newyreal/newyres);
    if (resize == GWY_ROTATE_RESIZE_SAME_SIZE) {
        newxreal /= sqrt(q);
        newyreal *= sqrt(q);
    }
    else if (q > 1.0) {
        /* X pixel size is larger.  So reduce xreal when cutting, enlarge yreal
         * when expanding. */
        if (resize == GWY_ROTATE_RESIZE_CUT)
            xreal /= q;
        else if (resize == GWY_ROTATE_RESIZE_EXPAND)
            yreal *= q;
    }
    else if (q < 1.0) {
        /* Y pixel size is larger.  So reduce yreal when cutting, enlarge xreal
         * when expanding. */
        if (resize == GWY_ROTATE_RESIZE_CUT)
            yreal *= q;
        else if (resize == GWY_ROTATE_RESIZE_EXPAND)
            xreal /= q;
    }
    h = sqrt(newxreal/newxres * newyreal/newyres);

    cphi = cos(angle);
    sphi = sin(angle);
    axx = h/dx*cphi;
    axy = -h/dx*sphi;
    ayx = h/dy*sphi;
    ayy = h/dy*cphi;
    bx = 0.5*xres + 0.5*h/dx*(-(newxres-1)*cphi + (newyres-1)*sphi);
    by = 0.5*yres - 0.5*h/dy*((newxres-1)*sphi + (newyres-1)*cphi);

    result = gwy_data_field_new(newxres, newyres, newxreal, newyreal, FALSE);
    result->xoff = dfield->yoff + 0.5*(yreal - newxreal);
    result->yoff = dfield->xoff + 0.5*(xreal - newyreal);
    gwy_data_field_copy_units(dfield, result);

    if (exterior_mask) {
        gwy_serializable_clone(G_OBJECT(result), G_OBJECT(exterior_mask));
        gwy_data_field_clear(exterior_mask);
        g_object_ref(exterior_mask);
    }
    else
        exterior_mask = gwy_data_field_new_alike(result, TRUE);

    dest = result->data;
    m = exterior_mask->data;

    sf = -((suplen - 1)/2);
    st = suplen/2;

    avg = 0.0;
    n = 0;
#ifdef _OPENMP
#pragma omp parallel if(gwy_threads_are_enabled()) default(none) \
            reduction(+:avg,n) \
            shared(src,dest,m,xres,yres,newxres,newyres,suplen,sf,st,axx,axy,ayx,ayy,bx,by,interp)
#endif
    {
        gdouble *coeff = g_new(gdouble, suplen*suplen);
        gint newifrom = gwy_omp_chunk_start(newyres);
        gint newito = gwy_omp_chunk_end(newyres);
        gint newi, newj;

        for (newi = newifrom; newi < newito; newi++) {
            for (newj = 0; newj < newxres; newj++) {
                gdouble x = axx*newj + axy*newi + bx;
                gdouble y = ayx*newj + ayy*newi + by;
                gdouble v = 0.0;
                gint i, j, ii, jj, oldi, oldj;

                if (x >= 0.0 && y >= 0.0 && x < xres && y < yres) {
                    oldi = (gint)floor(y);
                    y -= oldi;
                    oldj = (gint)floor(x);
                    x -= oldj;
                    for (i = sf; i <= st; i++) {
                        ii = (oldi + i + 2*st*yres) % (2*yres);
                        if (G_UNLIKELY(ii >= yres))
                            ii = 2*yres-1 - ii;
                        for (j = sf; j <= st; j++) {
                            jj = (oldj + j + 2*st*xres) % (2*xres);
                            if (G_UNLIKELY(jj >= xres))
                                jj = 2*xres-1 - jj;
                            coeff[(i - sf)*suplen + j - sf] = src[ii*xres + jj];
                        }
                    }
                    v = gwy_interpolation_interpolate_2d(x, y, suplen, coeff, interp);
                    avg += v;
                    n++;
                }
                else
                    m[newxres*newi + newj] = 1.0;
                dest[newxres*newi + newj] = v;
            }
        }
        g_free(coeff);
    }

    if (n != newxres*newyres) {
        avg /= n;
        for (n = 0; n < newxres*newyres; n++) {
            if (m[n])
                dest[n] = avg;
        }
    }

    g_object_unref(coeffield);
    gwy_data_field_invalidate(exterior_mask);
    g_object_unref(exterior_mask);

    return result;
}

static void
invert_array_in_place(gdouble *d, guint n)
{
    gdouble *e = d + n-1;

    n /= 2;
    while (n--) {
        GWY_SWAP(gdouble, *d, *e);
        d++;
        e--;
    }
}

/**
 * gwy_data_field_invert:
 * @data_field: A data field.
 * @x: %TRUE to reflect Y, i.e. rows within the XY plane.
 * @y: %TRUE to reflect X, i.e. columns within the XY plane.
 * @z: %TRUE to invert values.
 *
 * Reflects and/or inverts a data field.
 *
 * In the case of value reflection, it's inverted about the mean value.
 *
 * Note that the axis parameter convention is confusing and different from gwy_brick_invert() and
 * gwy_data_line_invert().  Parameters @x an @y correspond the axes around which to flip (which themselves stay
 * unchanged). You may need to swap @x and @y arguments compared what you would pass naturally.
 **/
void
gwy_data_field_invert(GwyDataField *data_field,
                      gboolean x,
                      gboolean y,
                      gboolean z)
{
    gint xres, yres, i, j, n;
    gdouble avg;
    gdouble *data, *flip;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    xres = data_field->xres;
    yres = data_field->yres;
    n = xres*yres;

    if (z) {
        avg = gwy_data_field_get_avg(data_field);
        data = data_field->data;
        for (i = 0; i < n; i++)
            data[i] = 2.0 * avg - data[i];

        /* We can transform stats */
        data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED) | CBIT(ARF) | CBIT(ART)
                              | CBIT(ARE) | CBIT(VAR);
        CVAL(data_field, MIN) = 2.0 * avg - CVAL(data_field, MIN);
        CVAL(data_field, MAX) = 2.0 * avg - CVAL(data_field, MAX);
        GWY_SWAP(gdouble, CVAL(data_field, MIN), CVAL(data_field, MAX));
        CVAL(data_field, SUM) = 2.0 * n * avg - CVAL(data_field, SUM);
        /* RMS doesn't change */
        CVAL(data_field, MED) = 2.0 * avg - CVAL(data_field, MED);
        CVAL(data_field, ARF) = 2.0 * avg - CVAL(data_field, ARF);
        CVAL(data_field, ART) = 2.0 * avg - CVAL(data_field, ART);
        GWY_SWAP(gdouble, CVAL(data_field, ARF), CVAL(data_field, ART));
        /* Area doesn't change */
    }

    if (x && y) {
        invert_array_in_place(data_field->data, n);
    }
    else if (y) {
        for (i = 0; i < yres; i++)
            invert_array_in_place(data_field->data + i*xres, xres);
    }
    else if (x) {
        for (i = 0; i < yres/2; i++) {
            data = data_field->data + i*xres;
            flip = data_field->data + (yres-1 - i)*xres;
            for (j = 0; j < xres; j++, data++, flip++)
                GWY_SWAP(gdouble, *data, *flip);
        }
    }
    else
        return;

    /* No cached value changes */
    data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED) | CBIT(ARF) | CBIT(ART) | CBIT(ARE)
                          | CBIT(VAR);
}

/* Block sizes are measured in destination, in source, the dims are swapped. */
static inline void
swap_block(const gdouble *sb, gdouble *db,
           guint xblocksize, guint yblocksize,
           guint dxres, guint sxres)
{
    guint i, j;

    for (i = 0; i < yblocksize; i++) {
        const gdouble *s = sb + i;
        gdouble *d = db + i*dxres;
        for (j = xblocksize; j; j--, d++, s += sxres)
            *d = *s;
    }
}

static void
transpose_to(const GwyDataField *source,
             guint col, guint row,
             guint width, guint height,
             GwyDataField *dest,
             guint destcol, guint destrow)
{
    guint dxres = dest->xres, sxres = source->xres;
    guint jmax = height/BLOCK_SIZE * BLOCK_SIZE;
    guint imax = width/BLOCK_SIZE * BLOCK_SIZE;
    const gdouble *sbase = source->data + sxres*row + col;
    gdouble *dbase = dest->data + dxres*destrow + destcol;
    guint ib, jb;

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(ib, jb) \
            shared(sbase,dbase,sxres,dxres,imax,jmax,height)
#endif
    for (ib = 0; ib < imax; ib += BLOCK_SIZE) {
        for (jb = 0; jb < jmax; jb += BLOCK_SIZE)
            swap_block(sbase + (jb*sxres + ib), dbase + (ib*dxres + jb), BLOCK_SIZE, BLOCK_SIZE, dxres, sxres);
        if (jmax != height)
            swap_block(sbase + (jmax*sxres + ib), dbase + (ib*dxres + jmax), height - jmax, BLOCK_SIZE, dxres, sxres);
    }
    if (imax != width) {
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            private(jb) \
            shared(sbase,dbase,sxres,dxres,imax,jmax,width)
#endif
        for (jb = 0; jb < jmax; jb += BLOCK_SIZE)
            swap_block(sbase + (jb*sxres + imax), dbase + (imax*dxres + jb), BLOCK_SIZE, width - imax, dxres, sxres);
        if (jmax != height) {
            swap_block(sbase + (jmax*sxres + imax), dbase + (imax*dxres + jmax),
                       height - jmax, width - imax, dxres, sxres);
        }
    }
}

/**
 * gwy_data_field_flip_xy:
 * @src: Source data field.
 * @dest: Destination data field.
 * @minor: %TRUE to mirror about the minor diagonal; %FALSE to mirror about
 *         major diagonal.
 *
 * Copies data from one data field to another with transposition.
 *
 * The destination data field is resized as necessary, its real dimensions set to transposed @src dimensions and its
 * offsets are reset.  Units are not updated.
 *
 * Since: 2.49
 **/
void
gwy_data_field_flip_xy(GwyDataField *src, GwyDataField *dest,
                       gboolean minor)
{
    gint xres, yres;

    g_return_if_fail(GWY_IS_DATA_FIELD(src));
    g_return_if_fail(GWY_IS_DATA_FIELD(dest));
    xres = src->xres;
    yres = src->yres;
    gwy_data_field_resample(dest, yres, xres, GWY_INTERPOLATION_NONE);
    transpose_to(src, 0, 0, xres, yres, dest, 0, 0);
    if (minor)
        invert_array_in_place(dest->data, xres*yres);
    dest->yreal = src->xreal;
    dest->xreal = src->yreal;
    dest->xoff = dest->yoff = 0.0;
}

/**
 * gwy_data_field_area_flip_xy:
 * @src: Source data field.
 * @col: Upper-left column coordinate in @src.
 * @row: Upper-left row coordinate in @src.
 * @width: Area width (number of columns) in @src.
 * @height: Area height (number of rows) in @src.
 * @dest: Destination data field.
 * @minor: %TRUE to mirror about the minor diagonal; %FALSE to mirror about
 *         major diagonal.
 *
 * Copies data from a rectangular part of one data field to another with transposition.
 *
 * The destination data field is resized as necessary, its real dimensions set to transposed @src area dimensions and
 * its offsets are reset.  Units are not updated.
 *
 * Since: 2.49
 **/
void
gwy_data_field_area_flip_xy(GwyDataField *src,
                            gint col, gint row, gint width, gint height,
                            GwyDataField *dest,
                            gboolean minor)
{
    if (!_gwy_data_field_check_area(src, col, row, width, height))
        return;
    g_return_if_fail(GWY_IS_DATA_FIELD(dest));

    gwy_data_field_resample(dest, height, width, GWY_INTERPOLATION_NONE);
    transpose_to(src, col, row, width, height, dest, 0, 0);
    if (minor)
        invert_array_in_place(dest->data, width*height);

    dest->yreal = dest->yres * gwy_data_field_get_dx(src);
    dest->xreal = dest->xres * gwy_data_field_get_dy(src);
    dest->xoff = dest->yoff = 0.0;
}

/**
 * gwy_data_field_fill:
 * @data_field: A data field.
 * @value: Value to be entered.
 *
 * Fills a data field with given value.
 **/
void
gwy_data_field_fill(GwyDataField *data_field, gdouble value)
{
    gint i;
    gdouble *p = data_field->data;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    for (i = data_field->xres * data_field->yres; i; i--, p++)
        *p = value;

    /* We can precompute stats */
    set_cache_for_constant_field(data_field, value);
}

static void
set_cache_for_constant_field(GwyDataField *data_field, gdouble value)
{
    data_field->cached = CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED) | CBIT(ARF) | CBIT(ART)
                         | CBIT(ARE) | CBIT(VAR) | CBIT(MSQ);
    CVAL(data_field, MIN) = value;
    CVAL(data_field, MAX) = value;
    CVAL(data_field, SUM) = data_field->xres * data_field->yres * value;
    CVAL(data_field, RMS) = 0.0;
    CVAL(data_field, MED) = value;
    CVAL(data_field, ARF) = value;
    CVAL(data_field, ART) = value;
    CVAL(data_field, ARE) = data_field->xreal * data_field->yreal;
    CVAL(data_field, VAR) = 0.0;
    CVAL(data_field, MSQ) = value*value;
}

/**
 * gwy_data_field_area_fill:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @value: Value to be entered
 *
 * Fills a rectangular part of a data field with given value.
 **/
void
gwy_data_field_area_fill(GwyDataField *data_field,
                         gint col, gint row, gint width, gint height,
                         gdouble value)
{
    gint i, j;
    gdouble *drow;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;

    for (i = 0; i < height; i++) {
        drow = data_field->data + (row + i)*data_field->xres + col;

        for (j = 0; j < width; j++)
            *(drow++) = value;
    }
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_area_fill_mask:
 * @data_field: A data field.
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @mode: Masking mode to use.  See the introduction for description of masking modes.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @value: Value to be entered
 *
 * Fills a masked rectangular part of a data field with given value.
 *
 * Since: 2.44
 **/
void
gwy_data_field_area_fill_mask(GwyDataField *data_field,
                              GwyDataField *mask,
                              GwyMaskingType mode,
                              gint col, gint row, gint width, gint height,
                              gdouble value)
{
    gint i, j;
    gdouble *drow;
    const gdouble *mrow;

    if (!mask || mode == GWY_MASK_IGNORE) {
        gwy_data_field_area_fill(data_field, col, row, width, height, value);
        return;
    }
    if (!_gwy_data_field_check_area(data_field, col, row, width, height)
        || !_gwy_data_field_check_mask(data_field, &mask, &mode))
        return;

    for (i = 0; i < height; i++) {
        drow = data_field->data + (row + i)*data_field->xres + col;
        mrow = mask->data + (row + i)*mask->xres + col;

        for (j = 0; j < width; j++) {
            if ((mode == GWY_MASK_INCLUDE && *mrow > 0.0) || (mode == GWY_MASK_EXCLUDE && *mrow < 1.0))
                *drow = value;
            drow++;
            mrow++;
        }
    }
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_clear:
 * @data_field: A data field.
 *
 * Fills a data field with zeroes.
 **/
void
gwy_data_field_clear(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_clear(data_field->data, data_field->xres*data_field->yres);

    /* We can precompute stats */
    set_cache_for_constant_field(data_field, 0.0);
}

/**
 * gwy_data_field_area_clear:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Fills a rectangular part of a data field with zeroes.
 **/
void
gwy_data_field_area_clear(GwyDataField *data_field,
                          gint col, gint row, gint width, gint height)
{
    gint i;
    gdouble *drow;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;

    gwy_data_field_invalidate(data_field);
    if (height == 1 || (col == 0 && width == data_field->xres)) {
        gwy_clear(data_field->data + data_field->xres*row + col, width*height);
        return;
    }

    for (i = 0; i < height; i++) {
        drow = data_field->data + (row + i)*data_field->xres + col;
        gwy_clear(drow, width);
    }
}

/**
 * gwy_data_field_multiply:
 * @data_field: A data field.
 * @value: Value to multiply @data_field with.
 *
 * Multiplies all values in a data field by given value.
 **/
void
gwy_data_field_multiply(GwyDataField *data_field, gdouble value)
{
    gint i;
    gdouble *p;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    p = data_field->data;
    for (i = data_field->xres * data_field->yres; i; i--, p++)
        *p *= value;

    /* We can transform stats */
    data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED) | CBIT(ARF) | CBIT(ART) | CBIT(MSQ);
    CVAL(data_field, MIN) *= value;
    CVAL(data_field, MAX) *= value;
    CVAL(data_field, SUM) *= value;
    CVAL(data_field, RMS) *= fabs(value);
    CVAL(data_field, MED) *= value;
    CVAL(data_field, ARF) *= value;
    CVAL(data_field, ART) *= value;
    CVAL(data_field, MSQ) *= value*value;
    if (value < 0) {
        GWY_SWAP(gdouble, CVAL(data_field, MIN), CVAL(data_field, MAX));
        GWY_SWAP(gdouble, CVAL(data_field, ARF), CVAL(data_field, ART));
    }
}

/**
 * gwy_data_field_area_multiply:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @value: Value to multiply area with.
 *
 * Multiplies values in a rectangular part of a data field by given value
 **/
void
gwy_data_field_area_multiply(GwyDataField *data_field,
                             gint col, gint row, gint width, gint height,
                             gdouble value)
{
    gint i, j;
    gdouble *drow;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(col >= 0 && row >= 0 && width >= 0 && height >= 0
                     && col + width <= data_field->xres && row + height <= data_field->yres);

    for (i = 0; i < height; i++) {
        drow = data_field->data + (row + i)*data_field->xres + col;

        for (j = 0; j < width; j++)
            *(drow++) *= value;
    }
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_add:
 * @data_field: A data field.
 * @value: Value to be added to data field values.
 *
 * Adds given value to all values in a data field.
 **/
void
gwy_data_field_add(GwyDataField *data_field, gdouble value)
{
    gint i;
    gdouble *p;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    p = data_field->data;
    for (i = data_field->xres * data_field->yres; i; i--, p++)
        *p += value;

    /* We can transform stats */
    data_field->cached &= CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS) | CBIT(MED) | CBIT(ARF) | CBIT(ART) | CBIT(ARE)
                          | CBIT(VAR);
    CVAL(data_field, MIN) += value;
    CVAL(data_field, MAX) += value;
    CVAL(data_field, SUM) += data_field->xres * data_field->yres * value;
    /* RMS doesn't change */
    CVAL(data_field, MED) += value;
    CVAL(data_field, ARF) += value;
    CVAL(data_field, ART) += value;
    /* Area doesn't change */
    /* There is transformation formula for MSQ, but it can be prone to ugly cancellation errors. */
}

/**
 * gwy_data_field_abs:
 * @data_field: A data field.
 *
 * Takes absolute value of all values in a data field.
 *
 * Since: 2.52
 **/
void
gwy_data_field_abs(GwyDataField *data_field)
{
    gint i;
    gdouble *p;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    p = data_field->data;
    for (i = data_field->xres * data_field->yres; i; i--, p++)
        *p = fabs(*p);

    /* The only cached stat we could transform is the maximum.  That's probably not worth the fuss. */
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_area_add:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 * @value: Value to be added to area values.
 *
 * Adds given value to all values in a rectangular part of a data field.
 **/
void
gwy_data_field_area_add(GwyDataField *data_field,
                        gint col, gint row, gint width, gint height,
                        gdouble value)
{
    gint i, j;
    gdouble *drow;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(col >= 0 && row >= 0 && width >= 0 && height >= 0
                     && col + width <= data_field->xres && row + height <= data_field->yres);

    for (i = 0; i < height; i++) {
        drow = data_field->data + (row + i)*data_field->xres + col;

        for (j = 0; j < width; j++)
            *(drow++) += value;
    }
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_area_abs:
 * @data_field: A data field.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Takes absolute value of values in a rectangular part of a data field.
 *
 * Since: 2.52
 **/
void
gwy_data_field_area_abs(GwyDataField *data_field,
                        gint col, gint row, gint width, gint height)
{
    gint i, j;
    gdouble *drow;

    if (!_gwy_data_field_check_area(data_field, col, row, width, height))
        return;

    for (i = 0; i < height; i++) {
        drow = data_field->data + (row + i)*data_field->xres + col;

        for (j = 0; j < width; j++, drow++)
            *drow = fabs(*drow);
    }
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_get_row:
 * @data_field: A data field.
 * @data_line: A data line.  It will be resized to width ot @data_field.
 * @row: Row index.
 *
 * Extracts a data field row into a data line.
 **/
void
gwy_data_field_get_row(GwyDataField *data_field,
                       GwyDataLine* data_line,
                       gint row)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(row >= 0 && row < data_field->yres);

    gwy_data_line_resample(data_line, data_field->xres, GWY_INTERPOLATION_NONE);
    data_line->real = data_field->xreal;
    gwy_assign(data_line->data, data_field->data + row*data_field->xres, data_field->xres);
    gwy_data_field_copy_units_to_data_line(data_field, data_line);
}


/**
 * gwy_data_field_get_column:
 * @data_field: A data field
 * @data_line: A data line.  It will be resized to height of @data_field.
 * @col: Column index.
 *
 * Extracts a data field column into a data line.
 **/
void
gwy_data_field_get_column(GwyDataField *data_field,
                          GwyDataLine* data_line,
                          gint col)
{
    gint k;
    gdouble *p;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(col >= 0 && col < data_field->xres);

    gwy_data_line_resample(data_line, data_field->yres, GWY_INTERPOLATION_NONE);
    data_line->real = data_field->yreal;
    p = data_field->data + col;
    for (k = 0; k < data_field->yres; k++)
        data_line->data[k] = p[k*data_field->xres];
    gwy_data_field_copy_units_to_data_line(data_field, data_line);
}

/**
 * gwy_data_field_get_row_part:
 * @data_field: A data field.
 * @data_line: A data line.  It will be resized to the row part width.
 * @row: Row index.
 * @from: Start column index.
 * @to: End column index + 1.
 *
 * Extracts part of a data field row into a data line.
 **/
void
gwy_data_field_get_row_part(GwyDataField *data_field,
                            GwyDataLine *data_line,
                            gint row,
                            gint from,
                            gint to)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(row >= 0 && row < data_field->yres);
    GWY_ORDER(gint, from, to);

    if (data_line->res != (to - from))
        gwy_data_line_resample(data_line, to - from, GWY_INTERPOLATION_NONE);

    data_line->real = data_field->xreal*(to - from)/data_field->xres;
    gwy_assign(data_line->data, data_field->data + row*data_field->xres + from, to - from);
    gwy_data_field_copy_units_to_data_line(data_field, data_line);
}

/**
 * gwy_data_field_get_column_part:
 * @data_field: A data field.
 * @data_line: A data line.  It will be resized to the column part height.
 * @col: Column index.
 * @from: Start row index.
 * @to: End row index + 1.
 *
 * Extracts part of a data field column into a data line.
 **/
void
gwy_data_field_get_column_part(GwyDataField *data_field,
                               GwyDataLine *data_line,
                               gint col,
                               gint from,
                               gint to)
{
    gint k;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(col >= 0 && col < data_field->xres);
    GWY_ORDER(gint, from, to);

    if (data_line->res != (to - from))
        gwy_data_line_resample(data_line, to-from, GWY_INTERPOLATION_NONE);

    data_line->real = data_field->yreal*(to - from)/data_field->yres;
    for (k = 0; k < to - from; k++)
        data_line->data[k] = data_field->data[(k+from)*data_field->xres + col];
    gwy_data_field_copy_units_to_data_line(data_field, data_line);
}

/**
 * gwy_data_field_set_row_part:
 * @data_field: A data field.
 * @data_line: A data line.
 * @row: Row index.
 * @from: Start row index.
 * @to: End row index + 1.
 *
 * Puts a data line into a data field row.
 *
 * If data line length differs from @to-@from, it is resampled to this length.
 **/
void
gwy_data_field_set_row_part(GwyDataField *data_field,
                            GwyDataLine *data_line,
                            gint row,
                            gint from,
                            gint to)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(row >= 0 && row < data_field->yres);
    GWY_ORDER(gint, from, to);

    if (data_line->res != (to - from))
        gwy_data_line_resample(data_line, to-from, GWY_INTERPOLATION_LINEAR);

    gwy_assign(data_field->data + row*data_field->xres + from, data_line->data, to - from);
    gwy_data_field_invalidate(data_field);
}


/**
 * gwy_data_field_set_column_part:
 * @data_field: A data field.
 * @data_line: A data line.
 * @col: Column index.
 * @from: Start row index.
 * @to: End row index + 1.
 *
 * Puts a data line into data field column.
 *
 * If data line length differs from @to-@from, it is resampled to this length.
 **/
void
gwy_data_field_set_column_part(GwyDataField *data_field,
                               GwyDataLine* data_line,
                               gint col,
                               gint from,
                               gint to)
{
    gint k;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(col >= 0 && col < data_field->xres);
    GWY_ORDER(gint, from, to);

    if (data_line->res != (to-from))
        gwy_data_line_resample(data_line, to-from, GWY_INTERPOLATION_LINEAR);

    for (k = 0; k < to-from; k++)
        data_field->data[(k+from)*data_field->xres + col] = data_line->data[k];
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_set_row:
 * @data_field: A data field.
 * @data_line: A data line.
 * @row: Row index.
 *
 * Sets a row in the data field to values of a data line.
 *
 * Data line length must be equal to width of data field.
 **/
void
gwy_data_field_set_row(GwyDataField *data_field,
                       GwyDataLine* data_line,
                       gint row)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(row >= 0 && row < data_field->yres);
    g_return_if_fail(data_field->xres == data_line->res);

    gwy_assign(data_field->data + row*data_field->xres, data_line->data, data_field->xres);
    gwy_data_field_invalidate(data_field);
}


/**
 * gwy_data_field_set_column:
 * @data_field: A data field.
 * @data_line: A data line.
 * @col: Column index.
 *
 * Sets a column in the data field to values of a data line.
 *
 * Data line length must be equal to height of data field.
 **/
void
gwy_data_field_set_column(GwyDataField *data_field,
                          GwyDataLine* data_line,
                          gint col)
{
    gint k;
    gdouble *p;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(col >= 0 && col < data_field->xres);
    g_return_if_fail(data_field->yres == data_line->res);

    p = data_field->data + col;
    for (k = 0; k < data_field->yres; k++)
        p[k*data_field->xres] = data_line->data[k];
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_get_profile:
 * @data_field: A data field.
 * @data_line: A data line.  It will be resized to @res samples.  It is possible to pass %NULL to instantiate and
 *             return a new #GwyDataLine.
 * @scol: The column the line starts at (inclusive).
 * @srow: The row the line starts at (inclusive).
 * @ecol: The column the line ends at (inclusive).
 * @erow: The row the line ends at (inclusive).
 * @res: Requested resolution of data line (the number of samples to take). If nonpositive, data line resolution is
 *       chosen to match @data_field's.
 * @thickness: Thickness of line to be averaged.
 * @interpolation: Interpolation type to use.
 *
 * Extracts a possibly averaged profile from data field to a data line.
 *
 * Returns: @data_line itself if it was not %NULL, otherwise a newly created data line.
 **/
GwyDataLine*
gwy_data_field_get_profile(GwyDataField *data_field,
                           GwyDataLine *data_line,
                           gint scol, gint srow,
                           gint ecol, gint erow,
                           gint res,
                           gint thickness,
                           GwyInterpolationType interpolation)
{
    gint k, j;
    gdouble cosa, sina, size, mid, sum;
    gdouble col, row, srcol, srrow;
    gint xres, yres;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    g_return_val_if_fail(!data_line || GWY_IS_DATA_LINE(data_line), NULL);
    xres = data_field->xres;
    yres = data_field->yres;
    g_return_val_if_fail(scol >= 0 && srow >= 0 && ecol >= 0 && erow >= 0
                         && srow < yres && scol < xres && erow < yres && ecol < xres,
                         NULL);

    size = hypot(abs(scol - ecol) + 1, abs(srow - erow) + 1);
    size = MAX(size, 1.0);
    if (res <= 0)
        res = GWY_ROUND(size);

    cosa = (ecol - scol)/(res - 1.0);
    sina = (erow - srow)/(res - 1.0);

    /* Extract regular one-pixel line */
    if (data_line)
        gwy_data_line_resample(data_line, res, GWY_INTERPOLATION_NONE);
    else
        data_line = gwy_data_line_new(res, 1.0, FALSE);

    for (k = 0; k < res; k++)
        data_line->data[k] = gwy_data_field_get_dval(data_field, scol + 0.5 + k*cosa, srow + 0.5 + k*sina,
                                                     interpolation);
    data_line->real = hypot(abs(scol - ecol)*data_field->xreal/xres, abs(srow - erow)*data_field->yreal/yres);
    data_line->real *= res/(res - 1.0);
    gwy_data_field_copy_units_to_data_line(data_field, data_line);

    if (thickness <= 1)
        return data_line;

    /*add neighbour values to the line*/
    for (k = 0; k < res; k++) {
        mid = data_line->data[k];
        sum = 0;
        srcol = scol + 0.5 + k*cosa;
        srrow = srow + 0.5 + k*sina;
        for (j = -thickness/2; j < thickness - thickness/2; j++) {
            col = srcol + j*sina;
            row = srrow - j*cosa;
            if (col >= 0 && col < (xres-1) && row >= 0 && row < (yres-1))
                sum += gwy_data_field_get_dval(data_field, col, row, interpolation);
            else
                sum += mid;
        }
        data_line->data[k] = sum/(gdouble)thickness;
    }

    return data_line;
}

/**
 * gwy_data_field_get_profile_mask:
 * @data_field: A data field.
 * @ndata: Location where to store the actual number of extracted points, which may differ from @res.
 * @mask: Mask specifying which values to take into account/exclude, or %NULL.
 * @masking: Masking mode to use.
 * @xfrom: The real @x-coordinate where the line starts.
 * @yfrom: The real @y-coordinate where line starts.
 * @xto: The real @x-coordinate where the line ends.
 * @yto: The real @y-coordinate where line ends.
 * @res: Requested resolution, i.e. the number of samples to take. If nonpositive, sampling is chosen to match
 *       @data_field's.
 * @thickness: Thickness of line to be averaged.
 * @interpolation: Interpolation type to use.
 *
 * Extracts a possibly averaged profile from data field, with masking.
 *
 * The extracted profile can contain holes due to masking.  It can also contain no points at all if the all data
 * values along the profile were excluded due to masking – in this case %NULL is returned.
 *
 * Unlike gwy_data_field_get_profile(), this function takes real coordinates (without offsets), not row and column
 * indices.
 *
 * Returns: A newly allocated array of #GwyXY coordinare pairs, or %NULL. The caller must free the returned array with
 *          g_free().
 *
 * Since: 2.49
 **/
GwyXY*
gwy_data_field_get_profile_mask(GwyDataField *dfield,
                                gint *ndata,
                                GwyDataField *mask,
                                GwyMaskingType masking,
                                gdouble xfrom, gdouble yfrom,
                                gdouble xto, gdouble yto,
                                gint res,
                                gint thickness,
                                GwyInterpolationType interpolation)
{
    gint k, i, j, kk;
    gdouble xreal, yreal, dx, dy, xstep, ystep, step, size, tx, ty, h;
    gint xres, yres, tres, n;
    GwyXY *xydata;
    const gdouble *m;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(dfield), NULL);
    g_return_val_if_fail(!mask || GWY_IS_DATA_FIELD(mask), NULL);
    g_return_val_if_fail(ndata, NULL);

    if (masking == GWY_MASK_IGNORE)
        mask = NULL;
    else if (!mask)
        masking = GWY_MASK_IGNORE;
    m = mask ? mask->data : NULL;

    xres = dfield->xres;
    yres = dfield->yres;
    xreal = dfield->xreal;
    yreal = dfield->yreal;
    dx = xreal/xres;
    dy = yreal/yres;

    size = hypot(fabs(xto - xfrom)/dx + 1, fabs(yto - yfrom)/dy + 1);
    size = MAX(size, 1.0);
    if (res <= 0)
        res = GWY_ROUND(size);

    gwy_debug("size: %g, res: %d", size, res);
    if (xto == xfrom && yto == yfrom) {
        xto += 0.2*dx;
        yto += 0.2*dy;
        xfrom -= 0.2*dx;
        yfrom -= 0.2*dy;
    }
    xstep = (xto - xfrom)/(res - 1.0);
    ystep = (yto - yfrom)/(res - 1.0);
    step = hypot(xstep, ystep);
    gwy_debug("step (%g, %g)", xstep, ystep);

    if (thickness <= 1) {
        tres = 0;
        tx = ty = 0.0;
    }
    else {
        tres = 2*(thickness - 1);
        tx = (yto - yfrom)/dy;
        ty = -(xto - xfrom)/dx;
        h = hypot(tx, ty);
        tx *= dx/h * 0.5*thickness/tres;
        ty *= dy/h * 0.5*thickness/tres;
    }
    gwy_debug("tres: %d, tstep (%g, %g)", tres, tx, ty);

    xydata = g_new0(GwyXY, res);

    n = 0;
    for (k = 0; k < res; k++) {
        gdouble xc = xfrom + xstep*k;
        gdouble yc = yfrom + ystep*k;
        gdouble z = 0.0;
        gint w = 0;

        for (kk = -tres; kk <= tres; kk++) {
            gdouble x = xc + kk*tx;
            gdouble y = yc + kk*ty;

            x = CLAMP(x, 0.0, 0.999999*xreal);
            y = CLAMP(y, 0.0, 0.999999*yreal);

            if (masking != GWY_MASK_IGNORE) {
                i = (gint)floor(y/dy);
                j = (gint)floor(x/dx);
                if ((masking == GWY_MASK_INCLUDE && m[i*xres + j] <= 0.0)
                    || (masking == GWY_MASK_EXCLUDE && m[i*xres + j] >= 1.0))
                    continue;
            }

            z += gwy_data_field_get_dval_real(dfield, x, y, interpolation);
            w++;
        }
        gwy_debug("[%d] %d", k, w);
        if (w) {
            xydata[n].x = step*k;
            xydata[n].y = z/w;
            n++;
        }
    }

    *ndata = n;
    if (!n)
        GWY_FREE(xydata);

    return xydata;
}

/**
 * gwy_data_field_new_binned:
 * @data_field: A data field.
 * @binw: Bin height (in pixels).
 * @binh: Bin width (in pixels).
 * @xoff: Horizontal offset of bins (in pixels).
 * @yoff: Vertical offset of bins (in pixels).
 * @trimlowest: Number of lowest values to discard.
 * @trimhighest: Number of highest values to discard.
 *
 * Creates a new data field by binning an existing one.
 *
 * The data field is divided into rectangles of dimensions @binw×@binh, offset by (@xoff, @yoff).  The values in each
 * complete rectangle are averaged and the average becomes the pixel value in the newly created, smaller data field.
 *
 * Note that the result is the average – not sum – of the individual values. Multiply the returned data field with
 * @binw×@binh if you want sum.
 *
 * By giving non-zero @trimlowest and @trimhighest you can change the plain average to a trimmed one (even turning it
 * to median in the extreme case). It must always hold that @trimlowest + @trimhighest is smaller than @binw×@binh.
 *
 * Returns: A newly created data field.
 *
 * Since: 2.50
 **/
GwyDataField*
gwy_data_field_new_binned(GwyDataField *data_field,
                          gint binw, gint binh,
                          gint xoff, gint yoff,
                          gint trimlowest, gint trimhighest)
{
    GwyDataField *result;
    gint binsize, xres, yres, newxres, newyres, i, j, k;
    gdouble xreal, yreal, xoffset, yoffset, z;
    gdouble *buf, *d, *r;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), NULL);
    g_return_val_if_fail(binw > 0, NULL);
    g_return_val_if_fail(binh > 0, NULL);
    binsize = binw*binh;
    g_return_val_if_fail(trimlowest >= 0 && trimlowest < binsize, NULL);
    g_return_val_if_fail(trimhighest >= 0 && trimhighest < binsize - trimlowest, NULL);
    if (binsize == 1)
        return gwy_data_field_duplicate(data_field);

    xoff = ((xoff % binw) + binw) % binw;
    yoff = ((yoff % binh) + binh) % binh;
    xres = data_field->xres;
    yres = data_field->yres;
    xreal = data_field->xreal;
    yreal = data_field->yreal;
    xoffset = data_field->xoff;
    yoffset = data_field->yoff;

    if (xres < binw + xoff || yres < binh + yoff) {
        g_warning("No complete bin can be formed.");
        result = gwy_data_field_new(1, 1, xreal, yreal, FALSE);
        gwy_data_field_copy_units(data_field, result);
        result->xoff = xoffset;
        result->yoff = yoffset;
        result->data[0] = gwy_data_field_get_avg(data_field);
        return result;
    }

    newxres = (xres - xoff)/binw;
    newyres = (yres - yoff)/binh;
    result = gwy_data_field_new(newxres, newyres, newxres*xreal*binw/xres, newyres*yreal*binh/yres, FALSE);
    result->xoff = xoffset + xoff*xreal/xres;
    result->yoff = yoffset + yoff*yreal/yres;
    gwy_data_field_copy_units(data_field, result);

    /* Prevent rounding errors from introducing different values in constants
     * field during resampling. */
    if (data_field_is_constant(data_field, &z)) {
        gwy_data_field_fill(result, z);
        return result;
    }

    d = data_field->data + yoff*xres + xoff;
    r = result->data;
    buf = g_new(gdouble, binsize);
    for (i = 0; i < newyres; i++) {
        for (j = 0; j < newxres; j++) {
            for (k = 0; k < binh; k++) {
                gwy_assign(buf + k*binw, d + (i*binh + k)*xres + j*binw, binw);
            }
            r[i*newxres + j] = gwy_math_trimmed_mean(binsize, buf, trimlowest, trimhighest);
        }
    }
    g_free(buf);

    return result;
}

/**
 * gwy_data_field_bin:
 * @data_field: A data field.
 * @target: Target data field.  It will be resized as necessary.
 * @binw: Bin height (in pixels).
 * @binh: Bin width (in pixels).
 * @xoff: Horizontal offset of bins (in pixels).
 * @yoff: Vertical offset of bins (in pixels).
 * @trimlowest: Number of lowest values to discard.
 * @trimhighest: Number of highest values to discard.
 *
 * Bins a data field into another data field.
 *
 * See gwy_data_field_new_binned() for a detailed description.
 *
 * Since: 2.55
 **/
void
gwy_data_field_bin(GwyDataField *data_field,
                   GwyDataField *target,
                   gint binw, gint binh,
                   gint xoff, gint yoff,
                   gint trimlowest, gint trimhighest)
{
    gint binsize, xres, yres, newxres, newyres, i, j, k;
    gdouble xreal, yreal, xoffset, yoffset, z;
    gdouble *buf, *d, *r;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(target));
    g_return_if_fail(binw > 0);
    g_return_if_fail(binh > 0);
    binsize = binw*binh;
    g_return_if_fail(trimlowest >= 0 && trimlowest < binsize);
    g_return_if_fail(trimhighest >= 0 && trimhighest < binsize - trimlowest);
    if (binsize == 1) {
        gwy_data_field_assign(target, data_field);
        return;
    }

    xoff = ((xoff % binw) + binw) % binw;
    yoff = ((yoff % binh) + binh) % binh;
    xres = data_field->xres;
    yres = data_field->yres;
    xreal = data_field->xreal;
    yreal = data_field->yreal;
    xoffset = data_field->xoff;
    yoffset = data_field->yoff;

    if (xres < binw + xoff || yres < binh + yoff) {
        g_warning("No complete bin can be formed.");
        gwy_data_field_resample(target, 1, 1, GWY_INTERPOLATION_NONE);
        gwy_data_field_copy_units(data_field, target);
        target->xoff = xoffset;
        target->yoff = yoffset;
        target->data[0] = gwy_data_field_get_avg(data_field);
        return;
    }

    newxres = (xres - xoff)/binw;
    newyres = (yres - yoff)/binh;
    gwy_data_field_resample(target, newxres, newyres, GWY_INTERPOLATION_NONE);
    target->xreal = newxres*xreal*binw/xres;
    target->yreal = newyres*yreal*binh/yres;
    target->xoff = xoffset + xoff*xreal/xres;
    target->yoff = yoffset + yoff*yreal/yres;
    gwy_data_field_copy_units(data_field, target);

    /* Prevent rounding errors from introducing different values in constants
     * field during resampling. */
    if (data_field_is_constant(data_field, &z)) {
        gwy_data_field_fill(target, z);
        return;
    }

    d = data_field->data + yoff*xres + xoff;
    r = target->data;
    buf = g_new(gdouble, binsize);
    for (i = 0; i < newyres; i++) {
        for (j = 0; j < newxres; j++) {
            for (k = 0; k < binh; k++) {
                gwy_assign(buf + k*binw, d + (i*binh + k)*xres + j*binw, binw);
            }
            r[i*newxres + j] = gwy_math_trimmed_mean(binsize, buf, trimlowest, trimhighest);
        }
    }
    g_free(buf);
    gwy_data_field_invalidate(target);
}

/**
 * gwy_data_field_get_xder:
 * @data_field: A data field.
 * @col: Column index.
 * @row: Row index.
 *
 * Computes central derivative in X direction.
 *
 * On border points, one-side derivative is returned.
 *
 * Returns: Derivative in X direction.
 **/
gdouble
gwy_data_field_get_xder(GwyDataField *data_field,
                        gint col, gint row)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    g_return_val_if_fail(gwy_data_field_inside(data_field, col, row), 0.0);
    return _gwy_data_field_xder(data_field, col, row);
}

/**
 * gwy_data_field_get_yder:
 * @data_field: A data field.
 * @col: Column index.
 * @row: Row index.
 *
 * Computes central derivative in Y direction.
 *
 * On border points, one-side derivative is returned.
 *
 * Note the derivative is for legacy reasons calulcated for the opposite y direction than is usual elsewhere in
 * Gwyddion, i.e. if values increase with increasing row number, the returned value is negative.
 *
 * Returns: Derivative in Y direction
 **/
gdouble
gwy_data_field_get_yder(GwyDataField *data_field,
                        gint col, gint row)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    g_return_val_if_fail(gwy_data_field_inside(data_field, col, row), 0.0);
    return _gwy_data_field_yder(data_field, col, row);
}

/**
 * gwy_data_field_get_angder:
 * @data_field: A data field.
 * @col: Column index.
 * @row: Row index.
 * @theta: Angle defining the direction (in radians, counterclockwise).
 *
 * Computes derivative in direction specified by given angle.
 *
 * Returns: Derivative in direction given by angle @theta.
 **/
gdouble
gwy_data_field_get_angder(GwyDataField *data_field,
                          gint col, gint row,
                          gdouble theta)
{
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0.0);
    g_return_val_if_fail(gwy_data_field_inside(data_field, col, row), 0.0);

    return _gwy_data_field_xder(data_field, col, row)*cos(theta)
           + _gwy_data_field_yder(data_field, col, row)*sin(theta);
}

/**
 * gwy_data_field_copy_units:
 * @data_field: A data field.
 * @target: Target data field.
 *
 * Sets lateral and value units of a data field to match another data field.
 *
 * Since: 2.49
 **/
void
gwy_data_field_copy_units(GwyDataField *data_field,
                          GwyDataField *target)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(target));

    _gwy_copy_si_unit(data_field->si_unit_xy, &target->si_unit_xy);
    _gwy_copy_si_unit(data_field->si_unit_z, &target->si_unit_z);
}

/**
 * gwy_data_field_copy_units_to_data_line:
 * @data_field: A data field to get units from.
 * @data_line: A data line to set units of.
 *
 * Sets lateral and value units of a data line to match a data field.
 **/
void
gwy_data_field_copy_units_to_data_line(GwyDataField *data_field,
                                       GwyDataLine *data_line)
{
    g_return_if_fail(GWY_IS_DATA_LINE(data_line));
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    _gwy_copy_si_unit(data_field->si_unit_xy, &data_line->si_unit_x);
    _gwy_copy_si_unit(data_field->si_unit_z, &data_line->si_unit_y);
}

gboolean
_gwy_data_field_check_area(GwyDataField *data_field,
                           gint col, gint row,
                           gint width, gint height)
{
    gint xres, yres;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), FALSE);
    xres = data_field->xres;
    yres = data_field->yres;
    g_return_val_if_fail(col >= 0 && col < xres, FALSE);
    g_return_val_if_fail(row >= 0 && row < yres, FALSE);
    g_return_val_if_fail(width > 0 && width <= xres - col, FALSE);
    g_return_val_if_fail(height > 0 && height <= yres - row, FALSE);

    return TRUE;
}

gboolean
_gwy_data_field_check_mask(GwyDataField *data_field,
                           GwyDataField **mask,
                           GwyMaskingType *masking)
{
    /* NULL @mask is a direct caller error.  We allow NULL @masking for old functions that do not have mode so masking
     * is implicitly done in INCLUDE mode. */
    g_assert(mask);
    if (!*mask) {
        if (masking)
            *masking = GWY_MASK_IGNORE;
        return TRUE;
    }
    if (masking) {
        if (*masking == GWY_MASK_IGNORE) {
            *mask = NULL;
            return TRUE;
        }
        g_return_val_if_fail(*masking == GWY_MASK_INCLUDE || *masking == GWY_MASK_EXCLUDE, FALSE);
    }
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), FALSE);
    g_return_val_if_fail((*mask)->xres == data_field->xres, FALSE);
    g_return_val_if_fail((*mask)->yres == data_field->yres, FALSE);
    return TRUE;
}

#undef gwy_data_field_invalidate
void
gwy_data_field_invalidate(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    data_field->cached = 0;
}

static void
fill_missing_points(GwyDataField *dfield, GwyDataField *mask)
{
    gint xres = gwy_data_field_get_xres(dfield);
    gint yres = gwy_data_field_get_yres(dfield);
    MaskedPoint *mpts;
    gdouble *d, *w;
    gint nmissing, k, kk, i, j;

    d = gwy_data_field_get_data(dfield);
    w = gwy_data_field_get_data(mask);

    nmissing = 0;
    for (kk = 0; kk < xres*yres; kk++) {
        if (w[kk])
            nmissing++;
    }
    if (!nmissing)
        return;

    if (nmissing == xres*yres) {
        gwy_data_field_clear(dfield);
        return;
    }

    /* This physically touches the mask data but does not change their
     * interpretation. */
    gwy_data_field_grain_simple_dist_trans(mask, GWY_DISTANCE_TRANSFORM_EUCLIDEAN, FALSE);

    mpts = g_new(MaskedPoint, nmissing);
    k = 0;
    for (kk = 0; kk < xres*yres; kk++) {
        if (w[kk]) {
            mpts[k].dist = w[kk];
            mpts[k].i = kk/xres;
            mpts[k].j = kk % xres;
            k++;
        }
    }
    qsort(mpts, nmissing, sizeof(MaskedPoint), gwy_compare_double);

    for (k = 0; k < nmissing; k++) {
        gdouble z = 0.0, dist = mpts[k].dist;
        gint n = 0;

        i = mpts[k].i;
        j = mpts[k].j;
        kk = i*xres + j;

        /* Cardinal. */
        if (i > 0 && w[kk - xres] < dist) {
            z += d[kk - xres];
            n++;
        }
        if (j > 0 && w[kk-1] < dist) {
            z += d[kk-1];
            n++;
        }
        if (j < xres-1 && w[kk+1] < dist) {
            z += d[kk+1];
            n++;
        }
        if (i < yres-1 && w[kk + xres] < dist) {
            z += d[kk + xres];
            n++;
        }
        z *= 2.0;
        n *= 2;

        /* Diagonal, half weight. */
        if (i > 0 && j > 0 && w[kk-1 - xres] < dist) {
            z += d[kk-1 - xres];
            n++;
        }
        if (i > 0 && j < xres-1 && w[kk+1 - xres] < dist) {
            z += d[kk+1 - xres];
            n++;
        }
        if (i < yres-1 && j > 0 && w[kk-1 + xres] < dist) {
            z += d[kk-1 + xres];
            n++;
        }
        if (i < yres-1 && j < xres-1 && w[kk+1 + xres] < dist) {
            z += d[kk+1 + xres];
            n++;
        }

        g_assert(n);
        d[kk] = z/n;
    }

    g_free(mpts);
}

static void
fill_missing_points_all(GwyDataField *dfield,
                        const GwyXYZ *points,
                        guint npoints)
{
    gdouble xc = dfield->xoff + 0.5*dfield->xreal;
    gdouble yc = dfield->yoff + 0.5*dfield->yreal;
    gdouble d2min = G_MAXDOUBLE;
    guint k, kmin = 0;

    g_return_if_fail(npoints);

    for (k = 0; k < npoints; k++) {
        const GwyXYZ *pt = points + k;
        gdouble x = pt->x - xc;
        gdouble y = pt->y - yc;
        gdouble d2 = x*x + y*y;
        if (d2 < d2min) {
            d2min = d2;
            kmin = k;
        }
    }

    gwy_data_field_fill(dfield, points[kmin].z);
}

/**
 * gwy_data_field_average_xyz:
 * @data_field: A data field to fill with regularised XYZ data.
 * @density_map: Optional data field to fill with XYZ point density map.  It can be %NULL.
 * @points: Array of XYZ points.  Coordinates X and Y represent positions in the plane; the Z-coordinate represents
 *          values.
 * @npoints: Number of points.
 *
 * Fills a data field with regularised XYZ data using a simple method.
 *
 * The real dimensions and offsets of @field determine the rectangle in the XY plane that will be regularised.  The
 * regularisation method is fast but simple and there are no absolute guarantees of quality, even though the result
 * will be usually quite acceptable.
 *
 * This especially applies to reasonable views of the XYZ data.  Unreasonable views can be rendered unreasonably.  In
 * particular if the rectangle does not contain any point from @points (either due to high zoom to an empty region or
 * by just being completely off) @data_field will be filled entirely with the value of the closest point or something
 * similar.
 *
 * Since: 2.44
 **/
void
gwy_data_field_average_xyz(GwyDataField *dfield,
                           GwyDataField *densitymap,
                           const GwyXYZ *points,
                           gint npoints)
{
    GwyDataField *extfield, *extweights;
    gdouble xoff, yoff, qx, qy;
    gint extxres, extyres, xres, yres, k;
    gint imin = G_MAXINT, imax = G_MININT, jmin = G_MAXINT, jmax = G_MININT;
    gdouble *d, *w;
    guint ninside;

    g_return_if_fail(GWY_IS_DATA_FIELD(dfield));
    if (densitymap) {
        g_return_if_fail(GWY_IS_DATA_FIELD(densitymap));
        g_return_if_fail(densitymap->xres == dfield->xres);
        g_return_if_fail(densitymap->yres == dfield->yres);
    }

    if (!points || !npoints) {
        gwy_data_field_clear(dfield);
        if (densitymap)
            gwy_data_field_clear(densitymap);
        return;
    }

    xres = dfield->xres;
    yres = dfield->yres;
    xoff = dfield->xoff;
    yoff = dfield->yoff;
    qx = dfield->xreal/xres;
    qy = dfield->yreal/yres;
    g_return_if_fail(qx > 0.0);
    g_return_if_fail(qy > 0.0);
    gwy_debug("dfield %dx%d", xres, yres);

#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(min:imin,jmin) reduction(max:imax,jmax) \
            private(k) \
            shared(points,npoints,qx,qy,xoff,yoff)
#endif
    for (k = 0; k < npoints; k++) {
        const GwyXYZ *pt = points + k;
        gdouble x = (pt->x - xoff)/qx;
        gdouble y = (pt->y - yoff)/qy;
        gint j = (gint)floor(x);
        gint i = (gint)floor(y);

        if (j < jmin)
            jmin = j;
        if (j > jmax)
            jmax = j;

        if (i < imin)
            imin = i;
        if (i > imax)
            imax = i;
    }

    /* Honour exterior if it is not too far away.  We do not want to construct useless huge data fields for zoom-in
     * scenarios. */
    gwy_debug("true extrange [%d,%d)x[%d,%d)", jmin, jmax, imin, imax);
    imin = CLAMP(imin, -(yres/4 + 16), 0);
    imax = CLAMP(imax, yres-1, yres + yres/4 + 15);
    jmin = CLAMP(jmin, -(xres/4 + 16), 0);
    jmax = CLAMP(jmax, xres-1, xres + xres/4 + 15);
    gwy_debug("extrange [%d,%d)x[%d,%d)", jmin, jmax, imin, imax);

    extxres = jmax+1 - jmin;
    extyres = imax+1 - imin;
    gwy_debug("extfield %dx%d", extxres, extyres);
    extfield = gwy_data_field_new(extxres, extyres, qx*extxres, qy*extyres, TRUE);
    extweights = gwy_data_field_new(extxres, extyres, qx*extxres, qy*extyres, TRUE);
    d = gwy_data_field_get_data(extfield);
    w = gwy_data_field_get_data(extweights);

    for (k = 0; k < npoints; k++) {
        const GwyXYZ *pt = points + k;
        gdouble x = (pt->x - xoff)/qx - jmin;
        gdouble y = (pt->y - yoff)/qy - imin;
        gdouble z = pt->z;
        gint j = (gint)floor(x);
        gint i = (gint)floor(y);
        gdouble xx = x - j;
        gdouble yy = y - i;
        gdouble ww;
        gint kk;

        /* Ensure we are always working in (j,j+1) x (i,i+1) rectangle. */
        if (xx < 0.5) {
            xx += 1.0;
            j--;
        }
        xx -= 0.5;
        if (yy < 0.5) {
            yy += 1.0;
            i--;
        }
        yy -= 0.5;

        kk = i*extxres + j;
        if (j >= 0 && j < extxres && i >= 0 && i < extyres) {
            ww = (1.0 - xx)*(1.0 - yy);
            d[kk] += ww*z;
            w[kk] += ww;
        }
        if (j+1 >= 0 && j+1 < extxres && i >= 0 && i < extyres) {
            ww = xx*(1.0 - yy);
            d[kk+1] += ww*z;
            w[kk+1] += ww;
        }
        if (j >= 0 && j < extxres && i+1 >= 0 && i+1 < extyres) {
            ww = (1.0 - xx)*yy;
            d[kk + extxres] += ww*z;
            w[kk + extxres] += ww;
        }
        if (j+1 >= 0 && j+1 < extxres && i+1 >= 0 && i+1 < extyres) {
            ww = xx*yy;
            d[kk + extxres+1] += ww*z;
            w[kk + extxres+1] += ww;
        }
    }

    if (densitymap)
        gwy_data_field_area_copy(extweights, densitymap, -jmin, -imin, xres, yres, 0, 0);

    ninside = 0;
#ifdef _OPENMP
#pragma omp parallel for if(gwy_threads_are_enabled()) default(none) \
            reduction(+:ninside) \
            private(k) \
            shared(d,w,extxres,extyres)
#endif
    for (k = 0; k < extyres; k++) {
        gint j;
        for (j = 0; j < extxres; j++) {
            gint kk = k*extxres + j;
            if (w[kk]) {
                d[kk] = d[kk]/w[kk];
                w[kk] = 0.0;
                ninside++;
            }
            else
                w[kk] = 1.0;
        }
    }

    gwy_debug("nfilled: %d, nmissing: %d", ninside, extxres*extyres - ninside);
    if (ninside) {
        fill_missing_points(extfield, extweights);
        gwy_data_field_area_copy(extfield, dfield, -jmin, -imin, xres, yres, 0, 0);
    }
    else {
        fill_missing_points_all(dfield, points, npoints);
    }

    g_object_unref(extfield);
    g_object_unref(extweights);
}

/************************** Documentation ****************************/

/**
 * SECTION:datafield
 * @title: GwyDataField
 * @short_description: Two-dimensional data representation
 *
 * #GwyDataField is an object that is used for representation of all two-dimensional data matrices. Most of the basic
 * data handling and processing functions in Gwyddion are declared here as they are connected with #GwyDataField.
 **/

/**
 * GwyDataField:
 *
 * The #GwyDataField struct contains private data only and should be accessed using the functions below.
 **/

/**
 * gwy_data_field_invalidate:
 * @data_field: A data field to invalidate.
 *
 * Invalidates cached data field stats.
 *
 * User code should rarely need this macro, as all #GwyDataField methods do proper invalidation when they change data,
 * as well as gwy_data_field_get_data() does.
 *
 * However, if you get raw data with gwy_data_field_get_data() and then mix direct changes to it with calls to methods
 * like gwy_data_field_get_max(), you may need to explicitely invalidate cached values to let gwy_data_field_get_max()
 * know it has to recompute the maximum.
 **/

/**
 * gwy_data_field_duplicate:
 * @data_field: A data field to duplicate.
 *
 * Convenience macro doing gwy_serializable_duplicate() with all the necessary typecasting.
 *
 * Use gwy_data_field_new_alike() if you don't want to copy data, only resolutions and units.
 **/

/**
 * gwy_data_field_assign:
 * @dest: Target data field.
 * @source: Source data field.
 *
 * Convenience macro making one data field identical to another.
 *
 * This is just a gwy_serializable_clone() wrapper with all the necessary typecasting.
 *
 * Since: 2.52
 **/

/**
 * gwy_data_field_get_xmeasure:
 * @data_field: A data field.
 *
 * Alias for gwy_data_field_get_dx().
 **/

/**
 * gwy_data_field_get_ymeasure:
 * @data_field: A data field.
 *
 * Alias for gwy_data_field_get_dy().
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
