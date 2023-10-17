/*
 *  $Id: gwysiunit.c 25120 2022-10-20 16:09:02Z yeti-dn $
 *  Copyright (C) 2004-2022 David Necas (Yeti), Petr Klapetek.
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
#include <libgwyddion/gwyutils.h>
#include <libgwyddion/gwymath.h>
#include <libgwyddion/gwyenum.h>
#include <libgwyddion/gwydebugobjects.h>
#include <libgwyddion/gwyserializable.h>
#include <libgwyddion/gwysiunit.h>

#define GWY_SI_UNIT_TYPE_NAME "GwySIUnit"

#define unit_index(u,i) g_array_index((u)->units,GwySimpleUnit,(i))

enum {
    VALUE_CHANGED,
    LAST_SIGNAL
};

enum {
    GWY_SI_UNIT_N_BASE = 7,
};

typedef struct {
    GQuark unit;
    gshort power;
    gshort traits;
} GwySimpleUnit;

typedef void (*FormatPowerFunc)(GString *str, gint power);

typedef struct {
    FormatPowerFunc format_power;
    const gchar *multiplier;
    const gchar *unit_times;
    const gchar *unit_division;
    const gchar *power_unit_separator;
} GwySIStyleSpec;

typedef struct {
    guint len;
    const gchar *symbol;
    const gchar *name;
} GwyUnitLongName;

typedef struct {
    const gchar *unit;
    gint powers[GWY_SI_UNIT_N_BASE];
    gdouble factor;
} GwySIUnitDecomposition;

static void         gwy_si_unit_finalize         (GObject *object);
static void         gwy_si_unit_serializable_init(GwySerializableIface *iface);
static GByteArray*  gwy_si_unit_serialize        (GObject *obj,
                                                  GByteArray *buffer);
static gsize        gwy_si_unit_get_size         (GObject *obj);
static GObject*     gwy_si_unit_deserialize      (const guchar *buffer,
                                                  gsize size,
                                                  gsize *position);
static GObject*     gwy_si_unit_duplicate_real   (GObject *object);
static void         gwy_si_unit_clone_real       (GObject *source,
                                                  GObject *copy);
static gboolean     gwy_si_unit_parse            (GwySIUnit *siunit,
                                                  const gchar *string);
static GwySIUnit*   gwy_si_unit_power_real       (GwySIUnit *siunit,
                                                  gint power,
                                                  GwySIUnit *result);
static GwySIUnit*   gwy_si_unit_canonicalize     (GwySIUnit *siunit);
static gboolean     gwy_si_unit_equal_do         (GwySIUnit *siunit1,
                                                  GwySIUnit *siunit2,
                                                  gboolean strict);
static gboolean     gwy_si_unit_equal_direct     (GwySIUnit *siunit1,
                                                  GwySIUnit *siunit2);
static const gchar* gwy_si_unit_prefix           (gint power);
static void         gwy_si_unit_format           (GwySIUnit *siunit,
                                                  const GwySIStyleSpec *fs,
                                                  GwySIValueFormat *vf,
                                                  gboolean magnitude_fixed);
static void         format_power_plain           (GString *string,
                                                  gint power);
static void         format_power_pango           (GString *string,
                                                  gint power);
static void         format_power_TeX             (GString *string,
                                                  gint power);
static void         format_power_unicode         (GString *string,
                                                  gint power);

/* Canonical form must be always first, because this table is used for reverse mapping too */
static const GwyEnum SI_prefixes[] = {
    { "k",     3  },
    { "c",    -2  },
    { "m",    -3  },
    { "M",     6  },
    { "µ",    -6  },
    /* People are extremely creative when it comes to µ replacements... */
    { "μ",    -6  },
    { "~",    -6  },
    { "u",    -6  },
    { "\265", -6  },
    { "G",     9  },
    { "n",    -9  },
    { "T",     12 },
    { "p",    -12 },
    { "P",     15 },
    { "f",    -15 },
    { "E",     18 },
    { "a",    -18 },
    { "Z",     21 },
    { "z",    -21 },
    { "Y",     24 },
    { "y",    -24 },
};

/* TODO: silly units we should probably support specially: kg */

/* Units that can conflict with prefixes */
static const gchar *known_units[] = {
    "deg", "Pa", "cd", "mol", "cal", "px", "pt", "cps", "cts", "Gy", "Gauss"
};

/* XXX: The base unit SI is kilogram.  But there is no way we can generally keep the kilo prefix there, so we
 * decompose to prefixable units -- and prefixable unit is gram.  Yes, this means we are not techncally decomposing to
 * *base* SI units.  */
static const gchar *base_si_units[GWY_SI_UNIT_N_BASE] = {
    "m", "g", "s", "A", "K", "mol", "cd",
};

static const GwySIUnitDecomposition derived_units[] = {
    /*           m   g   s   A   K mol  cd   */
    { "Hz",   {  0,  0, -1,  0,  0,  0,  0 }, 1.0            },
    { "N",    {  1,  1, -2,  0,  0,  0,  0 }, 1.0e3          },
    { "Pa",   { -1,  1, -2,  0,  0,  0,  0 }, 1.0e3          },
    { "J",    {  2,  1, -2,  0,  0,  0,  0 }, 1.0e3          },
    { "eV",   {  2,  1, -2,  0,  0,  0,  0 }, 1.60217653e-16 },
    { "W",    {  2,  1, -3,  0,  0,  0,  0 }, 1.0e3          },
    { "C",    {  0,  0,  1,  1,  0,  0,  0 }, 1.0            },
    { "V",    {  2,  1, -3, -1,  0,  0,  0 }, 1.0e3          },
    { "F",    { -2, -1,  4,  2,  0,  0,  0 }, 1.0e-3         },
    { "H",    {  2,  1, -2, -2,  0,  0,  0 }, 1.0e3          },
    { "Ω",    {  2,  1, -3, -2,  0,  0,  0 }, 1.0e3          },
    { "S",    { -2, -1,  3,  2,  0,  0,  0 }, 1.0e-3         },
    { "T",    {  0,  1, -2, -1,  0,  0,  0 }, 1.0e3          },
    { "Wb",   {  2,  1, -2, -1,  0,  0,  0 }, 1.0e3          },
    { "Bq",   {  0,  0, -1,  0,  0,  0,  0 }, 1.0            },
    { "Sv",   {  2,  0, -2,  0,  0,  0,  0 }, 1.0            },
    { "Gy",   {  2,  0, -2,  0,  0,  0,  0 }, 1.0            },
};

/* Unit formats */
static const GwySIStyleSpec format_style_plain = {
    &format_power_plain, NULL, " ", "/", " "
};
static const GwySIStyleSpec format_style_markup = {
    &format_power_pango, NULL, " ", "/", " "
};
static const GwySIStyleSpec format_style_vfmarkup = {
    &format_power_pango, "× ", " ", "/", " "
};
static const GwySIStyleSpec format_style_TeX = {
    &format_power_TeX, NULL, "\\,", "/", "\\,"
};
static const GwySIStyleSpec format_style_vfTeX = {
    &format_power_TeX, "\\times", "\\,", "/", "\\,"
};
static const GwySIStyleSpec format_style_unicode = {
    &format_power_unicode, NULL, " ", "/", " "
};
static const GwySIStyleSpec format_style_vfunicode = {
    &format_power_unicode, "× ", " ", "/", " "
};

static const GwySIStyleSpec *format_styles[] = {
    NULL,
    &format_style_plain,
    &format_style_markup,
    &format_style_vfmarkup,
    &format_style_TeX,
    &format_style_vfTeX,
    &format_style_unicode,
    &format_style_vfunicode,
};

static guint si_unit_signals[LAST_SIGNAL] = { 0 };

G_DEFINE_TYPE_EXTENDED
    (GwySIUnit, gwy_si_unit, G_TYPE_OBJECT, 0,
     GWY_IMPLEMENT_SERIALIZABLE(gwy_si_unit_serializable_init))

static void
gwy_si_unit_serializable_init(GwySerializableIface *iface)
{
    iface->serialize = gwy_si_unit_serialize;
    iface->deserialize = gwy_si_unit_deserialize;
    iface->get_size = gwy_si_unit_get_size;
    iface->duplicate = gwy_si_unit_duplicate_real;
    iface->clone = gwy_si_unit_clone_real;
}

static void
gwy_si_unit_class_init(GwySIUnitClass *klass)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(klass);

    gobject_class->finalize = gwy_si_unit_finalize;

/**
 * GwySIUnit::value-changed:
 * @gwysiunit: The #GwySIUnit which received the signal.
 *
 * The ::value-changed signal is emitted whenever SI unit changes.
 */
    si_unit_signals[VALUE_CHANGED]
        = g_signal_new("value-changed",
                       G_OBJECT_CLASS_TYPE(gobject_class),
                       G_SIGNAL_RUN_FIRST,
                       G_STRUCT_OFFSET(GwySIUnitClass, value_changed),
                       NULL, NULL,
                       g_cclosure_marshal_VOID__VOID,
                       G_TYPE_NONE, 0);
}

#if 0
static void
debug_print_unit(GArray *units, const gchar *name)
{
    guint i;

    g_printerr("%s: ", name);
    for (i = 0; i < units->len; i++) {
        g_printerr(" %s(%d)",
                   g_quark_to_string(g_array_index(units, GwySimpleUnit, i).unit),
                   g_array_index(units, GwySimpleUnit, i).power);
    }
    g_printerr("\n");
}
#endif

static void
gwy_si_unit_init(GwySIUnit *siunit)
{
    siunit->units = g_array_new(FALSE, FALSE, sizeof(GwySimpleUnit));
}

static void
gwy_si_unit_finalize(GObject *object)
{
    GwySIUnit *si_unit = (GwySIUnit*)object;

    if (si_unit->units)
        g_array_free(si_unit->units, TRUE);

    G_OBJECT_CLASS(gwy_si_unit_parent_class)->finalize(object);
}

static GByteArray*
gwy_si_unit_serialize(GObject *obj,
                      GByteArray *buffer)
{
    GwySIUnit *si_unit;
    GByteArray *retval;

    g_return_val_if_fail(GWY_IS_SI_UNIT(obj), NULL);

    si_unit = GWY_SI_UNIT(obj);
    {
        gchar *unitstr = gwy_si_unit_get_string(si_unit, GWY_SI_UNIT_FORMAT_PLAIN);
        GwySerializeSpec spec[] = {
            { 's', "unitstr", &unitstr, NULL, },
        };
        gwy_debug("unitstr = <%s>", unitstr);
        retval = gwy_serialize_pack_object_struct(buffer, GWY_SI_UNIT_TYPE_NAME, G_N_ELEMENTS(spec), spec);
        g_free(unitstr);
        return retval;
    }
}

static gsize
gwy_si_unit_get_size(GObject *obj)
{
    GwySIUnit *si_unit;
    gsize size;

    g_return_val_if_fail(GWY_IS_SI_UNIT(obj), 0);

    si_unit = GWY_SI_UNIT(obj);
    size = gwy_serialize_get_struct_size(GWY_SI_UNIT_TYPE_NAME, 0, NULL);
    /* Just estimate */
    size += 20*si_unit->units->len;

    return size;
}

static GObject*
gwy_si_unit_deserialize(const guchar *buffer,
                        gsize size,
                        gsize *position)
{
    gchar *unitstr = NULL;
    GwySerializeSpec spec[] = {
        { 's', "unitstr", &unitstr, NULL, },
    };
    GwySIUnit *si_unit;

    g_return_val_if_fail(buffer, NULL);

    if (!gwy_serialize_unpack_object_struct(buffer, size, position, GWY_SI_UNIT_TYPE_NAME, G_N_ELEMENTS(spec), spec))
        return NULL;

    if (unitstr && !*unitstr) {
        g_free(unitstr);
        unitstr = NULL;
    }
    si_unit = gwy_si_unit_new(unitstr);
    g_free(unitstr);

    return (GObject*)si_unit;
}


static GObject*
gwy_si_unit_duplicate_real(GObject *object)
{
    GwySIUnit *si_unit, *duplicate;

    g_return_val_if_fail(GWY_IS_SI_UNIT(object), NULL);
    si_unit = GWY_SI_UNIT(object);
    duplicate = gwy_si_unit_new_parse("", NULL);
    duplicate->power10 = si_unit->power10;
    g_array_append_vals(duplicate->units, si_unit->units->data, si_unit->units->len);

    return (GObject*)duplicate;
}

static void
gwy_si_unit_clone_real(GObject *source, GObject *copy)
{
    GwySIUnit *si_unit, *clone;

    g_return_if_fail(GWY_IS_SI_UNIT(source));
    g_return_if_fail(GWY_IS_SI_UNIT(copy));

    si_unit = GWY_SI_UNIT(source);
    clone = GWY_SI_UNIT(copy);
    if (gwy_si_unit_equal_direct(si_unit, clone))
        return;

    g_array_set_size(clone->units, 0);
    g_array_append_vals(clone->units, si_unit->units->data, si_unit->units->len);
    clone->power10 = si_unit->power10;
    g_signal_emit(copy, si_unit_signals[VALUE_CHANGED], 0);
}

/**
 * gwy_si_unit_new:
 * @unit_string: Unit string (it can be %NULL for an empty unit).
 *
 * Creates a new SI unit from string representation.
 *
 * Unit string represents unit with no prefixes (e. g. "m", "N", "A", etc.)
 *
 * Returns: A new SI unit.
 **/
GwySIUnit*
gwy_si_unit_new(const char *unit_string)
{
    return gwy_si_unit_new_parse(unit_string, NULL);
}

/**
 * gwy_si_unit_new_parse:
 * @unit_string: Unit string (it can be %NULL for an empty unit).
 * @power10: Where power of 10 should be stored (or %NULL).
 *
 * Creates a new SI unit from string representation.
 *
 * This is a more powerful version of gwy_si_unit_new(): @unit_string may be a relatively complex unit, with prefixes,
 * like "pA/s" or "km^2". Beside conversion to a base SI unit like "A/s" or "m^2" it also computes the power of 10 one
 * has to multiply the base unit with to get an equivalent of @unit_string.
 *
 * For example, for <literal>"pA/s"</literal> it will store -12 to @power10 because 1 pA/s is 1e-12 A/s, for
 * <literal>"km^2"</literal> it will store 6 to @power10 because 1 km^2 is 1e6 m^2.
 *
 * Returns: A new SI unit.
 **/
GwySIUnit*
gwy_si_unit_new_parse(const char *unit_string,
                      gint *power10)
{
    GwySIUnit *siunit;

    siunit = g_object_new(GWY_TYPE_SI_UNIT, NULL);
    gwy_si_unit_parse(siunit, unit_string);
    if (power10)
        *power10 = siunit->power10;

    return siunit;
}

/**
 * gwy_si_unit_set_from_string:
 * @siunit: An SI unit.
 * @unit_string: Unit string to set @siunit from (it can be %NULL for an empty unit).
 *
 * Sets string that represents unit.
 *
 * It must be base unit with no prefixes (e. g. "m", "N", "A", etc.).
 **/
void
gwy_si_unit_set_from_string(GwySIUnit *siunit,
                            const gchar *unit_string)
{
    g_return_if_fail(GWY_IS_SI_UNIT(siunit));
    gwy_si_unit_set_from_string_parse(siunit, unit_string, NULL);
}

/**
 * gwy_si_unit_set_from_string_parse:
 * @siunit: An SI unit.
 * @unit_string: Unit string to set @siunit from (it can be %NULL for an empty
 *               unit).
 * @power10: Where power of 10 should be stored (or %NULL).
 *
 * Changes an SI unit according to string representation.
 *
 * This is a more powerful version of gwy_si_unit_set_from_string(), please see gwy_si_unit_new_parse() for some
 * discussion.
 **/
void
gwy_si_unit_set_from_string_parse(GwySIUnit *siunit,
                                  const gchar *unit_string,
                                  gint *power10)
{
    g_return_if_fail(GWY_IS_SI_UNIT(siunit));

    gwy_si_unit_parse(siunit, unit_string);
    if (power10)
        *power10 = siunit->power10;

    g_signal_emit(siunit, si_unit_signals[VALUE_CHANGED], 0);
}

static inline const GwySIStyleSpec*
gwy_si_unit_find_style_spec(GwySIUnitFormatStyle style)
{
    if (style == GWY_SI_UNIT_FORMAT_NONE || (guint)style >= G_N_ELEMENTS(format_styles)) {
        g_critical("Invalid format style %d.", style);
        style = GWY_SI_UNIT_FORMAT_PLAIN;
    }

    return format_styles[style];
}

/**
 * gwy_si_unit_get_string:
 * @siunit: An SI unit.
 * @style: Unit format style.
 *
 * Obtains string representing a SI unit.
 *
 * Returns: A newly allocated string that represents the base unit (with no prefixes).
 **/
gchar*
gwy_si_unit_get_string(GwySIUnit *siunit,
                       GwySIUnitFormatStyle style)
{
    GwySIValueFormat vf;
    gchar *s;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);
    siunit->power10 = 0;
    vf.magnitude = 1.0;
    vf.precision = 3;
    vf.units_gstring = g_string_new(NULL);

    gwy_si_unit_format(siunit, gwy_si_unit_find_style_spec(style), &vf, TRUE);
    s = vf.units_gstring->str;
    g_string_free(vf.units_gstring, FALSE);

    return s;
}

/* Be more adverse to powers of 10 when we have plain numbers. */
static void
unpower_format_for_plain_numbers(GwySIValueFormat *vf, gdouble value,
                                 GwySIUnit *siunit)
{
    if (siunit->units->len || vf->magnitude == 1.0)
        return;

    if (vf->magnitude == 1e-3 && value >= 1e-3) {
        vf->magnitude = 1.0;
        vf->precision += 3;
    }
}

/**
 * gwy_si_unit_get_format_for_power10:
 * @siunit: An SI unit.
 * @style: Unit format style.
 * @power10: Power of 10, in the same sense as gwy_si_unit_new_parse() returns it.
 * @format: A value format to set-up, may be %NULL, a new value format is allocated then.
 *
 * Finds format for representing a specific power-of-10 multiple of a unit.
 *
 * The values should be then printed as value/@format->magnitude [@format->units] with @format->precision decimal
 * places.
 *
 * This function does not change the precision field of @format.
 *
 * Returns: The value format.  If @format was %NULL, a newly allocated format is returned, otherwise (modified)
 *          @format itself is returned.
 **/
GwySIValueFormat*
gwy_si_unit_get_format_for_power10(GwySIUnit *siunit,
                                   GwySIUnitFormatStyle style,
                                   gint power10,
                                   GwySIValueFormat *format)
{
    const GwySIStyleSpec *spec;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);

    spec = gwy_si_unit_find_style_spec(style);
    if (!format) {
        format = (GwySIValueFormat*)g_new0(GwySIValueFormat, 1);
        format->units_gstring = g_string_new(NULL);
    }

    siunit->power10 = power10;
    format->magnitude = pow10(power10);
    gwy_si_unit_format(siunit, spec, format, TRUE);

    return format;
}

/**
 * gwy_si_unit_get_format:
 * @siunit: An SI unit.
 * @style: Unit format style.
 * @value: Value the format should be suitable for.
 * @format: A value format to set-up, may be %NULL, a new value format is allocated then.
 *
 * Finds a good format for representing a value.
 *
 * The values should be then printed as value/@format->magnitude [@format->units] with @format->precision decimal
 * places.
 *
 * Returns: The value format.  If @format was %NULL, a newly allocated format is returned, otherwise (modified)
 *          @format itself is returned.
 **/
GwySIValueFormat*
gwy_si_unit_get_format(GwySIUnit *siunit,
                       GwySIUnitFormatStyle style,
                       gdouble value,
                       GwySIValueFormat *format)
{
    const GwySIStyleSpec *spec;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);

    spec = gwy_si_unit_find_style_spec(style);
    if (!format) {
        format = (GwySIValueFormat*)g_new0(GwySIValueFormat, 1);
        format->units_gstring = g_string_new(NULL);
    }

    value = fabs(value);
    if (!value) {
        format->magnitude = 1;
        format->precision = 2;
    }
    else {
        format->magnitude = gwy_math_humanize_numbers(value/36, value, &format->precision);
        unpower_format_for_plain_numbers(format, value, siunit);
    }
    siunit->power10 = GWY_ROUND(log10(format->magnitude));
    gwy_si_unit_format(siunit, spec, format, FALSE);

    return format;
}

/**
 * gwy_si_unit_get_format_with_resolution:
 * @siunit: An SI unit.
 * @style: Unit format style.
 * @maximum: The maximum value to be represented.
 * @resolution: The smallest step (approximately) that should make a visible difference in the representation.
 * @format: A value format to set-up, may be %NULL, a new value format is allocated then.
 *
 * Finds a good format for representing a range of values with given resolution.
 *
 * The values should be then printed as value/@format->magnitude [@format->units] with @format->precision decimal
 * places.
 *
 * Returns: The value format.  If @format was %NULL, a newly allocated format is returned, otherwise (modified)
 *          @format itself is returned.
 **/
GwySIValueFormat*
gwy_si_unit_get_format_with_resolution(GwySIUnit *siunit,
                                       GwySIUnitFormatStyle style,
                                       gdouble maximum,
                                       gdouble resolution,
                                       GwySIValueFormat *format)
{
    const GwySIStyleSpec *spec;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);

    spec = gwy_si_unit_find_style_spec(style);
    if (!format) {
        format = (GwySIValueFormat*)g_new0(GwySIValueFormat, 1);
        format->units_gstring = g_string_new(NULL);
    }

    maximum = fabs(maximum);
    resolution = fabs(resolution);
    if (!maximum) {
        format->magnitude = 1;
        format->precision = 2;
    }
    else {
        format->magnitude = gwy_math_humanize_numbers(resolution, maximum, &format->precision);
        unpower_format_for_plain_numbers(format, maximum, siunit);
    }
    siunit->power10 = GWY_ROUND(log10(format->magnitude));
    gwy_si_unit_format(siunit, spec, format, FALSE);

    return format;
}

/**
 * gwy_si_unit_get_format_with_digits:
 * @siunit: An SI unit.
 * @style: Unit format style.
 * @maximum: The maximum value to be represented.
 * @sdigits: The number of significant digits the value should have.
 * @format: A value format to set-up, may be %NULL, a new value format is allocated then.
 *
 * Finds a good format for representing a values with given number of significant digits.
 *
 * The values should be then printed as value/@format->magnitude [@format->units] with @format->precision decimal
 * places.
 *
 * Returns: The value format.  If @format was %NULL, a newly allocated format is returned, otherwise (modified)
 *          @format itself is returned.
 **/
GwySIValueFormat*
gwy_si_unit_get_format_with_digits(GwySIUnit *siunit,
                                   GwySIUnitFormatStyle style,
                                   gdouble maximum,
                                   gint sdigits,
                                   GwySIValueFormat *format)
{
    const GwySIStyleSpec *spec;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);

    spec = gwy_si_unit_find_style_spec(style);
    if (!format) {
        format = (GwySIValueFormat*)g_new0(GwySIValueFormat, 1);
        format->units_gstring = g_string_new(NULL);
    }

    maximum = fabs(maximum);
    if (!maximum) {
        format->magnitude = 1;
        format->precision = sdigits;
    }
    else {
        format->magnitude = gwy_math_humanize_numbers(maximum/pow10(sdigits), maximum, &format->precision);
        unpower_format_for_plain_numbers(format, maximum, siunit);
    }
    siunit->power10 = GWY_ROUND(log10(format->magnitude));
    gwy_si_unit_format(siunit, spec, format, FALSE);

    return format;
}

/**
 * gwy_si_unit_equal:
 * @siunit: First unit.
 * @siunit2: Second unit.
 *
 * Checks whether two SI units are equal.
 *
 * Returns: %TRUE if the units are equal.
 **/
gboolean
gwy_si_unit_equal(GwySIUnit *siunit1, GwySIUnit *siunit2)
{
    /* Here we do not want eV to be equal to J. */
    return gwy_si_unit_equal_do(siunit1, siunit2, TRUE);
}

/**
 * gwy_si_unit_equal_string:
 * @siunit: An SI unit.
 * @unit_string: Unit string (it can be %NULL for an empty unit).
 *
 * Checks whether an SI unit corresponds to given string.
 *
 * Any power-of-ten prefixes are ignored.  This function is mostly useful for quick commensurability checks with
 * simple units such as "m" and for checking whether a unit is non-empty (by comparing with %NULL or an empty string).
 *
 * Returns: %TRUE if the units is equivalent to the given string.
 *
 * Since: 2.49
 **/
gboolean
gwy_si_unit_equal_string(GwySIUnit *siunit,
                         const gchar *unit_string)
{
    GwySIUnit *tmpunit;
    gboolean retval;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), FALSE);

    if (!unit_string || !*unit_string)
        return !siunit->units->len;

    if (siunit->units->len == 1) {
        GwySimpleUnit *unit = &unit_index(siunit, 0);
        if (gwy_strequal(g_quark_to_string(unit->unit), unit_string))
            return TRUE;
        /* If they differ we cannot make any conclusion. */
    }

    tmpunit = gwy_si_unit_new_parse(unit_string, NULL);
    retval = gwy_si_unit_equal_do(siunit, tmpunit, FALSE);
    g_object_unref(tmpunit);

    return retval;
}

static const GwyUnitLongName*
check_for_long_name(const gchar *s, guint len,
                    const GwyUnitLongName *longnames, guint n,
                    gboolean exact_match)
{
    guint i;

    for (i = 0; i < n; i++) {
        const gchar *name = longnames[i].name;
        guint ll = longnames[i].len;

        if (ll > len)
            return NULL;

        if (exact_match) {
            if (ll != len)
                continue;
            if (memcmp(s, name, len) == 0)
                return longnames + i;
        }
        else {
            if (g_ascii_strncasecmp(s, name, ll) == 0)
                return longnames + i;
        }
    }

    return NULL;
}

/* fix all kinds of sloppy and strange notations */
static void
fix_unit_name(GString *str)
{
    /* Keep the lists sorted by length so that we can give up quickly and only attempt to translate long names. */

    /* Special silly names which cannot take prefixes and we only accept as exact matches. */
    static const GwyUnitLongName odd_names[] = {
        { 1, "Å",   "\xc5",    },
        { 1, "deg", "\xb0",    },
        { 1, "deg", "\xba",    },
        { 2, "deg", "°",       },
        { 2, "Å",   "AA",      },
        { 2, "Å",   "Å",       },
        { 3, "Å",   "Ang",     },
        { 3, "Å",   "ang",     },
        { 3, "µm",  "mum",     },
        { 4, "",    "a.u.",    },
        { 5, "",    "a. u.",   },
        { 6, "",    "counts",  },
        { 6, "µm",  "micron",  },
        { 7, "µm",  "microns", },
        { 7, "µm",  "micro m", },
    };

    /* Long names and their prefixes. */
    static const GwyUnitLongName long_names[] = {
        { 3, "s",   "sec",      },
        { 3, "deg", "Deg",      },
        { 3, "Ω",   "Ohm",      },
        { 4, "V",   "Volt",     },
        { 4, "W",   "Watt",     },
        { 5, "Hz",  "Hertz",    },
        { 5, "F",   "Farad",    },
        { 5, "H",   "Henry",    },
        { 5, "J",   "Joule",    },
        { 5, "m",   "meter",    },
        { 5, "m",   "metre",    },
        { 5, "px",  "pixel",    },
        { 5, "T",   "Tesla",    },
        { 5, "Wb",  "Weber",    },
        { 6, "A",   "Ampere",   },
        { 6, "deg", "degree",   },
        { 6, "K",   "Kelvin",   },
        { 6, "N",   "Newton",   },
        { 6, "Pa",  "Pascal",   },
        { 6, "rad", "radian",   },
        { 6, "s",   "second",   },
        { 7, "cd",  "candela",  },
        { 7, "S",   "Siemens",  },
        { 8, "Å",   "Angstrom", },
    };

    static const GwyUnitLongName long_pfxes[] = {
        { 4, "a", "atto", },
        { 4, "G", "giga", },
        { 4, "k", "kilo", },
        { 4, "M", "mega", },
        { 4, "n", "nano", },
        { 4, "P", "peta", },
        { 4, "p", "pico", },
        { 4, "T", "tera", },
        { 5, "c", "centi", },
        { 5, "f", "femto", },
        { 5, "µ", "micro", },
        { 5, "m", "milli", },
    };

    const GwyUnitLongName *xname;
    const gchar *s = str->str, *prefix = "";
    guint len = str->len;

    if (!len)
        return;

    /* b0 = degree sign
     * ba = masculine ordinal indicator (yes, people use it for degrees) */
    if (len > 1 && (s[0] == '\xb0' || s[0] == '\xba')) {
        g_string_erase(str, 0, 1);
        g_string_prepend(str, "°");
        return;
    }

    /* Fix silly names. */
    if ((xname = check_for_long_name(s, len, odd_names, G_N_ELEMENTS(odd_names), TRUE))) {
        g_string_assign(str, xname->symbol);
        return;
    }

    /* Fix long names, possibly prefixed. */
    if ((xname = check_for_long_name(s, len, long_pfxes, G_N_ELEMENTS(long_pfxes), FALSE))) {
        prefix = xname->symbol;
        s += xname->len;
        len -= xname->len;
    }
    if ((xname = check_for_long_name(s, len, long_names, G_N_ELEMENTS(long_names), FALSE))
        && (xname->len == len || (xname->len+1 == len && g_ascii_tolower(s[len-1]) == 's'))) {
        g_string_assign(str, prefix);
        g_string_append(str, xname->symbol);
    }
}

static gboolean
parse_exponent_unicode(const gchar *s, gint *n, const gchar **end)
{
    static const guint lengths[10] = { 3, 2, 2, 2, 3, 3, 3, 3, 3, 3 };
    static const gchar digits[] = "⁰¹²³⁴⁵⁶⁷⁸⁹";
    const gchar *p = s;
    guint k, pos, diglen;
    gint i = 0, sign = 1;

    if (strncmp(p, "⁻", 3) == 0) {
        sign = -1;
        p += 3;
        /* Do not count standalone minus as a number. */
        s = p;
    }

    while (*p) {
        for (k = pos = 0; k < 10; k++) {
            diglen = lengths[k];
            if (strncmp(p, digits + pos, diglen) == 0) {
                i = 10*i + k;
                p += diglen;
                break;
            }
            pos += diglen;
        }
        if (k == 10)
            break;
    }

    if (p == s)
        return FALSE;

    *n = i*sign;
    if (end)
        *end = p;
    return TRUE;
}

static const gchar*
parse_numerical_multiplier(GwySIUnit *siunit, const gchar *string, gint sign)
{
    const gchar *end;
    gdouble q;
    gint n, power10 = 0;

    q = g_ascii_strtod(string, (gchar**)&end);
    /* Unfortunately, g_ascii_strtod() can parses Nanometre as NaN + ometre.
     * Accept only known numbers as multipliers. */
    if (end != string && !gwy_isinf(q) && !gwy_isnan(q)) {
        string = end;
        power10 = GWY_ROUND(log10(q));
        if (q <= 0 || fabs(log(q/pow10(power10))) > 1e-13) {
            /* g_warning("Bad multiplier %g", q); */
            power10 = 0;
        }
        else if (g_str_has_prefix(string, "<sup>")) {
            string += strlen("<sup>");
            n = strtol(string, (gchar**)&end, 10);
            if (end == string) {
                /* g_warning("Bad exponent %s", string); */
            }
            else if (!g_str_has_prefix(end, "</sup>")) {
                /* g_warning("Expected </sup> after exponent"); */
            }
            else
                power10 *= n;
            string = end;
        }
        else if (string[0] == '^') {
            string++;
            n = strtol(string, (gchar**)&end, 10);
            if (end == string) {
                /* g_warning("Bad exponent %s", string); */
            }
            else
                power10 *= n;
            string = end;
        }
        else if (parse_exponent_unicode(string, &n, &end)) {
            power10 *= n;
            string = end;
        }
    }

    siunit->power10 += sign*power10;
    return string;
}

static gboolean
gwy_si_unit_parse(GwySIUnit *siunit,
                  const gchar *string)
{
    GwySimpleUnit unit;
    const gchar *end;
    gchar *p, *e;
    gint i, j, pfpower;
    GString *buf;
    gboolean dividing = FALSE, beginning = TRUE, may_split_prefix;

    g_array_set_size(siunit->units, 0);
    siunit->power10 = 0;

    if (!string || !*string)
        return TRUE;

    /* give up when it looks too wild */
    end = strpbrk(string,
                  "\177\001\002\003\004\005\006\007\010\011\012\013\014\015\016\017"
                  "\020\021\022\023\024\025\026\027\030\031\032\033\034\035\036\037"
                  "!#$&()*,:;=?@\\[]_`|{}");
    if (end) {
        /* g_warning("Invalid character 0x%02x", *end); */
        return FALSE;
    }

    /* may start with a multiplier, but it must be a power of 10 */
    while (g_ascii_isspace(string[0]))
        string++;

    if (string[0] == '*')
        string++;
    else if (strncmp(string, "×", sizeof("×")-1) == 0)
        string += sizeof("×")-1;

    buf = g_string_new(NULL);

    /* the rest are units */
    while (*string) {
        /* here we actually parse the number, also after a division sign */
        if (beginning || dividing) {
            string = parse_numerical_multiplier(siunit, string, dividing ? -1 : 1);
            while (g_ascii_isspace(*string))
                string++;
        }

        /* units are separated with whitespace and maybe a division sign */
        end = string;
        do {
            end = strpbrk(end, " /");
            if (!end || end == string || *end != '/' || *(end-1) != '<')
                break;
            end++;
        } while (TRUE);
        if (!end)
            end = string + strlen(string);

        g_string_set_size(buf, 0);
        g_string_append_len(buf, string, end - string);
        fix_unit_name(buf);

        /* get prefix, but be careful not to split mol to mili-ol */
        pfpower = 0;
        may_split_prefix = buf->len > 1;
        if (may_split_prefix) {
            for (i = 0; i < G_N_ELEMENTS(known_units); i++) {
                if (g_str_has_prefix(buf->str, known_units[i]) && !g_ascii_isalpha(buf->str[strlen(known_units[i])])) {
                    may_split_prefix = FALSE;
                    break;
                }
            }
        }
        /* also don't split prefixes of long words, they are unlikely to be symbols. */
        if (may_split_prefix && buf->len > 4) {
           for (i = 0; i < buf->len; i++) {
               if (!g_ascii_isalpha(buf->str[i]))
                   break;
           }
           if (i == buf->len)
               may_split_prefix = FALSE;
        }

        if (may_split_prefix) {
            for (i = 0; i < G_N_ELEMENTS(SI_prefixes); i++) {
                const gchar *pfx = SI_prefixes[i].name;

                if (g_str_has_prefix(buf->str, pfx) && g_ascii_isalpha(buf->str[strlen(pfx)])) {
                    pfpower = SI_prefixes[i].value;
                    g_string_erase(buf, 0, strlen(pfx));
                    break;
                }
            }
        }

        /* get unit power */
        unit.power = 1;
        if (buf->len && (p = strstr(buf->str + 1, "<sup>"))) {
            unit.power = strtol(p + strlen("<sup>"), &e, 10);
            if (e == p + strlen("<sup>") || !g_str_has_prefix(e, "</sup>")) {
                /* g_warning("Bad power %s", p); */
                unit.power = 1;
            }
            else if (!unit.power || abs(unit.power) > 12) {
                /* g_warning("Bad power %d", unit.power); */
                unit.power = 1;
            }
            g_string_truncate(buf, p - buf->str);
        }
        else if (buf->len && (p = strchr(buf->str + 1, '^'))) {
            unit.power = strtol(p + 1, &e, 10);
            if (e == p + 1 || *e) {
                /* g_warning("Bad power %s", p); */
                unit.power = 1;
            }
            else if (!unit.power || abs(unit.power) > 12) {
                /* g_warning("Bad power %d", unit.power); */
                unit.power = 1;
            }
            g_string_truncate(buf, p - buf->str);
        }
        else if (buf->len) {
            /* Try to find a Unicode exponent first by looking for a non-letter character. */
            i = 1;
            while (g_ascii_isalpha(buf->str[i]))
                i++;
            if (parse_exponent_unicode(buf->str + i, &j, (const gchar**)&e)) {
                unit.power = j;
                g_string_truncate(buf, i);
            }
            else {
                /* Are we really desperate?  Yes, we are! */
                i = buf->len;
                while (i && (g_ascii_isdigit(buf->str[i-1]) || buf->str[i-1] == '-'))
                    i--;
                if (i != buf->len) {
                    unit.power = strtol(buf->str + i, NULL, 10);
                    if (!unit.power || abs(unit.power) > 12) {
                        /* g_warning("Bad power %d", unit.power); */
                        unit.power = 1;
                    }
                    g_string_truncate(buf, i);
                }
            }
        }

        /* handle some ugly, but quite common units */
        if (gwy_strequal(buf->str, "Å")) {
            pfpower -= 10;
            g_string_assign(buf, "m");
        }
        else if (gwy_strequal(buf->str, "%")) {
            pfpower -= 2;
            g_string_assign(buf, "");
        }
        else if (gwy_strequal(buf->str, "м")) {
            g_string_assign(buf, "m");
        }

        /* elementary sanity */
        if (!g_utf8_validate(buf->str, -1, (const gchar**)&p)) {
            /* g_warning("Unit string is not valid UTF-8"); */
            g_string_truncate(buf, p - buf->str);
        }
        if (!buf->len) {
            /* maybe it's just percentage.  cross fingers and proceed. */
            if (dividing)
                unit.power = -unit.power;
            siunit->power10 += unit.power * pfpower;
        }
        else if (!g_ascii_isalpha(buf->str[0]) && (guchar)buf->str[0] < 128) {
            /* g_warning("Invalid base unit: %s", buf->str); */
        }
        else {
            /* append it */
            unit.unit = g_quark_from_string(buf->str);
            if (dividing)
                unit.power = -unit.power;
            gwy_debug("<%s:%u> %d\n", buf->str, unit.unit, unit.power);
            siunit->power10 += unit.power * pfpower;
            g_array_append_val(siunit->units, unit);
        }

        /* TODO: scan known obscure units */
        unit.traits = 0;

        /* get to the next token, looking for division */
        while (g_ascii_isspace(*end))
            end++;
        if (*end == '/') {
            if (dividing) {
                /* g_warning("Cannot group multiple divisions"); */
            }
            dividing = TRUE;
            end++;
            while (g_ascii_isspace(*end))
                end++;
        }
        string = end;
        beginning = FALSE;
    }

    g_string_free(buf, TRUE);

    gwy_si_unit_canonicalize(siunit);

    return TRUE;
}

/**
 * gwy_si_unit_multiply:
 * @siunit1: An SI unit.
 * @siunit2: An SI unit.
 * @result: An SI unit to set to product of @siunit1 and @siunit2.  It is safe to pass one of @siunit1, @siunit2. It
 *          can be %NULL too, a new SI unit is created then and returned.
 *
 * Multiplies two SI units.
 *
 * Returns: When @result is %NULL, a newly created SI unit that has to be dereferenced when no longer used later.
 *          Otherwise @result itself is simply returned, its reference count is NOT increased.
 **/
GwySIUnit*
gwy_si_unit_multiply(GwySIUnit *siunit1,
                     GwySIUnit *siunit2,
                     GwySIUnit *result)
{
    return gwy_si_unit_power_multiply(siunit1, 1, siunit2, 1, result);
}

/**
 * gwy_si_unit_divide:
 * @siunit1: An SI unit.
 * @siunit2: An SI unit.
 * @result: An SI unit to set to quotient of @siunit1 and @siunit2.  It is safe to pass one of @siunit1, @siunit2. It
 *          can be %NULL too, a new SI unit is created then and returned.
 *
 * Divides two SI units.
 *
 * Returns: When @result is %NULL, a newly created SI unit that has to be dereferenced when no longer used later.
 *          Otherwise @result itself is simply returned, its reference count is NOT increased.
 **/
GwySIUnit*
gwy_si_unit_divide(GwySIUnit *siunit1,
                   GwySIUnit *siunit2,
                   GwySIUnit *result)
{
    return gwy_si_unit_power_multiply(siunit1, 1, siunit2, -1, result);
}

/**
 * gwy_si_unit_power:
 * @siunit: An SI unit.
 * @power: Power to raise @siunit to.
 * @result: An SI unit to set to power of @siunit.  It is safe to pass @siunit itself.  It can be %NULL too, a new SI
 *          unit is created then and returned.
 *
 * Computes a power of an SI unit.
 *
 * Returns: When @result is %NULL, a newly created SI unit that has to be dereferenced when no longer used later.
 *          Otherwise @result itself is simply returned, its reference count is NOT increased.
 **/
GwySIUnit*
gwy_si_unit_power(GwySIUnit *siunit,
                  gint power,
                  GwySIUnit *result)
{
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);
    g_return_val_if_fail(!result || GWY_IS_SI_UNIT(result), NULL);

    if (!result)
        result = gwy_si_unit_new(NULL);

    gwy_si_unit_power_real(siunit, power, result);
    g_signal_emit(result, si_unit_signals[VALUE_CHANGED], 0);

    return result;
}

static GwySIUnit*
gwy_si_unit_power_real(GwySIUnit *siunit,
                       gint power,
                       GwySIUnit *result)
{
    GArray *units;
    GwySimpleUnit *unit;
    gint j;

    units = g_array_new(FALSE, FALSE, sizeof(GwySimpleUnit));
    result->power10 = power*siunit->power10;

    if (power) {
        g_array_append_vals(units, siunit->units->data, siunit->units->len);
        for (j = 0; j < units->len; j++) {
            unit = &g_array_index(units, GwySimpleUnit, j);
            unit->power *= power;
        }
    }

    g_array_set_size(result->units, 0);
    g_array_append_vals(result->units, units->data, units->len);
    g_array_free(units, TRUE);

    return result;
}

/**
 * gwy_si_unit_nth_root:
 * @siunit: An SI unit.
 * @ipower: The root to take: 2 means a quadratic root, 3 means cubic root, etc.
 * @result: An SI unit to set to power of @siunit.  It is safe to pass @siunit itself.  It can be %NULL too, a new SI
 *          unit is created then and returned.
 *
 * Calulates n-th root of an SI unit.
 *
 * This operation fails if the result would have fractional powers that are not representable by #GwySIUnit.
 *
 * Returns: On success: When @result is %NULL, a newly created SI unit that has to be dereferenced when no longer used
 *          later, otherwise @result itself is simply returned, its reference count is NOT increased. On failure %NULL
 *          is always returned.
 *
 * Since: 2.5
 **/
GwySIUnit*
gwy_si_unit_nth_root(GwySIUnit *siunit,
                     gint ipower,
                     GwySIUnit *result)
{
    GArray *units;
    GwySimpleUnit *unit;
    gint j;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);
    g_return_val_if_fail(!result || GWY_IS_SI_UNIT(result), NULL);
    g_return_val_if_fail(ipower > 0, NULL);

    /* Check applicability */
    for (j = 0; j < siunit->units->len; j++) {
        unit = &unit_index(siunit, j);
        if (unit->power % ipower != 0)
            return NULL;
    }

    if (!result)
        result = gwy_si_unit_new(NULL);

    /* XXX: Applicability not required */
    result->power10 = siunit->power10/ipower;

    units = g_array_new(FALSE, FALSE, sizeof(GwySimpleUnit));
    g_array_append_vals(units, siunit->units->data, siunit->units->len);
    for (j = 0; j < units->len; j++) {
        unit = &g_array_index(units, GwySimpleUnit, j);
        unit->power /= ipower;
    }

    g_array_set_size(result->units, 0);
    g_array_append_vals(result->units, units->data, units->len);
    g_array_free(units, TRUE);

    g_signal_emit(result, si_unit_signals[VALUE_CHANGED], 0);

    return result;
}

/**
 * gwy_si_unit_power_multiply:
 * @siunit1: An SI unit.
 * @power1: Power to raise @siunit1 to.
 * @siunit2: An SI unit.
 * @power2: Power to raise @siunit2 to.
 * @result: An SI unit to set to @siunit1^@power1*@siunit2^@power2. It is safe to pass @siunit1 or @siunit2.  It can
 *          be %NULL too, a new SI unit is created then and returned.
 *
 * Computes the product of two SI units raised to arbitrary powers.
 *
 * This is the most complex SI unit arithmetic function.  It can be easily chained when more than two units are to be
 * multiplied.
 *
 * Returns: When @result is %NULL, a newly created SI unit that has to be dereferenced when no longer used later.
 *          Otherwise @result itself is simply returned, its reference count is NOT increased.
 *
 * Since: 2.4
 **/
GwySIUnit*
gwy_si_unit_power_multiply(GwySIUnit *siunit1,
                           gint power1,
                           GwySIUnit *siunit2,
                           gint power2,
                           GwySIUnit *result)
{
    GwySIUnit *op2 = NULL;
    GwySimpleUnit *unit, *unit2;
    gint i, j;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit1), NULL);
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit2), NULL);
    g_return_val_if_fail(!result || GWY_IS_SI_UNIT(result), NULL);

    if (!result)
        result = gwy_si_unit_new(NULL);

    /* Try to avoid hard work by making siunit2 the simplier one */
    if ((!siunit1->units->len && siunit2->units->len)
        || (!power1 && power2)
        || (siunit2 == result && siunit1 != result)) {
        GWY_SWAP(GwySIUnit*, siunit1, siunit2);
        GWY_SWAP(gint, power1, power2);
    }
    if (!power2 || !siunit2->units->len) {
        /* We can do this only if we won't use siunit2 for anything as it
         * can be the same object as result. */
        gwy_si_unit_power_real(siunit1, power1, result);
        gwy_si_unit_canonicalize(result);
        return result;
    }

    /* Operate on a temporary copy in the general case to ensure siunit2 and
     * result are different objects.*/
    if (siunit2 == result) {
        op2 = gwy_si_unit_duplicate(siunit2);
        siunit2 = op2;
    }
    gwy_si_unit_power_real(siunit1, power1, result);

    result->power10 += power2*siunit2->power10;
    for (i = 0; i < siunit2->units->len; i++) {
        unit2 = &unit_index(siunit2, i);

        for (j = 0; j < result->units->len; j++) {
            unit = &unit_index(result, j);
            gwy_debug("[%d] %u == [%d] %u", i, unit2->unit, j, unit->unit);
            if (unit2->unit == unit->unit) {
                unit->power += power2*unit2->power;
                break;
            }
        }
        if (j == result->units->len) {
            g_array_append_val(result->units, *unit2);
            unit = &g_array_index(result->units, GwySimpleUnit, result->units->len - 1);
            unit->power *= power2;
        }
    }
    gwy_si_unit_canonicalize(result);
    GWY_OBJECT_UNREF(op2);
    g_signal_emit(result, si_unit_signals[VALUE_CHANGED], 0);

    return result;
}

/**
 * gwy_si_unit_factor_to_base:
 * @siunit: An SI unit.
 * @result: An SI unit to set to decomposed @siunit.  It is safe to pass @siunit itself.  It can be %NULL too, a new
 *          SI unit is created then and returned.
 * @mfactor: Location to store multiplicative factor between @siunit and @result.  For instance, for electronvolt the
 *           value 1.60217653e-16 would be stored (the factor 1000 comes from kilogram).  Pass %NULL if you are only
 *           interested dimension equality.
 *
 * Factors a possibly derived SI unit to base units.
 *
 * For instance, if @siunit was set to "N/m" the result will be "kg/s^2".
 *
 * Normally the result will consist only of the base seven SI units.  However, recognised non-SI units (and
 * pseudounits) in @siunit, such as "px" are left intact in the decomposition.
 *
 * Also note that the decomposition is done to prefixable units. Kilogram is not prefixable (gram is) and there is no
 * general way to keep the kilo- on the kilograms when deriving units for different powers of 10.  Therefore, the
 * calculated factor corresponds to decomposition to grams.
 *
 * You must multiply the corresponding data with @factor if you intend to use the @result for them instead of @unit!
 *
 * Returns: When @result is %NULL, a newly created SI unit that has to be dereferenced when no longer used later.
 *          Otherwise @result itself is simply returned, its reference count is NOT increased.
 *
 * Since: 2.51
 **/
GwySIUnit*
gwy_si_unit_factor_to_base(GwySIUnit *siunit,
                           GwySIUnit *result,
                           gdouble *mfactor)
{
    GwySIUnit *factored = NULL;
    GwySimpleUnit *unit;
    GwySimpleUnit bunit;
    const GwySIUnitDecomposition *decomp;
    const gchar *unitstr;
    gdouble mf = 1.0;
    guint i, j, k;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit), NULL);
    g_return_val_if_fail(!result || GWY_IS_SI_UNIT(result), NULL);

    if (!result)
        result = gwy_si_unit_new(NULL);
    else if (result == siunit)
        result = factored = gwy_si_unit_new(NULL);
    else
        g_array_set_size(result->units, 0);

    for (i = 0; i < siunit->units->len; i++) {
        unit = &unit_index(siunit, i);
        unitstr = g_quark_to_string(unit->unit);
        /* Try to find the simple unit in known derived_units[] */
        for (j = 0; j < G_N_ELEMENTS(derived_units); j++) {
            decomp = derived_units + j;
            if (gwy_strequal(unitstr, decomp->unit)) {
                for (k = 0; k < GWY_SI_UNIT_N_BASE; k++) {
                    if (!decomp->powers[k])
                        continue;

                    bunit.unit = g_quark_from_static_string(base_si_units[k]);
                    bunit.power = decomp->powers[k] * unit->power;
                    bunit.traits = 0;
                    g_array_append_val(result->units, bunit);
                    mf *= decomp->factor * unit->power;
                }
                break;
            }
        }
        /* If we cannot find it, just copy it over. */
        if (j == G_N_ELEMENTS(derived_units)) {
            bunit = *unit;
            g_array_append_val(result->units, bunit);
        }
    }

    gwy_si_unit_canonicalize(result);
    /* If caller passes @result == @siunit, we have to modify @siunit and throw away the allocated unit.  Otherwise
     * @result is either some other unit or newly allocated and we just return it.  */
    if (result == factored) {
        g_array_set_size(siunit->units, 0);
        g_array_append_vals(siunit->units, result->units->data, result->units->len);
        result = siunit;
        g_object_unref(factored);
    }

    if (mfactor)
        *mfactor = mf;

    return result;
}

/* FIXME: Consider sorting SI units g, m, s, A, cd, mol. */
static GwySIUnit*
gwy_si_unit_canonicalize(GwySIUnit *siunit)
{
    GwySimpleUnit *dst, *src;
    gint i, j;

    /* consolidate multiple occurences of the same unit */
    i = 0;
    while (i < siunit->units->len) {
        src = &unit_index(siunit, i);

        for (j = 0; j < i; j++) {
            dst = &unit_index(siunit, j);
            if (src->unit == dst->unit) {
                dst->power += src->power;
                g_array_remove_index(siunit->units, i);
                break;
            }
        }

        if (j == i)
            i++;
    }

    /* remove units with zero power */
    i = 0;
    while (i < siunit->units->len) {
        if (unit_index(siunit, i).power)
            i++;
        else {
            g_array_remove_index(siunit->units, i);
        }
    }

    return siunit;
}

static gboolean
gwy_si_unit_equal_do(GwySIUnit *siunit1, GwySIUnit *siunit2, gboolean strict)
{
    GwySIUnit *factored1, *factored2;
    gdouble mf1 = 0.0, mf2 = 0.0;
    gboolean result = FALSE;

    if (siunit2 == siunit1)
        return TRUE;

    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit1), FALSE);
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit2), FALSE);
    if (gwy_si_unit_equal_direct(siunit1, siunit2))
        return TRUE;

    /* For strict comparison only accept as equal units that are identical, including any factors.  Otherwise just
     * check dimensional equality. */
    factored1 = gwy_si_unit_factor_to_base(siunit1, NULL, &mf1);
    factored2 = gwy_si_unit_factor_to_base(siunit2, NULL, &mf2);
    if ((!strict || fabs(log(mf1/mf2)) < 1e-12) && gwy_si_unit_equal_direct(factored1, factored2))
        result = TRUE;

    g_object_unref(factored1);
    g_object_unref(factored2);

    return result;
}

static gboolean
gwy_si_unit_equal_direct(GwySIUnit *siunit1, GwySIUnit *siunit2)
{
    guint i, j;

    if (siunit2 == siunit1)
        return TRUE;

    if (siunit2->units->len != siunit1->units->len)
        return FALSE;

    for (i = 0; i < siunit1->units->len; i++) {
        GwySimpleUnit *unit = &unit_index(siunit1, i);

        for (j = 0; j < siunit2->units->len; j++) {
            if (unit_index(siunit2, j).unit == unit->unit) {
                if (unit_index(siunit2, j).power != unit->power)
                    return FALSE;
                break;
            }
        }
        if (j == siunit2->units->len)
            return FALSE;
    }

    return TRUE;
}

static void
gwy_si_unit_format(GwySIUnit *siunit,
                   const GwySIStyleSpec *fs,
                   GwySIValueFormat *format,
                   gboolean magnitude_fixed)
{
    GString *string = format->units_gstring;
    const gchar *prefix = "No GCC, this can't be used uninitialized";
    GwySimpleUnit *unit;
    gint i, prefix_bearer, move_me_to_end, moveby;

    g_assert(string);
    g_string_truncate(string, 0);

    /* If there is a single unit with negative exponent, move it to the end
     * TODO: we may want more sophistication here */
    move_me_to_end = -1;
    if (siunit->units->len > 1) {
        for (i = 0; i < siunit->units->len; i++) {
            unit = &unit_index(siunit, i);
            if (unit->power < 0) {
                if (move_me_to_end >= 0) {
                    move_me_to_end = -1;
                    break;
                }
                move_me_to_end = i;
            }
        }
    }

    /* Find a victim to prepend a prefix to.  Mwhahaha. */
    prefix_bearer = -1;
    if (siunit->power10) {
        for (i = 0; i < siunit->units->len; i++) {
            if (i == move_me_to_end)
                continue;
            unit = &unit_index(siunit, i);
            if (siunit->power10 % (3*abs(unit->power)) == 0) {
                prefix_bearer = i;
                break;
            }
        }
    }
    if (siunit->power10 && prefix_bearer < 0 && move_me_to_end >= 0) {
        unit = &unit_index(siunit, move_me_to_end);
        if (siunit->power10 % (3*abs(unit->power)) == 0)
            prefix_bearer = move_me_to_end;
    }
    /* When we did not find any suitable prefix bearer, try moving the magnitude while keeping significant digits. */
    if (!magnitude_fixed && siunit->power10 && prefix_bearer < 0) {
        /* First try moving decimal dot right. */
        for (moveby = 3; prefix_bearer < 0 && moveby <= format->precision + 1; moveby += 3) {
            for (i = 0; i < siunit->units->len; i++) {
                if (i == move_me_to_end)
                    continue;
                unit = &unit_index(siunit, i);
                if ((siunit->power10 - moveby) % (3*abs(unit->power)) == 0) {
                    siunit->power10 -= moveby;
                    format->magnitude /= pow10(moveby);
                    format->precision -= moveby;
                    /* We allow moving one digit to far there, but we must avoid making the precsion negative. */
                    format->precision = MAX(format->precision, 0);
                    prefix_bearer = i;
                    break;
                }
            }
        }
        /* Then left, just once because we do not have any sanity check here.
         * XXX: is this part a good idea? */
        if (prefix_bearer < 0 && format->precision < 3) {
            moveby = 3;
            for (i = 0; i < siunit->units->len; i++) {
                if (i == move_me_to_end)
                    continue;
                unit = &unit_index(siunit, i);
                if ((siunit->power10 + moveby) % (3*abs(unit->power)) == 0) {
                    siunit->power10 += moveby;
                    format->magnitude *= pow10(moveby);
                    format->precision += moveby;
                    prefix_bearer = i;
                    break;
                }
            }
        }
    }
    /* Check whether we are not out of prefix range. */
    if (prefix_bearer >= 0) {
        unit = &unit_index(siunit, prefix_bearer);
        prefix = gwy_si_unit_prefix(siunit->power10/unit->power);
        if (!prefix)
            prefix_bearer = -1;
    }

    /* If we were unable to place the prefix, we must add a power of 10. */
    if (siunit->power10 && prefix_bearer < 0) {
        if (fs->multiplier)
            g_string_append(string, fs->multiplier);
        switch (siunit->power10) {
            case -1:
            g_string_append(string, "0.1");
            break;

            case 1:
            g_string_append(string, "10");
            break;

            case 2:
            g_string_append(string, "100");
            break;

            default:
            g_string_append(string, "10");
            fs->format_power(string, siunit->power10);
            break;
        }
        if (fs->power_unit_separator && siunit->units->len)
            g_string_append(string, fs->power_unit_separator);
    }

    /* Append units. */
    for (i = 0; i < siunit->units->len; i++) {
        if (i == move_me_to_end)
            continue;
        if (i > 1 || (i && move_me_to_end)) {
            g_string_append(string, fs->unit_times);
        }
        unit = &unit_index(siunit, i);
        if (i == prefix_bearer)
            g_string_append(string, prefix);
        g_string_append(string, g_quark_to_string(unit->unit));
        if (unit->power != 1)
            fs->format_power(string, unit->power);
    }
    if (move_me_to_end >= 0) {
        g_string_append(string, fs->unit_division);
        unit = &unit_index(siunit, move_me_to_end);
        if (move_me_to_end == prefix_bearer)
            g_string_append(string, prefix);
        g_string_append(string, g_quark_to_string(unit->unit));
        if (unit->power != -1)
            fs->format_power(string, -unit->power);
    }

    format->units = format->units_gstring->str;
}

static const gchar*
gwy_si_unit_prefix(gint power)
{
    gint i;

    for (i = 0; i < G_N_ELEMENTS(SI_prefixes); i++) {
        if (SI_prefixes[i].value == power)
            return SI_prefixes[i].name;
    }
    return NULL;
}

static void
format_power_plain(GString *string, gint power)
{
    g_string_append_printf(string, "^%d", power);
}

static void
format_power_pango(GString *string, gint power)
{
    g_string_append_printf(string, "<sup>%d</sup>", power);
}

static void
format_power_TeX(GString *string, gint power)
{
    if (power >= 0 && power <= 9)
        g_string_append_printf(string, "^%d", power);
    else
        g_string_append_printf(string, "^{%d}", power);
}

static void
format_power_unicode(GString *string, gint power)
{
    gchar buf[16];
    guint i;

    g_snprintf(buf, sizeof(buf), "%d", power);
    for (i = 0; buf[i]; i++) {
        if (buf[i] == '0' || (buf[i] >= '4' && buf[i] <= '9'))
            g_string_append_unichar(string, 0x2070 + buf[i] - '0');
        else if (buf[i] == '1')
            g_string_append_len(string, "¹", sizeof("¹")-1);
        else if (buf[i] == '2')
            g_string_append_len(string, "²", sizeof("²")-1);
        else if (buf[i] == '3')
            g_string_append_len(string, "³", sizeof("³")-1);
        else if (buf[i] == '-')
            g_string_append_len(string, "⁻", sizeof("⁻")-1);
        else {
            g_warning("Weird digits in exponent %s\n", buf);
            g_string_append_c(string, buf[i]);
        }
    }
}

/************************** Documentation ****************************/

/**
 * SECTION:gwysiunit
 * @title: GwySIUnit
 * @short_description: SI unit representation, physical quantitiy formatting
 *
 * #GwySIUnit object represents a physical SI unit (or any other unit), it can be created from a unit string with
 * gwy_si_unit_new().
 *
 * GwySIUnit is also responsible for prefixes selection and generally formatting of physical quantities (see also
 * gwymath for pure number formatting functions).  There are several functions computing value format (as
 * a #GwySIValueFormat structure) with given resolution -- gwy_si_unit_get_format_with_resolution(), or number of
 * significant digits -- gwy_si_unit_get_format_with_digits().
 **/

/**
 * GwySIUnit:
 *
 * The #GwySIUnit struct contains private data only and should be accessed
 * using the functions below.
 **/

/**
 * GwySIUnitFormatStyle:
 * @GWY_SI_UNIT_FORMAT_NONE: No units.  This value is unused by #GwySIUnit itself and must not be requested as
 *                           a format style.
 * @GWY_SI_UNIT_FORMAT_PLAIN: Plain style, as one would use on a text terminal.
 * @GWY_SI_UNIT_FORMAT_MARKUP: Pango markup, for units usable standalone.
 * @GWY_SI_UNIT_FORMAT_VFMARKUP: Pango markup, for units directly after value.
 * @GWY_SI_UNIT_FORMAT_TEX: Representation that can be typeset by TeX, for units usable standalone.
 * @GWY_SI_UNIT_FORMAT_VFTEX: Representation that can be typeset by TeX, for units directly after value.  (Since 2.50)
 * @GWY_SI_UNIT_FORMAT_UNICODE: Representation in which exponents are rendered as Unicode characters, for units usable
 *                              standalone.  (Since 2.50)
 * @GWY_SI_UNIT_FORMAT_VFUNICODE: Representation in which exponents are rendered as Unicode characters, for units
 *                                directly after value.  (Since 2.50)
 *
 * Physical quantity formatting style.
 *
 * The VF-variants differ from tne non-VF ones by a multiplication sign at the start of units (where appropriate).
 **/

/**
 * GwySIValueFormat:
 * @magnitude: Number to divide a quantity by (a power of 1000).
 * @precision: Number of decimal places to format a quantity to.
 * @units: Units to put after quantity divided by @magnitude.  This is actually an alias to @units_gstring->str.
 * @units_gstring: #GString used to represent @units internally.
 *
 * A physical quantity formatting information.
 *
 * The @magnitude and @precision fields can be directly modified if necessary. Units must be always set with
 * gwy_si_unit_value_format_set_units() to update the internal representation properly.
 */

/**
 * gwy_si_unit_duplicate:
 * @siunit: An SI unit to duplicate.
 *
 * Convenience macro doing gwy_serializable_duplicate() with all the necessary typecasting.
 **/

/**
 * gwy_si_unit_assign:
 * @dest: Target SI unit.
 * @source: Source SI unit.
 *
 * Convenience macro making one SI unit equal to another.
 *
 * This is just a gwy_serializable_clone() wrapper with all the necessary typecasting.
 *
 * Since: 2.51
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
