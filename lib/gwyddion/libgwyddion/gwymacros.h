/*
 *  $Id: gwymacros.h 24951 2022-08-26 15:24:04Z yeti-dn $
 *  Copyright (C) 2003-2016 David Necas (Yeti), Petr Klapetek.
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

#ifndef __GWY_MACROS_H__
#define __GWY_MACROS_H__

#include <stdarg.h>
#include <string.h>
#include <glib.h>
#include <gwyconfig.h>

/* FIXME: move to gwyconfig.h? or just config.h? */
/* XXX: Most of this is available in gi18n.h */
#ifdef ENABLE_NLS
#include <libintl.h>
#else
#define gettext(x) (x)
#define ngettext(sing, plur, n) ((n) == 1 ? (sing) : (plur))
#endif
#define _(x) gettext(x)

#ifdef gettext_noop
#define N_(x) gettext_noop(x)
#else
#define N_(x) (x)
#endif

#define GWY_SWAP(t, x, y) \
    do { \
    t __gwy_unsafe_swap = x; \
    x = y; \
    y = __gwy_unsafe_swap; \
    } while (0)

#define GWY_ORDER(t, x, y) \
    do { \
        if ((y) < (x)) { \
            t __gwy_unsafe_swap = x; \
            x = y; \
            y = __gwy_unsafe_swap; \
        } \
    } while (0)

#define gwy_strequal(a, b) \
    (!strcmp((a), (b)))

#define GWY_CLAMP(x, low, hi) \
    (G_UNLIKELY((x) > (hi)) ? (hi) : (G_UNLIKELY((x) < (low)) ? (low) : (x)))

#define gwy_clear(array, n)\
    memset((array), 0, (n)*sizeof((array)[0]))

#define gwy_assign(dest, source, n) \
    memcpy((dest), (source), (n)*sizeof((dest)[0]))

#define GWY_OBJECT_UNREF(obj) \
    do { \
        if (obj) { \
            g_object_unref(obj); \
            (obj) = NULL; \
        } \
    } while (0)

#define gwy_object_unref GWY_OBJECT_UNREF

#define GWY_SIGNAL_HANDLER_DISCONNECT(obj, hid) \
    do { \
        if (hid && obj) \
            g_signal_handler_disconnect(obj, hid); \
        (hid) = 0; \
    } while (0)

#define gwy_signal_handler_disconnect GWY_SIGNAL_HANDLER_DISCONNECT

#define GWY_FIND_PSPEC(type, id, spectype) \
    G_PARAM_SPEC_##spectype(g_object_class_find_property (G_OBJECT_CLASS(g_type_class_peek(type)), id))

#define GWY_SI_VALUE_FORMAT_FREE(vf) \
    do { \
        if (vf) { \
            gwy_si_unit_value_format_free(vf); \
            (vf) = NULL; \
        } \
    } while (0)

#define GWY_FREE(ptr) \
    do { \
        if (ptr) { \
            g_free(ptr); \
            (ptr) = NULL; \
        } \
    } while (0)

G_BEGIN_DECLS

#ifdef G_HAVE_GNUC_VARARGS
#  ifdef DEBUG
#    define gwy_debug(format...) \
            gwy_debug_gnu(G_LOG_DOMAIN,\
                          __FILE__ ":" G_STRINGIFY (__LINE__), \
                          G_STRFUNC, \
                          format)
#  else
#    define gwy_debug(format...) /* */
#  endif
#  define gwy_info(format...) \
          g_log(G_LOG_DOMAIN, G_LOG_LEVEL_INFO, format)
#elif defined(G_HAVE_ISO_VARARGS)
#  ifdef DEBUG
#    define gwy_debug(...) \
            gwy_debug_gnu(G_LOG_DOMAIN, \
                          __FILE__ ":" G_STRINGIFY(__LINE__), \
                          G_STRFUNC, \
                          __VA_ARGS__)
#  else
#    define gwy_debug(...) /* */
#  endif
#  define gwy_info(...) \
          g_log(G_LOG_DOMAIN, G_LOG_LEVEL_INFO, __VA_ARGS__)
#else
/* no varargs macros */
#  ifdef DEBUG
static inline void
gwy_debug(const gchar *format, ...)
{
    va_list args;
    va_start(args, format);
    g_logv(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, format, args);
    va_end(args);
}

gwy_info(const gchar *format, ...)
{
    va_list args;
    va_start(args, format);
    g_logv(G_LOG_DOMAIN, G_LOG_LEVEL_INFO, format, args);
    va_end(args);
}
#  else
static inline void
gwy_debug(const gchar *format, ...)
{
}
#  endif
#endif /* varargs macros */

void gwy_debug_gnu(const gchar *domain,
                   const gchar *fileline,
                   const gchar *funcname,
                   const gchar *format,
                   ...) G_GNUC_PRINTF(4, 5);

G_END_DECLS

#endif /* __GWY_MACROS_H__ */

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
