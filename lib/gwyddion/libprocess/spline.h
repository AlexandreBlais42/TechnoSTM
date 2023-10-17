/*
 *  $Id: spline.h 20678 2017-12-18 18:26:55Z yeti-dn $
 *  Copyright (C) 2016-2017 David Necas (Yeti).
 *  E-mail: yeti@gwyddion.net.
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

#ifndef __GWY_SPLINE_H__
#define __GWY_SPLINE_H__

#include <glib.h>
#include <libgwyddion/gwymath.h>

G_BEGIN_DECLS

typedef struct _GwySpline GwySpline;

#define GWY_TYPE_SPLINE (gwy_spline_get_type())

GType        gwy_spline_get_type        (void)               G_GNUC_CONST;
GwySpline*   gwy_spline_new             (void)               G_GNUC_MALLOC;
void         gwy_spline_free            (GwySpline *spline);
GwySpline*   gwy_spline_copy            (GwySpline *spline)  G_GNUC_MALLOC;
GwySpline*   gwy_spline_new_from_points (const GwyXY *xy,
                                         guint n);
guint        gwy_spline_get_npoints     (GwySpline *spline);
const GwyXY* gwy_spline_get_points      (GwySpline *spline);
const GwyXY* gwy_spline_get_tangents    (GwySpline *spline);
gdouble      gwy_spline_get_slackness   (GwySpline *spline);
gboolean     gwy_spline_get_closed      (GwySpline *spline);
void         gwy_spline_set_points      (GwySpline *spline,
                                         const GwyXY *xy,
                                         guint n);
void         gwy_spline_set_slackness   (GwySpline *spline,
                                         gdouble slackness);
void         gwy_spline_set_closed      (GwySpline *spline,
                                         gboolean closed);
gdouble      gwy_spline_length          (GwySpline *spline);
const GwyXY* gwy_spline_sample_naturally(GwySpline *spline,
                                         guint *n);
gdouble      gwy_spline_sample_uniformly(GwySpline *spline,
                                         GwyXY *xy,
                                         GwyXY *t,
                                         guint n);

G_END_DECLS

#endif /* __GWY_SPLINE_H__ */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
