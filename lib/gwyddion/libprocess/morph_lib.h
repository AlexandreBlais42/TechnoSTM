/*
 *  $Id: morph_lib.h 24187 2021-09-23 15:40:24Z yeti-dn $
 *  Copyright (C) 2003 David Necas (Yeti), Petr Klapetek.
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

/*< private_header >*/

#ifndef __GWY_MORPH_LIB_H__
#define __GWY_MORPH_LIB_H__

#include <libprocess/tip.h>

G_BEGIN_DECLS

G_GNUC_INTERNAL
gint **_gwy_morph_lib_iallocmatrix(gint yres, gint xres);

G_GNUC_INTERNAL
void _gwy_morph_lib_ifreematrix(gint **mptr);

/*simple routines - integer arithmetics*/
G_GNUC_INTERNAL
gint **_gwy_morph_lib_ireflect(const gint *const *image,
                               gint xres, gint yres);

G_GNUC_INTERNAL
gint **_gwy_morph_lib_ierosion(const gint *const *image, gint xres, gint yres,
                               const gint *const *tip, gint txres, gint tyres,
                               gint xc, gint yc,
                               GwySetMessageFunc set_message,
                               GwySetFractionFunc set_fraction);

G_GNUC_INTERNAL
gint **_gwy_morph_lib_icmap(const gint *const *image, gint xres, gint yres,
                            const gint *const *tip, gint txres, gint tyres,
                            const gint *const *rsurf, gint xc, gint yc,
                            GwySetMessageFunc set_message,
                            GwySetFractionFunc set_fraction);

/*tip estimation routines - all in integer artithmetics*/
G_GNUC_INTERNAL
gint _gwy_morph_lib_itip_estimate(const gint *const *image,
                                  gint xres, gint yres,
                                  gint txres, gint tyres,
                                  gint xc, gint yc,
                                  gint **tip0,
                                  gint thresh,
                                  gboolean use_edges,
                                  GwySetMessageFunc set_message,
                                  GwySetFractionFunc set_fraction);

G_GNUC_INTERNAL
gint _gwy_morph_lib_itip_estimate0(const gint *const *image,
                                   gint xres, gint yres,
                                   gint txres, gint tyres,
                                   gint xc, gint yc,
                                   gint **tip0,
                                   gint thresh,
                                   gboolean use_edges,
                                   GwySetMessageFunc set_message,
                                   GwySetFractionFunc set_fraction);


G_END_DECLS

#endif

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
