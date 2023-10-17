/*
 *  $Id: linestats.h 20678 2017-12-18 18:26:55Z yeti-dn $
 *  Copyright (C) 2003-2017 David Necas (Yeti), Petr Klapetek.
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

#ifndef __GWY_PROCESS_LINESTATS_H__
#define __GWY_PROCESS_LINESTATS_H__

#include <libprocess/dataline.h>

G_BEGIN_DECLS

gdouble gwy_data_line_get_max           (GwyDataLine *data_line);
gdouble gwy_data_line_get_min           (GwyDataLine *data_line);
void    gwy_data_line_get_min_max       (GwyDataLine *data_line,
                                         gdouble *min,
                                         gdouble *max);
gdouble gwy_data_line_min_pos_i         (GwyDataLine *data_line);
gdouble gwy_data_line_max_pos_i         (GwyDataLine *data_line);
gdouble gwy_data_line_min_pos_r         (GwyDataLine *data_line);
gdouble gwy_data_line_max_pos_r         (GwyDataLine *data_line);
gdouble gwy_data_line_get_avg           (GwyDataLine *data_line);
gdouble gwy_data_line_get_rms           (GwyDataLine *data_line);
gdouble gwy_data_line_get_tan_beta0     (GwyDataLine *data_line);
gdouble gwy_data_line_get_variation     (GwyDataLine *data_line);
gdouble gwy_data_line_get_sum           (GwyDataLine *data_line);
gdouble gwy_data_line_get_ra            (GwyDataLine *data_line);
gdouble gwy_data_line_get_skew          (GwyDataLine *data_line);
gdouble gwy_data_line_get_kurtosis      (GwyDataLine *data_line);
gdouble gwy_data_line_part_get_max      (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_min      (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
void    gwy_data_line_part_get_min_max  (GwyDataLine *data_line,
                                         gint from,
                                         gint to,
                                         gdouble *min,
                                         gdouble *max);
gdouble gwy_data_line_part_get_avg      (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_rms      (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_tan_beta0(GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_variation(GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_sum      (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_ra       (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_skew     (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_part_get_kurtosis (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_get_modus         (GwyDataLine *data_line,
                                         gint histogram_steps);
gdouble gwy_data_line_part_get_modus    (GwyDataLine *data_line,
                                         gint from,
                                         gint to,
                                         gint histogram_steps);
gdouble gwy_data_line_get_median        (GwyDataLine *data_line);
gdouble gwy_data_line_part_get_median   (GwyDataLine *data_line,
                                         gint from,
                                         gint to);
gdouble gwy_data_line_get_length        (GwyDataLine *data_line);
gdouble gwy_data_line_get_xpm           (GwyDataLine *data_line,
                                         gint m,
                                         gint k);
gdouble gwy_data_line_get_xvm           (GwyDataLine *data_line,
                                         gint m,
                                         gint k);
gdouble gwy_data_line_get_xtm           (GwyDataLine *data_line,
                                         gint m,
                                         gint k);
gint    gwy_data_line_get_kth_peaks     (GwyDataLine *data_line,
                                         gint m,
                                         gint rank,
                                         gboolean peaks,
                                         gboolean average,
                                         gdouble pthreshold,
                                         gdouble vthreshold,
                                         gdouble *peakvalues);
gint    gwy_data_line_count_peaks       (GwyDataLine *data_line,
                                         gboolean peaks,
                                         gdouble pthreshold,
                                         gdouble vthreshold);
void    gwy_data_line_distribution      (GwyDataLine *data_line,
                                         GwyDataLine *distribution,
                                         gdouble ymin,
                                         gdouble ymax,
                                         gboolean normalize_to_unity,
                                         gint nstats);
void    gwy_data_line_dh                (GwyDataLine *data_line,
                                         GwyDataLine *target_line,
                                         gdouble ymin,
                                         gdouble ymax,
                                         gint nsteps);
void    gwy_data_line_cdh               (GwyDataLine *data_line,
                                         GwyDataLine *target_line,
                                         gdouble ymin,
                                         gdouble ymax,
                                         gint nsteps);
void    gwy_data_line_da                (GwyDataLine *data_line,
                                         GwyDataLine *target_line,
                                         gdouble ymin,
                                         gdouble ymax,
                                         gint nsteps);
void    gwy_data_line_cda               (GwyDataLine *data_line,
                                         GwyDataLine *target_line,
                                         gdouble ymin,
                                         gdouble ymax,
                                         gint nsteps);
void    gwy_data_line_acf               (GwyDataLine *data_line,
                                         GwyDataLine *target_line);
void    gwy_data_line_hhcf              (GwyDataLine *data_line,
                                         GwyDataLine *target_line);
/* XXX: Fix the crappy arguments to enums, cannot do it now as the type differs
 * so it would be an API breakage. */
void    gwy_data_line_psdf              (GwyDataLine *data_line,
                                         GwyDataLine *target_line,
                                         gint windowing,
                                         gint interpolation);

G_END_DECLS

#endif /* __GWY_PROCESS_LINESTATS_H__ */


/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

