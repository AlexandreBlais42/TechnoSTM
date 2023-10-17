/*
 *  $Id: gwynlfitpreset.c 24845 2022-06-27 10:16:28Z yeti-dn $
 *  Copyright (C) 2000-2003 Martin Siler.
 *  Copyright (C) 2004-2020 David Necas (Yeti), Petr Klapetek.
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

#include "config.h"
#include <string.h>
#include <libgwyddion/gwymacros.h>
#include <libgwyddion/gwymath.h>
#include <libgwyddion/gwydebugobjects.h>
#include <libgwyddion/gwynlfitpreset.h>
#include "gwyddioninternal.h"

static GwyNLFitPreset*
gwy_nlfit_preset_new_static(const GwyNLFitPresetBuiltin *data);

G_DEFINE_TYPE(GwyNLFitPreset, gwy_nlfit_preset, GWY_TYPE_RESOURCE)

G_GNUC_UNUSED
static gdouble
_gwy_gamma(gdouble x)
{
    static const gdouble p[] = {
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };
    gdouble xi = 0.99999999999980993, t;
    guint i;

    if (x < 0.5)
        return G_PI/(sin(G_PI*x)*_gwy_gamma(1.0 - x));

    x -= 1.0;
    for (i = 0; i < G_N_ELEMENTS(p); i++)
        xi += p[i]/(x + i + 1);

    t = x + G_N_ELEMENTS(p) - 0.5;
    return sqrt(2.0*G_PI) * pow(t, x + 0.5) * exp(-t) * xi;
}

/*********************** gaussian *****************************/
static gdouble
gauss_func(gdouble x,
           G_GNUC_UNUSED gint n_param,
           const gdouble *b,
           G_GNUC_UNUSED gpointer user_data,
           gboolean *fres)
{
    gdouble c;

    if (b[3] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = (x - b[0])/b[3];
    return b[2] * exp(-c*c) + b[1];
}

static void
gauss_guess(gint n_dat,
            const gdouble *x,
            const gdouble *y,
            gdouble *param,
            gboolean *fres)
{
    gint i, imin, imax, nv;
    gdouble c0, cm, vec[5];

    param[0] = 0.0;
    param[1] = 0.0;
    param[2] = 0.0;
    param[3] = 1.0;

    imin = imax = 0;
    for (i = 1; i < n_dat; i++) {
        if (y[imax] < y[i])
            imax = i;
        if (y[imin] > y[i])
            imin = i;
    }

    c0 = 0.0;
    for (i = 0; i < n_dat; i++)
        c0 += y[i];
    c0 /= n_dat;

    nv = G_N_ELEMENTS(vec);
    for (i = 0; i < nv; i++)
        vec[i] = y[i*(n_dat - 1)/(nv - 1)];
    cm = gwy_math_median(nv, vec);

    if (cm - y[imin] < y[imax] - cm) {
        param[1] = y[imin];
        param[0] = x[imax];
    }
    else {
        param[1] = y[imax];
        param[0] = x[imin];
    }

    c0 = c0 - param[1];
    c0 *= fabs(x[n_dat-1] - x[0]);
    c0 /= G_PI;
    param[2] = c0*fabs(c0)/(y[imax] - y[imin]);

    if (c0)
        param[3] = param[2]/c0;
    else
        param[3] = fabs(x[n_dat-1] - x[0])/4.0;

    param[2] /= param[3]*param[3];
    param[3] /= sqrt(G_LN2);

    *fres = TRUE;
}

/******************** gaussian PSDF ***************************/
static gdouble
gauss_psdf_func(gdouble x,
                G_GNUC_UNUSED gint n_param,
                const gdouble *b,
                G_GNUC_UNUSED gpointer user_data,
                gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x*b[1];

    return b[0]*b[0]*b[1]/(2.0*GWY_SQRT_PI) * exp(-c*c/4);
}

static void
gauss_psdf_guess(gint n_dat,
                 const gdouble *x,
                 const gdouble *y,
                 gdouble *param,
                 gboolean *fres)
{
    gdouble sw, w, xx, l, k2, k4, b0, b2, d, alpha, beta;
    gint i;

    param[0] = 0.0;
    param[1] = 1.0;

    sw = k2 = k4 = b0 = b2 = 0.0;
    for (i = 0; i < n_dat; i++) {
        if (x[i] <= 0.0 || y[i] <= 0.0)
            continue;

        w = y[i]*y[i];
        sw += w;
        xx = x[i]*x[i];
        k2 += w*xx;
        k4 += w*xx*xx;
        l = log(y[i]);
        b0 += w*l;
        b2 += w*xx*l;
    }
    if (!sw) {
        *fres = FALSE;
        return;
    }

    k2 /= sw;
    k4 /= sw;
    b0 /= sw;
    b2 /= sw;
    d = k4 - k2*k2;
    if (!d) {
        *fres = FALSE;
        return;
    }
    alpha = (b0*k4 - b2*k2)/d;
    beta = (b2 - k2*b0)/d;
    if (beta >= 0.0) {
        *fres = FALSE;
        return;
    }
    param[1] = 2*sqrt(-beta);
    param[0] = exp(alpha/2)*sqrt(2*GWY_SQRT_PI/param[1]);
    *fres = TRUE;
}

static GwySIUnit*
psdf_get_units(G_GNUC_UNUSED GwyNLFitPreset *preset,
               guint param,
               GwySIUnit *siunit_x,
               GwySIUnit *siunit_y)
{
    GwySIUnit *result;

    switch (param) {
        case 0:
        result = gwy_si_unit_multiply(siunit_x, siunit_y, NULL);
        if (!gwy_si_unit_nth_root(result, 2, result))
            gwy_si_unit_set_from_string(result, NULL);
        return result;
        break;

        case 1:
        return gwy_si_unit_power(siunit_x, -1, NULL);
        break;

        default:
        g_return_val_if_reached(NULL);
        break;
    }
}

/******************** two-gaussian PSDF ***************************/
static gdouble
gauss2_psdf_func(gdouble x,
                 G_GNUC_UNUSED gint n_param,
                 const gdouble *b,
                 G_GNUC_UNUSED gpointer user_data,
                 gboolean *fres)
{
    gdouble f1, f2;
    gboolean fres1, fres2;
    f1 = gauss_psdf_func(x, 2, b, NULL, &fres1);
    f2 = gauss_psdf_func(x, 2, b + 2, NULL, &fres2);
    *fres = fres1 && fres2;
    return f1 + f2;
}

static void
gauss2_psdf_guess(gint n_dat,
                  const gdouble *x,
                  const gdouble *y,
                  gdouble *param,
                  gboolean *fres)
{
    gauss_psdf_guess(n_dat, x, y, param, fres);

    if (!*fres) {
        param[2] = 0.0;
        param[3] = 1.0;
        return;
    }

    /* Don't actually try to guess.  Use the one-Gauss estimate, split sigma
     * and make sure Ts are somewhat different. */
    param[2] = param[0] /= G_SQRT2;

    param[3] = GWY_SQRT3*param[1];
    param[1] /= GWY_SQRT3;
}

static GwySIUnit*
psdf2_get_units(GwyNLFitPreset *preset,
                guint param,
                GwySIUnit *siunit_x,
                GwySIUnit *siunit_y)
{
    return psdf_get_units(preset, param % 2, siunit_x, siunit_y);
}

/******************* gaussian HHCF ********************************/
static gdouble
gauss_hhcf_func(gdouble x,
                G_GNUC_UNUSED gint n_param,
                const gdouble *b,
                G_GNUC_UNUSED gpointer user_data,
                gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x/b[1];

    return 2*b[0]*b[0] * (1 - exp(-c*c));
}

static void
gauss_hhcf_guess(gint n_dat,
                 const gdouble *x,
                 const gdouble *y,
                 gdouble *param,
                 gboolean *fres)
{
    gint i;

    param[0] = 0;
    for (i = (n_dat/2); i < n_dat; i++) {
        param[0] += y[i]/(n_dat/2);
    }
    param[0] = sqrt(param[0]);
    param[1] = x[n_dat-1]/50;

    *fres = TRUE;
}

static GwySIUnit*
cf_get_units(G_GNUC_UNUSED GwyNLFitPreset *preset,
             guint param,
             GwySIUnit *siunit_x,
             GwySIUnit *siunit_y)
{
    GwySIUnit *result;

    switch (param) {
        case 0:
        result = gwy_si_unit_duplicate(siunit_y);
        if (!gwy_si_unit_nth_root(result, 2, result))
            gwy_si_unit_set_from_string(result, NULL);
        return result;
        break;

        case 1:
        return gwy_si_unit_duplicate(siunit_x);
        break;

        default:
        g_return_val_if_reached(NULL);
        break;
    }
}

/****************    gaussian ACF  *****************************************/
static gdouble
gauss_acf_func(gdouble x,
               G_GNUC_UNUSED gint n_param,
               const gdouble *b,
               G_GNUC_UNUSED gpointer user_data,
               gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x/b[1];

    return b[0]*b[0] * exp(-c*c);
}

static void
gauss_acf_guess(gint n_dat,
                const gdouble *x,
                const gdouble *y,
                gdouble *param,
                gboolean *fres)
{
    gdouble sw, w, xx, l, k2, k4, b0, b2, d, alpha, beta;
    gint i;

    param[0] = 0.0;
    param[1] = 1.0;

    sw = k2 = k4 = b0 = b2 = 0.0;
    for (i = 0; i < n_dat; i++) {
        if (x[i] <= 0.0 || y[i] <= 0.0)
            continue;

        w = y[i]*y[i];
        sw += w;
        xx = x[i]*x[i];
        k2 += w*xx;
        k4 += w*xx*xx;
        l = log(y[i]);
        b0 += w*l;
        b2 += w*xx*l;
    }
    if (!sw) {
        *fres = FALSE;
        return;
    }

    k2 /= sw;
    k4 /= sw;
    b0 /= sw;
    b2 /= sw;
    d = k4 - k2*k2;
    if (!d) {
        *fres = FALSE;
        return;
    }
    alpha = (b0*k4 - b2*k2)/d;
    beta = (b2 - k2*b0)/d;
    if (beta >= 0.0) {
        *fres = FALSE;
        return;
    }
    param[0] = exp(alpha/2);
    param[1] = sqrt(-1.0/beta);
    *fres = TRUE;
}

/******************** gaussian IPSDF ***************************/
static gdouble
gauss_rpsdf_func(gdouble x,
                 G_GNUC_UNUSED gint n_param,
                 const gdouble *b,
                 G_GNUC_UNUSED gpointer user_data,
                 gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x*b[1];

    return b[0]*b[0]*b[1]*b[1]/2.0 * x*exp(-c*c/4);
}

static void
gauss_rpsdf_guess(gint n_dat,
                  const gdouble *x,
                  const gdouble *y,
                  gdouble *param,
                  gboolean *fres)
{
    gdouble sw, w, xx, l, k2, k4, b0, b2, d, alpha, beta;
    gint i;

    param[0] = 0.0;
    param[1] = 1.0;

    sw = k2 = k4 = b0 = b2 = 0.0;
    for (i = 0; i < n_dat; i++) {
        if (x[i] <= 0.0 || y[i] <= 0.0)
            continue;

        w = y[i]*y[i];
        sw += w;
        xx = x[i]*x[i];
        k2 += w*xx;
        k4 += w*xx*xx;
        l = log(y[i]/x[i]);
        b0 += w*l;
        b2 += w*xx*l;
    }
    if (!sw) {
        *fres = FALSE;
        return;
    }

    k2 /= sw;
    k4 /= sw;
    b0 /= sw;
    b2 /= sw;
    d = k4 - k2*k2;
    if (!d) {
        *fres = FALSE;
        return;
    }
    alpha = (b0*k4 - b2*k2)/d;
    beta = (b2 - k2*b0)/d;
    if (beta >= 0.0) {
        *fres = FALSE;
        return;
    }
    param[1] = 2*sqrt(-beta);
    param[0] = exp(alpha/2)*G_SQRT2/param[1];
    *fres = TRUE;
}

/**************** exponential ************************************/
static gdouble
exp_func(gdouble x,
         G_GNUC_UNUSED gint n_param,
         const gdouble *b,
         G_GNUC_UNUSED gpointer user_data,
         gboolean *fres)
{
    if (b[2] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;

    return b[1]*exp(x/b[2]) + b[0];
}

static void
exp_guess(gint n_dat,
          const gdouble *x,
          const gdouble *y,
          gdouble *param,
          gboolean *fres)
{
    gint i, imin, imax;
    gdouble s;

    imin = imax = 0;
    s = y[0];
    for (i = 1; i < n_dat; i++) {
        if (y[i] > y[imax])
            imax = i;
        if (y[i] < y[imin])
            imin = i;
        s += y[i];
    }

    if (y[imax] == y[imin]) {
        param[0] = 0.0;
        param[1] = y[imin];
        param[2] = 10*(x[n_dat-1] - x[0]);
        return;
    }

    s /= n_dat;
    if (2.0*s < y[imax] + y[imin])
        s -= y[imin];
    else
        s -= y[imax];
    s *= x[imax] - x[imin];

    param[2] = s/(y[imax] - y[imin]);
    param[1] = (y[imax] - y[imin])/(exp(x[imax]/param[2])
                                    - exp(x[imin]/param[2]));
    param[0] = y[imin] - param[1]*exp(x[imin]/param[2]);

    *fres = TRUE;
}

/**************** exponential PSDF **************************/
static gdouble
exp_psdf_func(gdouble x,
              G_GNUC_UNUSED gint n_param,
              const gdouble *b,
              G_GNUC_UNUSED gpointer user_data,
              gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x*b[1];

    return b[0]*b[0]*b[1]/G_PI/(1.0 + c*c);
}

static void
exp_psdf_guess(gint n_dat,
               const gdouble *x,
               const gdouble *y,
               gdouble *param,
               gboolean *fres)
{
    gdouble s0, s1, max;
    gint i, n;

    param[0] = 0.0;
    param[1] = 0.0;
    param[2] = 0.0;

    s0 = s1 = 0.0;
    max = -G_MAXDOUBLE;
    n = 0;
    for (i = 0; i < n_dat; i++) {
        if (x[i] < 0.0)
            continue;

        n++;
        s0 += y[i];
        s1 += sqrt(x[i])*y[i];
        if (x[i] > max)
            max = x[i];
    }
    if (!s0 || !s1) {
        *fres = FALSE;
        return;
    }

    max /= n;
    s0 *= max;
    s1 *= max;

    param[0] = sqrt(s0);
    param[1] = 2.0*s0*s0/(s1*s1);
    *fres = TRUE;
}

/***************** exponential HHCF ********************************/
static gdouble
exp_hhcf_func(gdouble x,
              G_GNUC_UNUSED gint n_param,
              const gdouble *b,
              G_GNUC_UNUSED gpointer user_data,
              gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x/b[1];

    return 2*b[0]*b[0] * (1 - exp(-c));
}

static void
exp_hhcf_guess(gint n_dat,
               const gdouble *x,
               const gdouble *y,
               gdouble *param,
               gboolean *fres)
{
    gint i;

    param[0] = 0;
    for (i = (n_dat/2); i < n_dat; i++)
        param[0] += y[i]/(n_dat/2);

    *fres = param[0] >= 0;
    param[0] = sqrt(param[0]);
    param[1] = x[n_dat-1]/50;
}

/*************** exponential ACF ************************************/
static gdouble
exp_acf_func(gdouble x,
             G_GNUC_UNUSED gint n_param,
             const gdouble *b,
             G_GNUC_UNUSED gpointer user_data,
             gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x/b[1];

    return b[0]*b[0] * exp(-c);
}

static void
exp_acf_guess(gint n_dat,
              const gdouble *x,
              const gdouble *y,
              gdouble *param,
              gboolean *fres)
{
    gdouble sw, w, xx, l, k2, k4, b0, b2, d, alpha, beta;
    gint i;

    param[0] = 0.0;
    param[1] = 0.0;
    param[2] = 0.0;

    sw = k2 = k4 = b0 = b2 = 0.0;
    for (i = 0; i < n_dat; i++) {
        if (x[i] <= 0.0 || y[i] <= 0.0)
            continue;

        w = y[i]*y[i];
        sw += w;
        xx = x[i];
        k2 += w*xx;
        k4 += w*xx*xx;
        l = log(y[i]);
        b0 += w*l;
        b2 += w*xx*l;
    }
    if (!sw) {
        *fres = FALSE;
        return;
    }

    k2 /= sw;
    k4 /= sw;
    b0 /= sw;
    b2 /= sw;
    d = k4 - k2*k2;
    if (!d) {
        *fres = FALSE;
        return;
    }
    alpha = (b0*k4 - b2*k2)/d;
    beta = (b2 - k2*b0)/d;
    if (beta >= 0.0) {
        *fres = FALSE;
        return;
    }
    param[0] = exp(alpha/2);
    param[1] = -1.0/beta;
    *fres = TRUE;
}

/******************** exponential IPSDF ***************************/
static gdouble
exp_rpsdf_func(gdouble x,
               G_GNUC_UNUSED gint n_param,
               const gdouble *b,
               G_GNUC_UNUSED gpointer user_data,
               gboolean *fres)
{
    gdouble c, d;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = x*b[1];
    d = 1 + c*c;
    d *= d*d;

    return b[0]*b[0]*b[1] * c/sqrt(d);
}

static void
exp_rpsdf_guess(gint n_dat,
                const gdouble *x,
                const gdouble *y,
                gdouble *param,
                gboolean *fres)
{
    gdouble s0, sm1, max;
    gint i, n;

    param[0] = 0.0;
    param[1] = 0.0;
    param[2] = 0.0;

    s0 = sm1 = 0.0;
    max = -G_MAXDOUBLE;
    n = 0;
    for (i = 0; i < n_dat; i++) {
        if (x[i] <= 0.0)
            continue;

        n++;
        s0 += y[i];
        sm1 += y[i]/x[i];
        if (x[i] > max)
            max = x[i];
    }
    if (!s0 || !sm1) {
        *fres = FALSE;
        return;
    }

    max /= n;
    s0 *= max;
    sm1 *= max;

    param[0] = sqrt(s0);
    param[1] = sm1/s0;
    *fres = TRUE;
}

/******************** K-correlated PSDF ***************************/

/* XXX: σ² = s² β(1/2, (C-1)/2)
 * But we cannot normalise to σ² because C and σ are terribly correlated. */
static gdouble
kcorr_psdf_func(gdouble x,
                G_GNUC_UNUSED gint n_param,
                const gdouble *b,
                G_GNUC_UNUSED gpointer user_data,
                gboolean *fres)
{
    gdouble s = b[0], T = b[1], C = b[2], t;

    if (T == 0.0 || C <= 0.0) {
        *fres = FALSE;
        return 0.0;
    }
    *fres = TRUE;
    t = T*x;
    return T*s*s/pow(1.0 + t*t, 0.5*(C + 1.0));
}

static void
kcorr_psdf_guess(gint n_dat,
                 const gdouble *x,
                 const gdouble *y,
                 gdouble *param,
                 gboolean *fres)
{
    exp_psdf_guess(n_dat, x, y, param, fres);
    /* Usually C is around 1.  We probably cannot guess it better... */
    param[2] = 0.2;
}

static GwySIUnit*
kcorr_psdf_get_units(GwyNLFitPreset *preset,
                     guint param,
                     GwySIUnit *siunit_x,
                     GwySIUnit *siunit_y)
{
    if (param == 2)
        return gwy_si_unit_new(NULL);

    return psdf_get_units(preset, param, siunit_x, siunit_y);
}

/******************** K-correlated PSDF ***************************/

static gdouble
powerexp_acf_func(gdouble x,
                  G_GNUC_UNUSED gint n_param,
                  const gdouble *b,
                  G_GNUC_UNUSED gpointer user_data,
                  gboolean *fres)
{
    gdouble c;

    if (b[1] == 0) {
        *fres = FALSE;
        return 0;
    }
    *fres = TRUE;
    c = fabs(x/b[1]);
    if (c > 0.0)
        c = pow(c, b[2]);

    return b[0]*b[0] * exp(-c);
}

static void
powerexp_acf_guess(gint n_dat,
                   const gdouble *x,
                   const gdouble *y,
                   gdouble *param,
                   gboolean *fres)
{
    exp_acf_guess(n_dat, x, y, param, fres);
    param[2] = 1.5;
}

static GwySIUnit*
powerexp_acf_get_units(GwyNLFitPreset *preset,
                       guint param,
                       GwySIUnit *siunit_x,
                       GwySIUnit *siunit_y)
{
    if (param == 2)
        return gwy_si_unit_new(NULL);

    return cf_get_units(preset, param, siunit_x, siunit_y);
}

/**************   polynomial 0th order ********************************/
static gdouble
poly_0_func(G_GNUC_UNUSED gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    *fres = TRUE;

    return b[0];
}

static void
poly_0_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gwy_math_fit_polynom(n_dat, x, y, 0, param);
    *fres = TRUE;
}

/*************** polynomial 1st order ********************************/
static gdouble
poly_1_func(gdouble x,
           G_GNUC_UNUSED gint n_param,
           const gdouble *b,
           G_GNUC_UNUSED gpointer user_data,
           gboolean *fres)
{
    *fres = TRUE;

    return b[0] + x*b[1];
}

static void
poly_1_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gwy_math_fit_polynom(n_dat, x, y, 1, param);
    *fres = TRUE;
}

/************* polynomial 2nd order **********************************/
static gdouble
poly_2_func(gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    *fres = TRUE;

    return b[0] + x*(b[1] + x*b[2]);
}

static void
poly_2_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gwy_math_fit_polynom(n_dat, x, y, 2, param);
    *fres = TRUE;
}

/************** polynomial 3rd order *****************************/
static gdouble
poly_3_func(gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    *fres = TRUE;

    return b[0] + x*(b[1] + x*(b[2] + x*b[3]));
}

static void
poly_3_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gwy_math_fit_polynom(n_dat, x, y, 3, param);
    *fres = TRUE;
}

/************** polynomial 4th order *****************************/
static gdouble
poly_4_func(gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    *fres = TRUE;

    return b[0] + x*(b[1] + x*(b[2] + x*(b[3] + x*b[4])));
}

static void
poly_4_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gwy_math_fit_polynom(n_dat, x, y, 4, param);
    *fres = TRUE;
}

/************** polynomial 4th order *****************************/
static gdouble
poly_5_func(gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    *fres = TRUE;

    return b[0] + x*(b[1] + x*(b[2] + x*(b[3] + x*(b[4] + x*b[5]))));
}

static void
poly_5_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gwy_math_fit_polynom(n_dat, x, y, 5, param);
    *fres = TRUE;
}

/******************* square signal ********************************/
static gdouble
square_func(gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    gdouble s, t;

    x = (x - b[1])/b[0];

    if (b[0] == 0.0) {
        *fres = FALSE;
        return 0.0;
    }

    /* This would be an exactly square wave but we cannot fit a step-like
     * function because the sum of squares has not meaningful derivatives
     * by parameters.  So fake an almost-square function using sine. */
    /* return (x - floor(x)) > 0.5 ? b[2] : b[3]; */

    *fres = TRUE;
    s = sin(2.0*G_PI*x);
    t = pow(fabs(s), 0.001) * (s > 0.0 ? -1.0 : 1.0);
    return 0.5*(b[2] + b[3]) + 0.5*t*(b[2] - b[3]);
}

static void
square_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gint i;
    gdouble min, max;

    param[0] = fabs(x[n_dat - 1] - x[0])/10.0;
    param[1] = 0;

    min = G_MAXDOUBLE;
    max = -G_MAXDOUBLE;
    for (i = 0; i < n_dat; i++) {
        if (min > y[i])
            min = y[i];
        if (max < y[i])
            max = y[i];
    }
    param[2] = min;
    param[3] = max;

    *fres = TRUE;
}

/******************* parabolic step ********************************/
static gdouble
parstep_func(gdouble x,
             G_GNUC_UNUSED gint n_param,
             const gdouble *b,
             G_GNUC_UNUSED gpointer user_data,
             gboolean *fres)
{
    gdouble xc = b[0], yc = b[1], w = b[2], h = b[3], s = b[4];
    gdouble xleft, xright;

    if (h == 0.0 || w == 0.0) {
        *fres = FALSE;
        return 0.0;
    }

    *fres = TRUE;
    x -= xc;
    xleft = MAX(1.0 + 0.5*s/h, 0.0);
    xright = MAX(1.0 - 0.5*s/h, 0.0);
    xleft = -0.5*fabs(w)*sqrt(xleft);
    xright = 0.5*fabs(w)*sqrt(xright);
    if (x <= xleft)
        return yc - 0.5*s;
    if (x >= xright)
        return yc + 0.5*s;

    x *= 2.0/w;
    return yc + h*(1.0 - x*x);
}

static void
parstep_guess(gint n_dat,
              const gdouble *x,
              const gdouble *y,
              gdouble *param,
              gboolean *fres)
{
    gint i;
    gdouble xmin, xmax, yxmin, yxmax, ymin, ymax, yminx, ymaxx;

    xmin = G_MAXDOUBLE;
    xmax = -G_MAXDOUBLE;
    ymin = G_MAXDOUBLE;
    ymax = -G_MAXDOUBLE;
    yxmax = yxmin = yminx = ymaxx = 0.0;
    for (i = 0; i < n_dat; i++) {
        if (ymin > y[i]) {
            ymin = y[i];
            yminx = x[i];
        }
        if (ymax < y[i]) {
            ymax = y[i];
            ymaxx = x[i];
        }
        if (xmin > x[i]) {
            xmin = x[i];
            yxmin = y[i];
        }
        if (xmax < x[i]) {
            xmax = x[i];
            yxmax = y[i];
        }
    }
    param[0] = param[1] = param[4] = 0.0;
    param[2] = param[3] = 1.0;
    *fres = FALSE;
    if (xmax <= xmin)
        return;

    param[1] = 0.5*(yxmin + yxmax);
    param[4] = yxmax - yxmin;
    param[2] = (xmax - xmin)/6.0;
    if (fabs(ymax - param[1]) >= fabs(ymin - param[1])) {
        /* Positive parabola. */
        param[0] = ymaxx;
        param[3] = ymax - param[1];
    }
    else {
        /* Negative parabola. */
        param[0] = yminx;
        param[3] = ymin - param[1];
    }

    *fres = TRUE;
}

/******************* smooth step ********************************/
static gdouble
smstep_func(gdouble x,
            G_GNUC_UNUSED gint n_param,
            const gdouble *b,
            G_GNUC_UNUSED gpointer user_data,
            gboolean *fres)
{
    gdouble xc = b[0], yc = b[1], w = b[2], h = b[3], alpha = b[4], beta = b[5];
    gdouble ax, f;

    *fres = TRUE;
    x -= xc;
    w = fabs(w);
    ax = fabs(x);
    f = yc + alpha*x + 0.5*beta*ax;
    if (w == 0.0)
        f += (x > 0.0) ? 0.5*h : -0.5*h;
    else
        f += h*erf(x/w)/2.0 + 0.5*beta*w*exp(-ax/w);

    return f;
}

static void
smstep_guess(gint n_dat,
             const gdouble *x,
             const gdouble *y,
             gdouble *param,
             gboolean *fres)
{
    gint i, nleft = 0, nright = 0;
    gdouble xmin, xmax, ymin, ymax, xymin, xymax, yleft, yright;

    xmin = G_MAXDOUBLE;
    xmax = -G_MAXDOUBLE;
    ymin = G_MAXDOUBLE;
    ymax = -G_MAXDOUBLE;
    xymin = xymax = 0.0;
    for (i = 0; i < n_dat; i++) {
        if (ymin > y[i]) {
            xymin = x[i];
            ymin = y[i];
        }
        if (ymax < y[i]) {
            ymax = y[i];
            xymax = x[i];
        }
        if (xmin > x[i])
            xmin = x[i];
        if (xmax < x[i])
            xmax = x[i];
    }
    param[0] = param[1] = param[4] = param[5] = 0.0;
    param[2] = param[3] = 1.0;
    *fres = FALSE;

    if (xmax <= xmin)
        return;

    yleft = yright = 0.0;
    for (i = 0; i < n_dat; i++) {
        gdouble q = (x[i] - xmin)/(xmax - xmin);
        if (q >= 0.2 && q <= 0.42) {
            yleft += y[i];
            nleft++;
        }
        if (q >= 0.58 && q <= 0.8) {
            yright += y[i];
            nright++;
        }
    }

    param[0] = 0.5*(xmin + xmax);
    param[1] = 0.5*(ymin + ymax);
    param[2] = (xmax - xmin)/8.0;
    if (nleft && nright)
        param[3] = (yleft/nleft < yright/nright) ? ymax - ymin : ymin - ymax;
    else
        param[3] = (xymax > xymin) ? ymax - ymin : ymin - ymax;
    /* Using zero here prevents L-M from making progress. */
    param[4] = 0.05*param[3]/(xmax - xmin);
    param[5] = param[4]/(xmax - xmin);

    *fres = TRUE;
}

/******************* erf bent step ********************************/
static gdouble
bentstep_func(gdouble x,
              G_GNUC_UNUSED gint n_param,
              const gdouble *b,
              G_GNUC_UNUSED gpointer user_data,
              gboolean *fres)
{
    gdouble xc = b[0], yc = b[1], w = b[2], h = b[3], alpha = b[4], beta = b[5];
    gdouble f;

    *fres = TRUE;
    x -= xc;
    w = fabs(w);
    f = yc + alpha*x + beta*x*x;
    if (w == 0.0)
        f += (x > 0.0) ? 0.5*h : -0.5*h;
    else
        f += h*erf(x/w)/2.0;

    return f;
}

/******************* tanh bent step ********************************/
static gdouble
boltzmannstep_func(gdouble x,
                   G_GNUC_UNUSED gint n_param,
                   const gdouble *b,
                   G_GNUC_UNUSED gpointer user_data,
                   gboolean *fres)
{
    gdouble xc = b[0], yc = b[1], w = b[2], h = b[3], alpha = b[4], beta = b[5];
    gdouble f;

    *fres = TRUE;
    x -= xc;
    w = fabs(w);
    f = yc + alpha*x + beta*x*x;
    if (w == 0.0)
        f += (x > 0.0) ? 0.5*h : -0.5*h;
    else
        f += h*tanh(x/w)/2.0;

    return f;
}

/******************* power function ********************************/
static gdouble
power_func(gdouble x,
           G_GNUC_UNUSED gint n_param,
           const gdouble *b,
           G_GNUC_UNUSED gpointer user_data,
           gboolean *fres)
{
    if (x < 0.0) {
        *fres = FALSE;
        return 0;
    }
    if (x == 0.0 && b[3] <= 0.0) {
        *fres = FALSE;
        return 0;
    }

    *fres = TRUE;
    return b[0] + b[1]*pow(x, b[2]);
}

static void
power_guess(gint n_dat,
            const gdouble *x,
            const gdouble *y,
            gdouble *param,
            gboolean *fres)
{
    gint i;
    gdouble q, la, lb, c1, c2;

    param[0] = 0.0;
    param[1] = 0.0;
    param[2] = 1.0;

    if (n_dat < 4 || x[0] <= 0.0) {
        *fres = FALSE;
        return;
    }

    i = GWY_ROUND(sqrt(n_dat-1));
    i = CLAMP(i, 2, n_dat-2);
    la = log(x[1]/x[i]);
    lb = log(x[n_dat-1]/x[i]);
    q = (y[n_dat-1] - y[i])/(y[i] - y[1]);
    if (fabs(q) < 0.4)
        param[2] = -2.0*(q*la + lb)/(q*la*la + lb*lb);
    else {
        c1 = log(q)/lb;
        c2 = -log(q)/la;
        if (c1 + c2 <= 0.0)
            param[2] = c2;
        else
            param[2] = c1;
    }

    q = (y[n_dat-1] - y[i])*(y[i] - y[1]);
    q /= (pow(x[n_dat - 1], param[2]) - pow(x[i], param[2]))
         *(pow(x[i], param[2]) - pow(x[1], param[2]));
    param[1] = sqrt(fabs(q));

    param[0] = (y[0] + y[i] + y[n_dat-1]
                - param[1]*(pow(x[0], param[2]) + pow(x[i], param[2])
                            + pow(x[n_dat-1], param[2])))/3.0;

    *fres = TRUE;
}

static GwySIUnit*
power_get_units(G_GNUC_UNUSED GwyNLFitPreset *preset,
                guint param,
                G_GNUC_UNUSED GwySIUnit *siunit_x,
                GwySIUnit *siunit_y)
{
    switch (param) {
        case 0:
        return gwy_si_unit_duplicate(siunit_y);
        break;

        case 1:
        /* XXX: Ugly, arbitrary real power */
        return gwy_si_unit_new(NULL);
        break;

        case 2:
        return gwy_si_unit_new(NULL);
        break;

        default:
        g_return_val_if_reached(NULL);
        break;
    }
}

/******************* lorentzian ********************************/
static gdouble
lorentz_func(gdouble x,
             G_GNUC_UNUSED gint n_param,
             const gdouble *b,
             G_GNUC_UNUSED gpointer user_data,
             gboolean *fres)
{
    *fres = TRUE;
    return b[0] + b[2]/(b[3]*b[3] + (x - b[1])*(x - b[1]));
}

static void
lorentz_guess(gint n_dat,
              const gdouble *x,
              const gdouble *y,
              gdouble *param,
              gboolean *fres)
{
    gint i, imin, imax, nv;
    gdouble c0, cm, vec[5];

    param[0] = 0.0;
    param[1] = 0.0;
    param[2] = 0.0;
    param[3] = 1.0;

    imin = imax = 0;
    for (i = 1; i < n_dat; i++) {
        if (y[imax] < y[i])
            imax = i;
        if (y[imin] > y[i])
            imin = i;
    }

    c0 = 0.0;
    for (i = 0; i < n_dat; i++)
        c0 += y[i];
    c0 /= n_dat;

    nv = G_N_ELEMENTS(vec);
    for (i = 0; i < nv; i++)
        vec[i] = y[i*(n_dat - 1)/(nv - 1)];
    cm = gwy_math_median(nv, vec);

    if (cm - y[imin] < y[imax] - cm) {
        param[0] = y[imin];
        param[1] = x[imax];
    }
    else {
        param[0] = y[imax];
        param[1] = x[imin];
    }

    c0 = c0 - param[0];
    c0 *= fabs(x[n_dat-1] - x[0]);
    c0 /= G_PI;
    param[2] = c0*fabs(c0)/(y[imax] - y[imin]);

    if (c0)
        param[3] = param[2]/c0;
    else
        param[3] = fabs(x[n_dat-1] - x[0])/4.0;

    *fres = TRUE;
}

/******************* sinc ********************************/
static gdouble
sinc_func(gdouble x,
          G_GNUC_UNUSED gint n_param,
          const gdouble *b,
          G_GNUC_UNUSED gpointer user_data,
          gboolean *fres)
{
    *fres = TRUE;
    x /= b[1];
    if (x == 0.0)
        return b[0];
    return b[0]*sin(x)/x;
}

static void
sinc_guess(gint n_dat,
           const gdouble *x,
           const gdouble *y,
           gdouble *param,
           gboolean *fres)
{
    gint i;
    gdouble max, min, xmin;

    max = min = y[0];
    xmin = x[0];
    for (i = 1; i < n_dat; i++) {
        if (G_UNLIKELY(max < y[i]))
            max = y[i];
        if (G_UNLIKELY(min > y[i])) {
            min = y[i];
            xmin = x[i];
        }
    }

    param[1] = xmin/4.493409457909064175307880927276;
    param[0] = max;
    *fres = TRUE;
}


/******************** preset default weights *************************/
static void
weights_linear_decrease(gint n_dat,
                        G_GNUC_UNUSED const gdouble *x,
                        G_GNUC_UNUSED const gdouble *y,
                        gdouble *weight)
{
    gint i;

    for (i = 0; i < n_dat; i++)
        weight[i] = 1 - (gdouble)i/(gdouble)n_dat;
}

/************************** generic scalers ****************************/

static GwySIUnit*
builtin_get_unit(GwyNLFitPreset *preset,
                 guint param,
                 GwySIUnit *siunit_x,
                 GwySIUnit *siunit_y)
{
    const GwyNLFitParam *p;

    g_return_val_if_fail(preset->builtin, NULL);
    g_return_val_if_fail(param < preset->builtin->nparams, NULL);

    if (preset->builtin->get_unit)
        return preset->builtin->get_unit(preset, param, siunit_x, siunit_y);

    /* Generic unit construction when there are no peculiarities */
    p = preset->builtin->param + param;
    return gwy_si_unit_power_multiply(siunit_x, p->power_x,
                                      siunit_y, p->power_y,
                                      NULL);
}

/************************** presets ****************************/

static const GwyNLFitParam gauss_params[] = {
   { "x<sub>0</sub>", 1, 0, },
   { "y<sub>0</sub>", 0, 1, },
   { "a",             0, 1, },
   { "b",             1, 0, },
};

static const GwyNLFitParam exp_params[] = {
   { "y<sub>0</sub>", 0, 1, },
   { "a",             0, 1, },
   { "b",             1, 0, },
};

/* XXX: σ power is fractional, attempted to fix case-by-case */
static const GwyNLFitParam gauss_two_params[] = {
   { "σ", 0, 0, },
   { "T", 1, 0, },
};

/* XXX: σ power is fractional, attempted to fix case-by-case */
static const GwyNLFitParam gauss2_four_params[] = {
   { "σ<sub>1</sub>", 0, 0, },
   { "T<sub>1</sub>", 1, 0, },
   { "σ<sub>2</sub>", 0, 0, },
   { "T<sub>2</sub>", 1, 0, },
};

/* XXX: s power is fractional, attempted to fix case-by-case */
static const GwyNLFitParam kcorr_psdf_params[] = {
   { "s", 0, 0, },
   { "T", 1, 0, },
   { "C", 0, 0, },
};

/* XXX: σ power is fractional, attempted to fix case-by-case */
static const GwyNLFitParam powerexp_acf_params[] = {
   { "σ", 0, 0, },
   { "T", 1, 0, },
   { "p", 0, 0, },
};

static const GwyNLFitParam poly5_params[] = {
   { "a",  0, 1, },
   { "b", -1, 1, },
   { "c", -2, 1, },
   { "d", -3, 1, },
   { "e", -4, 1, },
   { "f", -5, 1, },
};

static const GwyNLFitParam square_params[] = {
    { "T", 1, 0, },
    { "s", 1, 0, },
    { "t", 0, 1, },
    { "b", 0, 1, },
};

static const GwyNLFitParam parstep_params[] = {
    { "x<sub>0</sub>", 1, 0, },
    { "y<sub>0</sub>", 0, 1, },
    { "w", 1, 0, },
    { "h", 0, 1, },
    { "s", 0, 1, },
};

static const GwyNLFitParam smstep_params[] = {
    { "x<sub>0</sub>", 1, 0, },
    { "y<sub>0</sub>", 0, 1, },
    { "w", 1, 0, },
    { "h", 0, 1, },
    { "α", -1, 1, },
    { "β", -1, 1, },
};

static const GwyNLFitParam bentstep_params[] = {
    { "x<sub>0</sub>", 1, 0, },
    { "y<sub>0</sub>", 0, 1, },
    { "w", 1, 0, },
    { "h", 0, 1, },
    { "α", -1, 1, },
    { "β", -2, 1, },
};

static const GwyNLFitParam power_params[] = {
    { "a", 0, 1, },
    { "b", 0, 1, },
    { "c", 0, 0, },
};

static const GwyNLFitParam lorentz_params[] = {
    { "y<sub>0</sub>", 0, 1, },
    { "x<sub>0</sub>", 1, 0, },
    { "a",             2, 1, },
    { "b",             1, 0, },
};

static const GwyNLFitParam sinc_params[] = {
    { "a", 0, 1, },
    { "b", 1, 0, },
};

static const GwyNLFitPresetBuiltin fitting_presets[] = {
    {
        N_("Gaussian"),
        "<i>f</i>(<i>x</i>) "
            "= <i>y</i><sub>0</sub> "
            "+ <i>a</i> exp[−(<i>x</i> "
            "− <i>x</i><sub>0</sub>)<sup>2</sup>"
            "/b<sup>2</sup>]",
        &gauss_func,
        NULL,
        &gauss_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(gauss_params),
        gauss_params,
    },
    {
        N_("Gaussian (PSDF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sup>2</sup><i>T</i>/(2√π) "
            "exp[−(<i>x</i><i>T</i>/2)<sup>2</sup>]",
        &gauss_psdf_func,
        NULL,
        &gauss_psdf_guess,
        &psdf_get_units,
        NULL,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Gaussian (ACF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sup>2</sup> "
            "exp[−(<i>x</i>/<i>T</i>)<sup>2</sup>]",
        &gauss_acf_func,
        NULL,
        &gauss_acf_guess,
        &cf_get_units,
        &weights_linear_decrease,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Gaussian (HHCF)"),
        "<i>f</i>(<i>x</i>) "
            "= 2σ<sup>2</sup>"
            "[1 − exp(−(<i>x</i>/<i>T</i>)<sup>2</sup>)]",
        &gauss_hhcf_func,
        NULL,
        &gauss_hhcf_guess,
        &cf_get_units,
        &weights_linear_decrease,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Gaussian (RPSDF)"),
        "<i>f</i>(<i>x</i>) "
            "= (σ<i>T</i>)<sup>2</sup>/2 "
            "<i>x</i> exp[−(<i>x</i><i>T</i>/2)<sup>2</sup>]",
        &gauss_rpsdf_func,
        NULL,
        &gauss_rpsdf_guess,
        &psdf_get_units,
        NULL,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Two Gaussians (PSDF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sub>1</sub><sup>2</sup><i>T</i><sub>1</sub>/(2√π) "
            "exp[−(<i>x</i><i>T</i><sub>1</sub>/2)<sup>2</sup>] "
            "+ σ<sub>2</sub><sup>2</sup><i>T</i><sub>2</sub>/(2√π) "
            "exp[−(<i>x</i><i>T</i><sub>2</sub>/2)<sup>2</sup>]",
        &gauss2_psdf_func,
        NULL,
        &gauss2_psdf_guess,
        &psdf2_get_units,
        NULL,
        G_N_ELEMENTS(gauss2_four_params),
        gauss2_four_params,
    },
    {
        N_("Exponential"),
        "<i>f</i>(<i>x</i>) "
            "= <i>y</i><sub>0</sub> + <i>a</i> exp(<i>x</i>/<i>b</i>)",
        &exp_func,
        NULL,
        &exp_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(exp_params),
        exp_params,
    },
    {
        N_("Exponential (PSDF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sup>2</sup><i>T</i>/π "
            "1/[1 + <i>x</i><sup>2</sup><i>T</i><sup>2</sup>]",
        &exp_psdf_func,
        NULL,
        &exp_psdf_guess,
        &psdf_get_units,
        NULL,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Exponential (ACF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sup>2</sup> exp[−<i>x</i>/<i>T</i>]",
        &exp_acf_func,
        NULL,
        &exp_acf_guess,
        &cf_get_units,
        &weights_linear_decrease,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Exponential (HHCF)"),
        "<i>f</i>(<i>x</i>) "
            "= 2σ<sup>2</sup>"
            "[1 − exp(−<i>x</i>/<i>T</i>)]",
        &exp_hhcf_func,
        NULL,
        &exp_hhcf_guess,
        &cf_get_units,
        &weights_linear_decrease,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("Exponential (RPSDF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sup>2</sup><i>T</i> "
            "<i>xT</i>/[1 + <i>x</i><sup>2</sup><i>T</i><sup>2</sup>]"
            "<sup>3/2</sup>",
        &exp_rpsdf_func,
        NULL,
        &exp_rpsdf_guess,
        &psdf_get_units,
        NULL,
        G_N_ELEMENTS(gauss_two_params),
        gauss_two_params,
    },
    {
        N_("K-correlated (PSDF)"),
        "<i>f</i>(<i>x</i>) "
            "= s<sup>2</sup><i>T</i>"
            "/[1 + <i>x</i><sup>2</sup><i>T</i><sup>2</sup>]"
            "<sup><i>C</i>/2</sup>",
        &kcorr_psdf_func,
        NULL,
        &kcorr_psdf_guess,
        &kcorr_psdf_get_units,
        NULL,
        G_N_ELEMENTS(kcorr_psdf_params),
        kcorr_psdf_params,
    },
    {
        N_("Power-exponential (ACF)"),
        "<i>f</i>(<i>x</i>) "
            "= σ<sup>2</sup> "
            "exp[−(<i>x</i>/<i>T</i>)<sup><i>p</i></sup>]",
        &powerexp_acf_func,
        NULL,
        &powerexp_acf_guess,
        &powerexp_acf_get_units,
        &weights_linear_decrease,
        G_N_ELEMENTS(powerexp_acf_params),
        powerexp_acf_params,
    },
    {
        N_("Polynomial (order 0)"),
        "<i>f</i>(<i>x</i>) = <i>a</i>",
        &poly_0_func,
        NULL,
        &poly_0_guess,
        NULL,
        NULL,
        1,
        poly5_params,
    },
    {
        N_("Polynomial (order 1)"),
        "<i>f</i>(<i>x</i>) = <i>a</i> + <i>b</i><i>x</i>",
        &poly_1_func,
        NULL,
        &poly_1_guess,
        NULL,
        NULL,
        2,
        poly5_params,
    },
    {
        N_("Polynomial (order 2)"),
        "<i>f</i>(<i>x</i>) "
            "= <i>a</i> + <i>b</i><i>x</i> + <i>c</i><i>x</i><sup>2</sup>",
        &poly_2_func,
        NULL,
        &poly_2_guess,
        NULL,
        NULL,
        3,
        poly5_params,
    },
    {
        N_("Polynomial (order 3)"),
        "<i>f</i>(<i>x</i>) "
            "= <i>a</i> + <i>b</i><i>x</i> + <i>c</i><i>x</i><sup>2</sup> "
            "+ <i>d</i><i>x</i><sup>3</sup>",
        &poly_3_func,
        NULL,
        &poly_3_guess,
        NULL,
        NULL,
        4,
        poly5_params,
    },
    {
        N_("Polynomial (order 4)"),
        "<i>f</i>(<i>x</i>) "
            "= <i>a</i> + <i>b</i><i>x</i> + <i>c</i><i>x</i><sup>2</sup> "
            "+ <i>d</i><i>x</i><sup>3</sup> + <i>e</i><i>x</i><sup>4</sup>",
        &poly_4_func,
        NULL,
        &poly_4_guess,
        NULL,
        NULL,
        5,
        poly5_params,
    },
    {
        N_("Polynomial (order 5)"),
        "<i>f</i>(<i>x</i>) "
            "= <i>a</i> + <i>b</i><i>x</i> + <i>c</i><i>x</i><sup>2</sup> "
            "+ <i>d</i><i>x</i><sup>3</sup> + <i>e</i><i>x</i><sup>4</sup> "
            "+ <i>f</i><i>x</i><sup>5</sup>",
        &poly_5_func,
        NULL,
        &poly_5_guess,
        NULL,
        NULL,
        6,
        poly5_params,
    },
    {
        N_("Square wave"),
        "<i>f</i>(<i>x</i>) "
            "= (<i>t</i> + <i>b</i>)/2 "
            "+ 2(<i>t</i> − <i>b</i>)"
            "/π "
            "∑<sub>k</sub> "
            "sin(2π<i>k</i>(<i>x</i> − <i>s</i>)/<i>T</i>)/<i>k</i>",
        &square_func,
        NULL,
        &square_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(square_params),
        square_params,
    },
    {
        N_("Parabolic step"),
        "<i>f</i>(<i>x</i>) "
            "= y<sub>0</sub> "
            "+ {−<i>s</i>/2 "
            "or <i>h</i>[1 − 4(<i>x</i> − <i>x</i><sub>0</sub>)<sup>2</sup>"
            "/<i>w</i><sup>2</sup>] "
            "or <i>s</i>/2}",
        &parstep_func,
        NULL,
        &parstep_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(parstep_params),
        parstep_params,
    },
    {
        N_("Smooth slanted step"),
        "<i>f</i>(<i>x</i>) "
            "= y<sub>0</sub> "
            "+ <i>h</i>/2 erf(<i>ξ</i>/<i>w</i>) "
            "+ <i>α</i><i>ξ</i> "
            "+ <i>β</i>/2 [|<i>ξ</i>| "
            "+ <i>w</i> exp(-|<i>ξ</i>|/<i>w</i>)], "
            "<i>ξ</i> = <i>x</i> − <i>x</i><sub>0</sub>",
        &smstep_func,
        NULL,
        &smstep_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(smstep_params),
        smstep_params,
    },
    {
        N_("Smooth bent step"),
        "<i>f</i>(<i>x</i>) "
            "= y<sub>0</sub> "
            "+ <i>h</i>/2 erf(<i>ξ</i>/<i>w</i>) "
            "+ <i>α</i><i>ξ</i> "
            "+ <i>β</i><i>ξ</i><sup>2</sup>, "
            "<i>ξ</i> = <i>x</i> − <i>x</i><sub>0</sub>",
        &bentstep_func,
        NULL,
        &smstep_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(bentstep_params),
        bentstep_params,
    },
    {
        N_("Boltzmann bent step"),
        "<i>f</i>(<i>x</i>) "
            "= y<sub>0</sub> "
            "+ <i>h</i>/2 tanh(<i>ξ</i>/<i>w</i>) "
            "+ <i>α</i><i>ξ</i> "
            "+ <i>β</i><i>ξ</i><sup>2</sup>, "
            "<i>ξ</i> = <i>x</i> − <i>x</i><sub>0</sub>",
        &boltzmannstep_func,
        NULL,
        &smstep_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(bentstep_params),
        bentstep_params,
    },
    {
        N_("Power"),
        "<i>f</i>(<i>x</i>) "
            "= <i>a</i> + <i>b</i><i>x</i><sup><i>c</i></sup>",
        &power_func,
        NULL,
        &power_guess,
        &power_get_units,
        NULL,
        G_N_ELEMENTS(power_params),
        power_params,
    },
    {
        N_("Lorentzian"),
        "<i>f</i>(<i>x</i>) "
            "= <i>y</i><sub>0</sub> "
            "+ <i>a</i>/[<i>b</i><sup>2</sup> "
            "+ (<i>x</i> − <i>x</i><sub>0</sub>)<sup>2</sup>]",
        &lorentz_func,
        NULL,
        &lorentz_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(lorentz_params),
        lorentz_params,
    },
    {
        N_("Sinc"),
        "<i>f</i>(<i>x</i>) "
            "= <i>a</i> sinc(<i>x</i>/<i>b</i>)",
        &sinc_func,
        NULL,
        &sinc_guess,
        NULL,
        NULL,
        G_N_ELEMENTS(sinc_params),
        sinc_params,
    },
};

static void
gwy_nlfit_preset_class_init(GwyNLFitPresetClass *klass)
{
    GwyResourceClass *parent_class, *res_class = GWY_RESOURCE_CLASS(klass);

    parent_class = GWY_RESOURCE_CLASS(gwy_nlfit_preset_parent_class);
    res_class->item_type = *gwy_resource_class_get_item_type(parent_class);

    res_class->item_type.type = G_TYPE_FROM_CLASS(klass);
    res_class->name = "nlfitpresets";
    res_class->inventory = gwy_inventory_new(&res_class->item_type);
    gwy_inventory_forget_order(res_class->inventory);
}

static void
gwy_nlfit_preset_init(G_GNUC_UNUSED GwyNLFitPreset *preset)
{
}

static GwyNLFitPreset*
gwy_nlfit_preset_new_static(const GwyNLFitPresetBuiltin *data)
{
    GwyNLFitPreset *preset;

    preset = g_object_new(GWY_TYPE_NLFIT_PRESET, "is-const", TRUE, NULL);
    preset->builtin = data;
    g_string_assign(GWY_RESOURCE(preset)->name, data->name);

    return preset;
}

void
_gwy_nlfit_preset_class_setup_presets(void)
{
    GwyResourceClass *klass;
    GwyNLFitPreset *preset;
    guint i;

    /* Force class instantiation, this function is called before it's first
     * referenced. */
    klass = g_type_class_ref(GWY_TYPE_NLFIT_PRESET);

    for (i = 0; i < G_N_ELEMENTS(fitting_presets); i++) {
        preset = gwy_nlfit_preset_new_static(fitting_presets + i);
        gwy_inventory_insert_item(klass->inventory, preset);
        g_object_unref(preset);
    }
    gwy_inventory_restore_order(klass->inventory);

    /* The presets added a reference so we can safely unref it again */
    g_type_class_unref(klass);
}

/**
 * gwy_nlfit_preset_get_value:
 * @x: The point to compute value at.
 * @preset: A NL fitter function preset.
 * @params: Preset parameter values.
 * @fres: Set to %TRUE if succeeds, %FALSE on failure.
 *
 * Calculates preset function value in a single point with given parameters.
 *
 * Returns: The function value.
 **/
gdouble
gwy_nlfit_preset_get_value(GwyNLFitPreset *preset,
                           gdouble x,
                           const gdouble *params,
                           gboolean *fres)
{
    /* FIXME: builtin */
    return preset->builtin->function(x, preset->builtin->nparams, params,
                                     NULL, fres);
}

/**
 * gwy_nlfit_preset_get_formula:
 * @preset: A NL fitter function preset.
 *
 * Gets function formula of @preset (with Pango markup).
 *
 * Returns: The preset function formula.
 **/
const gchar*
gwy_nlfit_preset_get_formula(GwyNLFitPreset* preset)
{
    g_return_val_if_fail(GWY_IS_NLFIT_PRESET(preset), "");
    /* FIXME: builtin */
    return preset->builtin->formula;
}

/**
 * gwy_nlfit_preset_get_param_name:
 * @preset: A NL fitter function preset.
 * @param: A parameter number.
 *
 * Gets the name of a fitting parameter of a fitter preset.
 *
 * The name may contain Pango markup.
 *
 * Returns: The name of parameter @param.
 **/
const gchar*
gwy_nlfit_preset_get_param_name(GwyNLFitPreset* preset,
                                gint param)
{
    const GwyNLFitParam *par;

    g_return_val_if_fail(GWY_IS_NLFIT_PRESET(preset), "");
    /* FIXME: builtin */
    g_return_val_if_fail(param >= 0 && param < preset->builtin->nparams, NULL);
    par = preset->builtin->param + param;

    return par->name;
}

/**
 * gwy_nlfit_preset_get_param_units:
 * @preset: A NL fitter function preset.
 * @param: A parameter number.
 * @siunit_x: SI unit of abscissa.
 * @siunit_y: SI unit of ordinate.
 *
 * Derives the SI unit of a fitting parameter from the units of abscissa and
 * ordinate.
 *
 * Returns: A newly created #GwySIUnit with the units of the parameter @param.
 *          If the units of @param are not representable as #GwySIUnit,
 *          the result is unitless (i.e. it will be presented as a mere
 *          number).
 *
 * Since: 2.5
 **/
GwySIUnit*
gwy_nlfit_preset_get_param_units(GwyNLFitPreset *preset,
                                 gint param,
                                 GwySIUnit *siunit_x,
                                 GwySIUnit *siunit_y)
{
    g_return_val_if_fail(GWY_IS_NLFIT_PRESET(preset), NULL);
    /* FIXME: builtin */
    return builtin_get_unit(preset, param, siunit_x, siunit_y);
}

/**
 * gwy_nlfit_preset_get_nparams:
 * @preset: A NL fitter function preset.
 *
 * Return the number of parameters of @preset.
 *
 * Returns: The number of function parameters.
 **/
gint
gwy_nlfit_preset_get_nparams(GwyNLFitPreset* preset)
{
    g_return_val_if_fail(GWY_IS_NLFIT_PRESET(preset), 0);
    /* FIXME: builtin */
    return preset->builtin->nparams;
}

/**
 * gwy_nlfit_preset_guess:
 * @preset: A NL fitter function preset.
 * @n_dat: The number of data points (number of items in @x and @y).
 * @x: Abscissa points.
 * @y: Ordinate points.
 * @params: The array to fill with estimated parameter values.  It has to be
 *          at least gwy_nlfit_preset_get_nparams() long.
 * @fres: Set to %TRUE if succeeds, %FALSE on failure.
 *
 * Performs initial parameter estimate for a NL fitter.
 *
 * The initial estimate method depends on the function used.  There is no
 * absolute guarantee of quality, however if the data points approximately
 * match the fitted function the fit will typically converge from the returned
 * estimate.
 *
 * The parameters are filled also on failure, though just with some neutral
 * values that should not give raise to NaNs and infinities.
 **/
void
gwy_nlfit_preset_guess(GwyNLFitPreset *preset,
                       gint n_dat,
                       const gdouble *x,
                       const gdouble *y,
                       gdouble *params,
                       gboolean *fres)
{
    g_return_if_fail(GWY_IS_NLFIT_PRESET(preset));
    /* FIXME: builtin */
    preset->builtin->guess(n_dat, x, y, params, fres);
}

/**
 * gwy_nlfit_preset_create_fitter:
 * @preset: A NL fitter function preset.
 *
 * Creates a non-linear least-squares fitter for a preset function.
 *
 * The created fitter will be of the function data type, as created with
 * gwy_math_nlfit_new().
 *
 * If you do not need to modify the fitter settings you can use
 * gwy_nlfit_preset_fit() directly with %NULL fitter.
 *
 * Returns: A newly created fitter for @preset.
 *
 * Since: 2.47
 **/
GwyNLFitter*
gwy_nlfit_preset_create_fitter(GwyNLFitPreset *preset)
{
    g_return_val_if_fail(GWY_IS_NLFIT_PRESET(preset), NULL);
    return gwy_math_nlfit_new(preset->builtin->function,
                              preset->builtin->derive
                              ? preset->builtin->derive
                              : gwy_math_nlfit_derive);
}

/**
 * gwy_nlfit_preset_fit:
 * @preset: A NL fitter function preset.
 * @fitter: A Marquardt-Levenberg nonlinear fitter already initialized for
 *          @preset's function, or %NULL.
 * @n_dat: The number of data points (number of items in @x and @y).
 * @x: Abscissa points.
 * @y: Ordinate points.
 * @params: Initial parameter estimate (the number of parameters depends on
 *          the fitted preset and it can be obtained with
 *          gwy_nlfit_preset_get_nparams()).
 * @err: Array to store parameter errros to, may be %NULL.
 * @fixed_param: Which parameters should be treated as fixed (set
 *               corresponding element to %TRUE for them).  May be %NULL if
 *               all parameters are variable.
 *
 * Performs a nonlinear fit with a preset.
 *
 * See gwy_math_nlfit_fit_full() for details.
 *
 * Returns: Either @fitter itself, or a newly created fitter if it was %NULL.
 **/
GwyNLFitter*
gwy_nlfit_preset_fit(GwyNLFitPreset *preset,
                     GwyNLFitter *fitter,
                     gint n_dat,
                     const gdouble *x,
                     const gdouble *y,
                     gdouble *param,
                     gdouble *err,
                     const gboolean *fixed_param)
{
    gdouble *weight = NULL;
    gboolean ok;
    gint i;

    g_return_val_if_fail(GWY_IS_NLFIT_PRESET(preset), NULL);
    /* FIXME: builtin */
    /*use numerical derivation if necessary*/
    if (fitter) {
        /* XXX */
        g_return_val_if_fail(fitter->fmarq == preset->builtin->function, NULL);
    }
    else
        fitter = gwy_math_nlfit_new(preset->builtin->function,
                                    preset->builtin->derive
                                    ? preset->builtin->derive
                                    : gwy_math_nlfit_derive);

    /*load default weights for given function type*/
    if (preset->builtin->set_default_weights) {
        weight = g_new(gdouble, n_dat);
        preset->builtin->set_default_weights(n_dat, x, y, weight);
    }

    /* FIXME: builtin */
    ok = gwy_math_nlfit_fit_full(fitter, n_dat, x, y, weight,
                                 preset->builtin->nparams, param,
                                 fixed_param, NULL, preset) >= 0.0;

    if (ok && err && gwy_math_nlfit_get_covar(fitter)) {
    /* FIXME: builtin */
        for (i = 0; i < preset->builtin->nparams; i++)
            err[i] = gwy_math_nlfit_get_sigma(fitter, i);
    }

    g_free(weight);

    return fitter;
}

/**
 * gwy_nlfit_presets:
 *
 * Gets inventory with all the NLFit presets.
 *
 * Returns: NLFit preset inventory.
 **/
GwyInventory*
gwy_nlfit_presets(void)
{
    return
        GWY_RESOURCE_CLASS(g_type_class_peek(GWY_TYPE_NLFIT_PRESET))->inventory;
}

/************************** Documentation ****************************/

/**
 * SECTION:gwynlfitpreset
 * @title: GwyNLFitPreset
 * @short_description: NL fitter preset functions
 * @see_also: #GwyNLFitter
 *
 * <link linkend="GwyNLFitter">Non-linear fitter</link> presets are predefined
 * fitting functions, with guessing, weighting, etc. to ease fitting of
 * common functions.  Most of those currently implemented are related to
 * one-dimensional statistical characteristics of randomly rough surfaces,
 * however, several general-purpose fitting functions are also available.
 *
 * The presets are identified by name and can be obtained from the
 * corresponding #GwyInventory, see gwy_nlfit_presets().  All their properties
 * can be queried with methods such as gwy_nlfit_preset_get_nparams() or
 * gwy_nlfit_preset_get_formula().  A generic method to derive units of
 * fitting parameters from units of fitted data is also available:
 * gwy_nlfit_preset_get_param_units().
 *
 * As of version 2.7 the defined functions include:
 * <simplelist type='vert'>
 * <member><literal>"Gaussian"</literal></member>
 * <member><literal>"Gaussian (PSDF)"</literal></member>
 * <member><literal>"Gaussian (ACF)"</literal></member>
 * <member><literal>"Gaussian (HHCF)"</literal></member>
 * <member><literal>"Gaussian (RPSDF)"</literal></member>
 * <member><literal>"Exponential"</literal></member>
 * <member><literal>"Exponential (PSDF)"</literal></member>
 * <member><literal>"Exponential (ACF)"</literal></member>
 * <member><literal>"Exponential (HHCF)"</literal></member>
 * <member><literal>"Exponential (RPSDF)"</literal></member>
 * <member><literal>"Polynomial (order 0)"</literal></member>
 * <member><literal>"Polynomial (order 1)"</literal></member>
 * <member><literal>"Polynomial (order 2)"</literal></member>
 * <member><literal>"Polynomial (order 3)"</literal></member>
 * <member><literal>"Square wave"</literal></member>
 * <member><literal>"Power"</literal></member>
 * <member><literal>"Lorentzian"</literal></member>
 * <member><literal>"Sinc"</literal></member>
 * </simplelist>
 *
 * The following presets are available since version 2.40:
 * <simplelist type='vert'>
 * <member><literal>"Parabolic step"</literal></member>
 * <member><literal>"Two Gaussians (PSDF)"</literal></member>
 * </simplelist>
 *
 * The following presets are available since version 2.41:
 * <simplelist type='vert'>
 * <member><literal>"Smooth bent step"</literal></member>
 * <member><literal>"Smooth slanted step"</literal></member>
 * </simplelist>
 *
 * The following presets are available since version 2.48:
 * <simplelist type='vert'>
 * <member><literal>"K-correlated (PSDF)"</literal></member>
 * </simplelist>
 *
 * The following presets are available since version 2.49:
 * <simplelist type='vert'>
 * <member><literal>"Boltzmann bent step"</literal></member>
 * </simplelist>
 *
 * The result of the fitting is stored in a normal #GwyNLFitter, therefore the
 * typical use is:
 *
 * |[
 * GwyNLFitPreset *preset;
 * GwyNLFitter *fitter;
 * gdouble *params, *errors;
 *
 * preset = gwy_inventory_get_item(gwy_nlfit_presets(), "Gaussian");
 * params = g_new(gdouble, gwy_nlfit_preset_get_nparams(preset));
 * errors = g_new(gdouble, gwy_nlfit_preset_get_nparams(preset));
 * /&ast; Use a priori parameter value information, or: &ast;/
 * gwy_nlfit_preset_guess(preset, ndata, xdata, ydata, params, &ok);
 * fitter = gwy_nlfit_preset_fit(preset, NULL,
 *                               ndata, xdata, ydata, params,
 *                               params, errors, NULL);
 *
 * /&ast; Examine fitting status, parameters and errors... &ast;/
 *
 * g_free(params);
 * g_free(errors);
 * gwy_math_nlfit_free(fitter);
 * ]|
 **/

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
