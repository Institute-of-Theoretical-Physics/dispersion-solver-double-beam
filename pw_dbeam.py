# -*- coding: utf-8 -*-
"""
pw_dbeam.py

Solver for double-beam instability modes in collisionless plasmas

Yasuhito Narita, Uwe Motschmann, Horia Comisel, and Daniel Schmid
y.narita@tu-braunschweig.de

May 2026
License: MIT License

Usage:
    Refer to the README for detailed instructions.
    If you use this code, please cite it as follows:

    Narita, Y., Motschmann, U., Comi\c{s}el, H., and Schmid, D.:
    Double beam instability for the Mercury upstream waves,
    Astrophys. J., 983, 125 (2026.
    https://doi.org/10.3847/1538-4357/adc1bc
"""

# MIT License (Short form)
# Copyright (c) 2026 Y. Narita, U. Motschmann, H. Comisel, D. Schmid
#
# This software is released under the MIT License.
# http://opensource.org


import math as math
import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


# =================== #
# function disp_dbeam #
# =================== #

def disp_eibd():

    # ========================================= #
    # input parameters                          #
    # e: electrons  i: ions  b: beam1  d: beam2 #
    # ========================================= #
    
    # --------------------------------------------#
    # frequency parameters (in units of omega_ci) #
    # --------------------------------------------#
    omega_ce = -1836
    omega_ci = 1
    omega_cb = 1
    omega_cd = 1

    omega_pi = 3000
    omega_pb = 100
    omega_pd = 100

    # ---------------------------------------------#
    # omega_pe determined by the charge neutrality #
    # ---------------------------------------------#
    omega_pe=np.sqrt(
                 (-omega_pi**2/omega_ci
                  -omega_pb**2/omega_cb
                  -omega_pd**2/omega_cd) * omega_ce
                )

    
    # --------------------------------------------------#
    # Alfven velocity Va (in units of speed of light c) #
    # --------------------------------------------------#
    v_a=1/np.sqrt(
              (omega_pe/omega_ce)**2
             +(omega_pi/omega_ci)**2
             +(omega_pb/omega_cb)**2
             +(omega_pd/omega_cd)**2
             )

    # ---------------------------------------------#
    # beam velocities (in units of speed of light) #
    # ---------------------------------------------#
    v_b = 7.0 * v_a
    v_d = 5.0 * v_a

    # -----------------------------------------------------------------------#
    # wavenumber array (in units of ion inertial wavenumber as k*c/omega_pi) #
    # -----------------------------------------------------------------------#
    kcw_min  = 0.00
    kcw_max  = 0.50
    kcw_tics = 1000
    kcw_arr = np.linspace(kcw_min, kcw_max, num=kcw_tics, endpoint=False)

    karr = kcw_arr
    nrk = karr.size

    # ============================= #
    # further variables, file names #
    # ============================= #

    #----------------------------------------------------------#
    # electron and ion bulk velocities from current neutrality #
    #----------------------------------------------------------#
    v_e = ((omega_ce/omega_pe)**2/(omega_ce-omega_ci)*
          ((omega_pb/omega_cb)**2*(omega_ci-omega_cb)*v_b
          + (omega_pd/omega_cd)**2*(omega_ci-omega_cd)*v_d)
          )

    v_i = ((omega_ci/omega_pi)**2/(omega_ci-omega_ce)*
           ((omega_pb/omega_cb)**2*(omega_ce-omega_cb)*v_b
            + (omega_pd/omega_cd)**2*(omega_ce-omega_cd)*v_d)
           )

    # ---------------------------#
    # output file specifications #
    # ---------------------------#
    v_b_out = v_b / v_a
    v_d_out = v_d / v_a
    
    outfile = f"fig_vb_{v_b_out:04.1f}_vd_{v_d_out:04.1f}"
    pdffile = outfile + ".pdf"

    # ------------------------#
    # charge neutrality check #
    # ------------------------#
#    laneu=(omega_pe**2/omega_ce
#           + omega_pi**2/omega_ci
#           + omega_pb**2/omega_cb
#           + omega_pd**2/omega_cd)

    # ----------------------#
    # current density check #
    # ----------------------#
#    strom=(omega_pe**2/omega_ce*v_e
#           + omega_pi**2/omega_ci*v_i
#           + omega_pb**2/omega_cb*v_b
#           + omega_pd**2/omega_cd*v_d)

    # ---------------#
    # total momentum #
    # ---------------#
#    impuls=(omega_pe**2/omega_ce**2*v_e
#            + omega_pi**2/omega_ci**2*v_i
#            + omega_pb**2/omega_cb**2*v_b
#            + omega_pd**2/omega_cd**2*v_d)

    # -------------------#
    # firehose parameter #
    # -------------------#
#    fhp=((v_e*omega_pe/omega_ce)**2
#         + (v_i*omega_pi/omega_ci)**2
#         + (v_b*omega_pb/omega_cb)**2
#         + (v_d*omega_pd/omega_cd)**2)

    # -----------------#
    # console messages #
    # -----------------#
    print(' ')
    print('omega_pe = ',omega_pe)
    print('omega_pi = ',omega_pi)
    print('omega_pb = ',omega_pb)
    print('omega_pd = ',omega_pd)
    print('omega_ce = ',omega_ce)
    print('omega_ci = ',omega_ci)
    print('omega_cb = ',omega_cb)
    print('omega_cd = ',omega_cd)
    print(' ')
    print('v_e      = ',v_e)
    print('v_i      = ',v_i)
    print('v_b      = ',v_b)
    print('v_d      = ',v_d)
    print(' ')
    print('Alfven speed v_a     = ',v_a)
#    print('parameter firehose   = ', fhp)
#    print('parameter neutrality = ', laneu)
#    print('parameter current    = ', strom)
#    print('parameter momentum  = ', impuls)


    # ================ #
    # plot preparation #
    # ================ #
    cm = 1 / 2.54
    fig = plt.figure(figsize=(8.5 * cm, 6.0 * cm), facecolor="white")
    ax1 = fig.add_subplot(211, xlabel="")
    ax2 = fig.add_subplot(212, xlabel="")
    fig.subplots_adjust(bottom=0.2, left=0.2, top=0.88, right=0.9, hspace=0.1)

    # Storage arrays for unstable mode tracks
    k_keep, w_keep, g_keep = [], [], []
    all_roots_k, all_roots_real = [], []

    omega_max=2
    omega_min=-omega_max
    omega_min=-2
    gamma_max=0
#    kcw_max = 0.5
#    kcw_min = 0.0
#    kcw_tics = 300

#    kcw_delta = (kcw_max - kcw_min) / kcw_tics
#    kcw = kcw_min
#    n = 0

    # ------------------------------------------------ #
    # loop over wavenumbers, root-finding and plotting #
    # ------------------------------------------------ #
    for n, k_norm in enumerate(karr):

        #----------------------------------------------#
        # k is changed into kc by multiplying omega_pi #
        #----------------------------------------------#
        k = k_norm * omega_pi

        # ----------------------------------------- #
        # frequency factors for dispersion equation #
        # ----------------------------------------- #
        pk = np.poly1d([-k**2])

        pde  = np.poly1d([1,-k*v_e])
        pdec = pde+omega_ce

        pdi  = np.poly1d([1,-k*v_i])
        pdic = pdi+omega_ci

        pdb  = np.poly1d([1,-k*v_b])
        pdbc = pdb+omega_cb

        pdd  = np.poly1d([1,-k*v_d])
        pddc = pdd+omega_cd

        poe = np.poly1d([omega_pe**2])
        poi = np.poly1d([omega_pi**2])
        pob = np.poly1d([omega_pb**2])
        pod = np.poly1d([omega_pd**2])

        # --------------------------------------------------- #
        # determinant of dispersion matrix in polynomial form #
        # --------------------------------------------------- #
        p = (pk *    pdec*pdic*pdbc*pddc
             - poe*pde*     pdic*pdbc*pddc
             - poi*pdi*pdec*     pdbc*pddc
             - pob*pdb*pdec*pdic*     pddc
             - pod*pdd*pdec*pdic*pdbc     )

        # ------------------ #
        # numpy root-finding #
        # ------------------ #
        nstellen = np.roots(p)

         # -------------------------#
        # map current step outputs #
        #--------------------------#
#        for z in nstellen:
#            all_roots_k.append(k / omega_pi)
#            all_roots_real.append(z.real)
#            if z.imag > 1.0e-4:
#                k_keep.append(k / omega_pi)
#                w_keep.append(z.real)
#                g_keep.append(z.imag)

        for z in nstellen:
            all_roots_k.append(k / omega_pi)
            all_roots_real.append(z.real)
            if z.imag > 1.0e-4:
                k_keep.append(k / omega_pi)
                w_keep.append(z.real)
                g_keep.append(z.imag)


        # ------------- #
        # end_of_k_loop #
        # ------------- #

    # -------------------#
    # generate plot rows #
    # -------------------#
    ax1.plot(all_roots_k, all_roots_real, marker=".", color="0.6", ls="", ms=0.8)
    ax1.plot(k_keep, w_keep, marker=".", color="k", ls="", ms=1.2)
    ax2.plot(k_keep, g_keep, marker=".", color="k", ls="", ms=1.2)

    # --------------------#
    # axes bounds control #
    # --------------------#
    xmin, xmax = kcw_min, kcw_max
    ax1.set_xlim(xmin, xmax)
#    ax1.set_ylim(-1.5, 2.5)
    ax2.set_xlim(xmin, xmax)
#    ax2.set_ylim(0.0, 0.50)

#    xticks = [-0.2, 0, 0.2]
#    ax1.set_xticks(xticks)
#    ax1.set_yticks([-1, 0, 1, 2])
#    ax2.set_xticks(xticks)
#    ax2.set_yticks([0.0, 0.2, 0.4])
    ax1.set_xticklabels([])

    # Ticks density intervals
    ax1.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.05))

    for ax in (ax1, ax2):
        ax.tick_params(
            axis="x",
            which="both",
            length=4 if ax == ax2 else 2,
            labelbottom=ax == ax2,
            bottom=True,
            top=True,
            direction="in",
        )
        ax.tick_params(
            axis="y", which="both", left=True, right=True, direction="in"
        )

    # -----------------------#
    # labels and designations#
    # -----------------------#
#    title_text = (
#        r"$v_{\|\mathrm{b}}=20\,V_\mathrm{A}, v_{\perp\mathrm{b}}=20\,V_\mathrm{A}$"
#    )

    title_text = (
        rf"$v_{{\mathrm{{b}}}}={v_b_out:04.1f}\,V_\mathrm{{A}},  "
        rf"v_{{\mathrm{{d}}}}={v_d_out:04.1f}\,V_\mathrm{{A}}$"
    )

    ax1.set_title(title_text, fontsize=10.5)
    ax1.set_ylabel(r"$\omega/\Omega_\mathrm{i}$")
    ax2.set_xlabel(r"$kc/\omega_\mathrm{pi}$")
    ax2.set_ylabel(r"$\gamma/\Omega_\mathrm{i}$")

    # -----------------------------#
    # numpy save arrays in npz format #
    # -----------------------------#
    np.savez(
        outfile,
        wavenum=np.array(all_roots_k),
        freq_real=np.array(all_roots_real),
        unstable_k=np.array(k_keep),
        unstable_real=np.array(w_keep),
        unstable_imag=np.array(g_keep),
    )

    plt.tight_layout()
    plt.savefig(pdffile)
    print(f"Data saved to {outfile}.npz and plot generated as {pdffile}")


if __name__ == "__main__":
    disp_eibd()


