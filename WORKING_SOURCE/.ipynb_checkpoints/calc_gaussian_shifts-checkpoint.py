#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 14:55:44 2024

@author: coletamburri
"""

first_step_err_1 = np.sqrt(np.diag(fits_2g[0][1]))[1]
first_step_err_2 = np.sqrt(np.diag(fits_2g[0][1]))[4]

second_step_err_1 = np.sqrt(np.diag(fits_2g[1][1]))[1]
second_step_err_2 = np.sqrt(np.diag(fits_2g[1][1]))[4]

third_step_err_1 = np.sqrt(np.diag(fits_2g[2][1]))[1]
third_step_err_2 = np.sqrt(np.diag(fits_2g[2][1]))[4]

fourth_step_err_1 = np.sqrt(np.diag(fits_2g[3][1]))[1]
fourth_step_err_2 = np.sqrt(np.diag(fits_2g[3][1]))[4]

fifth_step_err_1 = np.sqrt(np.diag(fits_2g[4][1]))[1]
fifth_step_err_2 = np.sqrt(np.diag(fits_2g[4][1]))[4]

sixth_step_err_1 = np.sqrt(np.diag(fits_2g[5][1]))[1]
sixth_step_err_2 = np.sqrt(np.diag(fits_2g[5][1]))[4]

seventh_step_err_1 = np.sqrt(np.diag(fits_2g[6][1]))[1]
seventh_step_err_2 = np.sqrt(np.diag(fits_2g[6][1]))[4]

first_step_fit_1 = fits_2g[0][0][1]
first_step_fit_2 = fits_2g[0][0][4]

second_step_fit_1 = fits_2g[1][0][1]
second_step_fit_2 = fits_2g[1][0][4]

third_step_fit_1 = fits_2g[2][0][1]
third_step_fit_2 = fits_2g[2][0][4]

fourth_step_fit_1 = fits_2g[3][0][1]
fourth_step_fit_2 = fits_2g[3][0][4]

fifth_step_fit_1 = fits_2g[4][0][1]
fifth_step_fit_2 = fits_2g[4][0][4]

sixth_step_fit_1 = fits_2g[5][0][1]
sixth_step_fit_2 = fits_2g[5][0][4]

seventh_step_fit_1 = fits_2g[6][0][1]
seventh_step_fit_2 = fits_2g[6][0][4]

shift1s = [first_step_fit_1,second_step_fit_1,third_step_fit_1,fourth_step_fit_1,
         fifth_step_fit_1,sixth_step_fit_1,seventh_step_fit_1]

shift2s = [first_step_fit_2,second_step_fit_2,third_step_fit_2,fourth_step_fit_2,
         fifth_step_fit_2,sixth_step_fit_2,seventh_step_fit_2]

err1s = [first_step_err_1,second_step_err_1,third_step_err_1,fourth_step_err_1,
         fifth_step_err_1,sixth_step_err_1,seventh_step_err_1]

err2s = [first_step_err_2,second_step_err_2,third_step_err_2,fourth_step_err_2,
         fifth_step_err_2,sixth_step_err_2,seventh_step_err_2]

fig,ax = plt.subplots()
ax.errorbar(np.linspace(1,7,7),shift1s,yerr=err1s)
ax1=ax.twinx()
ax1.errorbar(np.linspace(1,7,7),shift2s,yerr=err2s)
fig.show








