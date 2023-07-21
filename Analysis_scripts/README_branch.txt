Differences to the main branch:

--------------------------------------

New files:
labelstoplot.m : add labels to the plot from TDTR_Bidirectional_MAIN_FIT.m

voltage_check.m : Plots for Vin, Vout, ratio, different averagings

----------------------------------------------------

Changes to files:
TDTR_Bidirectional_SUB_C.m : Implemented a penalty - if any parameters are negative, the RMSE is set to
				infinity. Effectively prevents negative parameter values to converge.

TDTR_Bidirectional_MAIN_SIM.m : Basically Tarmo's version of the script. Is only used for sensitivity
				analysis.

TDTR_Bidirectional_MAIN_FIT.m : Added a possibility to use different averagings.