#! /usr/bin/python


def printKick(element, L, B, Bhro):
    
    kick  = L * B / Brho # rad
    kick /= 2 # considering the powering scheme the current is split in ring 1-2 and ring 3-4
    kick *= 1000. # mrad
    
    print '%s, kick = %f mrad' % (element,kick)

# ===============================================================================

momentum = 2.141756 # GeV/c for a proton with E_kin=1.4
Brho = 3.3356 * momentum # T * m

L_MDBAA = 0.470 # m
L_MDBAB = 0.366 # m

# From LUC Sermeus measured during LS1
Calibration_MDBAA = 1.20055E-4 # T * m / A
Calibration_MDBAB = 1.56247E-4 # T * m / A

Current_BSW14L4 = 524.1 # A
Current_BSW15L1 = 519.1 # A
Current_BSW15L4 = 525.1 # A

printKick("BSW14L4", L_MDBAA, (Calibration_MDBAA*Current_BSW14L4/L_MDBAA), Brho)

printKick("BSW15L1", L_MDBAB, (Calibration_MDBAB*Current_BSW15L1/L_MDBAB), Brho)

printKick("BSW15L4", L_MDBAA, (Calibration_MDBAA*Current_BSW15L4/L_MDBAA), Brho)

# ===============================================================================


