/******************************************************************************************
 * MATCH KSW bump 
 *
 ******************************************************************************************/

call, file='../cmd/set_KSW.cmd';

/*******
 MATCH, USE_MACRO;
	VARY, NAME=kLBIKSW1L4,  STEP=1.E-8;
	VARY, NAME=kLBIKSW2L1,  STEP=1.E-8;
	VARY, NAME=kLBIKSW16L1,  STEP=1.E-8;
	VARY, NAME=kLBIKSW16L4,  STEP=1.E-8;
	USE_MACRO, name=assign_KSW_strength;
  USE_MACRO, name=ptc_twiss_macro;
	CONSTRAINT, EXPR=Table(ptc_twiss, PSB$END, X) = -0.035;
	CONSTRAINT, EXPR=Table(ptc_twiss, PSB$END, PX) = 0.00;
	CONSTRAINT, EXPR=Table(ptc_twiss, MYMARKER, X) = 0.00;
	CONSTRAINT, EXPR=Table(ptc_twiss, MYMARKER, PX) = 0.00;
	LMDIF, CALLS=1000, TOLERANCE=1.E-22;
 ENDMATCH;
 value, kLBIKSW1L4, kLBIKSW2L1, kLBIKSW16L1, kLBIKSW16L4;
 stop;
*******/


 // for 4.22/4.45
 ! klbiksw1l4  := 1.51009e-03 * ksw_factor; //max allowed 1.85e-3
 ! klbiksw2l1  := 5.02462e-03 * ksw_factor; //max allowed 6.60e-3
 ! klbiksw16l1 := 5.63305e-03 * ksw_factor; //max allowed 6.70e-3
 ! klbiksw16l4 := 1.09381e-03 * ksw_factor; //max allowed 2.05e-3

 // for 4.30/4.45
 klbiksw1l4  := 1.82636e-03 * ksw_factor ; //max allowed 1.85e-3
 klbiksw2l1  := 4.92291e-03 * ksw_factor ; //max allowed 6.60e-3
 klbiksw16l1 := 5.66090e-03 * ksw_factor ; //max allowed 6.70e-3
 klbiksw16l4 := 1.33963e-03 * ksw_factor ; //max allowed 2.05e-3

 // for 4.40/4.45
 ! klbiksw1l4  := 2.22258e-03 * ksw_factor; //max allowed 1.85e-3
 ! klbiksw2l1  := 4.79289e-03 * ksw_factor; //max allowed 6.60e-3
 ! klbiksw16l1 := 5.69721e-03 * ksw_factor; //max allowed 6.70e-3
 ! klbiksw16l4 := 1.65185e-03 * ksw_factor; //max allowed 2.05e-3

 MATCH, USE_MACRO;
  VARY, NAME=ksw_factor,  STEP=1.E-8;
  USE_MACRO, name=assign_KSW_strength;
  USE_MACRO, name=ptc_twiss_macro;
  CONSTRAINT, EXPR=Table(ptc_twiss, PSB$END, X) = -0.035;
  LMDIF, CALLS=1000, TOLERANCE=1.E-22;
 ENDMATCH;

 PLOT, table=ptc_twiss, HAXIS=S, VAXIS1=BETX,BETY, VAXIS2=X, COLOUR=100, VMIN=0,-0.08, VMAX=25,0.02, file='../output/psb_injection', interpolate;
 value, kLBIKSW1L4, kLBIKSW2L1, kLBIKSW16L1, kLBIKSW16L4;
 kBIKSW1L4  = kLBIKSW1L4  / BI1.KSW1L4->L;
 kBIKSW2L1  = kLBIKSW2L1  / BI1.KSW2L1->L;
 kBIKSW16L1 = kLBIKSW16L1 / BI1.KSW16L1->L;
 kBIKSW16L4 = kLBIKSW16L4 / BI1.KSW16L4->L;
 x_bump = table(ptc_twiss, PSB$END, X);

 create, table=KSW, column = x_bump, kBIKSW1L4, kBIKSW2L1, kBIKSW16L1, kBIKSW16L4;
 fill, table=KSW;
 write, table=KSW, file='../output/KSW.tfs';

 kLBIKSW1L4 = 0;
 kLBIKSW2L1 = 0;
 kLBIKSW16L1 = 0;
 kLBIKSW16L4 = 0;
 use, sequence = PSB;
