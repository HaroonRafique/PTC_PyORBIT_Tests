/******************************************************************************************
 * CALCULATE COMPENSATION OF INJECTION CHICANE BETA BEATING
 *
 ******************************************************************************************/

EOPTION, ADD=false;

if(correction_option == 0){
  kd3  := kkd + kd_trm;
  kd14 := kkd + kd_trm;
  BR.QDE3,  K1 :=  kd3;
  BR.QDE14, K1 :=  kd14;
  };
if(correction_option == 1){ 
  kd3  := kkd + kd_trm + kd_bal;
  kd14 := kkd + kd_trm - kd_bal;
  BR.QDE3,  K1 :=  kd3;
  BR.QDE14, K1 :=  kd14;
  };



! install a marker (needed for matching)
use, sequence = psb;
twiss;
seqedit, sequence=psb;
  mymarker: marker;
  flatten;
  ! install, element= mymarker, at = table(summ,Length)/2., from=BI1.TSTR1L1;
  install, element= mymarker, at = 1.32695+BR.BHZ82->L/2, from=BR.BHZ82;
  flatten;
endedit;
use, sequence = psb;

SELECT, FLAG=PTC_TWISS, COLUMN=NAME, KEYWORD, L, S, BETX, BETY, ALFX, ALFY, Dispx, mu1, mu2, x, y;
exec, ptc_twiss_macro;
PLOT, table=ptc_twiss, HAXIS=S, VAXIS1=BETX,BETY, VAXIS2=X, COLOUR=100, VMIN=0,-0.08, VMAX=25,0.02, file='../output/psb_injection', interpolate;
write, table=ptc_twiss, file='reference_twiss.tfs';
readmytable, file='reference_twiss.tfs', table=reference_twiss; 
system,'rm reference_twiss.tfs';
kkf_ref = kkf;
kkd_ref = kkd;


Qx   := Table(ptc_twiss, PSB$end, mu1);
Qy   := Table(ptc_twiss, PSB$end, mu2);
alfy := Table(ptc_twiss, mymarker, ALFY);

BSW1_K0 := + BSW_K0L/BI1.BSW1L1.1->L;
BSW2_K0 := - BSW_K0L/BI1.BSW1L1.2->L;
BSW3_K0 := - BSW_K0L/BI1.BSW1L1.3->L;
BSW4_K0 := + BSW_K0L/BI1.BSW1L1.4->L;
BSW1_K2 := + BSW_K2L/BI1.BSW1L1.1->L;
BSW2_K2 := - BSW_K2L/BI1.BSW1L1.2->L;
BSW3_K2 := - BSW_K2L/BI1.BSW1L1.3->L;
BSW4_K2 := + BSW_K2L/BI1.BSW1L1.4->L;
create, table=mytable, column=BSW_t, BSW_K0L, BSW1_K0, BSW2_K0, BSW3_K0, BSW4_K0, 
								BSW_K2L, BSW1_K2, BSW2_K2, BSW3_K2, BSW4_K2, kkf, 
								kkd, kd_trm, kd_bal, kd3, kd14, Qx, Qy, alfy, betabeaty_max;


readtable, file='../strength/BSW_collapse.tfs';
nrows = table(BSWTABLE, tablelength);

SETVARS,TABLE=BSWTABLE,ROW=1;
exec, assign_BSW_strength;
if (use_BSW_alignment==1){exec, assign_BSW_alignment;};
exec, ptc_twiss_macro;
PLOT, table=ptc_twiss, HAXIS=S, VAXIS1=BETX,BETY, VAXIS2=X, COLOUR=100, VMIN=0,-0.08, VMAX=25,0.02, file='../output/psb_injection', interpolate;


row = 0;
while(row<nrows){
! while(row<24){
  
   row = row + 1;
   iteration = 0;

   kd_bal = 0;

   while(iteration<2){ 
     value, iteration;

     ! take reference of bety at marker
     option, -info;
     BSW_K0L = 0;
     BSW_K2L = 0;
     kd_trm_tmp = kd_trm;
     kd_bal_tmp = kd_bal;
     kd_trm = 0;
     kd_bal = 0;
     use, sequence=PSB;
     ! kkf = kkf_ref;
     ! kkd = kkd_ref;

     exec, ptc_twiss_macro;
     bety_target = Table(ptc_twiss, mymarker, BETY);
     alfy_target = Table(ptc_twiss, mymarker, ALFY);

     kd_trm = kd_trm_tmp;
     kd_bal = kd_bal_tmp;
     SETVARS,TABLE=BSWTABLE,ROW=row;
     exec, assign_BSW_strength;
     if (use_BSW_alignment==1){exec, assign_BSW_alignment;};

     value, BSW_T, BSW_K0L, BSW_K2L;
     option, info;

     option, -warn;

    
     MATCH, USE_MACRO;
      VARY, NAME = kd_trm, STEP = 1e-6;
      USE_MACRO, name=ptc_twiss_macro;
      Constraint, EXPR=Table(ptc_twiss, mymarker, BETY) = bety_target;
      value, bety_target;
      LMDIF, TOLERANCE = 1.0E-18;
      ! JACOBIAN,CALLS=1000,TOLERANCE=1.0E-15,STRATEGY=3; !,COOL=real,BALANCE=real, random=real;
     ENDMATCH;
     

     if(correction_option > 0){
         MATCH, USE_MACRO;
          VARY, NAME = kd_bal, STEP = 1e-6;
          USE_MACRO, name=ptc_twiss_macro;
          Constraint, EXPR=Table(ptc_twiss, mymarker, ALFY) = alfy_target;
          ! LMDIF, TOLERANCE = 1.0E-18;
          JACOBIAN,CALLS=1000,TOLERANCE=1.0E-18,STRATEGY=3; !,COOL=real,BALANCE=real, random=real;
         ENDMATCH;
       };
      
      /*******
      kd_trm_max = -0.004;
      kd_trm = kd_trm_max - kd_trm_max/0.005*BSW_T;
      ******/

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! activate to achieve only partial compensation
     ! kd_trm = kd_trm * 0.5;
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     MATCH, USE_MACRO;
      VARY, NAME = kkf, STEP = 1e-8;
      VARY, NAME = kkd, STEP = 1e-8;
      USE_MACRO, name=ptc_twiss_macro;
      CONSTRAINT, EXPR=Table(ptc_twiss, PSB$end, mu1) = QH;
      CONSTRAINT, EXPR=Table(ptc_twiss, PSB$end, mu2) = QV;
      LMDIF, TOLERANCE = 1.0E-18;
     ENDMATCH;  

     option, warn;

     iteration = iteration + 1;
   }

   option, -info;
   exec, calculate_betabeat(ptc_twiss,reference_twiss);
   option, info;
   value, betabeaty_max;
   fill, table=mytable;
   ! write, table=mytable, file='../output/BSW_betabeat_correction.tfs';


   if(plot_corrected_optics ==1){
      PLOT, table=ptc_twiss, HAXIS=S, VAXIS1=BETX,BETY, VAXIS2=X, COLOUR=100, VMIN=0,-0.08, VMAX=25,0.02, file='../output/psb_injection', interpolate;
   	};
};

write, table=mytable, file='../output/BSW_betabeat_correction.tfs';


plot, table=mytable, HAXIS=BSW_t, VAXIS=KKD, KD3, KD14;
plot, table=mytable, HAXIS=BSW_t, VAXIS=betabeaty_max;
plot, table=mytable, HAXIS=BSW_t, VAXIS=BSW_K0L, BSW_K2L;

