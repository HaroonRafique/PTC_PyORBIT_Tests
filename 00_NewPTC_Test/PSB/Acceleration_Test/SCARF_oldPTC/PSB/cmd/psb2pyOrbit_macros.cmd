! Miscellanous macros to for converting the psb madx lattice to a PTC
! flat file and  related plotting etc. 

! 2019-04-17

! elisabeth.renner@cern.ch
! partly from /afs/cern.ch/user/p/pyorbit/public/pyorbit_examples_repository/PyORBIT_Examples/Create_A_Lattice/Flat_file.madx
!----------------------------------------------------------------------
! List of macros

! remove_elements(rm_option)                 : Macro to remove elements required for PTC-pyOrbit simulations
! edit_cavities()                            : Macro to remove def. cavities and insert 1 cavity element for pyOrbit 
! ptc_twiss_tune_macro_false_split()         : PTC macro to match tune
! myTWISS(filename)                          : Macro to generate twiss table with corrected dispersion dispx
! install_fixed_mask()                       : Macro to install fixed absorber downstream of movable absorber
!
!
!

!**********************************************************************
!**********************************************************************
!**********************************************************************
! Macro to remove elements required for PTC-pyOrbit simulations
! 	- select the desired elements to be removed:
!	- rm_option=
!		- 1: all elements of class instrument, monitor, marker*
!		- 2: user specified element selection: manually deselect/select
!		     classes to be removed (e.g. instruments, marker, monitor
!		     apart from wirescanner)
!		- 3: all elements of class marker
!		- 4: all elements of class marker, monitor, instrument, and
!		     stscrap (dipole scrapers). Subsequent redefinition of dipole
!		     aperture to stscrap aperture*
!	  *BI1.TSTR1L1 (charge exchange foil, class instrument) is reinstalled. 
!
remove_elements(rm_option): macro = {
  ! Case 1: Instruments, marker, monitor			   
  if (rm_option == 1) {
    use, sequence=psb;
    select, flag = seqedit, clear;
    select, flag = seqedit, class = marker;
    select, flag = seqedit, class = instrument;		       
    select, flag = seqedit, class = monitor;
    select, flag = seqedit, class = MDBAA; ! extraction bumpers
    select, flag = seqedit, class = MDBAB; ! extraction bumpers
    seqedit, sequence = psb;
      flatten;
      remove, element = selected;
      install, element = BI1.TSTR1L1,at=-0.226,from=BI1.BSW1L1.3;
      flatten;
    endedit;
  }
  ! Case 2:  user specified, here: instruments, marker, monitor APART FROM FORE SCANNERS  
  elseif (rm_option == 2) {
    use, sequence=psb;
    select, flag = seqedit, clear;
    //---------------------- MARKER        ---------------------------------------------
    
    select, flag = seqedit, class = marker;
    //---------------------- MONITOR        ---------------------------------------------
    select, flag = seqedit, class = BCSEA;    ! Beam Current Stripping-foil Efficiency, type A
    select, flag = seqedit, class = BCSMA;    ! Beam Current Stripping foil Monitor, type A
    select, flag = seqedit, class = BTVBP016  ! Beam TV Booster Type, Pneumatic, variant 016
    select, flag = seqedit, class = BTVTS;    ! Beam observation TV for Stripping foil (Booster)
    ! Do not remove TKSTR002 (foil)!
    //---------------------- MONITOR        ---------------------------------------------
    select, flag = seqedit, class = BPMTA;    ! Beam Position Monitor for Tune, type A
    select, flag = seqedit, class = BPPBA001; ! Beam Position Phase Pick-up, Booster type A
    select, flag = seqedit, class = BPPBB001; ! Beam Position Phase Pick-up, Booster type B
    ! select, flag = seqedit, class = BWSH;     ! Wire Scanner Profile Monitor, Horizontal
    ! select, flag = seqedit, class = BWSRA002; ! Beam Wire Scanner, Rotational, type A (vertical)
    ! select, flag = seqedit, class = BWSV;     ! Wire Scanner Profile Monitor, Vertical
    select, flag = seqedit, class = KQM__002; ! Q kicker tune measurement
    select, flag = seqedit, class = TDC;      ! DC Transformer
    ! The following three beam transformers (TFA, TMD and TSW) are inside the same casing. It is the casing that is 0.5293 m long.
    select, flag = seqedit, class = TFA;     ! Beam transformer slow, transfo fast (+transfo slow + transfo intermediaire, but obsolete electronics)
    select, flag = seqedit, class = UES;     ! UES Pick-up Monitor
    select, flag = seqedit, class = UES__003 ! UES Pick-up MonitorU
    select, flag = seqedit, class = QD;      ! UQD Pick-up Monitor quadrupolar
    select, flag = seqedit, class = UTFB;    ! PU Transverse feedback, electrostatic PU; type F
    select, flag = seqedit, class = UWBL;    ! Wide band PU, PU wide band (BSM, tomoscope)
    //---------------------- Extraction Bumper  ------------------------------------------
    select, flag = seqedit, class = MDBAA; ! extraction bumpers
    select, flag = seqedit, class = MDBAB; ! extraction bumpers
    seqedit, sequence = psb;
      flatten;
      remove, element = selected;
      flatten;
    endedit;
  }
  ! Case 3:  markers
  elseif (rm_option == 3){
    use, sequence=psb;
    select, flag = seqedit, clear;
    select, flag = seqedit, class = marker;
    seqedit, sequence = psb;
      flatten;
      remove, element = selected;
      flatten;
    endedit;
  }
  ! Case 4: markers, instruments, monitors and stscrap with dipole aperture
  ! redefinition
  elseif (rm_option == 4){
    use, sequence=psb;
    select, flag = seqedit, clear;
    select, flag = seqedit, class = marker;
    select, flag = seqedit, class = instrument;		       
    select, flag = seqedit, class = monitor;
    select, flag = seqedit, class = STSCRAP;
    //---------------------- Extraction Bumper  ------------------------------------------
    select, flag = seqedit, class = MDBAA; ! extraction bumpers
    select, flag = seqedit, class = MDBAB; ! extraction bumpers
    call, file = '../aperture/reducedBHZap.dbx';
    seqedit, sequence = psb;
      flatten;
      remove, element = selected;
      install, element = BI1.TSTR1L1,at=-0.226,from=BI1.BSW1L1.3;
      flatten;
    endedit;
  }
};
!**********************************************************************
! Macro to
! 	1) remove defined cavities in lattice and
! 	2) insert cavity element for pyOrbit simulations 
!
! Note:  
! 	 - position not relevant -> installed in period 5, downstream of
!   	 BR.BHZ51 (position in reference to BR.BHZ51)
! 	 - length 1 m
! 	 - Name: BR.C02 to be consistent with previous simulations. 
! 	  - no_cavity_totalpath ! 
! 	  (reference: H. Bartosik, 2019/03)
! 
edit_cavities(): macro = {
    use, sequence=psb;
    select, flag = seqedit, clear;
    select, flag = seqedit, class = ACWFB;
    seqedit, sequence = psb;
      flatten;
      remove, element = selected;
      flatten;
    endedit;
    
    use, sequence=psb;
    CAVFINEM	: RFCAVITY, L:= 1.0, VOLT:= 0.000, HARMON = 1, LAG:= 0.0, no_cavity_totalpath;
    ! Exact position not important. Install downstream of BHZ51 
    seqedit, sequence = psb;
      flatten;
      install, element = BR.C02, class=CAVFINEM, at=-2.0 ,from=BR.BHZ51;
      flatten;
    endedit;
};
!*********************************************************************
! Macro to generate twiss table with corrected dispersion dispx
! 	Delta_x=dispx*Delta_p/p;
!
myTWISS(filename): macro = {
  use, sequence = psb; 
  select, flag=twiss, clear;
  beta=sqrt(1-1/beam->gamma^2);
  twiss, sequence = psb,centre=true, table;
  dispx:=beta*table(twiss,dx); ! Delta_x=dispx*Delta_p/p
  select, flag=twiss, column=name,keyword,parent, s,l,angle,x,alfx,alfy,
  	  	      betx,bety,dx,dispx,dy,apertype,aper_1,aper_2,aper_3,
		      aper_4,mux,muy;	    
  twiss, sequence = psb,centre=true, file = filename, table;
  value, table(twiss,PSB$END,mux),table(twiss,PSB$END,muy);
};
!*********************************************************************
! PTC macro to match tune, taken from 
!
ptc_twiss_macro: macro={
  beta=sqrt(1-1/beam->gamma^2);
  dispx:=beta*table(ptc_twiss,disp1);
  ptc_create_universe;
  ptc_create_layout, model=2, method=2, nst=5, exact=true;
  ptc_setswitch, debuglevel=0, nocavity=false, fringe=true,
             exact_mis=true, time=true, totalpath=true;
  PTC_ALIGN;
  ptc_twiss, closed_orbit, table=ptc_twiss, icase=4, no=2, summary_table=ptc_twiss_summary;
  qx0=table(ptc_twiss_summary,Q1);
  qy0=table(ptc_twiss_summary,Q2);
  value, qx0, qy0;
  ptc_end;
};
!*********************************************************************
! Macro to install fixed absorber downstream of small mask
!   - positon and aperture taken from EDMS 1578463
!
install_fixed_mask(): macro{
  use, sequence = psb;
  seqedit, sequence = psb;
  flatten;
  install, element=BR1.TSAB8L4, class=TSAMA, AT = 0.47355 + BR.QDE8->L/2, from = BR.QDE8; ! AT = 0.1765, from = BR1.TSA8L4; 
  flatten;
  endedit;
  BR1.TSAB8L4, APERTYPE=RECTELLIPSE, APERTURE={0.0480, 0.0370, 1.0, 1.0}; ! RECTANGLE
  use, sequence = psb;
};
!*********************************************************************
! Macro to remove small mask
!   - positon and aperture taken from EDMS 1578463
!
remove_small_mask(): macro{
  use, sequence = psb;
  seqedit, sequence = psb;
  flatten;
  remove, element=BR1.TSA8L4; 
  flatten;
  endedit;
  use, sequence = psb;
};
!*********************************************************************
! Macro to generate flatfile 
!
generate_flatfile: macro= {
  ptc_create_universe;
  ptc_create_layout, model=2, method=2, nst=5,exact=true, resplit;
  ptc_setswitch, debuglevel=0, nocavity=false, fringe=true,
       exact_mis=true, time=true, totalpath=true;
  ptc_align;
  ptc_script, file="../cmd/resplit.ptc";
  ptc_script, file="../cmd/print_flat_file.ptc";
  ptc_twiss, icase=5, no=4, deltap_dependency, closed_orbit, table=ptc_twiss;
  ptc_end;
};

