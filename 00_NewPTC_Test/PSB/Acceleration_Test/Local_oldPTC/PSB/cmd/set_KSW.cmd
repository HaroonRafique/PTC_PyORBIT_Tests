assign_KSW_strength: macro = {
     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.KSW1L4;
     EFCOMP, DKN:={kLBIKSW1L4};   

     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.KSW2L1;
     EFCOMP, DKN:={kLBIKSW2L1};   

     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.KSW16L1;
     EFCOMP, DKN:={kLBIKSW16L1};   

     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.KSW16L4;
     EFCOMP, DKN:={kLBIKSW16L4};   
 };
