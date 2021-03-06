assign_BSW_strength: macro = {
     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.1;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.4;
     EFCOMP, DKN:={+BSW_K0L, 0, +BSW_K2L};   

     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.2;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.3;
     EFCOMP, DKN:={-BSW_K0L, 0, -BSW_K2L}; 
 };

 assign_BSW_alignment: macro = {
     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.1;
     EALIGN, DX=-0.0057;

     SELECT,FLAG=ERROR,CLEAR;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.2;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.3;
     SELECT,FLAG=ERROR,PATTERN=BI1.BSW1L1.4;
     EALIGN, DX=-0.0442;
 };
