 get_betax_value(tablename,row): macro = {betax=table(tablename,betx,row);}
 get_betay_value(tablename,row): macro = {betay=table(tablename,bety,row);}

 calculate_betabeat(twisstable_dis,twisstable_ref) : macro = {
    nrows_2check = table(twisstable_dis,tablelength);
    betabeatx_max = 0;
    betabeaty_max = 0;
    n = 0;
    row_ = 0;
    while(row_<nrows_2check){
      row_ = row_+1;
      exec, get_betax_value(twisstable_ref,$row_);
      exec, get_betay_value(twisstable_ref,$row_);
      bx_reference = betax;
      by_reference = betay;
      exec, get_betax_value(twisstable_ref,$row_);
      exec, get_betay_value(twisstable_dis,$row_);
      bx_distorted = betax;
      by_distorted = betay;
      betabeatx = (bx_distorted-bx_reference)/bx_reference;
      betabeaty = (by_distorted-by_reference)/by_reference;
      if (abs(betabeatx)>betabeatx_max) {betabeatx_max=betabeatx;}
      if (abs(betabeaty)>betabeaty_max) {betabeaty_max=betabeaty;}
    }
 }
