#! /bin/csh
#

set f2test = (psb_injection.madx \
   psb_orbit.madx \
   psb_extraction_QH4.17_QV4.23.madx)

# echo $f2test
# echo $f2test[1]
# echo $f2test[20]

echo "Executing MAD-X for the input files"
foreach f ($f2test)
  echo ""
  echo ""
  echo "   Input file :" $f
  echo ""
  echo ""
  if (-f $f.log) then
    rm -rf $f.log
  endif
  ~/MAD-X/madx < ../madx/$f > $f.log
  grep -i warning $f.log
  grep -i error $f.log
  echo ""
  echo -n "Enter return to continue with next file or 0 to exit  "
  set word = $<
  if ( $word == 0 ) break
end
