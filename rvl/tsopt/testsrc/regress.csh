#!/bin/csh
echo " "
echo "*************************"
echo "* TSOPT REGRESSION TEST *"
echo "*************************"
echo " "

echo "                *************************" >  regress.out
echo "                * TSOPT REGRESSION TEST *" >> regress.out
echo "                *************************" >> regress.out
echo " "                                            >> regress.out
foreach i ( `/bin/ls -1 *.x` )
    echo " " >> regress.out
    ./$i >>& regress.out
    if ( $? == 1 ) then
        echo "$i FAILED"
    else 
        echo "$i succeeded"
    endif
end
