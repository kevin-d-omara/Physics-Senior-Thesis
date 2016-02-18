DATE1=`date +%D`
DATE2=`date +%H:%M:%S`
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo '                     'beginning $type SHERPA on $nuclide
echo '                                ' on $DATE1
echo '                                ' at $DATE2
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sherpa14.x << INPUT
$shell
$pnum
$nnum
n
1
$interaction
1. 1. 1. 1.
n
$seed
$nsd
10
-1
w
$nuclide
o
x
x
x
INPUT
