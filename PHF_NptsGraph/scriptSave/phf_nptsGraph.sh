DATE1=`date +%D`
DATE2=`date +%H:%M:%S`
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo '                    'beginning $type phf_nptsGraph on $nuclide
echo '                                ' on $DATE1
echo '                                ' at $DATE2
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
outputType='_phf_nptsGraph_'
phf_nptsGraph.x << INPUT
j
$shell
$interaction
1. 1. 1. 1.
end
$pnum
$nnum
1
$nuclide
$JMIN
$JMAX
$tolerance
$npts
y
$nuclide
o
x
INPUT


#$nuclide$outputType$npts
