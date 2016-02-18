#path='Batch/'		# this is default path [note it is from where the .x is executed (i.e. in Runs/PHF)]

# ". test.sh" -> '. ' means 'source' which means run in same shell (i.e. retain access to all variables)

############################################# --- README --- #############################################
# Check your shell and interaction files, both here and in auto.sh!  Check auto.sh to make sure the order
#   is correct (i.e. One, Full, Random, etc. AND that variables are assigned appropriately).
##########################################################################################################

nsd=400
npts=25
tolerance=0.01
#start=15
#end=40
#step=5
path='Batch/'
output='Output/'
type='varNptsGraph'

#
# Begin EVEN A
#

JMIN=0
JMAX=10
decimal=''

# NUCLIDE Al28

	# FULL
	shell='sd'
	interaction='usdb'
	nuclide='Al28'
	pnum=5
	nnum=7
	seed=223
#	. $path'auto.sh'

#
# Begin ODD A
#

JMIN=0.5
JMAX=10.5
decimal='.5'

# NUCLIDE Na23

	# FULL
	shell='sd'
	interaction='usdb'
	nuclide='Na23'
	pnum=3
	nnum=4
	seed=153
	. $path'auto.sh'

clear
echo .
echo ..
echo ...
DATE1=`date +%D`
DATE2=`date +%H:%M:%S`
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo '                                 'Completed!
echo '                                ' on $DATE1
echo '                                ' at $DATE2
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo All finished!


# Backup + email
echo
echo Now backing up and emailing ...
echo 

mv *.dat Batch/data/

#touch EnergyHartreeFock.dat

echo Done!

# run me by typing: bash Batch/main.sh | tee Batch/terminal_output.txt

















