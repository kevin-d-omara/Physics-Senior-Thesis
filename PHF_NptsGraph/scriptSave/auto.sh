# Sherpa then phf.x for many npts
#

# SHERPA
. $path'sherpa.sh' | tee $path$output$nuclide'_Sherpa.txt'

# create blank files to be populated
directory=''
score='_'
extension='.dat'
for (( temp=0; temp<=10; temp=temp+1 ))
do
    touch $directory$nuclide$score$temp$decimal$extension
done

# PHF npts
start=10
end=50
step=1
for (( npts=start; npts<=end; npts=npts+step ))
do
	type=$npts
	. $path'phf_nptsGraph.sh' | tee $path$output$nuclide'_'$type'.txt.'
done
