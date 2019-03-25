#!/usr/bin/bash
SDIR='.' #Source code directory
EXE='streamlinesnt3d' #Executable name
DATADIR='.' #Data directory
#ddd
SID='TEST' #Simulation ID (name of the specific directory)
OPTION='SE'
MODE=0
if [ $# -eq 1 ]
then
	OPTION=$1
elif [ $# -eq 2 ]
then
	OPTION=$1
	SID=$2
elif [ $# -eq 3 ]
then
	OPTION=$1
	SID=$2
	MODE=$3
fi

#SIMDIR='../simulations/'$SID #Simulation directory
SIMDIR='.'
if [ ! -d $SIMDIR ]
then
	mkdir $SIMDIR
fi

####SIMULATION PARAMETERS####
#For clusters
RUNINGTIME='00:30'
NCORES=1
export OMP_NUM_THREADS=$NCORES
BOPTIONS=$BOPTIONS"-n ${NCORES} -W ${RUNINGTIME}"

MEM=2048 #[MB per cores]
SCRATCHMEM=128 #[MB per cores]
#BOPTIONS=$BOPTIONS" -R "'"'"rusage[mem=${MEM},scratch=${SCRATCHMEM}]"'"'

JOBNAME="AH0$SID"
BOPTIONS=$BOPTIONS" -J ${JOBNAME}"
echo boption : $BOPTIONS
#For program
Q1FILE=$DATADIR'/Ux.dat'
Q2FILE=$DATADIR'/Uy.dat'
Q3FILE=$DATADIR'/Uz.dat'
GFILE=$DATADIR'/voxelMesh.dat'
DM=0.0E-9
UREF=1.0
TREF=1.0
NX=300
NY=300
NZ=300
DX=3.0E-6 #$(awk "BEGIN {print 1.0/${NX}"})
RFLAG=0
#One stream line computation (mode=0)
SLSMAX=200
TMAX=0.0E2
CSF=1
#Particle plume transport (mode=1,2)
NP=100
NT=4
TL=(0.0 10.0 20.0 30.0 40.0 50.0)

#Automatic naming
if [ $MODE -eq 0 ]
then
	SNAME='M'$MODE'_DM'$DM'_SLSM'$SLSMAX'_TMAX'$TMAX'_CSF'$CSF
else
	SNAME='M'$MODE'_DM'$DM'_NP'$NP'_NT'$NT
fi
#IFILE=$SIMDIR'/input_'$SNAME'.in' #Input file
IFILE='input.in'
#OFILE=$SIMDIR'/sl_out_'$SNAME'.dat' #Output file
OFILE='streamlines.txt'
BSUBFILE=$SIMDIR'/bsub_out_'$SNAME'.txt' #Bsub output file

####ANSWERS TO INTERACTIVE PROGRAM####
cat > $IFILE << EOF
$MODE
$Q1FILE
$Q2FILE
$Q3FILE
$GFILE
$DM
$UREF $TREF
$NX $NY $NZ
$DX
$RFLAG
EOF

####MODE VARIATION####
if [ $MODE -eq 0 ]
then
	printf "%s\n" "$OFILE" >> $IFILE
	printf "%s\n" "$SLSMAX" >> $IFILE
	printf "%s\n" "$TMAX" >> $IFILE
	printf "%s\n" "$CSF" >> $IFILE
else
	printf "%s\n" "$NP" >> $IFILE
	printf "%s\n" "$NT" >> $IFILE
	printf "%s\n" "${TL[@]}" >> $IFILE
fi

####COMPILATION####
#cd $SDIR
make clean
make all
#cd ../workspace

####EXECUTION####
#STRAIGHTFORWARD EXECUTION
if [ $OPTION == 'SE' ]
then
	./$EXE < $IFILE -n $NCORES
#INTERACTIVE
elif [ $OPTION == 'I' ]
then
	bsub $BOPTIONS -I ./$EXE < $IFILE
#OUTPUT FILE OVERWRITE VERSION
elif [ $OPTION == 'OO' ]
then
	bsub $BOPTIONS -oo $BSUBFILE ./$EXE < $IFILE
#OUPUT FILE DEFAULT VERSION
elif [ $OPTION == 'O' ]
then
	bsub $BOPTIONS -o $BSUBFILE ./$EXE < $IFILE
fi

