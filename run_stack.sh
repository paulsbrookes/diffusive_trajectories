#!/bin/sh
STACK='stack'
N_TEMP_FOLDERS=(./$STACK'_temp'*)
N_TEMP_FOLDERS=${#N_TEMP_FOLDERS[@]}
NEXT_TEMP_FOLDER=$((N_TEMP_FOLDERS+1))
TEMP=$STACK'_temp'_$NEXT_TEMP_FOLDER
find ./$STACK -name "*.cfg~" -type f -delete
rm -r $TEMP
cp -r sourcecode $TEMP
N_FILES=(./$STACK/*)
N_FILES=${#N_FILES[@]} 
while [ $N_FILES -gt 0 ]
do
    shopt -s nullglob
    FILES=( ./$STACK/* )
    source ${FILES[0]}
    rm $TEMP'/parameters.h'
    echo "string folder = \"$folder\";
    const double endtime = $endtime;
    const int cavity_levels = $cavity_levels;
    const int transmon_levels = $transmon_levels;
    const int snap_to_snap = $snap_to_snap;
    const int bins = $bins;
    const int snapshots_number = $snapshots_number;
    double kappa = $kappa;
    double gamma_c = $gamma_c;
    double g = $g;
    double dt = endtime/bins;
    double epsilon = $epsilon;
    double omega_c = $omega_c;
    double omega_q = $omega_q;
    double omega_d = $omega_d;
    double dc = omega_c - omega_d;
    double dq = omega_q - omega_d;
    const int samplesize = $samplesize;
    const double theta = $theta;
    const double phi = $phi;
    const int photon_number = $photon_number;
    const double chi = $chi;" > $TEMP'/parameters.h'
    rm ${FILES[0]}
    rm $TEMP'/sim_9'
    g++ -std=c++0x -fopenmp -mcmodel=large -g -O -lgsl -lgslcblas 'D_files'/$transmon_levels't_'$cavity_levels'c'/functions.o $TEMP'/diffusive_sim.cpp' -o $TEMP'/sim'
    ulimit -s 10000000
    $TEMP/sim
    find ./$STACK -name "*.cfg~" -type f -delete
    N_FILES=(./$STACK/*)
    N_FILES=${#N_FILES[@]}
done
rm -r $TEMP
