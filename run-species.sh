echo "cd `pwd`; ./run_sd.sh -t 128 -s all" | qsub -N reverb-species -l nodes=node02:ppn=128 -v WORKDIR=`pwd` -j oe -o ~/log/reverb-species.log
