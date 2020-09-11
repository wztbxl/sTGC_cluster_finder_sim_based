
n_Evts=1000
i=0

while [ $i -le $n_Evts ]
do
    echo "./output/Evts_10_$i.root"
    echo "./log/Evts_10_$i.log"
    root -l -b -q 'finder2_v1.cxx("../stgc-cluster-sim/output/Evts_10_'$i'.root"," ./out/Evts_10/Cluster_output_Evts10_'$i'_v1.root")' > ./log/Cluster_output_Evts10_${i}_v1.log 
    echo "Run $i ended !"
    echo " "
    mkdir ./plots/10Evt/$i/
    mv ./plots/*.png ./plots/10Evt/$i/
    let i=i+1
done


