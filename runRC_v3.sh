
n_Evts=1000
i=0

while [ $i -le $n_Evts ]
do
    echo "./output/Evts_10_$i.root"
    echo "./log/Evts_10_$i.log"
    root -l -b -q 'finder2_v3.cxx("../stgc-cluster-sim/output/Evts_10_'$i'.root"," ./out/Evts_10_v3/Cluster_output_Evts10_'$i'_v1.root")' > ./log/v3/Cluster_output_Evts10_${i}_v3.log 
    echo "Run $i ended !"
    echo " "
    mkdir ./plots/10Evt_v3/$i/
    mv ./plots/*.png ./plots/10Evt_v3/$i/
    let i=i+1
done


