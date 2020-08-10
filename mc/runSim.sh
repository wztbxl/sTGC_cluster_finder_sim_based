
n_Evts=100
i=0

while [ $i -le $n_Evts ]
do
    echo "./output/Evts_10_$i.root"
    echo "./log/Evts_10_$i.log"
    ./sim ./output/Evts_10_$i.root ./log/Evts_10_$i.log
    echo "Run $i ended !"
    echo " "
    let i=i+1
done


