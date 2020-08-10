# sTGC_cluster_finder_sim_based

finder2_v1.cxx is the version I used before 2020/07/21 STAR FWD f2f meeting. I try divide one module into 8 pieces and combine 1D hits to 2D hits.
but this verision did not work now.
finder2_v2.cxx is the version I try handle high multiplicity events
finder2_v3.cxx is the version I try with Zhangbu's method

root -l -b -q 'finder2_v1.cxx("../stgc-cluster-sim/output/Evts_10_'$i'.root"," ./out/Evts_10/Cluster_output_Evts10_'$i'_v1.root")' > ./log/Cluster_output_Evts10_${i}_v1.log
this is how to run the finder2_v1.cxx

root -l -b -q 'finder2.cxx("simulation_file_name","output name)'

1. run the simulation code 
    ./runSim.sh
2. reconstruct cluster
    ./runRC.sh
3. check the reconstruct result
    root -l -q -b CheckReconstruct.cxx


