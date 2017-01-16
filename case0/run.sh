#dataset
dat=0
#algorithm
al=5_0_1_2_3_6
#trj=test2
#fl=2_0_1
for trj in test1 test2
do
for fl in 1_0 1_1 2_0_1
do
step=$(wc -l ../traj$dat/$trj"_loc" | awk '{print $1}')
key=D$dat"_T"$trj"_F"$fl"_A"$al
echo $key
echo $step
sh genCase.sh $dat $trj $fl $al $step $key > $key.cfg
nohup ./gpbf $key.cfg &
done
done
