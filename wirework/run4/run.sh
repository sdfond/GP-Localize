#dataset
dat=1
#algorithm
al=1_3
#trj=test2
#fl=2_0_1
for trj in work
do
for fl in 1_0 1_1 1_2 1_3 1_4 1_5 6_0_1_2_3_4_5
do
s=$(wc -l ../traj$dat/$trj"_loc" | awk '{print $1}')
step=$(expr $s - 1)
key=D$dat"_T"$trj"_F"$fl"_A"$al
echo $key
echo $step
sh genCase.sh $dat $trj $fl $al $step $key > $key.cfg
nohup ../gpbf $key.cfg &
done
done
