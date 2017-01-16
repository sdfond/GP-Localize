#dataset
dat=2
#algorithm
al=4_0_1_2_3
#trj=test2
#fl=2_0_1
for trj in real
do
for fl in 1_0 1_1 1_2 2_0_1 3_0_1_2
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
