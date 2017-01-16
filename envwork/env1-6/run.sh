#dataset
dat=0
#algorithm
al=4_0_1_2_3
#trj=test2
#fl=2_0_1
for trj in working
do
for fl in 1_0 1_1 2_0_1
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
