for i=1:6
    d = strcat('size100-', num2str(i))
    cd(d)
    pitc = load('D2_Treal_F1_1_A2_2_3_OPITC_err');
    sum(pitc)/length(pitc)
    cd ..
end