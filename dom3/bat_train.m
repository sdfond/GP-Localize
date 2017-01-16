for i=0:8
in=['input_',num2str(i)];
out=['output_',num2str(i)];
hypf=['hyper_',num2str(i)];
p.x=load(in);
p.y=load(out);
hyp=getHypARD(p,600);
fp=fopen(hypf,'w');
fprintf(fp,'%.4f %.4f %.4f %.4f %.4f', exp(hyp.lik), exp(hyp.cov'), hyp.mean);
fclose(fp);
end