r=1e-08
m=1e-08
alpha=1.9
fragment=10000000
for nrep in 0 1 2 3 4 5 6 7 8 9; do
	echo "arg-sample --sites ./sites/beta_sim_"$alpha"_r"$r"m"$m"_"$nrep".sites --region 1-"$fragment"  --mutrate "$m" --recombrate "$r" --popsize 1000000 --output ./weaved_args/beta_"$alpha"_arg-sample_r"$r"m"$m"_1-"$fragment"_"$nrep" --iters 1000 --ntimes 40"
done

r=1e-08
m=1e-08
alpha=1.7
fragment=10000000
for nrep in 0 1 2 3 4 5 6 7 8 9; do
	echo "arg-sample --sites ./sites/beta_sim_"$alpha"_r"$r"m"$m"_"$nrep".sites --region 1-"$fragment"  --mutrate "$m" --recombrate "$r" --popsize 1000000 --output ./weaved_args/beta_"$alpha"_arg-sample_r"$r"m"$m"_1-"$fragment"_"$nrep" --iters 1000 --ntimes 40"
done

r=1e-07
m=1e-07
alpha=1.5
fragment=10000000
for nrep in 0 1 2 3 4 5 6 7 8 9; do
	echo "arg-sample --sites ./sites/beta_sim_"$alpha"_r"$r"m"$m"_"$nrep".sites --region 1-"$fragment"  --mutrate "$m" --recombrate "$r" --popsize 1000000 --output ./weaved_args/beta_"$alpha"_arg-sample_r"$r"m"$m"_1-"$fragment"_"$nrep" --iters 1000 --ntimes 40" 
done

r=1e-06
m=1e-06
alpha=1.3
fragment=10000000
for nrep in 0 1 2 3 4 5 6 7 8 9; do
	echo "arg-sample --sites ./sites/beta_sim_"$alpha"_r"$r"m"$m"_"$nrep".sites --region 1-"$fragment"  --mutrate "$m" --recombrate "$r" --popsize 1000000 --output ./weaved_args/beta_"$alpha"_arg-sample_r"$r"m"$m"_1-"$fragment"_"$nrep" --iters 1000 --ntimes 40" 
done


