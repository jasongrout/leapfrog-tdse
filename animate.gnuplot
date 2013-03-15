reset
set term gif animate
set output "animate.gif"
set xrange[1:2]
set yrange[0:1]
n=50 # Number of frames
dn=5000/n # Assuming 5000 total points (1e4 time steps)
do for [i=0:n]{
   plot "./leapfrog_norm.out" u 1:dn*i w l
}
set output