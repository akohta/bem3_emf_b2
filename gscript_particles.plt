# Gnuplot script file
# Usage : gnuplot -e 'file="particles_filename"' gscript_particles.plt

if ( !exists("file") ) file="ex.particles"

set terminal postscript eps color enhanced "Arial" 20 size 6in,8in
set output "particles.eps"

#set palette model RGB rgbformulae 35,13,10
set palette model RGB rgbformulae 22,13,-31
set view equal xyz
set view 75,32,1,1
set ticslevel 0
unset colorbox
set grid x y z vertical
set pointsize 0.4
#set title "nodes for surface integral"

set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}" 
set zlabel "{/Arial-Italic z}"

splot for[i=0:*] file using 1:2:3:(($4==i)? $4 : 1/0) pt 7 palette title sprintf("object %d",i)
