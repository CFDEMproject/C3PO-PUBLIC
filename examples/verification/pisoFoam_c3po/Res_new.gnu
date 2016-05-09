reset

set term x11 0 noraise



set title "   "
set ylabel 'Residuen'
set xlabel 'Iterationen'

set logscale y
#set xrange [0:28500]
set yrange [0.000001:1]

plot "< cat log | grep 'Solving for Uz' | cut -d' ' -f9 | tr -d ','" title 'Res Uz_initial' with lines,\
"< cat log | grep 'Solving for Ux' | cut -d' ' -f9 | tr -d ','" every 1 title 'Res Ux_initial' with lines,\
"< cat log | grep 'Solving for Uy' | cut -d' ' -f9 | tr -d ','" every 1 title 'Res Uy_initial' with lines,\
"< cat log | grep 'Solving for k' | cut -d' ' -f9 | tr -d ','" every 1 title 'Res k_initial' with lines,\
"< cat log | grep 'Solving for epsilon' | cut -d' ' -f9 | tr -d ','" every 1 title 'Res epsilon_initial' with lines,\
"< cat log | grep 'Solving for T' | cut -d' ' -f9 | tr -d ','" every 1 title 'Res T_initial' with lines,\
"< cat log | grep 'Solving for p_rgh' | cut -d' ' -f9 | tr -d ','" every 3 title 'Res p_rgh_initial' with lines,\
"< cat log | grep 'Solving for p' | cut -d' ' -f9 | tr -d ','" every 2 title 'Res p_initial' with lines
#"< cat log | grep 'Solving for Ux' | cut -d' ' -f13 | tr -d ','" every 2 title 'Res Ux_final' with lines,\
#"< cat log.reactingParcelFoam | grep 'Solving for Uy' | cut -d' ' -f13 | tr -d ','" every 1 title 'Res Uy_final' with lines,\
#"< cat log.reactingParcelFoam | grep 'Solving for Uz' | cut -d' ' -f13 | tr -d ','" title 'Res Uz_final' with lines,\
#"< cat log.reactingParcelFoam | grep 'Solving for p' | cut -d' ' -f13 | tr -d ','" every 1 title 'Res p_final' with lines






#"< cat log.SRFSimpleFoam | grep 'time step continuity errors' | cut -d' ' -f15 | tr -d ','" title 'sum global' with lines,\
#"< cat log.SRFSimpleFoam | grep 'Solving for Uy' | cut -d' ' -f9 | tr -d ','" title 'Uy' with lines,\
#"< cat log.SRFSimpleFoam | grep 'Solving for Uz' | cut -d' ' -f9 | tr -d ','" title 'Uz' with lines,\
#"< cat log.SRFSimpleFoam | grep 'Solving for epsilon' | cut -d' ' -f9 | tr -d ','" title 'epsilon' with lines,\
#"< cat log.SRFSimpleFoam | grep 'time step continuity errors' | cut -d' ' -f15 | tr -d ','" title 'sum global' with lines\
#"< cat log.SRFSimpleFoam | grep 'Solving for k' | cut -d' ' -f9 | tr -d ','" title 'k' with lines,\



pause 5
reread
