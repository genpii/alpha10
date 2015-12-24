set multiplot layout 96,1
unset key
unset tics
unset border
set bmargin 0

do for[i=0:95]{
	p 'rf55-54.dat' every :::i::i w l linewidth 0.1
}
