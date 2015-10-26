set multiplot layout 96,1
unset key
unset tics
unset border
set bmargin 0

do for[i=0:95]{
	p 'out.dat' every :::i::i w l
}
