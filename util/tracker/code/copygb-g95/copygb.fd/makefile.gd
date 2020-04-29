w3=../w3lib/w3lib.a

libs=../bacio/libbacio.a  ../iplib/iplib.a ../splib/splib.a  ${w3}


copygb:	copygb.f  ${libs}
	g95 -o copygb -fbounds-check copygb.f ${libs}

${libs}:
	cd $(*D) ; make
