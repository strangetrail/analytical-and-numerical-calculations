#!/usr/bin/bash
#
#<filename>.<cas>.tex
LSARGC=$#;
for LSARGV in $*; do
	LSARGVFILE=${LSARGV%tex};
	#sed -i.bak s/\\$//g $LSARGVFILE'tex';
	latex $LSARGVFILE'tex';
	dvipng -D 175 $LSARGVFILE'dvi';
done;
#
