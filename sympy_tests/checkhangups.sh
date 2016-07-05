#!/usr/bin/bash
#
ps -AF | grep -i "tex\|test\|make\|symp\|pyth\|pypy\|maxima\|dvipng" | grep -i -v "vim\|emacs\|make clean\|guake\|grep\|sbcl\|firewall\-applet\|evince";
printf "\n";
#
