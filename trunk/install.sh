#!/bin/bash

#DIALOG=gdialog
DIALOG=dialog

# Graphical install to appear soon
#$DIALOG --title "Step 1" --backtitle "LIA_RAL Installation" --yesno "Welcome to LIARAL installation, please proceed to continue" 10 50

perl -pe 'my $dir=`pwd`;s/(LIA_RAL_TOPDIR=)/\1$dir/g' config.txt.ex > config.txt

