#!/bin/bash

#DIALOG=gdialog
DIALOG=dialog

##### Copier les données de config.txt.ex dans le config.txt
perl -pe 'my $dir=`pwd`;s/(LIA_RAL_TOPDIR=)/\1$dir/g' config.txt.ex > config.txt

#### créer la librairie
echo ----- création de la librairie LIA_SpkTools -------------
cd LIA_SpkTools
make
cd ..

### créer les outils de LIA RAL
echo ----- compilation des outils LIA_RAL --------------------
make


