#!/bin/bash
# synthesizer_dir="/Users/stephenwilkins/Dropbox/Research/data/synthesizer/"
synthesizer_dir="/research/astrodata/highz/synthesizer/" # apollo
machine="apollo"
# sps="bpass-2.2.1-bin_chabrier03-0.1,100.0 "
# sps="bpass-2.2.1-bin_bpl-0.1,1.0,100.0-1.3,2.0 bpass-2.2.1-bin_bpl-0.1,1.0,100.0-1.3,2.35 bpass-2.2.1-bin_bpl-0.1,1.0,100.0-1.3,2.7 bpass-2.2.1-bin_bpl-0.1,1.0,300.0-1.3,2.0 bpass-2.2.1-bin_bpl-0.1,1.0,300.0-1.3,2.35 bpass-2.2.1-bin_chabrier03-0.1,100.0 bpass-2.2.1-bin_chabrier03-0.1,300.0 bpass-2.2.1-sin_bpl-0.1,1.0,100.0-1.3,2.0 bpass-2.2.1-sin_bpl-0.1,1.0,100.0-1.3,2.35 bpass-2.2.1-sin_bpl-0.1,1.0,100.0-1.3,2.7 bpass-2.2.1-sin_bpl-0.1,1.0,300.0-1.3,2.0 bpass-2.2.1-sin_bpl-0.1,1.0,300.0-1.3,2.35 bpass-2.2.1-sin_bpl-0.1,1.0,300.0-1.3,2.7 bpass-2.2.1-sin_chabrier03-0.1,100.0 bpass-2.2.1-sin_chabrier03-0.1,300.0"
p="default_param.yaml"
c=$CLOUDY17

cd ..
python3 make_cloudy_input_grid.py -dir $synthesizer_dir -m $machine -sps $sps  -p $p  -c $c
