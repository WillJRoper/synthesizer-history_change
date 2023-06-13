

MACHINE="apollo"
SYNTHESIZER_DATA_DIR="/research/astrodata/highz/synthesizer/"
GRID_NAME="agn_cloudy_T_alpha_Z_U_nH" 
# GRID_NAME="agn_feltre16_alpha_Z_U_nH" 
CLOUDY=$CLOUDY17

python create_cloudy_input_grid.py -machine $MACHINE -synthesizer_data_dir $SYNTHESIZER_DATA_DIR -grid_name $GRID_NAME -cloudy $CLOUDY