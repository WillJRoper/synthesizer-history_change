
SYNTHESIZER_DATA_DIR="/research/astrodata/highz/synthesizer/"
GRID_NAME="agn_cloudy_T_alpha_Z_U_nH" 
# GRID_NAME="agn_feltre16_alpha_Z_U_nH" 

python create_synthesizer_grid.py -grid_name $GRID_NAME -synthesizer_data_dir $SYNTHESIZER_DATA_DIR

 