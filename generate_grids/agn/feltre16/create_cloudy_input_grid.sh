

MACHINE="apollo"
SYNTHESIZER_DATA_DIR="/research/astrodata/highz/synthesizer/"
GRID_NAME="feltre16"
CLOUDY_PARAMS="default.yaml"
CLOUDY=$CLOUDY17

# python create_cloudy_input_grid.py -machine $MACHINE -synthesizer_data_dir $SYNTHESIZER_DATA_DIR -grid_name $GRID_NAME -cloudy $CLOUDY -cloudy_params=$CLOUDY_PARAMS -dry_run True 
python create_cloudy_input_grid.py -machine $MACHINE -synthesizer_data_dir $SYNTHESIZER_DATA_DIR -grid_name $GRID_NAME -cloudy $CLOUDY -cloudy_params=$CLOUDY_PARAMS