

MACHINE="apollo"
SYNTHESIZER_DATA_DIR="/research/astrodata/highz/synthesizer/"
CLOUDY=$CLOUDY17

while read GRID_NAME; do
  echo "$GRID_NAME"
  python create_cloudy_input_grid.py -machine $MACHINE -synthesizer_data_dir $SYNTHESIZER_DATA_DIR -grid_name $GRID_NAME -cloudy $CLOUDY
done < grids.txt



