
SYNTHESIZER_DATA_DIR="/research/astrodata/highz/synthesizer/"
while read GRID_NAME; do
  echo "$GRID_NAME"
  python create_synthesizer_grid.py -grid_name $GRID_NAME -synthesizer_data_dir $SYNTHESIZER_DATA_DIR
done < grids.txt




 