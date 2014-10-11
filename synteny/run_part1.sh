#!/bin/sh

# run the first part of the pipeline for all four species

./050_sliding_window.py barley
./050_sliding_window.py brachy
./050_sliding_window.py rice
./050_sliding_window.py sorghum


./055_sliding_window_oat.py barley
./055_sliding_window_oat.py brachy
./055_sliding_window_oat.py rice
./055_sliding_window_oat.py sorghum


./060_make_synteny_plot_data_filtered.py barley
./060_make_synteny_plot_data_filtered.py brachy
./060_make_synteny_plot_data_filtered.py rice
./060_make_synteny_plot_data_filtered.py sorghum
