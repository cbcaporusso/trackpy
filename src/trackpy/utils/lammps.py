# This module contains functions for parsing lammps output files.
# WARNING - It include some magic constants that are specific to the
# LAMMPS ouput files that I have been using. It may not work for
# other files.

def lammps_header_parser(file, column_name: str):
    with open(file, 'r') as file:
        for line in file:
            if line.startswith('ITEM: ATOMS'):
                headers = line.split()[2:]
                column_index = headers.index(column_name)
                return column_index

def lammps_boxsize_parser(file):
    with open(file, 'r') as file:
        for line in file:
            if line.startswith('ITEM: BOX BOUNDS'):
                line_x = next(file)
                line_y = next(file)

                lx = float(line_x.split()[1]) 
                ly = float(line_y.split()[1])

                return lx, ly