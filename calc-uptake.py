#!/usr/bin/env python

import re
import argparse
import collections

import numpy as np
from numpy import linalg
import periodic


def get_init_symbs(lines):

    """Gets the symbols of the initial set of atoms

    :param lines: The lines in the output file

    """

    init_coord_beg_sentinel = re.compile(
        'Fractional coordinates of asymmetric unit'
        )
    init_coord_end_sentinel = re.compile(
        'Molecule list generated from bond lengths'
        )
    coord_re = re.compile(r'\d+\s+(\w+)\s+[cs]\s+\d+\.\d+\s+\d+\.\d+')

    beg_sentinel_idx = next(i for i, v in enumerate(lines) 
                            if init_coord_beg_sentinel.search(v) != None)
    end_sentinel_idx = next(i for i, v in enumerate(lines) 
                            if init_coord_end_sentinel.search(v) != None)

    coord_search = [
        coord_re.search(l) for l in lines[beg_sentinel_idx:end_sentinel_idx]
        ]
    symbols = [i.group(1) for i in coord_search if i != None]

    return symbols

def get_volume(lines):

    """Gets the volume of the unit cell

    :param lines: The lines in the output file

    """

    begin_sentinel = re.compile(r'^ *Cartesian lattice vectors \(Angstroms\)')
    end_sentinel = re.compile(r'^ *Cell parameters')
    lattice_re = re.compile(r'(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')

    begin_idx = next(i for i, v in enumerate(lines) 
		     if begin_sentinel.search(v) != None)
    end_idx = next(i for i, v in enumerate(lines)
		   if end_sentinel.search(v) != None)

    lattice_search = [lattice_re.search(v) for v in lines[begin_idx:end_idx]]
    lattices = [ i.group(1, 2, 3) for i in lattice_search if i != None]

    lattice_matrix = np.array(lattices, dtype=np.float64)

    return linalg.det(lattice_matrix)

def get_atm_weight(symb):

    """Gets the weight of a given element symbol

    :param symb: The element symbol, which is able to be suffixed by some
        numerals.

    """

    elem_symb = re.search(r'([a-zA-Z]+)\d*', symb).group(1)
    return periodic.element(elem_symb).mass

def get_subs_weight(symbs, ads_symb):

    """Gets the weight of the substrate

    :param symbs: A list of symbols for the system
    :param ads_symb: The element symbol for the adsorbant
    :returns: The weight of the substrate and the number of atoms

    """

    subs = [i for i in symbs if i != ads_symb]
    cnt = collections.Counter(subs)

    wgt = 0.0
    for i, v in cnt.iteritems():
	wgt += get_atm_weight(i) * v
	continue

    return wgt, len(subs)

def get_mean_n_atms(lines):
    
    """Gets the mean number of atoms

    :param lines: The lines in the output file

    """

    trials_re = re.compile('Trials:.*')
    match_res = [trials_re.search(i) for i in lines]
    last_trial = [i for i in match_res if i != None][-1]

    return float(
	last_trial.group(0).split()[-1]
	)

def main():

    """The main function"""

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--ads', default='H99', metavar='SYMBOL',
			help='The element symbol for the adsorbant')
    parser.add_argument('FILE', help='The GULP output file to process')
    args = parser.parse_args()
    ads_symb = args.ads
    file_name = args.FILE

    output_file = open(file_name)
    lines = [i for i in output_file]
    output_file.close()

    symbs = get_init_symbs(lines)
    vol = get_volume(lines)
    subs_wgt, subs_n_atm = get_subs_weight(symbs, ads_symb)
    ads_wgt = (get_mean_n_atms(lines) - subs_n_atm) * get_atm_weight(ads_symb)

    wgt_uptake = ads_wgt / subs_wgt * 100.0
    vol_uptake = (ads_wgt * 1.672622e-24) / (vol * 1e-27) # gram per L

    print "Weight uptake: %f %%; Volume uptake: %f g/L" % (wgt_uptake,
							   vol_uptake)

    return 0

if __name__ == '__main__':
    main()



















