import sys
import h5py
import numpy as np

from .galaxy import Galaxy


def load_FLARES(f, region, tag):
    with h5py.File(f, 'r') as hf:
        lens = hf[f'{region}/{tag}/Galaxy/S_Length'][:]
        ages = hf[f'{region}/{tag}/Particle/S_Age'][:]  # Gyr
        # coods = hf[f'{tag}/Particle/S_Coordinates'][:]
        mass = hf[f'{region}/{tag}/Particle/S_Mass'][:]  # 1e10 Msol
        # imass = hf[f'{tag}/Particle/S_MassInitial'][:]
        metals = hf[f'{region}/{tag}/Particle/S_Z'][:]
        # ids = hf[f'{tag}/Particle/S_ID'][:]
        # index = hf[f'{tag}/Particle/S_Index'][:]
        # hf[f'{pre}/S_Vel']
        # hf[f'{pre}/S_Z_smooth']

    
    # ages = np.log10(ages * 1e9)  # log10(yr)
    ages = (ages * 1e9)  # log10(yr)
    mass = mass * 1e10  # Msol

    begin, end = get_len(lens)

    galaxies = [None] * len(begin)
    for i, (b, e) in enumerate(zip(begin, end)):
        galaxies[i] = Galaxy()
        galaxies[i].load_stars(mass[b:e], ages[b:e], metals[b:e])

    return galaxies


def get_len(Length):
    begin = np.zeros(len(Length), dtype=np.int64)
    end = np.zeros(len(Length), dtype=np.int64)
    begin[1:] = np.cumsum(Length)[:-1]
    end = np.cumsum(Length)
    return begin, end


def main():
    # e.g.
    # '/cosma7/data/dp004/dc-love2/codes/flares/data/FLARES_30_sp_info.hdf5'
    _f = sys.argv[1]
    tag = sys.argv[2]  # '/010_z005p000/'
    load_FLARES(_f, tag)
    return None


if __name__ == "__main__":
    main()
