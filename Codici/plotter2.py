#!/usr/bin/env python3
# ——— monkey-patch per KDTree di pynbody ———
# import diretto del modulo che definisce KDTree
import pynbody.sph.kdtree as _kdtree

# 1) patch __init__ per ricast mass a pos.dtype prima di chiamare kdmain.init()
_orig_KDTree_init = _kdtree.KDTree.__init__
def _patched_KDTree_init(self, pos, mass, *args, **kwargs):
    if pos.dtype != mass.dtype:
        mass = mass.astype(pos.dtype)
    _orig_KDTree_init(self, pos, mass, *args, **kwargs)
_kdtree.KDTree.__init__ = _patched_KDTree_init

# 2) patch set_array_ref per riallineare mass al dtype di pos
_orig_set_array_ref = _kdtree.KDTree.set_array_ref
def _patched_set_array_ref(self, name, ar):
    if hasattr(self, '_pos') and ar.dtype != self._pos.dtype:
        ar = ar.astype(self._pos.dtype)
    return _orig_set_array_ref(self, name, ar)
_kdtree.KDTree.set_array_ref = _patched_set_array_ref
# ————————————————————————————————

import numpy as np
import pynbody
import argparse
from matplotlib import pyplot as plt
from pynbody.sph import kdtree


def gsoft(z, e=0.740):
    """
    Compute the softening length e (TNG100 default) as in mordor2.py.
    """
    return e if z <= 1 else (2*e)/(1+z)

def get_hmr(stars):
    """
    Compute half-mass radius of star particles.
    """
    prof = pynbody.analysis.profile.Profile(stars, ndim=3, type='log')
    return np.min(prof['rbins'][prof['mass_enc'] > 0.5 * prof['mass_enc'][-1]])

def center_and_rotate(gal):
    """
    Centre galaxy and align angular momentum exactly as in mordor2.py,
    ma con fallback su vel=False se mancano troppe particelle.
    Returns half-mass radius (hmr) in kpc.
    """

    # 1) centre on entire galaxy, con fallback per velocità
    try:
        pynbody.analysis.halo.center(gal, wrap=True, mode='hyb')
    except ValueError as e:
        if "Insufficient particles around center to get velocity" in str(e):
            pynbody.analysis.halo.center(gal, wrap=True, mode='hyb', vel=False)
        else:
            raise

    # 2) compute half-mass radius
    hmr = get_hmr(gal.s)

    # 3) check center-of-mass offset, con lo stesso fallback
    try:
        sc = pynbody.analysis.halo.center(gal.s, retcen=True, mode='hyb')
    except ValueError as e:
        if "Insufficient particles around center to get velocity" in str(e):
            sc = pynbody.analysis.halo.center(gal.s, retcen=True, mode='hyb', vel=False)
        else:
            raise

    if np.sqrt(np.sum(sc*sc)) > max(0.5 * hmr, 2.8 * gsoft(gal.properties['z'])):
        try:
            pynbody.analysis.halo.center(gal.s, wrap=True, mode='hyb')
        except:
            pynbody.analysis.halo.center(gal.s, wrap=True, mode='hyb', cen_size='3 kpc')
        hmr = get_hmr(gal.s)

    # 4) align face-on
    size_ang = max(3 * hmr, 2.8 * gsoft(gal.properties['z']))
    pynbody.analysis.angmom.faceon(gal.s, disk_size='%g kpc' % size_ang, already_centered=True)

    return hmr

def calc_k_rot(gal, hmr):
    """
    Calculate k_rot = sum(m * v_tan^2) / sum(m * v_tot^2) for particles within 3*hmr.
    """
    # positions in kpc
    x = gal.s['x'].in_units('kpc')
    y = gal.s['y'].in_units('kpc')
    z = gal.s['z'].in_units('kpc')
    # velocities in km/s
    vel = gal.s['vel'].in_units('km s**-1')
    vx = vel[:, 0]
    vy = vel[:, 1]
    vz = vel[:, 2]
    # spherical radius
    r = np.sqrt(x**2 + y**2)
    # tangential unit vectors
    phi_x = -y / r
    phi_y = x / r
    # tangential velocity
    v_tan = vx * phi_x + vy * phi_y
    # total velocity squared
    v_tot2 = vx**2 + vy**2 + vz**2
    # spherical mask within 3*hmr
    dist = np.sqrt(x**2 + y**2 + z**2)
    mask = dist <=  hmr / 4.5         # Mantenere 3*hmr
    m = gal.s['mass'][mask]
    # compute ratio
    k_rot = np.sum(m * v_tan[mask]**2) / np.sum(m * v_tot2[mask])
    return float(k_rot)

def main(index, v):
    # galaxy filename
    name = f'gal_{index:06d}.hdf5'
    print(f"Loading {name}...")
    gal = pynbody.load(name)
    gal.physical_units()
    # define softening length identically to mordor2.py
    gal['eps'] = pynbody.array.SimArray(
        2.8 * gsoft(gal.properties['z']) * np.ones_like(gal['x']), 'kpc'
    )
    # centre and align
    hmr = center_and_rotate(gal)
    print(f"Half-mass radius: {hmr:.3f} kpc")
    # compute k_rot
    krot = calc_k_rot(gal, hmr)
    title_k = f"$k_{{\\rm rot}} = {krot:.4f}$"
    print(f"k_rot: {krot:.4f}")
    # density map face-on
    width = v[2] * hmr
    pynbody.plot.sph.image(
        gal.s, 'rho', units='Msol pc^-2', width=f'{width:g} kpc',
        log=True, cmap='jet', vmin=v[0], vmax=v[1]
    )
    plt.suptitle(title_k, y=0.95)
    fname = f'gal_{index:06d}_faceon.png'
    plt.savefig(fname, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Salvata mappa face-on: {fname}")
    plt.close()
    # density map edge-on
    gal.rotate_x(90)
    pynbody.plot.sph.image(
        gal.s, 'rho', units='Msol pc^-2', width=f'{width:g} kpc',
        log=True, cmap='jet', vmin=v[0], vmax=v[1]
    )
    plt.suptitle(title_k, y=0.95)
    fname = f'gal_{index:06d}_edgeon.png'
    plt.savefig(fname, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Salvata mappa edge-on: {fname}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate density maps and compute k_rot for galaxy index'
    )
    parser.add_argument('index', type=int, help='Numeric index of galaxy file')
    parser.add_argument(
        '--v', nargs=3, type=float, default=[1e2, 1e4, 4],
        metavar=('VMIN', 'VMAX', 'FACTOR'),
        help='vmin vmax factor for image width = factor * hmr'
    )
    args = parser.parse_args()
    main(args.index, args.v)

