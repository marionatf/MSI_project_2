"""given a set of molecular dynamics trajectories, computes the volumetric map of 
the chemical potential for water."""

from htmd import *
from htmd.molecule.util import maxDistance
from htmd.molecule.util import writeVoxels
import time

def prepare_traject(s,ref_traj=None):
    traj = (Molecule(s.molfile))
    traj.read(s.trajectory)
    traj.wrap('protein')
    if ref_traj is None:
        traj.align('protein')
        ref_traj=traj.copy()
    else:
        traj.align('protein',ref_traj)
    coo = traj.get('coords')
    c = np.mean(coo, axis=0)
    traj.moveBy(-c)
    traj_wat= traj.copy()
    traj_wat.filter("water and name OH2")    
    return(traj,traj_wat,ref_traj)

def define_extremes(traj_wat):
    D = maxDistance(traj_wat)
    min_d = 0
    max_d = 2*int(D+1)    
    return(min_d,max_d)

def create_count_dict(min_d,max_d):
    w_count={}
    for x in range(min_d,max_d+1):
        for y in range(min_d,max_d+1):
            for z in range(min_d,max_d+1):
                w_count[(x,y,z)]=0
    return(w_count)

def calculate_occupancy(traj_wat,w_count,min_d,max_d):
    print("Calculating water occupancy....")
    for atom_c in traj_wat.coords:
        for i in range(traj_wat.numFrames):
            x=int(atom_c[0,i])
            y=int(atom_c[1,i])
            z=int(atom_c[2,i])
            if min_d < x < max_d and min_d < y < max_d and min_d < z < max_d:
                w_count[(x,y,z)]+=1
                w_count[(x-1,y,z)]+=0.5
                w_count[(x+1,y,z)]+=0.5
                w_count[(x,y-1,z)]+=0.5
                w_count[(x,y+1,z)]+=0.5
                w_count[(x,y,z-1)]+=0.5
                w_count[(x,y,z+1)]+=0.5   
    w_count_list= np.array([w_count[key] for key in sorted(w_count)]).reshape(max_d+1,max_d+1,max_d+1)   
    avg_wat = w_count_list / (traj_wat.numFrames*traj_wat.numAtoms)   
    print("Done.")    
    return(avg_wat)

def chemical_pot(occupancy):
    kB=0.001987191 
    T=298 
    pot = -kB*T*np.log(occupancy)        
    return(pot)

def view_vmd(traj, max_d):    
    traj.moveBy([max_d/2,max_d/2,max_d/2])
    htmd.config(viewer='vmd')
    traj.reps.remove()   
    traj.reps.add(sel='protein', style='NewCartoon', color='2')
    traj.view()
    time.sleep(6000)


def volumetric_map():
    sims = simlist(glob('data/*/'), glob('data/*/structure.pdb'), glob('input/*/'))
    sim_num=len(sims)
    for s in sims:
        print("\nSimulation", s.simid +1, "of", sim_num)
        if s.simid == 0:
            (traj,traj_wat,ref_traj)=prepare_traject(s)
            (min_d,max_d)=define_extremes(traj_wat)
        else:
            (traj,traj_wat,ref_traj)=prepare_traject(s,ref_traj)
        traj_wat.moveBy([max_d/2,max_d/2,max_d/2])
        w_count=create_count_dict(min_d,max_d)
        avg_wat = calculate_occupancy(traj_wat,w_count,min_d,max_d)
        try:
            occup_all += avg_wat
        except NameError:
            occup_all = avg_wat
    occup_all=occup_all/sim_num
    occup_all += 1e-40
    pot=chemical_pot(occup_all)
    writeVoxels(pot,"cubefile.cube",np.array([min_d,min_d,min_d]),np.array([max_d,max_d,max_d]),np.array([1,1,1]))
    print("\nVolumetric file created.\n")
    view_vmd(traj,max_d)

if __name__ == "__main__":
    volumetric_map()