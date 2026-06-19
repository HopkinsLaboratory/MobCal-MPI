import numpy as np

from atom_info import chemical_symbols, atomic_masses_legacy, covalent_radii

def find_connectivity(labels,xyz):
    '''Determine connectivity and return com of main group.'''
    covalent = {}
    for index, label in enumerate(chemical_symbols):
        covalent[label] = {'Radii':covalent_radii[index],'Mass':atomic_masses_legacy[index]}

    masses = [covalent[x]['Mass'] for x in labels] #get the masses of each atom

    bond_lengths = {}
    maximum = 0.0
    for label1 in sorted(labels):
        for label2 in sorted(labels):
            if ''.join(sorted(label1+label2)) in bond_lengths.keys(): continue
            bl = round((covalent[label1]['Radii']+covalent[label2]['Radii'])*1.1,2)
            bond_lengths[label1+label2] = bl
            if bl > maximum: maximum = bl

    def Distance(Group,mcom):
        '''Find the distance between a group of coordinates and a center.'''
        dist = (Group-mcom)**2
        dist = dist.sum(axis=-1)
        dist = np.sqrt(dist)
        return dist

    connectivity = {}
    for index, atom in enumerate(xyz):
        Dists = Distance(xyz[index+1:],atom)
        PB = np.argwhere(Dists<maximum) #Get distances shorter than maximum
        con = []
        for di in PB: #For distance index in passing bonds
            key = ''.join(sorted(labels[index]+labels[index+di[0]+1]))
            bl = bond_lengths[key]
            if Dists[di] < bl: con.append(index+di[0]+2) #add one because of skip and another to match gaussian labels
        if len(con) > 0: connectivity[index+1] = con

    molecules = []

    def search(connectivity,atom,molecule):
        for x in connectivity[atom]:
            if x in connectivity.keys(): search(connectivity,x,molecule)
            if x not in molecule: molecule.append(x)
        return molecule

    molecule = []
    for x in connectivity.keys():
        if x in molecule: continue
        if x not in molecule: molecule = []
        molecule = search(connectivity,x,molecule)
        molecule.append(x)
        molecules.append(sorted(molecule))
    #Find the core molecule
    if len(molecules) > 1: core_index = np.argmax([len(x) for x in molecules])
    else: core_index = 0
    core = molecules[core_index]
    core_xyz = np.array([xyz[x-1] for x in core])
    core_masses = np.array([masses[x-1] for x in core])
    core_weighted = core_xyz*np.reshape(core_masses,(len(core_masses),1))
    mcom = np.sum(core_weighted,axis=0)/sum(core_masses)
    #Subtract the com of the core from xyz
    xyz = xyz-mcom
    if len(molecules) > 1: 
        xyz = np.array([xyz[x] for x in range(len(labels)) if x+1 not in core])
        masses = np.array([masses[x] for x in range(len(labels)) if x+1 not in core])
    #Calculate distance of atom to com
    distances = np.array([Distance(x,mcom) for x in xyz])
    #Calculate mass weighted distances
    mw = sorted(list(distances*masses))
    return mw
