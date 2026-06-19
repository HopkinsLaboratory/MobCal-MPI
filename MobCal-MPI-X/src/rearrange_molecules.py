import numpy as np

import atom_info as AI

class GEOM:
    def __init__(self,geom):
        self.interperate(geom) #Read in molecule data
        self.comCentre() #Shift Molecule to Center
        self.getCOMDist() #Get distances to center of mass
        self.geom[:,1:] 
        self.index = None
        self.sort()

    def interperate(self,geom):
        atoms = geom.strip().split("\n")
        self.natoms = len(atoms)
        self.geom = np.zeros((len(atoms),4))
        for i,a in enumerate(atoms):
            a = a.strip().split()
            self.geom[i][0] = AI.atomic_masses_legacy[AI.chemical_symbols.index(a[0])]
            self.geom[i][1] = float(a[1])
            self.geom[i][2] = float(a[2])
            self.geom[i][3] = float(a[3])
        self.masses = self.geom[:,0]
        self.totalMass = np.sum(self.masses)

    def comCentre(self):
        CentreOfMass = np.zeros(3)
        for i in range(self.natoms):
            m = self.masses[i]
            if m >0.001:
                #geom[i][0] is mass and [i][1] is pc
                CentreOfMass[0] += self.masses[i]*self.geom[i][1] #x
                CentreOfMass[1] += self.masses[i]*self.geom[i][2] #y
                CentreOfMass[2] += self.masses[i]*self.geom[i][3] #z
        self.com = CentreOfMass/self.totalMass
        self.geom[:,1:] -= self.com

    def getCOMDist(self):
        self.comDist = np.zeros(self.natoms)
        for i in range(self.natoms):
            self.comDist[i] = R12(self.geom[i],self.com,com=True)

    def sort(self):         
        idx = np.argsort(self.comDist) 
        try: self.index = np.asarray(self.index)[idx]
        except IndexError: pass
        self.geom = self.geom[idx] #order atoms 
        self.order = np.array([x for x in range(self.natoms)])[idx] #order indexes

def R12(u,v,com=False,cc=False):
    if com:
        r = u[1:]-v
    elif cc:
        r = u-v
    else:
        r = u[1:]-v[1:]
    return np.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])         

if __name__ == "__main__":
    mol1 = '''
    C    -0.495540    -0.251825    0.285975
    H    -1.453131    0.250860    0.527272
    H    -0.399707    -1.172167    0.888992
    O    0.545824    0.595961    0.799789
    H    0.557857    1.458577    0.300488
    H    2.734917    0.036967    1.261451
    C    -0.291422    -0.460758    -1.197307
    H    -1.141858    -0.993909    -1.650570
    H    -0.194417    0.493961    -1.735613
    H    0.609998    -1.042799    -1.427561
    C    5.566822    1.233954    0.512592
    H    5.469949    1.495923    1.586236
    H    5.899599    2.127046    -0.046593
    O    4.220213    0.945751    0.174671
    H    1.955988    0.012905    0.710940
    C    6.463570    0.030969    0.286293
    H    6.613890    -0.185405    -0.777989
    H    7.462491    0.199064    0.718289
    H    6.062990    -0.874940    0.755885
    C    1.299373    3.923798    0.180446
    H    1.091219    3.420206    1.144452
    H    0.621017    4.789605    0.063472
    O    0.951198    2.929747    -0.783515
    H    0.699135    3.325730    -1.639296
    C    2.760930    4.286858    0.011052
    H    2.962561    4.723812    -0.975516
    H    3.077331    5.027101    0.761154
    H    3.411037    3.409028    0.128383
    C    3.914429    0.653743    -1.187202
    H    4.231874    -0.389189    -1.405117
    H    2.523150    -0.708213    0.911929
    O    2.511515    0.585159    -1.149307
    H    2.077541    1.478907    -1.140851
    C    4.395598    1.707129    -2.155714
    H    4.147966    2.723801    -1.822838
    H    5.483173    1.666199    -2.308632
    H    3.930440    1.566446    -3.145254
    '''

    mol2 = '''
    O    -0.036511    0.035060    -0.041379
    H    0.923870    -0.002623    0.015572
    H    -0.344311    -0.553895    0.641256
    O    2.871893    -0.079049    0.131112
    H    3.300511    -0.617725    -0.534012
    H    3.326479    0.762073    0.109662
    '''

    print(GEOM(mol1).geom)
    print(GEOM(mol2).geom)
