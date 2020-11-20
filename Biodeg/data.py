import re
import numpy as np 



def get_features(data : str, feature: str):
    ids_features = re.findall( r'interpretation\((.*),{}\((.*)\)\)'.format(feature), data)
    return ids_features

def molecule_feature(data : str, molecule : str, feature: str):
    return float(re.findall( r'interpretation\({},{}\((.*)\)\)'.format(molecule,feature), data)[0])
def molecule_features (data : str, molecule_id: str ):
    return [molecule_feature(data, molecule_id, feature) for feature in FEATURES]
def molecule_atoms(data: str, molecule_id : str):
    return  re.findall( r'interpretation\({},atm\((.*),(.*)\)\)'.format(molecule_id), data)
def get_relation(data : str, molecule_id : str, relation :str):
    return re.findall( r'interpretation\({},{}\((.*),(.*)\)\)'.format(molecule_id, relation), data)

class Molecule(object):
    def __init__(self,data,  molecule_id):
        self.features = FEATURES
        self.id = molecule_id
        self.data = data
        self.atoms_names = molecule_atoms(data,molecule_id)
        self.atoms = [atom for atom, _ in self.atoms_names ]
        self.x_g = molecule_features(self.data, self.id)
    
    @property
    def x_a(self):
        if hasattr(self, '_x_a'):
            return self._x_a 
        fgroup = get_relation(self.data, self.id, 'fgroup')
        mapping = {id : func for id, func in fgroup}
        fgmember = get_relation(self.data, self.id, "fgmember")
        members = {atm :mapping[id] for id, atm in fgmember  }
        self.x_a = [[name,members.get(atom, None)] for atom, name in self.atoms_names ]
        return self.x_a
    @x_a.setter
    def x_a(self, x_a):
        self._x_a = x_a
    @property
    def edge_index(self):
        if hasattr(self, '_edge_index'):
            return self._edge_index 
        col = []
        row = []
        attr = []
        bonds = self.get_bonds()
        for i,j,val in bonds :
            col.append(self.atoms.index(i))
            row.append(self.atoms.index(j))
            attr.append(val)
        self.edge_index = col, row, attr
        return  self.edge_index
    @edge_index.setter
    def edge_index(self, edge_index):
        self._edge_index = edge_index
    def get_bonds(self):
        return re.findall( r'interpretation\({},bond\((.*),(.*),(.*)\)\)'.format(self.id), self.data)

    def __repr__(self):
        return 'Molecule(atoms = [{}], edge_index = [ 2,{}], id = {})'.format(len(self.atoms), 
        len(self.edge_index[0]),self.id)
if __name__ == "__main__":

    with open("biodeg_ext.txt", 'r') as f: dataset = f.readlines()
    with open("biodeg_ext.txt", 'r') as f: data = f.read()
    molecules = list(set(re.findall( r'interpretation\((.*),logp', data)))
    FEATURES = ["mweight", "activity", "logPnorm", "activitynorm", "logp", "mweightnorm"]
    relations = ["fgroup", "fgmember", "bond"]
    
    graphs = []
    for molecule_id in molecules:
        graphs.append(Molecule(data, molecule_id))
   
    print(graphs[5])
    print(graphs[100])
   

