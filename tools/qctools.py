# Tools for molecular modelling with Jupyter Notebooks
# 2021-2026
# antti.karttunen@aalto.fi

############## Functions to facilitate printing ##############

def print_info(info):
    # info: string
    banner = "-----------------------------------------------------------"
    print(f"{banner}\n{info}\n{banner}")
          
def print_error(error):
    # error: string
    err = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"      
    print(f"{err}\n{error}\n{err}")

############## Functions for loading molecules ##############

def load_xyz(xyzfile, silent = False):
    """Loads XYZ file
    xyzfile: file name or path as string
    silent: no printing
    
    Returns: ase.Atoms or None if fails
    """
    from ase.io import read    
    try:
        atoms = read(xyzfile, format = 'xyz')
        if not silent:
            print_info(f"Molecule was loaded from file {xyzfile}\n"
                       f"Atoms: {len(atoms)}\n"
                       f"Formula: {str(atoms.symbols)}")
    except OSError:
        if not silent:
            print_error(f"Failed to load file {xyzfile}")
        return None
    else:
        atoms.info['id'] = xyzfile
        return atoms

def load_xyz_as_traj(xyztraj, silent = False):
    """ Converts multi-XYZ file to ASE trajectory (to visualize with nglview)
    xyztraj: file name or path as a string (multi-XYZ file)    
    silent: no printing

    Returns: ase.io.Trajectory or None if fails
    """
    from ase.io import read, write, Trajectory

    # Read the XYZ file (index=':' reads all frames/images)
    try:
        frames = read(xyztraj, index=':')
        if not silent:
            print_info(f"XYZ trajectory was loaded from file {xyztraj}\n"
                       f"Frames: {len(frames)}")
    except OSError:
        if not silent:
            print_error(f"Failed to load file {xyzfile}")
        return None        
    
    # Write the frames to a temporary .traj file and load as Trajectory
    try:
        trajfile = xyztraj + '.traj'
        with Trajectory(trajfile, mode='w') as traj:
            for frame in frames:
                traj.write(frame)            
    except OSError:
        if not silent:
            print_error(f"Failed to save file {trajfile}")
        return None
        
    try:
        traj = Trajectory(trajfile)
    except OSError:
        if not silent:
            print_error(f"Failed to load trajectory file {trajfile}")
    else:
        return traj

def save_xyz(atoms, xyzfile, silent = False):
    """Saves XYZ file
    atoms: ase.Atoms
    xyzfile: file name or path as string
    silent: no printing
    
    no return value
    """
    from ase.io import write    
    try:
        write(xyzfile, images = atoms, format = 'xyz')
        if not silent:
            print_info(f"Molecule was saved to file {xyzfile}\n")
    except OSError:
        if not silent:
            print_error(f"Failed to save file {xyzfile}")

def load_molecule_pubchem(name=None, cid=None, xyzfile=None):
    # Returns: ase.Atoms, None if fails
    # if xyzfile is not None, creates also an xyz-file
    
    from ase.data.pubchem import pubchem_atoms_search
    # Catch UserWarning about conformers
    import warnings    

    if name is not None and cid is not None:        
        print_error("Only give name or cid, not both.")
        return None
    elif name is None and cid is None:
        print_error("Give name or cid.")
        return None
    elif name is not None:
        usename = True
        id = name
    else:
        usename = False
        id = cid
        
    atoms = None
    print_info(f"Loading molecule {id} from PubChem...")
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            if usename:
                atoms = pubchem_atoms_search(name=name)                
            else:
                atoms = pubchem_atoms_search(cid=cid)  
            if len(w) > 0:
                print("ATTENTION: Molecule has several conformers, PubChem returned the first of them.")
    except:
        print_error("Loading the molecule failed, check the name or cid from https://pubchem.ncbi.nlm.nih.gov/")       
    
    if atoms is not None:
        print_info(f"Molecule {id} was loaded from PubChem.\n"
                   f"Atoms: {len(atoms)}\n"
                   f"Formula: {str(atoms.symbols)}")
        atoms.info['id'] = str(id)
        if xyzfile is not None:
            from ase.io import write
            write(filename = xyzfile, images = atoms, format = 'xyz')
            print_info(f"Molecular structure was saved in file {xyzfile}")
          
    return atoms
    
############## Functions for visualizing molecules with nglview

def show_molecule(molecule, size = (500, 400), style = 'ball+stick', unitcell = None, labels = None, bg = 'black', gui = False):
    """ Shows a molecule using NGLWidget
    molecule: ase.Atoms or ase.io.Trajectory
    size: tuple of two integers (x_pixels, y_pixels)
    style: 'ball+stick', 'spacefill', 'licorice', 'line' (http://nglviewer.org/ngl/api/manual/molecular-representations.html)
    unitcell: None or 'white', 'orange', 'red', ...
    labels: None or 'atomname', 'atomindex, 'element'
    bg: 'black', 'white', ...
    gui: set to True to show GUI
    Returns: NGLWidget
    """
    import nglview
    from ase import Atoms
    from ase.io.trajectory import TrajectoryReader

    if isinstance(molecule, Atoms):
        nv = nglview.show_ase(molecule)
    elif isinstance(molecule, TrajectoryReader):
        nv = nglview.show_asetraj(molecule)
    else:
        print_error("Invalid molecule!")
        return None
    
    nv._set_size(f"{size[0]}px", f"{size[1]}px")
    # Calls: nv._remote_call('setSize', target='Widget', args=[w, h])
    
    nv.clear_representations()    
    nv.add_representation(style)
    
    if style == 'spacefill':
        # radiusType = 'covalent' looks better for bulk materials. Can be changed to 'vdw' via GUI
        nv.update_representation(component = 0, repr_index = 0, radiusType='covalent')

    if unitcell is not None:
        nv.add_representation('unitcell')        
        # Changing unit cell color does not work, not even via NGL GUI?
        # nv.add_representation('unitcell', colorValue = '#ffffff')
    if labels is not None:
        nv.add_representation('label', labelType = labels)
        
    nv.parameters = dict(backgroundColor = bg, clipDist = -100)
    nv.camera = 'orthographic'
    nv.display(gui = gui)
    return nv  

