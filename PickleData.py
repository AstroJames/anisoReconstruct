from header import *


def ReadFiles():
    """
    This will read in all of the files from my MHD folder.

    INPUTS:
    ----------

    OUTPUTS:
    ----------

    """

    files    = {}
    densities   = ["M2MA0.1","M2MA0.5", "M2MA1", "M2MA10", "M2MA2",
                    "M4MA0.1","M4MA0.5", "M4MA1", "M4MA10", "M4MA2",
                    "M10MA0.1","M10MA0.5", "M10MA1", "M10MA10", "M10MA2",
                    "M20MA0.1","M20MA0.5", "M20MA1", "M20MA10", "M20MA2"]

    iterCount = 0

    # Read in all files
    for subDir in densities:
        os.system("ls /Volumes/JamesBe/MHD/{}/Proj > tmpFiles.txt".format(subDir))
        f = open("tmpFiles.txt")
        files[densities[iterCount]] = f.readlines()
        f.close()
        iterCount +=1

    return files

def LoadPickles(files,iter):
    """
    Load all pickles into a dictionary for optimal data handling.

    INPUTS:
    ----------

    OUTPUTS:
    ----------
    """

    dens    = {}
    dNames  = files.keys()

    for dir in dNames:
        dens[dir] = load_obj("/Volumes/JamesBe/MHD/" + dir + "/Proj/" + files[dir][iter].split("\n")[0] )

    return dens

def save_obj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch "+ name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
