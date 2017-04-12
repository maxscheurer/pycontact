from PyContact.exampleData.datafiles import DCD, PSF
from PyContact.cy_modules import cy_sasa_cuda, cy_gridsearch
import MDAnalysis
from PyContact.core.Biochemistry import vdwRadius
import numpy as np
from PyContact.core.LogPool import *
from PyContact.core.multi_accumulation import chunks

def calculate_sasa_parallel(input_coords, natoms, pairdist, nprad,
                            surfacePoints, probeRadius, pointstyle,
                            restricted, restrictedList, rank):
    """Computes the SASA in parallel."""

    temp_sasa = []
    frames_processed = 0
    # sasaProgressDict[rank] = frames_processed
    # print(len(input_coords))
    for c in input_coords:
        coords = np.reshape(c, (1, natoms * 3))
        npcoords = np.array(coords, dtype=np.float32)
        # print("start C")
        # startC = time.time()
        asa = cy_gridsearch.cy_sasa(npcoords, natoms, pairdist, 0, -1, nprad, surfacePoints, probeRadius,
                                    pointstyle, restricted, restrictedList)
        # stopC = time.time()
        # print("time for grid search: ", (stopC - startC))
        # print("asa:", asa)
        temp_sasa.append(asa)
        frames_processed += 1
        # sasaProgressDict[rank] = frames_processed
    return temp_sasa

def calculate_sasa_cuda(input_coords, natoms, pairdist, nprad,
                        surfacePoints, probeRadius, pointstyle,
                        restricted, restrictedList):

    temp_sasa = []
    # print(len(input_coords))
    coords = np.reshape(input_coords, (1, natoms * 3))
    npcoords = np.array(coords, dtype=np.float32)
    npcoords = npcoords[0]
    # print("start C")
    # startC = time.time()
    asa = cy_sasa_cuda.cy_sasa_cuda(npcoords, natoms, pairdist, nprad, surfacePoints, probeRadius,
                                    pointstyle, restricted, restrictedList)
    # stopC = time.time()
    # print("time for grid search: ", (stopC - startC))
    # print("asa:", asa)
    return asa

psf = PSF
dcd = DCD

allSasas = []

# load psf and trajectory, make lists with radii and coordinates
if psf == "" or dcd == "":
    pass

try:
    u = MDAnalysis.Universe(psf, dcd)
except IOError:
    pass

probeRadius = 1.4

# seltext = "segid UBQ"
# resseltext = "segid UBQ and same residue as around 5.0 (segid RN11)"
# seltext = sasaSelection1TextField.text()
# seltext2 = sasaSelection2TextField.text()
# resseltext = sasaRestrictionTextField.text()
seltext = "segid RN11"
resseltext = ""

# 0=spiral, 1=random (VMD)
pointstyle = 0
# number of points to approximate the sphere
surfacePoints = 100
# pair distance
pairdist = 2.0 * (2.0 + 1.4)

if resseltext != "":
    restricted = 1
else:
    restricted = 0

selection = u.select_atoms(seltext)

radius = []
restrictedList = []
if restricted:
    ressel = u.select_atoms(resseltext)
    for s in selection.atoms:
        if s in ressel.atoms:
            restrictedList.append(1)
        else:
            restrictedList.append(0)
        radius.append(vdwRadius(s.name[0]))
else:
    restrictedList = [0]
    for s in selection.atoms:
        radius.append(vdwRadius(s.name[0]))
natoms = len(selection)
nprad = np.array(radius, dtype=np.float32)
restrictedList = np.array(restrictedList, dtype=np.int32)

# TODO: bug if selection is not static for all frames
# TODO: dynamic allocation of positions in every frame!
input_coords = []
for ts in u.trajectory:
    input_coords.append(selection.positions)

results = []
trajLength = len(u.trajectory)
totalFramesToProcess = trajLength
allSasasTrue = []
print(radius, nprad)
for frame in range(totalFramesToProcess):
    allSasas.append(calculate_sasa_cuda(input_coords[frame], natoms, pairdist, nprad,
                                       surfacePoints, probeRadius, pointstyle,
                                       restricted, restrictedList))

# print("SASA Results: ", allSasas)

input_coords = []
for ts in u.trajectory:
    # ressel = u.select_atoms(resseltext)
    # print("restricted: ", len(ressel.atoms))
    input_coords.append(selection.positions)

nprocs = 2
input_chunks = chunks(input_coords, nprocs)
pool = LoggingPool(nprocs)
results = []
rank = 0
trajLength = len(u.trajectory)
totalFramesToProcess = trajLength
for input_coords_chunk in input_chunks:
    results.append(pool.apply_async(calculate_sasa_parallel, args=(input_coords_chunk, natoms, pairdist, nprad,
                                                                   surfacePoints, probeRadius, pointstyle,
                                                                   restricted, restrictedList, rank)))
    rank += 1
print("ranks", rank)
pool.close()
pool.join()

for r in results:
    allSasasTrue.extend(r.get())
print("SASA Results: ", allSasasTrue)
