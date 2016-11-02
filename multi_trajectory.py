def loop_trajectory(sel1c,sel2c,config,suppl):
    print len(sel1c) , len(sel2c)
    indices1 = suppl[0]
    indices2 = suppl[1]
    cutoff, hbondcutoff, hbondcutangle = config
    # resname_array = comm.bcast(resname_array, root=0)
    # resid_array = comm.bcast(resid_array, root=0)
    # name_array = comm.bcast(name_array, root=0)
    type_array = suppl[2]
    bonds = suppl[3]
    # segids = comm.bcast(segids, root=0)
    # backbone = comm.bcast(backbone, root=0)
    heavyatoms = suppl[4]

    allRankContacts = []
    start = time.time()
    for s1, s2 in zip(sel1c, sel2c):
        frame = 0
        currentFrameContacts = []
        result = np.ndarray(shape=(len(s1), len(s2)), dtype=float)
        distarray = distances.distance_array(s1, s2, box=None, result=result)
        contacts = np.where(distarray <= cutoff)
        for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
            convindex1 = indices1[idx1]  # idx1 converted to global atom indexing
            convindex2 = indices2[idx2]  # idx2 converted to global atom indexing
            # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered, hydrogen bonds can still be detected!
            # if re.match("H(.*)", self.name_array[convindex1]) or re.match("H(.*)", self.name_array[convindex2]):
            #     continue
            # distance between atom1 and atom2
            distance = distarray[idx1, idx2]
            weight = weight_function(distance)
            hydrogenBonds = []
            type1 = next((x.htype for x in heavyatoms if x.name == type_array[convindex1]), AtomHBondType.none)
            type2 = next((x.htype for x in heavyatoms if x.name == type_array[convindex2]), AtomHBondType.none)
            if type1 != AtomHBondType.none and type2 != AtomHBondType.none:
                if (type1 == AtomHBondType.both and type2 == AtomHBondType.both) or \
                        (type1 == AtomHBondType.acc and type2 == AtomHBondType.don) or \
                        (type1 == AtomHBondType.don and type2 == AtomHBondType.acc) or \
                        (type1 == AtomHBondType.both and type2 == AtomHBondType.acc) or \
                        (type1 == AtomHBondType.acc and type2 == AtomHBondType.both) or \
                        (type1 == AtomHBondType.don and type2 == AtomHBondType.both) or \
                        (type1 == AtomHBondType.both and type2 == AtomHBondType.don):
                    # print "hbond? %s - %s" % (type_array[convindex1], type_array[convindex2])
                    # search for hatom, check numbering in bond!!!!!!!!!!
                    b1 = bonds[convindex1]
                    b2 = bonds[convindex2]
                    # search for hydrogen atoms bound to atom 1
                    bondcount1 = 0
                    hydrogenAtomsBoundToAtom1 = []
                    for b in b1.types:
                        hydrogen = next((x for x in b if x.startswith("H")), 0)
                        # print b
                        if hydrogen != 0:
                            # print "h bond to atom1"
                            bondindices1 = b1.to_indices()[bondcount1]
                            hydrogenidx = next(
                                ((j + 1) for j in bondindices1 if type_array[j + 1].startswith("H")), -1)
                            if hydrogenidx != -1:
                                # print type_array[hydrogenidx]
                                hydrogenAtomsBoundToAtom1.append(hydrogenidx)
                        bondcount1 += 1
                    # search for hydrogen atoms bound to atom 2
                    bondcount2 = 0
                    hydrogenAtomsBoundToAtom2 = []
                    for b in b2.types:
                        hydrogen = next((x for x in b if x.startswith("H")), 0)
                        # print b
                        if hydrogen != 0:
                            # print "h bond to atom2"
                            bondindices2 = b2.to_indices()[bondcount2]
                            hydrogenidx = next(
                                ((k + 1) for k in bondindices2 if type_array[k + 1].startswith("H")), -1)
                            if hydrogenidx != -1:
                                # print type_array[hydrogenidx]
                                hydrogenAtomsBoundToAtom2.append(hydrogenidx)
                        bondcount2 += 1
                    # check hbond criteria for hydrogen atoms bound to first atom
                    for global_hatom in hydrogenAtomsBoundToAtom1:
                        conv_hatom = indices1.index(global_hatom)
                        typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[convindex2]),
                                         AtomHBondType.none)
                        if typeHeavy == AtomHBondType.acc and (distarray[conv_hatom, idx2] <= hbondcutoff):
                            donorPosition = s1[idx1]
                            hydrogenPosition = s1[conv_hatom]
                            acceptorPosition = s2[idx2]
                            v1 = hydrogenPosition - acceptorPosition
                            v2 = hydrogenPosition - donorPosition
                            v1norm = np.linalg.norm(v1)
                            v2norm = np.linalg.norm(v2)
                            dot = np.dot(v1, v2)
                            angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                            if angle >= hbondcutangle:
                                dist = distarray[conv_hatom, idx2]
                                new_hbond = HydrogenBond(convindex1, convindex2, global_hatom, dist, angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)
                            # print str(convindex1) + " " + str(convindex2)
                            # print "hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2)
                            # print angle
                    for global_hatom in hydrogenAtomsBoundToAtom2:
                        conv_hatom = indices2.index(global_hatom)
                        typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[convindex1]),
                                         AtomHBondType.none)
                        if typeHeavy == AtomHBondType.acc and (distarray[idx1, conv_hatom] <= hbondcutoff):
                            donorPosition = s2[idx2]
                            hydrogenPosition = s2[conv_hatom]
                            acceptorPosition = s1[idx1]
                            v1 = hydrogenPosition - acceptorPosition
                            v2 = hydrogenPosition - donorPosition
                            v1norm = np.linalg.norm(v1)
                            v2norm = np.linalg.norm(v2)
                            dot = np.dot(v1, v2)
                            angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                            if angle >= hbondcutangle:
                                dist = distarray[idx1, conv_hatom]
                                new_hbond = HydrogenBond(convindex2, convindex1, global_hatom, dist, angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)
            newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1), int(convindex2),
                                 hydrogenBonds)
            currentFrameContacts.append(newAtomContact)
        allRankContacts.append(currentFrameContacts)
    return allRankContacts