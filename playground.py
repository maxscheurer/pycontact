from __future__ import print_function
#draft & sketch playground

# prototype writer for vmd visualization
# f = open('showHBondsInVMD.tcl', 'w')
# f.write('mol new %s \n'%psf)
# f.write('mol addfile %s \n'%dcd)
# f.write('mol delrep 0 top \n')
# f.write('mol representation NewCartoon \n')
# f.write('mol Color ColorID 3 \n')
# f.write('mol selection {all} \n')
# f.write('mol addrep top \n')
# for frame in contactResults:
# 	for contact in frame:
# 		for hbond in contact.hbondinfo:
# 			# hbond.toString()
# 			f.write('mol representation VDW \n')
# 			f.write('mol Color Name \n')
# 			f.write('mol selection {index %d %d %d} \n'%(hbond.donorIndex, hbond.acceptorIndex, hbond.hydrogenIndex))
# 			f.write('mol addrep top \n')
# f.close()

#draft, don't delete
# counter = 0
# 	for frame in hbond_key_frame_accumulate[key]:
# 		for item in frame:
# 			for hbond in item:
# 				print counter
# 				hbond.toString()
# 		counter += 1

# show memory information
	    if 0:
	        memory = []
	        for var, obj in locals().items():
	            memory.append([var, sys.getsizeof(obj)])
	        sorted_by_second = sorted(memory, key=lambda data: data[1], reverse=True)
	        for item in sorted_by_second:
	            print(item[0], humansize(item[1]))
