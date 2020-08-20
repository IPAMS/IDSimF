import json


def writeParticleJsonFile(particleCoordinates,filename):
	import json
	data = {"number_of_particles":len(particleCoordinates),"particle_coordinates":particleCoordinates}
	with open(filename, 'w') as outfile:
		json.dump(data, outfile,indent=2)


