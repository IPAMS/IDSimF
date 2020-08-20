import lib_pyGeometryGeneration as libGeo
import numpy as np


def generateRandomSphereShell(scale,nParticles):

	#http://mathworld.wolfram.com/SpherePointPicking.html

	z = 2 * np.random.rand(nParticles) - 1   # uniform in -1, 1
	t = 2 * np.pi * np.random.rand(nParticles)   # uniform in 0, 2*pi
	x = np.sqrt(1 - z**2) * np.cos(t)
	y = np.sqrt(1 - z**2) * np.sin(t)
	coords = np.transpose(np.vstack([x,y,z]))*scale
	return(coords)



r1= 0.004
nParticles = 80000

coords = generateRandomSphereShell(r1,nParticles)




#print(coords)
libGeo.writeParticleJsonFile(coords.tolist(),"/Users/dystopic/Library/Developer/Xcode/DerivedData/Simion-IonSimulation-cpp-fjdovxyhvcghupdpgpccybtmeybg/Build/Products/Debug/testSphereUniformHighDensity.json")