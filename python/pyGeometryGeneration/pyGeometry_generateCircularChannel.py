import lib_pyGeometryGeneration as libGeo
import numpy as np


def generateRandomCircularChannel(length,radius,nParticles):

	#http://mathworld.wolfram.com/SpherePointPicking.html

	z = np.random.uniform(-length/2.0,length/2.0,nParticles)
	phi = np.random.uniform(0,2.0*np.pi,nParticles)


	x = radius*np.sin(phi)
	y = radius*np.cos(phi)
	coords = np.transpose(np.vstack([x,y,z]))
	return(coords)



r1= 0.004
#nParticles = 800000
#nParticles =  30000000
nParticles =  3000000

coords = generateRandomCircularChannel(0.04,0.002,nParticles)




#print(coords)
libGeo.writeParticleJsonFile(coords.tolist(),"/Users/dystopic/Library/Developer/Xcode/DerivedData/Simion-IonSimulation-cpp-fjdovxyhvcghupdpgpccybtmeybg/Build/Products/Debug/testCircularChannel_fullParticleNumber.json")