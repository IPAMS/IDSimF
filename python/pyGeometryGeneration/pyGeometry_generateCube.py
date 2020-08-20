import lib_pyGeometryGeneration as libGeo
import numpy as np


def generateRandomParticleBlock(lowerCorner,upperCorner,nParticles):

	x = np.random.uniform(lowerCorner[0],upperCorner[0],nParticles)
	y = np.random.uniform(lowerCorner[1],upperCorner[1],nParticles)
	z = np.random.uniform(lowerCorner[2],upperCorner[2],nParticles)

	coords = np.transpose(np.vstack([x,y,z]))
	return(coords)



sideLength = 0.004
originDistance = 0.002
thickness = 0.0001
nParticles = 10000

coordsX1 = generateRandomParticleBlock(
			[-originDistance-thickness, -sideLength/2.0, -sideLength/2.0],
			[-originDistance, sideLength/2.0, sideLength/2.0],
			nParticles
)
coordsX2 = generateRandomParticleBlock(
			[originDistance, -sideLength/2.0, -sideLength/2.0],
			[originDistance+thickness, sideLength/2.0, sideLength/2.0],
			nParticles
)
coordsY1 = generateRandomParticleBlock(
			[ -sideLength/2.0, -originDistance-thickness, -sideLength/2.0],
			[ sideLength/2.0, -originDistance, sideLength/2.0],
			nParticles
)
coordsY2 = generateRandomParticleBlock(
			[ -sideLength/2.0, originDistance, -sideLength/2.0],
			[ sideLength/2.0, originDistance+thickness, sideLength/2.0],
			nParticles
)
coordsZ1 = generateRandomParticleBlock(
			[ -sideLength/2.0, -sideLength/2.0, -originDistance-thickness],
			[ sideLength/2.0, sideLength/2.0,-originDistance],
			nParticles
)
coordsZ2 = generateRandomParticleBlock(
			[ -sideLength/2.0, -sideLength/2.0, originDistance],
			[ sideLength/2.0, sideLength/2.0, originDistance+thickness],
			nParticles
)

coords = np.vstack([coordsX1,coordsX2,coordsY1,coordsY2,coordsZ1,coordsZ2])




#print(coords)
libGeo.writeParticleJsonFile(coords.tolist(),"/Users/dystopic/Library/Developer/Xcode/DerivedData/Simion-IonSimulation-cpp-fjdovxyhvcghupdpgpccybtmeybg/Build/Products/Debug/testCube.json")