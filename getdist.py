import Bio
import Bio.PDB
import sys



def main():
    # initialize parser
    parser = Bio.PDB.PDBParser()

    # initialize phosphate arrays
    # each entry: [chain,number]
    HJPhosphates = []
    PXPhosphates = []

    # load 2pfj
    HJStructure = parser.get_structure("2pfj","2pfj.pdb")
    # load PX model
    PXStructure = parser.get_structure("PX","PX.pdb")

    # store distances as a dictionary where the key is a tuple of the two
    # phosphates and the value is the distance between them in angstroms
    HJDict = getPDistances(HJPhosphates,HJStructure)
    PXDict = getPDistances(PXPhosphates,PXStructure)
    HJOut = open("HJDistances.txt","w+")
    PXOut = open("PXDistances.txt","w+")
    writeFile(HJDict)
    writeFile(PXDict)
    writeFile(HJDict,HJOut)
    writeFile(PXDict,PXOut)
        
def getPDistances(phosphateList,structure):
    '''Given a list of phosphates in the form of [chain,res#],outputs a 
    dictionary of interphosphate distances where the key is a tuple of the two 
    phosphates and the value is the distance between them'''
    outputDict = {}
    for phos1 in phosphateList:
        for phos2 in phosphateList:
            [chain1,res1] = phos1
            [chain2,res2] = phos2
            p1 = structure[0][chain1][res1]["P"]
            p2 = structure[0][chain2][res2]["P"]
            dist = p1 - p2
            outputDict[(phos1,phos2)] = dist
    return outputDict

def writeFile(distDict,outStream=sys.stdout):
    '''Outputs distDict as a text file'''
    keys = sorted(distDict.keys())
    for (phos1,phos2) in keys:
        [chain1,res1] = phos1
        [chain2,res2] = phos2
        outStream.write("%s_%i-%s_%i: "%(chain1,res1,chain2,res2))
        outStream.write(distDict[(phos1,phos2)])
        outStream.write("\n")
    outStream.close()

if __name__ == "__main__": main()        