import Bio
import Bio.PDB
import sys



def main():
    # initialize parser
    parser = Bio.PDB.PDBParser()

    # initialize phosphate arrays
    # each entry: [chain,number]
    HJPhosphates = [["Z",8],["Y",8]]
    PXPhosphates = [["W",8],["W",9],["W",19],["W",20],["X",19],["X",18],["X",8]]
    PXPhosphates += [["X",7],["Z",8],["Z",9],["Z",19],["Z",20],["Y",8]]
    PXPhosphates += [["Y",9],["Y",19],["Y",20]]

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
    HJOut.close()
    PXOut.close()
        
def getPDistances(phosphateList,structure):
    '''Given a list of phosphates in the form of [chain,res#],outputs a 
    dictionary of interphosphate distances where the key is a tuple of the two 
    phosphates and the value is the distance between them'''
    outputDict = {}
    for phos1 in phosphateList:
        for phos2 in phosphateList:
            if phos2 == phos1:
                continue
            [chain1,res1] = phos1
            [chain2,res2] = phos2
            p1 = structure[0][chain1][res1]["P"]
            p2 = structure[0][chain2][res2]["P"]
            dist = p1 - p2
            key = "%s_%i-%s_%i"%(chain1,res1,chain2,res2)
            outputDict[key] = dist
    return outputDict

def writeFile(distDict,outStream=sys.stdout):
    '''Outputs distDict as a text file'''
    keys = sorted(distDict.keys())
    for key in keys:
        outStream.write("%s: "%(key))
        outStream.write(str(distDict[key]))
        outStream.write("\n")

if __name__ == "__main__": main()        