import sys
import os


from FcfrpPDB import getStructure


def main():
    id = sys.argv[1]
    directory = sys.argv[2]
    
    structure = getStructure(id,modelChoose=None,path=directory)
    seqRes = structure.get_SeqRes()    
    #Show all residues and its chain
    s = ""
    for l in seqRes.getAllRes():
        s = s + l.getChainId() + " " + l.getIdl3() + "\n"
    print s
    #For each chain, show its residue
    for chainId in seqRes.getAllChain():
        s = ""
        NumRes = 0
        Num = seqRes.getNumResChain(chainId)
        for l in seqRes.getResChain(chainId):
            s = s + l.getIdl3() + " "
            NumRes += 1
        print chainId
        print s
        print "Residue Number read was %i and the number residue in getNumRes is %i. If they are match, all right. Otherwise, can be error occurred " % (NumRes, Num)

main()