import sys
from BirthDeathCython import Population, Lockdown

def ReadRates(fn):
    with open(fn) as f:
        line = next(f).rstrip()#header with version etc
        line = line.split(" ")

        line = next(f).rstrip()
        line = line.split(" ")
        hapFilled = False
        if line[0] == "H":
            hapFilled = True
        shift = int(hapFilled)
        dim = len(line) - shift
        hapNum = int( 4**(dim - 3) )
        bRate = []
        dRate = []
        sRate = []
        mRate = []
        if dim < 3:
            print("At least three rates (B, D, S) are expected")
            sys.exit(1)
        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            line = [el for el in line[shift:]]
            bRate.append(float(line[0]))
            dRate.append(float(line[1]))
            sRate.append(float(line[2]))
            mRate.append( [] )
            mutations = line[3:]
            for mut in mutations:
                a = mut.split(',')
                if len(a) == 1:
                    mRate[len(bRate)-1].append( [float(a[0]), 1.0/3.0, 1.0/3.0, 1.0/3.0] )
                elif len(a) == 4:
                    mRate[len(bRate)-1].append( [float(a[0]), float(a[1]), float(a[2]), float(a[3])] )
                else:
                    print("Error in mutations!!!")
                    sys.exit(1)
        return([bRate, dRate, sRate, mRate])

def ReadSusceptibility(fn):
    with open(fn) as f:
        line = next(f).rstrip()#header with version etc
        line = line.split(" ")

        line = next(f).rstrip()
        line = line.split(" ")
        hapFilled = False
        if line[0] == "H":
            hapFilled = True
        shift = int(hapFilled)

        susceptibility = []
        sType = []

        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            line = [float(el) for el in line[shift:]]
            susceptibility.append( line[1:] )
            sType.append( int( line[0] ) )
        return([susceptibility, sType])

def ReadPopulations(fn):
    with open(fn) as f:
        line = next(f).rstrip()#header with version etc
        line = line.split(" ")

        line = next(f).rstrip()
        line = line.split(" ")
        populations = []
        lockdown = []
        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            populations.append( Population(int(line[1]), float(line[2])) )
            if len(line) == 6:
                lockdown.append( Lockdown(float(line[3]), float(line[4]), float(line[5])) )
        return(populations, lockdown)

def ReadMigrationRates(fn):
    with open(fn) as f:
        line = next(f).rstrip()#header with version etc
        line = line.split(" ")
        migrationRates = []
        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            migrationRates.append( [float(v) for v in line] )
        for i in range(len(migrationRates)):
            migrationRates[i][i] = 0.0
        return(migrationRates)

def ReadSusceptibilityTransition(fn):
    with open(fn) as f:
        line = next(f).rstrip()#header with version etc
        line = line.split(" ")
        suscepTransition = []
        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            suscepTransition.append( [float(v) for v in line] )
        for i in range(len(suscepTransition)):
            suscepTransition[i][i] = 0.0
        return(suscepTransition)

def writeMutations(mut):
    #digits replacement
    alleles = ["A","T","C","G"]
    for i in [1,3]:
        for j in range(len(mut[i])):
            mut[i][j] = alleles[mut[i][j]]

    mutations_dict = {}
    for nodeId in mut[0]:
        if nodeId in mutations_dict: #adding mutation for existing node
            mutations_dict[nodeId] += str(mut[1][mut[0].index(nodeId)]) \
                                      + str(mut[2][mut[0].index(nodeId)]) \
                                      + str(mut[3][mut[0].index(nodeId)])+','
        else:
            mutations_dict[nodeId] = str(mut[1][mut[0].index(nodeId)]) \
                                     + str(mut[2][mut[0].index(nodeId)]) \
                                     + str(mut[3][mut[0].index(nodeId)])+','
    #removing extra comma
    for nodeId in mutations_dict:
        mutations_dict[nodeId] = mutations_dict[nodeId][:-1]

    f_mut = open('mutation_output.tsv', 'w')
    for i in range(len(pruferSeq)):
        if i in mutations_dict:
            f_mut.write(str(i)+'\t'+str(mutations_dict[i])+'\n')
        else:
            f_mut.write(str(i)+'\n')
    f_mut.close()

# count of childrens
def frequentCart(nodes, sequence):
    result = dict()
    #print("nodes=", nodes)
    for node in nodes:
        result[node] = 0
    for parent in sequence:
        if parent in nodes:
            result[parent] = result[parent] + 1
    return result

# place of childrens
def allChildrens(nodes, sequence):
    result = dict()
    for node in nodes:
        result[node] = []
    for index in range(len(sequence)):
        if sequence[index] in nodes:
            result[sequence[index]].append(index)
    return result

def getOutputDict(nodes, times):
    result = dict()
    for node in nodes:
        result[node] = "{0}:{1}".format(node, times[node])
    return result

def phase3_LookForParents(resultOutput, listOfLeefs, pruferSeq, allChildren):
    alreadyFinishedParent = []
    #print('phase 3')
    parentFutureLeeves = []
    futureLeeves = []
    alreadyFinishedLeeves = []
    for leef in listOfLeefs:
        if leef in alreadyFinishedLeeves:
            continue
        parent = int(pruferSeq[leef])
        # root is found
        if parent == -1:
            resultOutput[leef] = "(" + resultOutput[leef] + ")"
            continue
        if parent in alreadyFinishedParent:
            continue
        alreadyFinishedParent.append(parent)
        listOfChildren = allChildren[parent]
        isAllChildrenLeefs = True
        for child in listOfChildren:
            if not child in listOfLeefs:
                isAllChildrenLeefs = False
                break
        if isAllChildrenLeefs:
            parentFutureLeeves.append(parent)
            message = ""
            for child in listOfChildren:
                childSplit = resultOutput[child].split(':')
                absTimeChild = float(childSplit[-1])
                absTimeParent = float(resultOutput[parent].split(':')[-1])
                time = absTimeChild - absTimeParent
                childSplit[-1] = str(time)
                resultOutput[child] = ":".join(childSplit)

                message += resultOutput[child] + "," ##?????
                alreadyFinishedLeeves.append(child)
                resultOutput.pop(child)
            resultOutput[parent] = "(" + message[:-1] + ")" + resultOutput[parent]#!!!!!
        else:
            futureLeeves.append(leef)
    return parentFutureLeeves, futureLeeves


def writeGenomeNewick(pruferSeq, times):
    #pruferSeq = pruferSeq.astype(int)
    for i in range(len(pruferSeq)):
        if pruferSeq[i] == i:
            pruferSeq[i] = -1
    #number of nodes
    numberOfNodes = len(pruferSeq)
    listOfNodes = [i for i in range(numberOfNodes)]
    frequencyCart = frequentCart(listOfNodes, pruferSeq)
    allChildren = allChildrens(listOfNodes, pruferSeq)
    resultOutput = getOutputDict(listOfNodes, times)

    #phase 2: look for normal leefs and parents
    listOfLeefs = []
    for key in frequencyCart:
        if frequencyCart[key] == 0:
            # find leefs
            resultOutput[key] = "{0}:{1}".format(key, times[key])
            listOfNodes.remove(key)
            listOfLeefs.append(key)

    #phase 3: look for parents
    parentFutureLeeves, futureLeeves = phase3_LookForParents(resultOutput, listOfLeefs, pruferSeq, allChildren)

    #phase 4: union lists
    while(True):
        listsOfNextLeeves = parentFutureLeeves + futureLeeves
        noParents = True
        actualLeafList = []
        for leef in listsOfNextLeeves:
            if pruferSeq[int(leef)] != -1.0:
                actualLeafList.append(leef)
                noParents = False
        if noParents:
             break
        else:
            parentFutureLeeves, futureLeeves = phase3_LookForParents(resultOutput, actualLeafList, pruferSeq, allChildren)

    f_nwk = open('newick_output.nwk', 'w')
    for key in resultOutput:
        f_nwk.write(resultOutput[key])
    f_nwk.write(';')
    f_nwk.close()
    #print(len(times))
