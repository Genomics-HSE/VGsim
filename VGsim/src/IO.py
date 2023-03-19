import sys
import math

def read_rates(fn):
    hapFilled = False
    flag = False
    bRate = []
    dRate = []
    sRate = []
    mRate = []
    with open(fn) as f:
        line = next(f).rstrip().split(" ")#header with version etc
        line = next(f).rstrip().split(" ")
        if line[0] == "H":
            hapFilled = True
        shift = int(hapFilled)
        if line[2+shift] == "SP":
            flag = True
        samProbability = int(flag)
        dim = len(line) - shift
        hapNum = int( 4**(dim - 3) )
        if dim < 3:
            print("At least three rates (B, D, S) are expected")
            sys.exit(1)

        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip().split(" ")[shift:]
            bRate.append(float(line[0]))
            if flag == 0:
                dRate.append(float(line[1]))
                sRate.append(float(line[2]))
            else:
                dRate.append(float(line[1]) * (1 - float(line[2])))
                sRate.append(float(line[1]) * float(line[2]))

            mRate.append( [] )
            mutations = line[3:]
            for mut in mutations:
                a = mut.split(',')
                if len(a) == 1:
                    mRate[-1].append( [float(a[0]), 1.0/3.0, 1.0/3.0, 1.0/3.0] )
                elif len(a) == 4:
                    mRate[-1].append( [float(a[0]), float(a[1]), float(a[2]), float(a[3])] )
                else:
                    print("Error in mutations!!!")
                    sys.exit(1)
    mRate = update_mRate(mRate)
    return bRate, dRate, sRate, mRate

def update_mRate(mRate):
    if math.log(len(mRate), 4) != int(math.log(len(mRate), 4)):
        print("Error!")
        sys.exit(1)

    for i in range(len(mRate)):
        for j in range(len(mRate[0])):
            mRate[i][j].insert(calculate_allele(i, j, len(mRate[0]))+1, 0)
    return mRate

def calculate_allele(haplotype, site, sites):
    for _ in range(sites-site):
        allele = haplotype % 4
        haplotype = haplotype // 4
    return allele

def read_susceptibility(fn):
    hapFilled = False
    susceptibility = []
    sType = []
    with open(fn) as f:
        line = next(f).rstrip().split(" ")#header with version etc
        line = next(f).rstrip().split(" ")
        if line[0] == "H":
            hapFilled = True
        shift = int(hapFilled)

        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip().split(" ")[shift:]
            susceptibility.append( line[1:] )
            sType.append( int( line[0] ) )
        return susceptibility, sType

def read_populations(fn):
    sizes = []
    contactDensity = []
    contactAfter = []
    startLD = []
    endLD = []
    samplingMultiplier = []
    with open(fn) as f:
        line = next(f).rstrip().split(" ")#header with version etc
        line = next(f).rstrip().split(" ")

        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip().split(" ")
            sizes.append(int(line[1]))
            contactDensity.append(float(line[2]))
            if len(line) == 4:
                part_line = line[3].split(",")
                if len(part_line) == 1:
                    contactAfter.append(0)
                    startLD.append(1.0)
                    endLD.append(1.0)
                    samplingMultiplier.append(float(part_line[0]))
                elif len(part_line) == 3:
                    contactAfter.append(float(part_line[0]))
                    startLD.append(float(part_line[1]))
                    endLD.append(float(part_line[2]))
                    samplingMultiplier.append(1)
            elif len(line) == 5:
                part_line1 = line[3].split(",")
                part_line2 = line[4].split(",")
                if len(part_line1) == 1:
                    samplingMultiplier.append(float(part_line1[0]))
                    contactAfter.append(float(part_line2[0]))
                    startLD.append(float(part_line2[1]))
                    endLD.append(float(part_line2[2]))
                elif len(part_line1) == 3:
                    samplingMultiplier.append(float(part_line2[0]))
                    contactAfter.append(float(part_line1[0]))
                    startLD.append(float(part_line1[1]))
                    endLD.append(float(part_line1[2]))
    return sizes, contactDensity, contactAfter, startLD, endLD, samplingMultiplier

def read_matrix(fn):
    with open(fn) as f:
        line = next(f).rstrip()#header with version etc
        line = line.split(" ")
        matrix = []
        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            matrix.append( [float(v) for v in line] )
    return matrix

def writeMutations(mut, len_prufer, name_file, file_path):
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

    if file_path != None:
        f_mut = open(file_path + '/' + name_file + ".tsv", 'w')
    else:
        f_mut = open(name_file + ".tsv", 'w')

    for i in range(len_prufer):
        if i in mutations_dict:
            f_mut.write(str(i)+'\t'+str(mutations_dict[i])+'\n')
    f_mut.close()

class Vertex():
    def __init__(self, root, root_time, children, populations):
        self.__children = children
        self.__root = root
        self.__root_time = root_time
        self.__root_population_id = populations[self.__root_time]
        left_node = self.__children[root][0][0]
        left_time = self.__children[root][0][1]
        right_node = self.__children[root][1][0]
        right_time = self.__children[root][1][1]

        if left_node in self.__children:
            self.__left_child = Vertex(left_node, left_time, self.__children, populations)
        else:
            self.__left_child = Leaf(left_node, left_time, populations)
            
        if right_node in self.__children:
            self.__right_child = Vertex(right_node, right_time, self.__children, populations)
        else:
            self.__right_child = Leaf(right_node, right_time, populations)

    def get_children(self, base_time):
        return '({0},{1}){2}:{3}'.format(self.__left_child.get_children(self.__root_time), self.__right_child.get_children(self.__root_time), self.__root, self.__root_time - base_time)

    def write_population(self):
        return '{0}\t{1}\n'.format(self.__root, self.__root_population_id) + self.__left_child.write_population() + self.__right_child.write_population()

class Leaf(Vertex):
    def __init__(self, leaf, times, populations):
        self.__leaf = leaf
        self.__times = times
        self.__leaf_population_id = populations[times]

    def get_children(self, base_time):
        return '{0}:{1}'.format(self.__leaf, self.__times - base_time)

    def write_population(self):
        return '{0}\t{1}\n'.format(self.__leaf, self.__leaf_population_id)

#find list with childrens
def find_children(pruferSeq, times):
    children = {}
    for index in range(len(pruferSeq)):
        add_list = []
        add_list = [index, times[index]]
        if pruferSeq[index] in children:
            children[pruferSeq[index]].append(add_list)
        else:
            children[pruferSeq[index]] = []
            children[pruferSeq[index]].append(add_list)
    return children

def get_last(output_string):
    try:
        return output_string[-1]
    except:
        return "notDigit"

def writeGenomeNewick(pruferSeq, times, populations, name_file, file_path):
    children = find_children(pruferSeq, times)
    root = children[-1][0][0]
    root_time = children[-1][0][1]

    result = Vertex(root, root_time, children, populations)

    if file_path != None:
        f_nwk = open(file_path + '/' + name_file + '_tree.nwk', 'w')
        f_pop = open(file_path + '/' + name_file + '_sample_population.tsv', 'w')
    elif name_file != None:
        f_nwk = open(name_file + '_tree.nwk', 'w')
        f_pop = open(name_file + '_sample_population.tsv', 'w')
    else:
        f_nwk = open('tree.nwk', 'w')
        f_pop = open('sample_population.tsv', 'w')

    f_nwk.write(result.get_children(root_time))
    f_nwk.write(';')
    f_pop.write(result.write_population())

    f_nwk.close()
    f_pop.close()
