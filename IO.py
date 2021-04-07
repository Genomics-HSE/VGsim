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
