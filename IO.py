import sys
from BirthDeath import Population

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
            line = [float(el) for el in line[shift:]]
            bRate.append(line[0])
            dRate.append(line[1])
            sRate.append(line[2])
            mRate.append( line[3:] )
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
        for line in f:
            if line[0] == "#":
                next
            line = line.rstrip()
            line = line.split(" ")
            populations.append( Population(int(line[1]), float(line[2])) )
        return(populations)

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
            migrationRates[i,i] = 0.0
        return(migrationRates)
