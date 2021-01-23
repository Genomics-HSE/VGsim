import sys

def ReadRates(fn):
    with open(fn) as f:
        line = next(f).rstrip()
        line = line.split(" ")

        line = next(f).rstrip()
        line = line.split(" ")
        print(line)
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

        print(bRate)
        print(dRate)
        print(sRate)
        print(mRate)

        return([bRate, dRate, sRate, mRate])
