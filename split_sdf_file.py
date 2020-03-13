import os,sys

def read_sdf_file(sdf,nCompound):
    counter = 0
    counterTotal = 0
    bunch = 0
    out = ''
    for line in open(sdf):
        out += line
        if line.startswith('$$$'):
            counter += 1
            counterTotal += 1
            sys.stdout.write("\r-> found %d ligands in SDF file" %counterTotal)
            sys.stdout.flush()
            if counter == nCompound:
                print ' -> writing file:',sdf.replace('.sdf','_'+str(bunch))+'.sdf'
                f = open(sdf.replace('.sdf','_'+str(bunch))+'.sdf','w')
                f.write(out)
                f.close()
                out = ''
                counter = 0
                bunch += 1
                break

if __name__ == '__main__':
    sdf = sys.argv[1]
    nCompound = int(sys.argv[2])
    read_sdf_file(sdf,nCompound)
