import meshio
from os import listdir
from multiprocessing import Pool
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Python tool um Vtus vom ASCII format to Binary zu convertieren')
    parser.add_argument('-if','--inputfolder', type=str, default="", help='folder where the data lies to run the script on ',required=True)
    parser.add_argument('-of','--outputfolder',type=str, default="", help='folder where the data is saved to')
    parser.add_argument('-c','--clean', type=str, default="", help='fields to remove from the vtu',nargs='+')
    parser.add_argument('-p','--processes',type=int,default=4, help='nr of processes to run on')
    #parser.add_argument('--concat', type=str, default="", ehlp = ' possibility to concate')
    args = parser.parse_args()
    files = [f for f in listdir(args.inputfolder) if ".vtu" in f and "BIN" not in f]
    removeFields = args.clean
    outputFolder = args.outputfolder
    def workOn(file):
        mesh = meshio.read(args.inputfolder + file)
        for field in removeFields:
            if field in vtu.array_names:
                mesh.cell_data.pop(field)
        meshio.write(outputFolder+"BIN"+ file,mesh,binary=True)
    with Pool(args.processes) as p:
        p.map(workOn, files)
