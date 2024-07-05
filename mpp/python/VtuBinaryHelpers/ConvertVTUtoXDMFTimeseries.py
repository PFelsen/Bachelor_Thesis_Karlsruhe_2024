import meshio
import argparse
from meshio.xdmf import TimeSeriesWriter
import os





if __name__ == '__main__':
    parser = argparse.ArgumentParser('Tool um mehrere VTUs in ein xdmf')
    parser.add_argument('-if','--inputfolder', type=str, help='folder where the data lies to run the script on ',required=True)
    parser.add_argument('-sf','--stringformatting',type=str, help='How the files are named with formatting which are about to be input in one file',required = True)
    parser.add_argument('-co','--counter',type=int, help='increment of the timeseries',required = True)
    parser.add_argument('-N','--timesteps',type=int, help='nr of timesteps',required = True)
    parser.add_argument('--compressionlevel',type=int, default=1, help='compressionlevel for gzip',required = False)
    parser.add_argument('-dt','--deltat',type=float, default=0.001, help='increment of the timeseries',required = False)
    parser.add_argument('-o','--output', type=str,default="result.xdmf", help = 'Name of the xdmf')
    parser.add_argument('-p','--processes',type=int,default=4, help='nr of processes to run on')
    args = parser.parse_args()
    print(meshio.__file__)
    def numpy_to_xml_string(self, data):
        if self.data_format == "XML":
            s = BytesIO()
            fmt = dtype_to_format_string[data.dtype.name]
            np.savetxt(s, data.flatten(), fmt)
            return s.getvalue().decode()
        elif self.data_format == "Binary":
            bin_filename = f"{self.filename.stem}{self.data_counter}.bin"
            self.data_counter += 1
            # write binary data to file
            with open(bin_filename, "wb") as f:
                data.tofile(f)
            return bin_filename

        if self.data_format != "HDF":
            raise WriteError()
        name = f"data{self.data_counter}"
        self.data_counter += 1
        #self.h5_file.create_dataset(name, data=data)
        self.h5_file.create_dataset(name, data=data,compression="gzip",compression_opts=args.compressionlevel)
        return os.path.basename(self.h5_filename) + ":/" + name

    TimeSeriesWriter.numpy_to_xml_string = numpy_to_xml_string
    with TimeSeriesWriter(args.output) as writer:
        for i in range(args.timesteps):
            fn = args.stringformatting.format(i*args.counter) + ".vtu"
            data= meshio.read(fn)
            if(i == 0):
                writer.write_points_cells(data.points, data.cells)
            else:
                data.cell_data.pop("Subdomain")
                data.cell_data.pop("ProcLoad")
            writer.write_data(i*args.deltat, cell_data=data.cell_data)
