import os

if __name__ == "__main__":
    # Displaying the parent directory of the script
    libpath= os.path.dirname(__file__) + "/../../build/lib/"
    #print(libpath)

    files = os.listdir(libpath)
    for filename in files:
        d_index = filename.index("d.")
        new_filename = filename[:d_index] + filename[d_index+1:]
        print("Renamed", filename, "to", new_filename)
        os.rename(libpath+filename, libpath+new_filename)