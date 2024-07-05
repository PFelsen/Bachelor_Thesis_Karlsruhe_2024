filename = "TensorTestData.dat"
testdataPath = "../tests/"



def write_test_data():
    print("Creating Test Data from file " + filename)

    tensorData = open(filename, 'r').readlines()

    with open(str(testdataPath + "TestTensorData.hpp"), 'w') as output_file:
        output_file.write("#ifndef TESTTENSORDATA_H\n")
        output_file.write("#define TESTTENSORDATA_H\n\n")

        output_file.write("#include <array>\n")
        output_file.write("#include <vector>\n\n")

        output_file.write("struct TensorData{\n"+
            "\tint dimension;\n"+
            "\tstd::array < double, 9 > tensorValues;\n"+
            "\tstd::array < double, 9 > inverseValues;\n"+
        "};\n\n")

        output_file.write("static std::vector <TensorData> tensorTestData{\n")

        first = True

        for data in tensorData:
            data = data.split(" ")[:-1]
            if not first:
                output_file.write(",\n")
            output_file.write("\t{" + data[0] + ", {" +data[1])
            for t in data[2:10]:
                output_file.write(", " + t)
            output_file.write("}, {" +data[10])
            for t in data[11:]:
                output_file.write(", " + t)
            output_file.write("}}")
            first = False

        output_file.write("\n};")


        output_file.write("\n#endif // TESTTENSORDATA_H")


if __name__ == "__main__":
    write_test_data()