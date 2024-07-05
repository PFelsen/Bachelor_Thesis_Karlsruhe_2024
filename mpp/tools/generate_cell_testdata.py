filename = "CellTestData.dat"
testdataPath = "../tests/cells/"

def point_to_vectorstring(pointdata):
    vectorstr = "{"
    for p in pointdata:
        vectorstr += p + ", "
    return vectorstr[:-2] + "}"

def points_to_vectorstring(points):
    N = int(len(points) / 3)
    vectorstr = "{"
    for i in range(N):
        vectorstr += point_to_vectorstring(points[i*3:(i+1)*3]) + ", "
    return vectorstr[:-2] + "}"



def write_test_data():
    print("Creating Test Data from file " + filename)

    cellData = open(filename, 'r').readlines()

    with open(str(testdataPath + "TestCellData.hpp"), 'w') as output_file:
        output_file.write("#ifndef TESTCELLDATA_H\n")
        output_file.write("#define TESTCELLDATA_H\n\n")

        output_file.write("#include \"geometry/Point.h\"\n")
        output_file.write("#include \"cells/Celltype.hpp\"\n\n")

        output_file.write("struct CellTestParameters{\n"+
            "\tCELLTYPE celltype;\n"+
            "\tbool isDegenerated;\n\n"+
            "\tstd::vector<Point> corners;\n"+
            "\tstd::vector<Point> edges;\n"+
            "\tstd::vector<Point> faces;\n"+
            "\tPoint center;\n"+
        "};\n\n")

        output_file.write("static std::vector <CellTestParameters> cellTestData{\n")

        first = True

        for data in cellData:
            data = data.split(" ")[:-1]
            if not first:
                output_file.write(",\n")

            type = data[0].upper()
            degen = "True" if data[1] == "1" else "False"
            output_file.write("\t{" + type + ", " +data[1]+", ")

            c, e, f = 0, 0, 0
            if type == "INTERVAL":
                c, e, f = 2, 2, 2
            elif type == "TRIANGLE":
                c, e, f = 3, 3, 3
            elif type == "QUADRILATERAL":
                c, e, f = 4, 4, 4
            elif type == "TETRAHEDRON":
                c, e, f = 4, 6, 4
            elif type == "HEXAHEDRON":
                c, e, f = 8, 12, 6
            else:
                continue

            m = max([c, e, f])

            output_file.write(points_to_vectorstring(data[2:2+(3*c)]) + ", ")
            output_file.write(points_to_vectorstring(data[2+(3*c):2+(3*(c+e))]) + ", ")
            output_file.write(points_to_vectorstring(data[2+(3*(c+e)):2+(3*(c+e+f))]) + ", ")
            output_file.write(point_to_vectorstring(data[2+(3*(c+e+f)):]) + "}")

            first = False

        output_file.write("\n};")


        output_file.write("\n#endif // TESTCELLDATA_H")


if __name__ == "__main__":
    write_test_data()