import math
import numpy as np

testfile_name = "CellTestData"

'''Creates a numpy point'''


def npp(*values):
    return np.array(values)


def midpoint(*points):
    n = len(points)

    midp = np.zeros(3)
    for p in points:
        midp += p / n

    return midp


def p_to_str(point):
    return " " + '{0:.10f}'.format(point[0]) + \
           " " + '{0:.10f}'.format(point[1]) + \
           " " + '{0:.10f}'.format(point[2])


def create_int(is_degen, *corners):
    if len(corners) != 2:
        raise Exception

    t_str = "Interval " + str(1 if is_degen else 0)
    # Corners
    for p in corners:
        t_str += p_to_str(p)

    # Edges
    for p in corners:
        t_str += p_to_str(p)

    # Faces
    for p in corners:
        t_str += p_to_str(p)

    # Center
    t_str += p_to_str(midpoint(corners[0], corners[1])) + " \n"

    return t_str


def create_triangle(is_degen, *corners):
    if len(corners) != 3:
        raise Exception

    t_str = "Triangle " + str(1 if is_degen else 0)
    # Corners
    for p in corners:
        t_str += p_to_str(p)

    # Edges
    t_str += p_to_str(midpoint(corners[0], corners[1]))
    t_str += p_to_str(midpoint(corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[2], corners[0]))

    # Faces
    t_str += p_to_str(midpoint(corners[0], corners[1]))
    t_str += p_to_str(midpoint(corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[2], corners[0]))

    # Center
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[2])) + " \n"

    return t_str


def create_quad(is_degen, *corners):
    if len(corners) != 4:
        raise Exception

    t_str = "Quadrilateral " + str(1 if is_degen else 0)
    # Corners
    for p in corners:
        t_str += p_to_str(p)

    # Edges
    t_str += p_to_str(midpoint(corners[0], corners[1]))
    t_str += p_to_str(midpoint(corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[2], corners[3]))
    t_str += p_to_str(midpoint(corners[3], corners[0]))

    # Faces
    t_str += p_to_str(midpoint(corners[0], corners[1]))
    t_str += p_to_str(midpoint(corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[2], corners[3]))
    t_str += p_to_str(midpoint(corners[3], corners[0]))

    # Center
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[2], corners[3])) + " \n"

    return t_str


def create_tetra(is_degen, *corners):
    if len(corners) != 4:
        raise Exception

    t_str = "Tetrahedron " + str(1 if is_degen else 0)
    # Corners
    for p in corners:
        t_str += p_to_str(p)

    # Edges
    t_str += p_to_str(midpoint(corners[0], corners[1]))
    t_str += p_to_str(midpoint(corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[2], corners[0]))
    t_str += p_to_str(midpoint(corners[0], corners[3]))
    t_str += p_to_str(midpoint(corners[1], corners[3]))
    t_str += p_to_str(midpoint(corners[2], corners[3]))

    # Faces
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[3]))
    t_str += p_to_str(midpoint(corners[1], corners[2], corners[3]))
    t_str += p_to_str(midpoint(corners[2], corners[0], corners[3]))

    # Center
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[2], corners[3])) + " \n"

    return t_str


def create_hexa(is_degen, *corners):
    if len(corners) != 8:
        raise Exception

    t_str = "Hexahedron " + str(1 if is_degen else 0)
    # Corners
    for p in corners:
        t_str += p_to_str(p)

    # Edges
    t_str += p_to_str(midpoint(corners[0], corners[1]))
    t_str += p_to_str(midpoint(corners[1], corners[2]))
    t_str += p_to_str(midpoint(corners[2], corners[3]))
    t_str += p_to_str(midpoint(corners[3], corners[0]))

    t_str += p_to_str(midpoint(corners[0], corners[4]))
    t_str += p_to_str(midpoint(corners[1], corners[5]))
    t_str += p_to_str(midpoint(corners[2], corners[6]))
    t_str += p_to_str(midpoint(corners[3], corners[7]))

    t_str += p_to_str(midpoint(corners[4], corners[5]))
    t_str += p_to_str(midpoint(corners[5], corners[6]))
    t_str += p_to_str(midpoint(corners[6], corners[7]))
    t_str += p_to_str(midpoint(corners[7], corners[4]))

    # Faces
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[2], corners[3]))
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[5], corners[4]))
    t_str += p_to_str(midpoint(corners[1], corners[2], corners[6], corners[5]))
    t_str += p_to_str(midpoint(corners[2], corners[3], corners[7], corners[6]))
    t_str += p_to_str(midpoint(corners[3], corners[0], corners[4], corners[7]))
    t_str += p_to_str(midpoint(corners[4], corners[5], corners[6], corners[7]))

    # Center
    t_str += p_to_str(midpoint(corners[0], corners[1], corners[2], corners[3],
                               corners[4], corners[5], corners[6], corners[7])) + " \n"

    return t_str


def write_interval_tests(ofile):
    ofile.write(create_int(False, npp(0, 0, 0), npp(1, 0, 0)))
    ofile.write(create_int(False, npp(0, 0, 0), npp(-1, 0, 0)))
    ofile.write(create_int(True, npp(0, 0, 0), npp(0, 0, 0)))
    for i in range(1, 4):
        ofile.write(create_int(False, npp(-i * 4, 0, 0), npp(i * 4, 0, 0)))
        for j in range(1, 4):
            if i != j:
                ofile.write(create_int(False, npp(i * 4, 0, 0), npp(j * 4, 0, 0)))


def write_triangle_tests(ofile):
    ofile.write(create_triangle(False, npp(0, 0, 0), npp(1, 0, 0), npp(0, 1, 0)))
    ofile.write(create_triangle(False, npp(0, 0, 0), npp(0, 1, 0), npp(1, 0, 0)))
    ofile.write(create_triangle(True, npp(0, 0, 0), npp(1, 0, 0), npp(1, 0, 0)))
    for i in range(1, 4):
        for j in range(1, 4):
            ofile.write(create_triangle(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp(0, j * 4, 0)))
            ofile.write(create_triangle(False, npp(0, 0, 0), npp(0, i * 4, 0), npp(j * 4, 0, 0)))


def write_quad_tests(ofile):
    ofile.write(create_quad(False, npp(0, 0, 0), npp(1, 0, 0), npp(1, 1, 0), npp(0, 1, 0)))
    ofile.write(create_quad(False, npp(0, 0, 0), npp(0, 1, 0), npp(1, 1, 0), npp(1, 0, 0)))
    ofile.write(create_quad(True, npp(0, 0, 0), npp(1, 1, 0), npp(0, 1, 0), npp(1, 0, 0)))
    ofile.write(create_quad(True, npp(0, 0, 0), npp(1, 1, 0), npp(1, 0, 0), npp(0, 1, 0)))
    ofile.write(create_quad(True, npp(0, 0, 0), npp(0, 1, 0), npp(1, 0, 0), npp(1, 0, 0)))
    for i in range(1, 4):
        for j in range(1, 4):
            ofile.write(create_quad(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp(i * 4, j * 4, 0), npp(0, j * 4, 0)))
            ofile.write(create_quad(False, npp(0, 0, 0), npp(0, i * 4, 0), npp(i * 4, j * 4, 0), npp(j * 4, 0, 0)))
            #ofile.write(create_quad(True, npp(0, 0, 0), npp(i * 4, j * 4, 0), npp(0, j * 4, 0), npp(i * 4, 0, 0)))
            #ofile.write(create_quad(True, npp(0, 0, 0), npp(i * 4, j * 4, 0), npp(i * 4, 0, 0), npp(0, j * 4, 0)))

            ofile.write(
                create_quad(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp((i + j) * 4, (i + j) * 4, 0), npp(0, j * 4, 0)))


def write_tet_tests(ofile):
    ofile.write(create_tetra(False, npp(0, 0, 0), npp(1, 0, 0), npp(0, 1, 0), npp(0, 0, 1)))
    ofile.write(create_tetra(False, npp(0, 0, 0), npp(1, 0, 0), npp(0, 1, 0), npp(0, 0, -1)))
    ofile.write(create_tetra(True, npp(0, 0, 0), npp(1, 0, 0), npp(0, 1, 0), npp(0, 0, 0)))
    ofile.write(create_tetra(True, npp(0, 0, 0), npp(1, 0, 0), npp(0, 0, 0), npp(0, 0, 0)))
    ofile.write(create_tetra(True, npp(0, 0, 0), npp(0, 0, 0), npp(0, 0, 0), npp(0, 0, 0)))
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                ofile.write(create_tetra(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp(0, j * 4, 0), npp(0, 0, k * 4)))
                ofile.write(create_tetra(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp(0, j * 4, 0), npp(0, 0, -k * 4)))


def write_hexa_tests(ofile):
    ofile.write(create_hexa(False, npp(0, 0, 0), npp(1, 0, 0), npp(1, 1, 0), npp(0, 1, 0),
                            npp(0, 0, 1), npp(1, 0, 1), npp(1, 1, 1), npp(0, 1, 1)))
    ofile.write(create_hexa(False, npp(0, 0, 0), npp(0, 1, 0), npp(1, 1, 0), npp(1, 0, 0),
                            npp(0, 0, 1), npp(0, 1, 1), npp(1, 1, 1), npp(1, 0, 1)))
    ofile.write(create_hexa(True, npp(0, 0, 0), npp(1, 1, 0), npp(1, 0, 0), npp(0, 1, 0),
                            npp(0, 0, 1), npp(0, 1, 1), npp(1, 1, 1), npp(1, 0, 1)))
    ofile.write(create_hexa(True, npp(0, 0, 0), npp(1, 0, 0), npp(1, 1, 0), npp(0, 1, 0),
                            npp(0, 0, 1), npp(1, 1, 1), npp(0, 1, 1), npp(1, 0, 1)))
    ofile.write(create_hexa(True, npp(0, 0, 0), npp(1, 0, 1), npp(1, 1, 0), npp(0, 1, 1),
                            npp(0, 0, 1), npp(1, 0, 0), npp(1, 1, 1), npp(0, 1, 0)))

    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                ofile.write(create_hexa(False, npp(0, 0, 0), npp(0, j * 4, 0), npp(i * 4, j * 4, 0), npp(i * 4, 0, 0),
                                        npp(0, 0, k * 4), npp(0, j * 4, k * 4), npp(i * 4, j * 4, k * 4),
                                        npp(i * 4, 0, k * 4)))
                ofile.write(create_hexa(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp(i * 4, j * 4, 0), npp(0, j * 4, 0)
                                        , npp(0, 0, k * 4), npp(i * 4, 0, k * 4), npp(i * 4, j * 4, k * 4),
                                        npp(0, j * 4, k * 4)))
                ofile.write(create_hexa(False, npp(0, 0, 0), npp(i * 4, 0, 0), npp((i+j) * 4, (i+j) * 4, 0), npp(0, j * 4, 0)
                                        , npp(0, 0, k * 4), npp(i * 4, 0, k * 5), npp(i * 4, j * 4, k * 6),
                                        npp(0, j * 4, k * 5)))


if __name__ == "__main__":
    with open(testfile_name + ".dat", 'w') as output_file:
        write_interval_tests(output_file)
        write_triangle_tests(output_file)
        write_quad_tests(output_file)
        write_tet_tests(output_file)
        write_hexa_tests(output_file)
