# Geometry Files

By default, M++ reads geometry data from .geo files in the designated geometry folder.
These files are always structured the following way:

* HEADER: Global geometry information
* POINTS: Verteces of the geometry
* CELLS: Combining point indices to cells
* FACES: Boundary faces and information
* VDATA: VertexData - For each vertex in POINTS, there has to be a list of doubles containing data
* CDATA: CellData - For each cell in CELLS, there has to be a list of doubles containing data

## Default format

### Header

The global header may contain the following information:

* Cellformat [= Vtu] : Defines which cell format is used and affects the read of the CELLS part.
* Celltype [= Interval / Triangle / Quadrilateral / Tetrahedron / Hexahedron] : Replaces the specification of the
  celltype in CELLS.
* Subdomain [= any short] : Replaces the subdomain of all cells by this value.
* Boundary [= any short] : Replaces the boundary of all faces by this value.
* VertexData[ = [x y ...]] : Inserts the given data container for all vertices.
* CellData[ = [x y ...]] : Inserts the given data container for all cells.

### Points

Each row defines a new vertex. Vertices are defined by a set of doubles.

```0.0 1.0```

defines a 2D-Point

```0.0 1.0 2.0```

defines a 3D-Point

### Cells

Each row defines a new cell with the following structure

```CELLTYPE SUBDOMAIN INDEX_1 INDEX_2 ... INDEX_N```

* CELLTYPE has to be the correct VTU type:
    * Interval: 3
    * Triangle: 5
    * Quadrilateral: 9
    * Tetrahedron: 10
    * Hexahedron: 12
* SUBDOMAIN is problem specific
* INDEX_1, ... INDEX_N are the vertex indices of the corners

If 'Celltype' is defined in the header, CELLTYPE has to be omitted.
If 'Subdomain' is defined in the header, SUBDOMAIN has to be omitted.

### Faces

Each row defines a new boundary face with the following structure

```BOUNDARY INDEX_1 INDEX_2 ... INDEX_N```

* BOUNDARY problem specific boundary condition
* INDEX_1, ... INDEX_N are the vertex indices of the corners

If 'Boundary' is defined in the header, BOUNDARY has to be omitted.

### VData / CData

Each row defines a set of data for the corresponding vertex / cell. The indexing is the same as in POINTS / CELLS

```DATA_1 ... DATA_N```

## Legacy Format

### Header

The legacy format has no header

### Points

Each row defines a new vertex. Vertices are defined by a set of doubles.

```0.0 1.0```

defines a 2D-Point

```0.0 1.0 2.0```

defines a 3D-Point

### Cells

Each row defines a new cell with the following structure

```CORNERS SUBDOMAIN INDEX_1 INDEX_2 ... INDEX_N```

* CORNERS are the number of corners of the cell
* SUBDOMAIN is problem specific
* INDEX_1, ... INDEX_N are the vertex indices of the corners

### Faces

Each row defines a new boundary face with the following structure

```CORNERS BOUNDARY INDEX_1 INDEX_2 ... INDEX_N```

* CORNERS are the number of corners of the face
* BOUNDARY problem specific boundary condition
* INDEX_1, ... INDEX_N are the vertex indices of the corners

### VData / CData

Each row defines a set of data for the corresponding vertex / cell. The indexing is the same as in POINTS / CELLS

```N DATA_1 ... DATA_N```

* N are the number of data
* DATA_1, ... DATA_N are the vertex/cell data