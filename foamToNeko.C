/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Timofey Mukha
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamToNeko

Description
    Convert OpenFOAM mesh to nmsh format native to Neko.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hexCell.H"
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

class NekoVertex
{
    public:

    int id;
    double x;
    double y;
    double z;
   

    NekoVertex(int idv, double xv, double yv, double zv)
    {
        id = idv; x = xv; y = yv; z = zv;
    }
};

class NekoHex
{
    public:

    int id;
    std::vector<NekoVertex> vertices;


    NekoHex(int idv, std::vector<NekoVertex> v)
    {
       id = idv; vertices = v;

    }

    void write(std::ofstream & fs)
    {
        fs.write(reinterpret_cast<const char*>(&id), sizeof id);
        
        for(auto & vi : vertices)
        {
            fs.write(reinterpret_cast<const char*>(&vi.id), sizeof vi.id);
            fs.write(reinterpret_cast<const char*>(&vi.x), sizeof vi.x);
            fs.write(reinterpret_cast<const char*>(&vi.y), sizeof vi.y);
            fs.write(reinterpret_cast<const char*>(&vi.z), sizeof vi.z);
        }
    

    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    

    const auto & points = mesh.points();
    const cellList & cells = mesh.cells();
    const List<face> & faces = mesh.faces();

    const label nCells = mesh.nCells();
    const label gDim = 3;

    std::ofstream fs("mesh.nmsh", std::ios::out | std::ios::binary | std::ios::trunc);
    fs.write(reinterpret_cast<const char*>(&nCells), sizeof nCells);
    fs.write(reinterpret_cast<const char*>(&gDim), sizeof gDim);

    std::vector<NekoHex> nekoMesh;

    int elementID = 1;
    for (auto & cell : cells)
    {
        if (cell.nFaces() != 6)
        {
            Info << "Found a cells that does not have 6 faces, so not a hex!";
            FatalErrorIn("main()")
               << abort(FatalError);
        }

        const pointField cellPoints = cell.points(faces, points);
        const labelList cellLabels = cell.labels(faces);

//        FixedList<label,8> list(cellLabels);

//       hexCell hex(list);

        const face faceMinZ = faces[cell[4]];
        //const face faceMaxZ = faces[cell.opposingFaceLabel(cell[4], faces)];
        const face faceMaxZ = cell.opposingFace(cell[4], faces);

//        Info << cell << " " << faceMinZ << " " << faceMaxZ << cell.opposingFace(cell[4], faces) <<  nl;
        
        std::vector<NekoVertex> NekoVertices;


        // For the min z face we need to reorder from OpenFOAM to Neko
        // but the max z face is actually correct
        //std::vector<label> ind = {0, 3, 2, 1};
        std::vector<label> ind = {0, 1, 2, 3};

        for (int i=0; i<4; i++)
        {
            const point p = points[faceMinZ[ind[i]]];
            Info << p << nl;
            NekoVertices.push_back(NekoVertex(faceMinZ[ind[i]] + 1, p[0], p[1], p[2]));
        }
        for (int i=0; i<4; i++)
        {
            const point p = points[faceMaxZ[i]];
            Info << p << nl;
            NekoVertices.push_back(NekoVertex(faceMaxZ[i] + 1, p[0], p[1], p[2]));
        }

        nekoMesh.push_back(NekoHex(elementID, NekoVertices));

        elementID++;

    }

    HashTable<int, int> hash(8*nCells);

    int vertexID = 1;

    for (int i=0; i<nCells; i++)
    {
        for(int j=0; j<8; j++)
        {
            int & id = nekoMesh[i].vertices[j].id;

            if (hash.found(id))
            {
                //nekoMesh[i].vertices[j].id = hash.at(id);
                id = hash.at(id);
            }
            else
            {
                hash.insert(id, vertexID);
                nekoMesh[i].vertices[j].id = vertexID;
                vertexID++;
            }

        }
    }

    // Write the mesh
    for(auto & element : nekoMesh)
    {
        element.write(fs);
    }





    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
