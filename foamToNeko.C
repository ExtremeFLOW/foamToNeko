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

#include "Time.H"
#include "fvMesh.H"
#include "columnFvMesh.H"
#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"
#include "cyclicPolyPatch.H"
 
using namespace Foam;

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include "nekoMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Try to find a face face withing the element. Returns face index in the element or -1 on fail
int find_hex_face(const NekoHex & element, const face & faceI, const pointField & points, bool debug=true)
{
    // Search each hex face
    for (int f=1; f<=6; f++)
    {
        //Info << "Face " << f << nl;
        std::vector<int> vind = element.face(f);

        // How many points we found in the face
        int found = 0;

        for(int pi=0; pi<4; pi++)
        {
            //Info << pi << nl;
            // One of the points was not found
            if (found < pi)
            {
                break;
            }

            //Info << "Found: " << found << nl;
            point pointI = points[faceI[pi]];
            for(int vi=0; vi<4; vi++)
            {
                const NekoVertex & vertexI = element.vertices[vind[vi] - 1];

                //vertexI.print();

                scalar dist = Foam::sqrt(sqr(pointI[0] - vertexI.x) +
                                         sqr(pointI[1] - vertexI.y) +
                                         sqr(pointI[2] - vertexI.z)
                                        );
                // Found the point
                if (dist < SMALL)
                {
                    found++;
                    //Info << "Found!" << nl;
                    break;
                }
                

            }
        }

        if (found == 4)
        {
            if (debug)
            {
                Info<< "Found face" << nl;
                Info<< "    " <<  element.vertices[vind[0] - 1].x << " "
                    << element.vertices[vind[0] - 1].y << " "
                    << element.vertices[vind[0] - 1].z
                    << nl;
                Info<< "    " <<  element.vertices[vind[1] - 1].x << " "
                    << element.vertices[vind[1] - 1].y << " "
                    << element.vertices[vind[1] - 1].z
                    << nl;
                Info<< "    " <<  element.vertices[vind[2] - 1].x << " "
                    << element.vertices[vind[2] - 1].y << " "
                    << element.vertices[vind[2] - 1].z
                    << nl;
                Info<< "    " <<  element.vertices[vind[3] - 1].x << " "
                    << element.vertices[vind[3] - 1].y << " "
                    << element.vertices[vind[3] - 1].z
                    << nl;
            }
            return f;
        }
    }
    return -1;
}

void print_face(const face & f, const pointField & p)
{
    Info << p[f[0]] << nl;
    Info << p[f[1]] << nl;
    Info << p[f[2]] << nl;
    Info << p[f[3]] << nl;

}


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

    Info << "Processing mesh cells" << nl;

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
            //Info << p << nl;
            NekoVertices.push_back(NekoVertex(faceMinZ[ind[i]] + 1, p[0], p[1], p[2]));
        }
        for (int i=0; i<4; i++)
        {
            const point p = points[faceMaxZ[i]];
//            Info << p << nl;
            NekoVertices.push_back(NekoVertex(faceMaxZ[i] + 1, p[0], p[1], p[2]));
        }

        nekoMesh.push_back(NekoHex(elementID, NekoVertices));

        elementID++;

    }

    HashTable<int, int> hash(nCells);

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

    // Write the mesh elements
    for(auto & element : nekoMesh)
    {
        element.write(fs);
    }


    Info << "Processing boundaries" << nl;


    // Boundary zones
    const fvBoundaryMesh& boundary = mesh.boundary();

    // Write the number of zone faces
    int nZones = faces.size() - mesh.faceNeighbour().size();
    fs.write(reinterpret_cast<const char*>(&nZones), sizeof nZones);

    Info << "Mesh has " << nZones << " boundary faces"  << nl;

    // We only need to go through one of the cyclic patches in a pair
    // to create the periodic zone. So here we will hold the ids of
    // patches that are already taken care of by going through the neigbour
    std::vector<int> processedCyclics;

    for (auto & patch : boundary)
    {

        Info<< "Boundary " <<  patch.name() << nl;
        const label patchID = boundary.findPatchID(patch.name());

        // Add neighbour patch
/*
        if (patch.type() == "cyclic")
        {
            if (std::binary_search(processedCyclics.begin(), processedCyclics.end(), patchID))
            {
                Info<< "Has already been processed by periodic neighbour"<<nl;
                continue;
            }

            const cyclicPolyPatch & cyclic =
                dynamic_cast<const cyclicPolyPatch &>(patch.patch());

            processedCyclics.push_back(cyclic.neighbPatchID());


        }
*/

        //indices of cells that own the patch faces
        const auto & faceCells = patch.faceCells();

        // We now need to find each patch face as a face index in a NekoHex
        for(int i=0; i<patch.size(); i++)
        {

            const NekoHex & element = nekoMesh[faceCells[i]];
            face faceI = faces[i + patch.start()];

            //print_face(faceI, points);

            int faceIndex = find_hex_face(element, faceI, points, false);

            // Face found
            if (faceIndex > 0)
            {


                // Periodic boundaries need special treatment
                if (patch.type() == "cyclic")
                {

                    const cyclicPolyPatch & cyclic =  dynamic_cast<const cyclicPolyPatch &>(patch.patch());
                    const labelUList & nbrCells = cyclic.nbrCells();
                    
                    // Neighbour element index, 1-based
                    int pE = nbrCells[i] + 1;

                    face nbrFace = faces[i + cyclic.neighbPatch().start()];
                    const NekoHex & nbrElement = nekoMesh[nbrCells[i]];

                    Info << "nbr face" << " " << faceI << nl;
                    //print_face(faceI, points);

                    int faceIndexNbr = find_hex_face(nbrElement, nbrFace, points, false);
                    Info << faceIndexNbr << nl;

                    std::array<int, 4> globalIds = {0, 0, 0, 0};

                    for(int ovi=0; ovi<4; ovi++)
                    {
                        const NekoVertex & oVert = element.vertices[element.face(faceIndex)[ovi] - 1];
                        Info <<"Original vertex" << nl;
                        oVert.print();


                        for(int vi=0; vi<4; vi++)
                        {
                            const NekoVertex & nbrVert = nbrElement.vertices[nbrElement.face(faceIndexNbr)[vi] - 1];

                            nbrVert.print();

                            scalar distx = mag(nbrVert.x - oVert.x);
                            scalar disty = mag(nbrVert.y - oVert.y);
                            scalar distz = mag(nbrVert.z - oVert.z);

                            // Found the point if at least 2 distances are 0
                            if ((distx < SMALL && disty < SMALL) ||
                                (distx < SMALL && distz < SMALL) ||
                                (disty < SMALL && distz < SMALL))
                            {
                                Info << "Found!" << nl;
                                nbrVert.print();

                                // We take the minumum of the two ids
                                // Note that since thesea are boundary points
                                // no other elements will hold them
                                globalIds[ovi] = std::min(oVert.id, nbrVert.id);

                                break;
                            }
                            

                        }
                    }

                    NekoZone zoneFace(faceCells[i] + 1, faceIndex, pE, faceIndexNbr, globalIds, 5);
                    zoneFace.write(fs);

                }
                else
                {
                    std::array<int, 4> dummy = {0, 0, 0, 0};
                    NekoZone zoneFace(faceCells[i] + 1, faceIndex, 0, patchID + 1, dummy, 7);
                    zoneFace.write(fs);
                }
            }
            else
            {
                FatalErrorIn("Could not find boundary face in the element.")
                   << abort(FatalError);
            }
            
        } // face loop
    } // patch loop



    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //


