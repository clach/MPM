#include <Eigen/Core>

#define TINYOBJLOADER_IMPLEMENTATION 
#include "tiny_obj_loader.h"

#include "mesh_query/mesh_query.h"

#include <sys/stat.h>
#include <iostream>

template<class T, int dim>
class PointCreation {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using iV = Eigen::Matrix<int,dim,1>;

    std::vector<TV> static createPoints(int numPointsPerCell, iV gridRes, TV min, TV max, bool jitter)
    {
        std::vector<TV> points;
        iV numPoints = numPointsPerCell * gridRes - iV(1, 1, 1);

        T pointsSpacingX = (max[0] - min[0]) / numPoints[0];
        T pointsSpacingY = (max[1] - min[1]) / numPoints[1];
        T pointsSpacingZ = (max[2] - min[2]) / numPoints[2];

        //std::cout << "numParticles1D " << numPoints[0] << std::endl;
        //std::cout << "particleSpacing " << pointsSpacingX << std::endl;

        for (float x = pointsSpacingX; x < pointsSpacingX * numPoints[0]; x += pointsSpacingY) 
        {
            for (float y = pointsSpacingY; y < pointsSpacingY * numPoints[1]; y += pointsSpacingY) 
            {
                for (float z = pointsSpacingZ; z < pointsSpacingZ * numPoints[2]; z += pointsSpacingZ) 
                {
                    TV point = TV(x, y, z);
                    if (jitter) 
                    {
                        point += TV(0, 0, 0);
                    }
                    points.push_back(point);
                }
            }
        }
        //std::cout << "startingPointsSize " << points.size() << std::endl;
        return points;
    }

    std::vector<TV> static selectInBox(std::vector<TV> points, TV min, TV max)
    {
        std::vector<TV> validPoints;
        for (size_t i = 0; i < points.size(); i++) 
        {
            TV p = points[i];
            if (p[0] > min[0] && p[0] < max[0] &&
                p[1] > min[1] && p[1] < max[1] &&
                p[2] > min[2] && p[2] < max[2]) 
                {
                    validPoints.push_back(p);
                }
        }

        return validPoints;
    }

    std::vector<TV> static selectInMesh(std::vector<TV> points, std::string filepath)
    {
        /*
        tinyobj::attrib_t attrib;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;

        std::string warn;
        std::string err;

        bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str());

        if (!warn.empty()) {
            std::cout << warn << std::endl;
        }

        if (!err.empty()) {
            std::cerr << err << std::endl;
        }

        if (!ret) {
            exit(1);
        }

        int numVertices = 0;
        int numTriangles = 0;
        std::vector<T> vertexPositions;
        std::vector<int> triangles;

        // Loop over shapes
        for (size_t s = 0; s < shapes.size(); s++) {
            // Loop over faces (polygon)

            //std::vector<float> positions = shapes[s].mesh.positions;
            std::vector<unsigned int> indices = shapes[s].mesh.indices;
            // assume triangles

        
            size_t index_offset = 0;
            for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
                int fv = shapes[s].mesh.num_face_vertices[f];

                // Loop over vertices in the face.
                for (size_t v = 0; v < fv; v++) {
                    // access to vertex
                    tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                    tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
                    tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
                    tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];

                    vertexPositions.push_back((T)vx);
                    vertexPositions.push_back((T)vy);
                    vertexPositions.push_back((T)vz);
                    numVertices++;
                }
                index_offset += fv;

                numTriangles++;
            }

            vertexPositions = positions;
            triangles = indices;
        }
        numVertices = vertexPositions.size() / 3;
        numTriangles = triangles.size() / 3;

*/

        std::vector<TV> validPoints;
        // for (size_t i = 0; i < points.size(); i++) 
        // {
        //     TV p = points[i];
        //     if (true) 
        //     {
        //         validPoints.push_back(p);
        //     }
        // }

        return validPoints;
    }
};