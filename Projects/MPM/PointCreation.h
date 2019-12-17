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



    /// Returns a pointer to a newly constructed mesh object, given
// the number of vertices, a pointer to the array of vertex coordinates
// (three doubles per vertex), the number of triangles, and a pointer to the
// array of triangle vertex indices (three ints per triangle).
// Note that the vertex and triangle data is not copied, so it is up to the
// user to preserve the arrays over the lifetime of the mesh object.
// MeshObject*
// construct_mesh_object(int num_vertices,
//                       const double *positions,
//                       int num_triangles,
//                       const int *triangles);

    std::vector<TV> static selectInMesh(std::vector<TV> points, std::string filepath)
    {

        tinyobj::attrib_t attrib;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;

        std::string warn;
        std::string err;

        // triangulate mesh
        bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str(), 0, true);

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
        std::vector<double> vertexPositions;
        std::vector<int> indices;

        // Loop over shapes
        for (size_t s = 0; s < shapes.size(); s++) {
            //std::vector<float> positions = shapes[s].mesh.positions;
            //std::vector<unsigned int> indices = shapes[s].mesh.indices;
        
            // Loop over faces (polygon)
            size_t index_offset = 0;
            for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
                int fv = shapes[s].mesh.num_face_vertices[f];

                // Loop over vertices in the face.
                for (size_t v = 0; v < fv; v++) {
                    // access to vertex
                    tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                    indices.push_back(idx.vertex_index);

                    tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
                    tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
                    tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
                    //positions.push_back(glm::vec3(vx, vy, vz));

                    vertexPositions.push_back((T)vx);
                    vertexPositions.push_back((T)vy);
                    vertexPositions.push_back((T)vz);
                    numVertices++;
                }
                index_offset += fv;

                numTriangles++;
            }
        }

        std::cout << "Num triangles: " << numTriangles << std::endl;
        std::cout << "Num vertices: " << numVertices << std::endl;

        std::cout << "index count: " << indices.size() / 3.f << std::endl;

        std::cout << "vertex posiiton count: " << vertexPositions.size() / 3.f << std::endl;

        MeshObject* mesh = construct_mesh_object(numVertices, vertexPositions.data(), numTriangles, indices.data());

        std::vector<TV> validPoints;
        for (size_t i = 0; i < points.size(); i++) 
        {
            TV p = points[i];
            double point[3] = { p(0), p(1), p(2) };
            if (point_inside_mesh(point, mesh))
            {
                validPoints.push_back(p);
            }
        }

        return validPoints;
    }
};