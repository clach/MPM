#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>

template<class T, int dim>
class TetToPoly {
public:

    using TV = Eigen::Matrix<T,dim,1>;
    using Tet = Eigen::Matrix<int,4,1>;
    using Segment = Eigen::Matrix<int,2,1>;

    struct segment_hash {
        inline std::size_t operator()(const Segment& s) const {
            return s(0) * 31 + s(1);
        }
    };

    void convert(std::string filename, std::vector<TV>& inX, 
        std::vector<Segment>& inSegments) {

        std::cout << "Converting " << filename << "..." << std::endl;

        std::vector<TV> x;
        std::vector<Tet> tets;
        std::vector<Segment> segments;

        // read in tet file
        std::ifstream inputFile(filename);
        if (inputFile.is_open()) {
            std::string line;
            
            while (std::getline(inputFile, line)) {
                //std::cout << "line = " << line << std::endl;
                std::istringstream issLine(line);
                std::vector<std::string> tokens{std::istream_iterator<std::string>{issLine},
                      std::istream_iterator<std::string>{}};
                if (tokens.size() == 0) {
                    continue;
                }

                if (tokens[0].compare("POINTS") == 0) {
                    int numPoints = std::stoi(tokens[1]);
                    for (int i = 0; i < numPoints; i++) {
                        std::getline(inputFile, line);
                        std::istringstream issLine2(line);
                        std::vector<std::string> positions{std::istream_iterator<std::string>{issLine2},
                            std::istream_iterator<std::string>{}};
                        TV pos;
                        for (int j = 0; j < dim; j++) {
                            std::istringstream issPos(positions[j]);
                            T posComponent;
                            issPos >> posComponent;
                            pos(j) = posComponent;
                        }
                        //std::cout << "pos:     " << pos << std::endl;
                        x.push_back(pos);     
                    }
                } else if (tokens[0].compare("CELLS") == 0) {
                    int numCells = std::stoi(tokens[1]);
                    for (int i = 0; i < numCells; i++) {
                        std::getline(inputFile, line);
                        std::istringstream issLine2(line);
                        std::vector<std::string> tetIndices{std::istream_iterator<std::string>{issLine2},
                            std::istream_iterator<std::string>{}};

                        if (tetIndices.size() == 0) {
                            continue;
                        }

                        Tet tet;
                        for (int j = 1; j < 5; j++) { // want to get 4 vertices of tet
                            int index = std::stoi(tetIndices[j]);
                            tet(j - 1) = index;
                        }
                        tets.push_back(tet);     
                    }
                }
            }
        }
        
        std::unordered_set<Segment, segment_hash> hashSet;

        // for each tet
        for (size_t i = 0; i < tets.size(); i++) {
            // for each edge of each tet
            for (int j = 0; j < 4; j++) {
                for (int k = j + 1; k < 4; k++) {
                    // create sorted pair for each edge
                    int p = std::min(tets[i](j), tets[i](k));
                    int q = std::max(tets[i](j), tets[i](k));
                    Segment s = Segment(p, q);

                    if (hashSet.find(s) != hashSet.end()) {
                        continue;
                    }

                    hashSet.insert(s);
                }
            }
        }

        // add unique edges to segments
        for (Segment s : hashSet) {
            segments.push_back(s);
        }

        // modify input in place
        inX = x;
        inSegments = segments;
    }

};
