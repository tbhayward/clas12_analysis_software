#include "load_trees.h"
#include <iostream>

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    // Containers for different tree categories
    std::map<std::string, TTree*> dataTrees;
    std::map<std::string, TTree*> genMcTrees;
    std::map<std::string, TTree*> recMcTrees;
    std::map<std::string, TTree*> eppi0DataTrees;
    std::map<std::string, TTree*> eppi0GenMcTrees;
    std::map<std::string, TTree*> eppi0RecMcTrees;

    // Load all trees from files
    loadTrees(dataTrees,
              genMcTrees,
              recMcTrees,
              eppi0DataTrees,
              eppi0GenMcTrees,
              eppi0RecMcTrees);

    std::cout << "All trees loaded successfully." << std::endl;

    // TODO: Call further analysis steps here

    return 0;
}