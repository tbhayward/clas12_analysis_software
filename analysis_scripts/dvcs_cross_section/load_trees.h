#ifndef LOAD_TREES_H
#define LOAD_TREES_H

#include <map>
#include <string>
class TTree;

// Load all required ROOT trees into the provided maps
void loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees);

#endif // LOAD_TREES_H