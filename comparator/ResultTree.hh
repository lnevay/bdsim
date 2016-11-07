#ifndef COMPRESULTTREE_H
#define COMPRESULTTREE_H

#include "Result.hh"

#include <sstream>
#include <string>

/** 
 * @brief Result of comparing 2 TTrees. 
 * 
 * @author Stewart Boogert.
 */

class ResultTree: public Result
{ 
public:
  int t1NEntries;
  int t2NEntries;
  std::vector<std::string> offendingBranches;

  virtual std::string print() const
  {
    std::stringstream ss; 
    ss << Result::print();
    ss << "File (1/2) entries (" << t1NEntries << "/" << t2NEntries << ")\n";
    ss << "Offending branches: ";
    for (const auto& branchName : offendingBranches)
      {ss << branchName << " ";}
    ss << "\n";
    return ss.str();
  }
};

#endif
