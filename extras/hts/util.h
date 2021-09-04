#pragma once

#include <iostream>
#include <fstream>
//#include <vector>
//#include <unordered_set>
#include <string>

/*
template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &elems) {
  os << "[";
  for (auto iter = elems.begin(); iter != elems.end(); iter++) {
    os << *iter;
    if (std::next(iter) != elems.end()) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::unordered_set<T> &elems) {
  os << "{";
  for (auto iter = elems.begin(); iter != elems.end(); iter++) {
    os << *iter;
    if (std::next(iter) != elems.end()) {
      os << ", ";
    }
  }
  os << "}";
  return os;
}

template<typename ForwardIterator>
void print_range(std::ostream &os, ForwardIterator begin, const ForwardIterator &end) {
  os << "[";
  for (auto iter = begin; iter != end; iter++) {
    os << *iter;
    if (std::next(iter) != end) {
      os << ", ";
    }
  }
  os << "]";
}
*/

constexpr char path_seperator =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif

std::string path_join(const std::string &base, const std::string &ext) {
  if (base.empty()) return ext;
  if (ext.empty()) return base;
  if (base.back() == path_seperator || ext.front() == path_seperator) {
    return base + ext;
  } else {
    return base + path_seperator + ext;
  }
}

bool does_file_exist(const std::string &filename) {
  std::ifstream infile(filename);
  return infile.good();
}