#ifndef PATHADAPTOR_HPP
#define PATHADAPTOR_HPP

#include <boost/filesystem.hpp>

namespace boost{

namespace filesystemEx{

class pathadaptor{
 protected:
  boost::filesystem::path _dir;
 public:
  pathadaptor(boost::filesystem::path path){
    if(boost::filesystem::is_directory(path))
      _dir = path;
    else{
      _dir = path.parent_path();
    }
  }
  boost::filesystem::path absolutepath(boost::filesystem::path path){
    if(boost::filesystem::exists(path)) return path;
    boost::filesystem::path ret(_dir);
    ret /= path;
    if(boost::filesystem::exists(ret)) return ret;
    return boost::filesystem::path();
  }
};

};

};

#endif
