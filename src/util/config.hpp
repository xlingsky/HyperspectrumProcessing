#ifndef XLINGSKY_CONFIG_HPP
#define XLINGSKY_CONFIG_HPP

#include <boost/property_tree/ptree.hpp>

#include "util/pathadaptor.hpp"

namespace xlingsky {
namespace util{

    class config : public boost::property_tree::ptree {
        private:
        boost::filesystemEx::pathadaptor _path_adaptor;
        typedef boost::property_tree::ptree base;
        typedef base::path_type path_type;
        public:
          config(const boost::property_tree::ptree &tree,
                 const boost::filesystem::path &path)
              : boost::property_tree::ptree(tree), _path_adaptor(path) {}
          boost::filesystem::path filepath(const path_type& path) const {
            return _path_adaptor.realpath(get<std::string>(path));
          }
          boost::filesystem::path filepath(const path_type& path, const std::string& default_value) const {
            return _path_adaptor.realpath(get(path, default_value));
          }
    };

};
};

#endif