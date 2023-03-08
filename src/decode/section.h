#ifndef XLINGSKY_ISTREAM_SECTION_H
#define XLINGSKY_ISTREAM_SECTION_H

namespace xlingsky{
namespace istream{

class Section{
 public:
  class Iterator{
   public:
    void operator++(){
    }
    friend bool operator==(){
    }
  };
 protected:
 public:
  Section() {}
  virtual ~Section() {}
  virtual bool operator()(){
    for(auto it=begin(); it!= end(); ++it){
    }
  }
  Iterator begin(){
  }
  Iterator end(){
  }
};

class extractor{
 public:
  void run(){
  }
};

};
};

#endif
