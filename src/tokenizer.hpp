#ifndef TOKENIZER_HPP
#define TOKENIZER_HPP

#include "kmp.hpp"
#include <assert.h>

// #include <boost/multi_index_container.hpp>
// #include <boost/multi_index/identity.hpp>
// #include <boost/multi_index/ordered_index.hpp>
// #include <boost/multi_index/hashed_index.hpp>
// #include <boost/multi_index/sequenced_index.hpp>
// #include <boost/multi_index/random_access_index.hpp>
// #include <boost/multi_index/member.hpp>

namespace tokenizer{

template<typename T, typename IndexT=unsigned char>
struct TokenT{
  typedef IndexT IndexType;
  typedef T DataType;
  std::vector<DataType> key;
  std::vector<IndexType>  lps;
  TokenT(){}
  template<typename Iterator>
  TokenT(Iterator first, Iterator last){
    while(first!=last){
      key.push_back(*first);
      ++first;
    }
    build();
  }
  TokenT(const TokenT<T,IndexT>& t){
    key = t.key;
    lps = t.lps;
  }
  size_t size() const { return key.size(); }
  void build(){
    lps.resize(key.size());
    BuildPartialMatchTable(key.begin(), (IndexT)key.size(), lps.begin());
  }
};
template<typename T, typename IndexT=unsigned char>
struct SectionT{
  typedef IndexT IndexType;
  typedef T     DataType;
  IndexType id;
  const DataType* data;
  SectionT(IndexType id_, const DataType* data_)
      : id(id_), data(data_) {}
};

typedef TokenT<char> Token;
typedef SectionT<char> Section;

// typedef typename boost::multi_index_container<
//   Section,
//   boost::multi_index::indexed_by<
//     boost::multi_index::random_access<>,
//     boost::multi_index::hashed_non_unique<boost::multi_index::member<Section, Section::Index, &Section::id> >
//     >
//   > Sections;

template<class DataT, class IndexT>
class QInfoHelper{
 public:
  typedef DataT     			DataType;
  typedef IndexT    	IndexType;

  template<class TokenIterator>
  static DataType& key(TokenIterator t, int i) { return t->key[i]; }
  template<class TokenIterator>
  static size_t key_length(TokenIterator t) { return t->key.size(); }
  template<class TokenIterator>
  static IndexType    	lps(TokenIterator t, int i) { return t->lps[i]; }
  template<class SectionInserter, class DataRandomIterator>
  static void insert(SectionInserter inserter, DataRandomIterator it, int id){
    *inserter = Section(id, it);
    ++inserter;
  }
};


template<
  class DataRandomIterator, class TokenRandomIterator, class SectionInserter,
  class HelperType = QInfoHelper<typename std::iterator_traits<DataRandomIterator>::value_type, typename std::iterator_traits<TokenRandomIterator>::value_type::IndexType>
  >
size_t apply(
    DataRandomIterator data_first, DataRandomIterator data_last,
    TokenRandomIterator token_first, TokenRandomIterator token_last,
    SectionInserter section_inserter
    ){
  std::vector<typename HelperType::IndexType> ks(std::distance(token_first, token_last), 0);
  size_t count = 0;
  for(; data_first!=data_last; ++data_first){
    auto it_token = token_first;
    auto it_ks = ks.begin();
    for(; it_token!=token_last; ++it_token, ++it_ks){
      assert(HelperType::key_length(it_token)!=0);
      auto& k = *it_ks;
      while(k!=(typename HelperType::IndexType)-1 && *data_first!=HelperType::key(it_token, k)){
        k = HelperType::lps(it_token, k);
      }
      ++k;
      if(k==HelperType::key_length(it_token)){
        HelperType::insert( section_inserter, data_first+1, std::distance(token_first, it_token));
        k=0;
        ++count;
      }
    }
  }
  return count;
}

};


#endif
